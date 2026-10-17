"""The config space as numbers: features for the surrogate model, unit positions, valid assembly, random sampling.

Every tunable parameter maps to one feature in [0, 1]: numeric ranges linearly, or on a log scale when they span
two orders of magnitude or more (rates, heterozygosities); bools 0/1; enums by position. Parameters with a
single allowed value are fixed and not features, and so are the values in PINNED: settings that can never score,
and arguments GATK does not read.

GATK parameters sit in cost tiers. The names are copied from evaluator.gatk so search does not import evaluator.
A batched variant pays the most expensive stage it changes. Activity (~0.59) and assembly (~0.55) are only a
partial discount and add latency, so they count as expensive, with every parameter that forces its own GATK run.
PairHMM (~0.12) and genotyping (~0.01) are the cheap axes a family can sweep inside one run.
"""
from __future__ import annotations

import dataclasses
import math
import random
from typing import Dict, List, Optional

from common.params import ParamRanges, ParamSpec
from common.types import Config
from search.bamgroups import BQST, CONT, MBQ, MMQ, MRPAS, snap


# emit_ref_confidence: GVCF / BP_RESOLUTION make GATK write reference blocks; hap.py rejects that output on the
# validator too, so such a config scores nothing (checked on 20260928T092900Z with the validator-exact path).
# The two pruning floats are GATK's adaptive-pruning arguments. HaplotypeCaller leaves adaptive pruning off
# (useAdaptivePruning defaults false) and this space has no switch to turn it on, so the assembler builds
# LowWeightChainPruner(min_pruning) and never reads them. Holding them at the GATK defaults drops two dimensions
# that cannot change a call.
# recover_all_dangling_branches stays an axis: on 20261002T115300Z none of the 333 configs with it reached the top
# 200 of 7,067 true-scored, but on 20261002T104100Z a climb gained 0.04 by flipping it - the axis book decides per
# round, and a bool costs one move per refill.
PINNED = {"gatk": {"emit_ref_confidence": "NONE",
                   "adaptive_pruning_initial_error_rate": 0.001,
                   "pruning_lod_threshold": 2.302585}}
# assembly_region_padding below 5 is never proposed: the batch path's GATK dies in assembly for heads with padding
# 0-1 (6 of 6 plain runs on 20261002T104100Z; 5 and 13 pass), and on 20261002T115300Z none of the 1,322 configs with
# padding 0-4 reached the top 200 of 7,067 true-scored (best 85.67 against 85.84 at padding 13).
FLOORS = {"gatk": {"assembly_region_padding": 5}}

# Live /scoring/parameter-ranges is tighter than the validator GATK template (templates.tool_params).
# Local GATK scores the template bounds. Family A (PCR NONE/CONSERVATIVE × large assembly) peaked ~81.1
# on 20261010T044200Z against king 82.78. These are the remaining unclimbed values: template-legal,
# outside that published box. Adaptive-pruning floats stay PINNED: HaplotypeCaller does not read them
# without a switch, and that switch is not a published parameter.
VALIDATOR_EDGES = {
    "gatk": {
        "min_pruning": {"min": 1},
        "min_dangling_branch_length": {"min": 1},
        "active_probability_threshold": {"min": 0.0001},
        "max_num_haplotypes_in_population": {"max": 512},
        "contamination_fraction_to_filter": {"max": 0.5},
        "standard_min_confidence_threshold_for_calling": {"min": 0.0},
        "max_assembly_region_size": {"max": 1000},
        "max_reads_per_alignment_start": {"min": 0},
        "sample_ploidy": {"min": 1, "max": 10, "allowed": ()},
        "dont_use_soft_clipped_bases": {"allowed": ()},
    }
}

# Sparse rungs outside the published box. Dense ladders stay on the API range so certify does not
# grow by hundreds of neighbors (haplotypes 129–512, contamination 0.05–0.5).
EDGE_RUNGS = {
    "min_pruning": (1,),
    "min_dangling_branch_length": (1,),
    "active_probability_threshold": (0.0001, 0.0002, 0.0005),
    "max_num_haplotypes_in_population": (160, 192, 256, 384, 512),
    "contamination_fraction_to_filter": (0.08, 0.10, 0.15, 0.25, 0.5),
    "standard_min_confidence_threshold_for_calling": (0.0, 10.0, 20.0),
    "max_assembly_region_size": (800, 900, 1000),
    "max_reads_per_alignment_start": (0,),
    "sample_ploidy": (1, 3),
    "dont_use_soft_clipped_bases": (True,),
}

# Assemble and seed these; leave them out of names so setting-distance stays on the published axes.
SEED_ONLY = frozenset({"sample_ploidy", "dont_use_soft_clipped_bases"})


def bam_caps(bam: str, chrom: str, start: int, end: int, reads: int = 20000) -> Dict[str, object]:
    """Bounds the round's data sets: min_base_quality_score above the BAM's highest base quality keeps no base and
    calls nothing (409 zero-scored configs in one fleet run), so the cap is that maximum (sampled from `reads`
    reads of the region; 20k reads take ~0.1 s). {} when the BAM cannot be read (tests, a missing file)."""
    try:
        import pysam
        with pysam.AlignmentFile(bam) as f:
            top, n = 0, 0
            for r in f.fetch(chrom, max(0, start - 1), end):
                q = r.query_qualities
                if q is None:
                    continue
                top = max(top, max(q))
                n += 1
                if n >= reads:
                    break
    except Exception:  # noqa: BLE001 - no BAM, no index, no pysam: no caps
        return {}
    return {"min_base_quality_score": int(top)} if n else {}

# Copied from evaluator.gatk.BRANCH groups. Search must not import evaluator (CONTRIBUTING.md).
GENOTYPING_PARAMS = ("standard_min_confidence_threshold_for_calling", "indel_heterozygosity",
                     "max_alternate_alleles", "contamination_fraction_to_filter")
PAIRHMM_PARAMS = ("base_quality_score_threshold", "pcr_indel_model", "pair_hmm_gap_continuation_penalty",
                  "phred_scaled_global_read_mismapping_rate")
ASSEMBLY_PARAMS = ("min_pruning", "min_dangling_branch_length", "recover_all_dangling_branches",
                   "max_num_haplotypes_in_population", "adaptive_pruning_initial_error_rate",
                   "pruning_lod_threshold")
ACTIVITY_PARAMS = ("active_probability_threshold", "min_assembly_region_size", "max_assembly_region_size")
# Not branchable: each distinct value is its own GATK run.
FULL_RUN_PARAMS = ("min_base_quality_score", "min_mapping_quality_score", "assembly_region_padding",
                   "heterozygosity", "max_reads_per_alignment_start")

# Designed genotyping sweep, inside the platform ranges. One axis at a time so a family stays one GATK stage.
# 90 and 100 sit in front: a family is only a few slots, and the large-indel hole is extra calls,
# so the higher thresholds have to be the ones that actually get queued.
CALL_CONF_GRID = (90.0, 100.0, 80.0, 60.0, 45.0, 30.0)
CONTAMINATION_GRID = (0.0, 0.01, 0.02, 0.04)

class Space:
    def __init__(self, ranges: ParamRanges, tool: str = "gatk", pinned: Optional[Dict[str, object]] = None,
                 caps: Optional[Dict[str, object]] = None, bounds=None):
        """`caps`: upper bounds the round's data sets (bam_caps), on top of FLOORS; values beyond them are never
        proposed and are clamped when a config comes from outside (a seed, a portfolio member).
        `bounds`: a BamBounds from this round's BAM. The settings it groups are proposed only at their
        representatives, and a config from outside is snapped onto those representatives. Contamination fractions
        that remove the same number of reads at every observed depth share one ladder rung."""
        self.ranges, self.tool = ranges, tool
        self.bounds = bounds
        specs = ranges.tools[tool]
        self.floors = {n: f for n, f in FLOORS.get(tool, {}).items() if n in specs}
        merged = dict(bounds.caps) if bounds is not None else {}
        merged.update(caps or {})
        self.caps = {n: c for n, c in merged.items() if n in specs and specs[n].max is not None and c < specs[n].max}
        self.groups: Dict[str, tuple] = dict(bounds.groups) if bounds is not None else {}
        mmq = specs.get(MMQ)
        self._mmq_floor = int(mmq.min) if mmq is not None and mmq.min is not None else 0
        self.specs: Dict[str, ParamSpec] = {}
        for n, sp in specs.items():
            if n in self.floors:
                sp = dataclasses.replace(sp, min=max(sp.min, self.floors[n]))
            if n in self.caps:
                sp = dataclasses.replace(sp, max=max(sp.min, self.caps[n]))
            reps = self.groups.get(n)
            if reps:
                sp = dataclasses.replace(sp, allowed=tuple(reps))
                if len(reps) == 1:
                    sp = dataclasses.replace(sp, default=reps[0], min=reps[0], max=reps[0])
            self.specs[n] = sp
        self.pinned = dict(PINNED.get(tool, {}) if pinned is None else pinned)
        self.reductions = None                       # set by the engine: caps and flat thresholds learned from runs
        # Published API box, before validator-edge widening. Ladders stay dense here; EDGE_RUNGS add the rest.
        self.api = {n: {"min": s.min, "max": s.max, "allowed": s.allowed, "type": s.type}
                    for n, s in self.specs.items()}
        self._widen(tool)
        self.names: List[str] = [n for n, s in sorted(self.specs.items())
                                 if not self._fixed(s) and n not in self.pinned and n not in SEED_ONLY]

    def _widen(self, tool: str) -> None:
        """Open validator-template bounds the published API left closed, and keep client/submit in sync."""
        extra = VALIDATOR_EDGES.get(tool, {})
        published = self.ranges.tools.get(tool, {})
        for n, over in extra.items():
            s = self.specs.get(n)
            if s is None:
                continue
            # A BAM whose starts never reach the published min makes max_reads one run. 0 is that same run.
            if n == "max_reads_per_alignment_start" and self._fixed(s):
                continue
            kw: Dict[str, object] = {}
            if "min" in over and (s.min is None or over["min"] < s.min):
                kw["min"] = over["min"]
            if "max" in over and (s.max is None or over["max"] > s.max):
                kw["max"] = over["max"]
            if "allowed" in over:
                kw["allowed"] = tuple(over["allowed"])
            elif s.allowed:
                merged = list(s.allowed)
                for v in EDGE_RUNGS.get(n, ()):
                    if v not in merged:
                        merged.append(v)
                if len(merged) != len(s.allowed):
                    kw["allowed"] = tuple(merged)
            if not kw:
                continue
            new = dataclasses.replace(s, **kw)
            self.specs[n] = new
            if n in published:
                published[n] = new

    @staticmethod
    def _fixed(s: ParamSpec) -> bool:
        return (s.allowed and len(s.allowed) == 1) or (s.type in ("int", "float") and s.min == s.max)

    def _log(self, s: ParamSpec) -> bool:
        pub = self.api.get(s.name) if getattr(self, "api", None) else None
        if pub and pub["min"] is not None and pub["max"] is not None:
            lo, hi = float(pub["min"]), float(pub["max"])
        elif s.min is not None and s.max is not None:
            lo, hi = float(s.min), float(s.max)
        else:
            return False
        return s.type == "float" and lo > 0 and hi / lo >= 100 - 1e-9

    # ---- config <-> features -------------------------------------------------------------------------------
    def encode(self, config: Config) -> List[float]:
        return [self._to_unit(self.specs[n], config.params.get(n, self.specs[n].default)) for n in self.names]

    def _snap(self, name: str, value, mmq: int):
        if self.bounds is None or name not in self.specs:
            return value
        spec = self.specs[name]
        if spec.min is None or spec.max is None:
            return value
        lo, hi = (float(spec.min), float(spec.max)) if name == CONT else (int(spec.min), int(spec.max))
        return snap(self.bounds.facts, name, value, mmq=mmq, lo=lo, hi=hi)

    def _to_unit(self, s: ParamSpec, v) -> float:
        pub = self.api.get(s.name) if getattr(self, "api", None) else None
        allowed = pub["allowed"] if pub else s.allowed
        typ = pub["type"] if pub else s.type
        lo = pub["min"] if pub else s.min
        hi = pub["max"] if pub else s.max
        if allowed:
            if v not in allowed and s.name in self.groups:
                v = self._snap(s.name, v, self._mmq_floor)
            return allowed.index(v) / max(len(allowed) - 1, 1) if v in allowed else 0.0
        if typ == "bool":
            return 1.0 if v else 0.0
        lo, hi, x = float(lo), float(hi), float(v)
        if hi == lo:
            return 0.0
        if self._log(s):
            return (math.log(x) - math.log(lo)) / (math.log(hi) - math.log(lo))
        return (x - lo) / (hi - lo)

    def _from_unit(self, s: ParamSpec, u: float):
        u = min(1.0, max(0.0, u))
        if s.allowed:
            return s.allowed[min(int(u * len(s.allowed)), len(s.allowed) - 1)]
        if s.type == "bool":
            return u >= 0.5
        lo, hi = float(s.min), float(s.max)
        x = math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo))) if self._log(s) else lo + u * (hi - lo)
        return int(round(x)) if s.type == "int" else round(x, 6)

    def _same(self, spec: ParamSpec, a, b) -> bool:
        if spec.type == "float":
            return abs(float(a) - float(b)) <= 1e-6 * max(1.0, abs(float(a)), abs(float(b)))
        return a == b

    def _assemble(self, params: Dict[str, object]) -> Optional[Config]:
        """A full valid config from `params`, with pinned values and the two GATK cross-constraints repaired."""
        full = {n: params.get(n, self.pinned.get(n, s.default)) for n, s in self.specs.items()}
        full.update(self.pinned)
        for n, f in self.floors.items():
            if isinstance(full.get(n), (int, float)) and full[n] < f:
                full[n] = f
        for n, c in self.caps.items():
            if isinstance(full.get(n), (int, float)) and full[n] > c:
                full[n] = c
        lo, hi = full.get("min_assembly_region_size"), full.get("max_assembly_region_size")
        if isinstance(lo, (int, float)) and isinstance(hi, (int, float)) and lo > hi:
            hi_max = self.specs["max_assembly_region_size"].max
            if lo <= hi_max:
                full["max_assembly_region_size"] = int(lo)
            else:
                full["min_assembly_region_size"] = int(hi)
        bq = full.get(BQST)
        if isinstance(bq, (int, float)) and not isinstance(bq, bool) and bq < 6:
            full[BQST] = 6
        if self.bounds is not None:
            mmq = full.get(MMQ, self._mmq_floor)
            mmq = int(mmq) if isinstance(mmq, (int, float)) and not isinstance(mmq, bool) else self._mmq_floor
            for n in (MBQ, BQST, MRPAS, CONT):
                if n in full and n in self.specs:
                    full[n] = self._snap(n, full[n], mmq)
        if self.reductions is not None:
            full = self.reductions.adjust(full)
        config = Config(self.tool, full)
        return None if self.ranges.errors(config) else config

    # ---- sampling ----------------------------------------------------------------------------------------------
    def random(self, rng: random.Random) -> Config:
        """A uniform random config that passes the platform's ranges and cross-parameter constraints."""
        for _ in range(1000):
            params = {n: self._from_unit(s, rng.random()) if n in self.names
                      else self.pinned.get(n, s.allowed[0] if s.allowed else s.default) for n, s in self.specs.items()}
            c = self._assemble(params)
            if c is not None:
                return c
        raise RuntimeError("could not sample a valid config")
