"""The config space as numbers: features for the surrogate model, random sampling, and local moves.

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
from typing import Dict, List, Optional, Sequence, Tuple

from common.params import ParamRanges, ParamSpec
from common.types import Config


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

# First full runs, in the order the budget is spent. Each tuple is one expensive parameter of the GATK default
# and the values next to that default: the default is already a strong config, and a one-step move is the
# cheapest place to learn whether that axis is flat. Round-robin, so every axis is tried once before any axis
# is tried twice. Cheap parameters stay at defaults and are swept inside the run.
_NEAR_STARTS = (
    ("min_pruning", (3, 4)),
    ("min_dangling_branch_length", (2, 8)),
    ("recover_all_dangling_branches", (True,)),
    ("max_num_haplotypes_in_population", (64, 32)),
    ("active_probability_threshold", (0.001, 0.005)),
    ("min_assembly_region_size", (25, 100)),
    ("max_assembly_region_size", (200, 400)),
    ("assembly_region_padding", (50, 150)),
    ("min_base_quality_score", (15, 20)),
    ("min_mapping_quality_score", (10, 30)),
    ("max_reads_per_alignment_start", (25, 100)),
    ("heterozygosity", (0.0003, 0.003)),
)

# Two starts outside that neighborhood, so the first wave is not a single hill. One is stricter than the
# default (fewer, cleaner reads and haplotypes). One keeps more indel paths and more reads.
_FAR_STARTS = (
    {"min_pruning": 6, "min_base_quality_score": 20, "min_mapping_quality_score": 40,
     "max_reads_per_alignment_start": 25, "active_probability_threshold": 0.01,
     "assembly_region_padding": 50, "max_num_haplotypes_in_population": 32},
    {"recover_all_dangling_branches": True, "min_dangling_branch_length": 2,
     "max_reads_per_alignment_start": 150, "assembly_region_padding": 200,
     "min_assembly_region_size": 25, "max_assembly_region_size": 500,
     "active_probability_threshold": 0.001},
)


class Space:
    def __init__(self, ranges: ParamRanges, tool: str = "gatk", pinned: Optional[Dict[str, object]] = None,
                 caps: Optional[Dict[str, object]] = None):
        """`caps`: upper bounds the round's data sets (bam_caps), on top of FLOORS; values beyond them are never
        proposed and are clamped when a config comes from outside (a seed, a portfolio member)."""
        self.ranges, self.tool = ranges, tool
        specs = ranges.tools[tool]
        self.floors = {n: f for n, f in FLOORS.get(tool, {}).items() if n in specs}
        self.caps = {n: c for n, c in (caps or {}).items() if n in specs and specs[n].max is not None and c < specs[n].max}
        self.specs: Dict[str, ParamSpec] = {}
        for n, sp in specs.items():
            if n in self.floors:
                sp = dataclasses.replace(sp, min=max(sp.min, self.floors[n]))
            if n in self.caps:
                sp = dataclasses.replace(sp, max=max(sp.min, self.caps[n]))
            self.specs[n] = sp
        self.pinned = dict(PINNED.get(tool, {}) if pinned is None else pinned)
        self.names: List[str] = [n for n, s in sorted(self.specs.items()) if not self._fixed(s) and n not in self.pinned]
        cheap = set(GENOTYPING_PARAMS) | set(PAIRHMM_PARAMS)
        self.genotyping: List[str] = [n for n in self.names if n in GENOTYPING_PARAMS]
        self.pairhmm: List[str] = [n for n in self.names if n in PAIRHMM_PARAMS]
        # Anything not known to be cheap is expensive, including a parameter the platform adds later.
        self.expensive: List[str] = [n for n in self.names if n not in cheap]
        self.expensive_ints: List[str] = [n for n in self.expensive if self.specs[n].type == "int"]

    @staticmethod
    def _fixed(s: ParamSpec) -> bool:
        return (s.allowed and len(s.allowed) == 1) or (s.type in ("int", "float") and s.min == s.max)

    def _log(self, s: ParamSpec) -> bool:
        return s.type == "float" and s.min is not None and s.min > 0 and s.max / s.min >= 100 - 1e-9

    # ---- config <-> features -------------------------------------------------------------------------------
    def encode(self, config: Config) -> List[float]:
        return [self._to_unit(self.specs[n], config.params.get(n, self.specs[n].default)) for n in self.names]

    def _to_unit(self, s: ParamSpec, v) -> float:
        if s.allowed:
            return s.allowed.index(v) / max(len(s.allowed) - 1, 1) if v in s.allowed else 0.0
        if s.type == "bool":
            return 1.0 if v else 0.0
        lo, hi, x = float(s.min), float(s.max), float(v)
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
        bq = full.get("base_quality_score_threshold")
        if isinstance(bq, (int, float)) and bq < 6:
            full["base_quality_score_threshold"] = 6
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

    def initial_design(self, rng: random.Random, n: int) -> List[Config]:
        """n configs to run first. The GATK default, then one step off it on each expensive axis, then two
        starts outside that neighborhood. Cheap parameters stay at defaults.

        A spread across the cube spends the only full runs we can afford on configs that are nowhere near the
        one point already known to score. Further starts, once this list is used up, are small moves of the
        default on expensive parameters only.
        """
        if n <= 0:
            return []
        out: List[Config] = []
        seen = set()

        def add(params: Dict[str, object]) -> bool:
            c = self._assemble(params)
            if c is None or c.key in seen:
                return False
            seen.add(c.key)
            out.append(c)
            return True

        add({})
        near: List[Tuple[str, object]] = []
        waves = max((len(values) for _, values in _NEAR_STARTS), default=0)
        for wave in range(waves):
            for name, values in _NEAR_STARTS:
                if wave < len(values) and name in self.names:
                    near.append((name, values[wave]))
        n_far = 2 if n >= 8 else (1 if n >= 4 else 0)
        n_far = min(n_far, max(0, n - len(out)))
        for name, value in near:
            if len(out) >= n - n_far:
                break
            add({name: value})
        for params in _FAR_STARTS:
            if len(out) >= n:
                break
            if n_far <= 0:
                break
            if add(params):
                n_far -= 1
        anchor = out[0] if out else None
        while len(out) < n and anchor is not None:
            c = self.neighbour(anchor, rng, n_params=1, scale=0.08, pool=self.expensive)
            if c is None or c.key in seen:
                name = rng.choice(self.expensive) if self.expensive else None
                if name is None or not add({name: self._from_unit(self.specs[name], rng.random())}):
                    break
                continue
            seen.add(c.key)
            out.append(c)
        return out

    def with_cheap_defaults(self, config: Config) -> Optional[Config]:
        """Keep `config`'s expensive parameters; put genotyping and PairHMM back at defaults for this round to retune."""
        params = {n: config.params.get(n, s.default) for n, s in self.specs.items()}
        for n in self.genotyping + self.pairhmm:
            params[n] = self.specs[n].default
        return self._assemble(params)

    def clamp_frozen(self, config: Optional[Config], incumbent: Optional[Config],
                     frozen: Sequence[str]) -> Optional[Config]:
        """Pin low-importance parameters at the incumbent so later samples spend their moves on the rest."""
        if config is None or incumbent is None or not frozen:
            return config
        params = dict(config.params)
        for n in frozen:
            if n in incumbent.params:
                params[n] = incumbent.params[n]
        return self._assemble(params)

    def neighbour(self, config: Config, rng: random.Random, n_params: int = 2, scale: float = 0.15,
                  pool: Optional[Sequence[str]] = None) -> Optional[Config]:
        """Move 1..n_params parameters: numeric by a Gaussian step of `scale` in unit space, bool/enum to a random
        value (possibly the same). `pool` restricts which parameters may move. None if no valid, different
        neighbour was found."""
        names = self.names if pool is None else [n for n in pool if n in self.names]
        if not names:
            return None
        for _ in range(50):
            params = dict(config.params)
            for n in rng.sample(names, min(len(names), rng.randint(1, max(1, n_params)))):
                s = self.specs[n]
                if s.allowed or s.type == "bool":
                    params[n] = self._from_unit(s, rng.random())
                else:
                    params[n] = self._from_unit(s, self._to_unit(s, params.get(n, s.default)) + rng.gauss(0, scale))
            c = self._assemble(params)
            if c is not None and c.key != config.key:
                return c
        return None

    def _axis_variants(self, config: Config, axes: Sequence[Tuple[str, Sequence]]) -> List[Config]:
        """One-parameter variants of `config` along each (name, values) axis, skipping the current value."""
        out: List[Config] = []
        base = dict(config.params)
        for name, values in axes:
            spec = self.specs.get(name)
            if spec is None or name not in self.names:
                continue
            current = base.get(name, spec.default)
            for value in values:
                if spec.allowed and value not in spec.allowed:
                    continue
                if self._same(spec, current, value):
                    continue
                params = dict(base)
                params[name] = value
                c = self._assemble(params)
                if c is not None and c.key != config.key:
                    out.append(c)
        return out

    def genotyping_variants(self, config: Config) -> List[Config]:
        """Calling-threshold and contamination values. Each differs from `config` in one genotyping parameter."""
        axes = []
        if "standard_min_confidence_threshold_for_calling" in self.names:
            axes.append(("standard_min_confidence_threshold_for_calling", CALL_CONF_GRID))
        if "contamination_fraction_to_filter" in self.names:
            axes.append(("contamination_fraction_to_filter", CONTAMINATION_GRID))
        return self._axis_variants(config, axes)

    def pairhmm_variants(self, config: Config) -> List[Config]:
        """The other `pcr_indel_model` settings. One enum, four values, and it prices as a PairHMM variant."""
        spec = self.specs.get("pcr_indel_model")
        if spec is None or "pcr_indel_model" not in self.names:
            return []
        return self._axis_variants(config, [("pcr_indel_model", tuple(spec.allowed))])
