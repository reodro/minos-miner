"""Configs that provably make the same calls on a round: run one GATK job for all of them. Owner: evaluator.

Some parameters only ever reach GATK 4.5.0.0 through comparisons against values present in the round's reads. Two
settings that no read value can tell apart make identical runs, so a config is mapped to a canonical one (same
everything, these parameters moved to their class's representative) and the executor's dedup reuses the result.
Proven from the 4.5.0.0 source (and the batching patch, which keeps these call sites unchanged):

- max_reads_per_alignment_start (AssemblyRegionWalker -> PositionalDownsampler -> ReservoirDownsampler.submit): reads
  are grouped by assigned start; a group draws a random number and discards reads only once it holds more reads than
  the cap. When no start position has more reads than the cap, nothing is discarded and no random number is drawn
  (so the QD-jitter stream is untouched too). Every cap >= the most reads sharing a start is the same run.
- min_base_quality_score (mbq): compared as qual >= mbq (AbstractReadThreadingGraph, assembly), qual <= mbq is skipped
  (ReferenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny, activity) and qual <= mbq - 1 is clipped
  (ReadClipper.hardClipLowQualEnds via finalizeRegion). Pileup detection, the only other use, is off. The values
  compared are the reads' base qualities, or after overlapping-pair correction (FragmentUtils) min(qual, 20) or 0.
  Two cutoffs are the same run when no such value lies in the closed interval between them.
- base_quality_score_threshold (bqst): compared as q < bqst on q = min(qual, read MAPQ) (PairHMMLikelihoodCalculation
  Engine.capMinimumReadQualities), reads with MAPQ < min_mapping_quality_score having been filtered out. Two
  thresholds are the same run when no such q lies in [low, high). Values below 6 are left alone (GATK rejects them).
- contamination_fraction_to_filter: AlleleBiasedDownsamplingUtils.selectAlleleBiasedEvidence removes
  (int)(totalEvidence * fraction) reads. Two fractions are the same run when that count agrees at every pileup
  depth the BAM contains. A fraction the 0.0025 ladder does not represent is left as written.

Everything is conservative: the value sets come from every read in the BAM (a superset of what any config reads),
always include 0 and 20, and every mapping quality the filter can let through. Coverage is the raw pileup, a
superset of the reads that reach genotyping.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass
from typing import Dict, FrozenSet, Optional

from common.types import Config

OVERLAP_QUALS = (0, 20)        # FragmentUtils: conflicting mates -> 0; agreeing -> min(q, HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL)
MIN_USABLE_Q = 6               # QualityUtils.MIN_USABLE_Q_SCORE: GATK refuses a lower base_quality_score_threshold
MRPAS, MBQ, BQST, MMQ = ("max_reads_per_alignment_start", "min_base_quality_score", "base_quality_score_threshold",
                         "min_mapping_quality_score")
CONT = "contamination_fraction_to_filter"
CONT_STEP = 0.0025                 # the search ladder; same constant as search/bamgroups.py


@dataclass(frozen=True)
class RoundFacts:
    """What the round's BAM lets configs see: the only inputs to canonical()."""
    quals: FrozenSet[int]          # every base quality in the BAM
    mapqs: FrozenSet[int]          # every mapping quality in the BAM
    max_reads_per_start: int       # most reads sharing one assigned start position
    reads: int = 0
    depths: FrozenSet[int] = frozenset()   # pileup depths; empty leaves contamination alone

    @property
    def compared_quals(self) -> list:
        """Every value min_base_quality_score can be compared with (sorted)."""
        return sorted(set(self.quals) | set(OVERLAP_QUALS))


def scan(bam: str) -> RoundFacts:
    """One pass over every read of the BAM (mapped, unmapped, any region)."""
    import pysam
    seen_q = bytearray(256)
    chunk = bytearray()
    mapqs = set()
    best = run = 0
    last = None
    n = 0
    events = []
    with pysam.AlignmentFile(bam) as f:
        for r in f.fetch(until_eof=True):
            n += 1
            q = r.query_qualities
            if q is not None:
                chunk.extend(q)
                if len(chunk) > 1 << 22:
                    for v in set(chunk):
                        seen_q[v] = 1
                    chunk.clear()
            mapqs.add(r.mapping_quality)
            if r.reference_id < 0 or r.reference_start < 0:   # no assigned position: never downsampled
                continue
            pos = (r.reference_id, r.reference_start)
            run = run + 1 if pos == last else 1
            last = pos
            best = max(best, run)
            end = r.reference_end
            if end is None or end <= r.reference_start:
                continue
            events.append((r.reference_start, 1))
            events.append((end, -1))
    events.sort()
    depth = 0
    prev = None
    seen_depth = set()
    for pos, delta in events:
        if prev is not None and pos > prev and depth > 0:
            seen_depth.add(depth)
        depth += delta
        prev = pos
    for v in set(chunk):
        seen_q[v] = 1
    return RoundFacts(frozenset(i for i in range(256) if seen_q[i]), frozenset(mapqs), best, n, frozenset(seen_depth))


def _class_low(value: int, values: list, closed: bool) -> int:
    """The smallest setting equivalent to `value`: no element of `values` (sorted) may lie between them, in the
    closed interval when `closed` (mbq) or in [low, high) otherwise (bqst)."""
    if closed and value in values:
        return value                                    # a compared value itself: its class is just itself
    i = bisect.bisect_left(values, value)               # values[i-1] < value <= values[i]
    return values[i - 1] + 1 if i > 0 else value


def canonical(config: Config, facts: Optional[RoundFacts]) -> Config:
    """The config every provably identical config maps to on this round (itself when nothing applies)."""
    if facts is None or config.tool != "gatk":
        return config
    p: Dict = dict(config.params)
    changed = False
    v = p.get(MRPAS)
    if isinstance(v, int) and not isinstance(v, bool) and v > 0 and v >= facts.max_reads_per_start > 0:
        p[MRPAS] = facts.max_reads_per_start
    v = p.get(MBQ)
    if isinstance(v, int) and not isinstance(v, bool) and 0 <= v <= 127:
        p[MBQ] = _class_low(v, facts.compared_quals, closed=True)
    v, mmq = p.get(BQST), p.get(MMQ)
    if isinstance(v, int) and not isinstance(v, bool) and MIN_USABLE_Q <= v <= 127 and isinstance(mmq, int):
        capped = {min(q, m) for q in facts.compared_quals for m in facts.mapqs if m >= mmq}
        low = _class_low(v, sorted(capped), closed=False)
        p[BQST] = max(low, MIN_USABLE_Q)
    v = p.get(CONT)
    if facts.depths and isinstance(v, (int, float)) and not isinstance(v, bool) and 0.0 <= float(v) <= 0.05:
        ordered = tuple(sorted(d for d in facts.depths if d > 0))
        sig = tuple(int(d * float(v)) for d in ordered)
        n = int(round(0.05 / CONT_STEP))
        for i in range(n + 1):
            rep = round(i * CONT_STEP, 6)
            if tuple(int(d * rep) for d in ordered) == sig:
                p[CONT] = rep
                break
    for k in (MRPAS, MBQ, BQST, CONT):
        changed |= p.get(k) != config.params.get(k)
    return Config(config.tool, p) if changed else config
