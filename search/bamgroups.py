"""Values of three GATK settings the round's BAM cannot tell apart. Owner: search.

The rules are the ones in evaluator/equivalence.py, copied here so search does not import evaluator.
min_base_quality_score and base_quality_score_threshold only ever get compared with base qualities (and, for the
second, mapping qualities) that occur in the BAM. max_reads_per_alignment_start does nothing once it is at least
the most reads sharing one start. Two settings in the same group are one GATK run.

base_quality_score_threshold's groups depend on min_mapping_quality_score. The list the search proposes is the
finest one, taken at the lowest mapping quality the platform allows, so a later move of that setting is not
collapsed into a coarser group. A config is still snapped to the group its own mapping quality defines.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass
from typing import Dict, FrozenSet, Optional, Tuple

from common.params import ParamRanges

MBQ = "min_base_quality_score"
BQST = "base_quality_score_threshold"
MRPAS = "max_reads_per_alignment_start"
MMQ = "min_mapping_quality_score"
OVERLAP_QUALS = (0, 20)          # FragmentUtils: conflicting mates -> 0; agreeing -> min(qual, 20)
MIN_USABLE_Q = 6                 # QualityUtils.MIN_USABLE_Q_SCORE: GATK refuses a lower base_quality_score_threshold


@dataclass(frozen=True)
class BamFacts:
    quals: FrozenSet[int]
    mapqs: FrozenSet[int]
    max_reads_per_start: int

    @property
    def compared_quals(self) -> list:
        return sorted(set(self.quals) | set(OVERLAP_QUALS))


@dataclass(frozen=True)
class BamBounds:
    """`caps` clamp a range. `groups` replaces a range with the representatives the search proposes.
    `count` is the number of distinct values the search is left with, against `full` published integers."""
    caps: Dict[str, int]
    groups: Dict[str, Tuple[int, ...]]
    full: Dict[str, int]
    count: Dict[str, int]
    facts: BamFacts

    @property
    def report(self) -> str:
        bits = []
        factor = 1.0
        for name in (MBQ, BQST, MRPAS):
            published, left = self.full.get(name), self.count.get(name)
            if not published or not left:
                continue
            factor *= published / left
            bits.append(f"{name} {published} -> {left}")
        return ", ".join(bits) + f"; these three axes are 1/{factor:.0f} of the integer grid"


def _class_low(value: int, values: list, closed: bool) -> int:
    """Smallest setting equivalent to `value`. Closed interval for min_base_quality_score, half-open for the
    PairHMM threshold. Same function as evaluator.equivalence._class_low."""
    if closed and value in values:
        return value
    i = bisect.bisect_left(values, value)
    return values[i - 1] + 1 if i > 0 else value


def _span(spec) -> Optional[Tuple[int, int, int]]:
    if spec is None or spec.type != "int" or spec.min is None or spec.max is None:
        return None
    lo, hi = int(spec.min), int(spec.max)
    return lo, hi, hi - lo + 1


def mbq_groups(compared: list, lo: int, hi: int) -> Tuple[int, ...]:
    seen = set()
    reps = []
    for v in range(lo, hi + 1):
        r = _class_low(v, compared, True)
        if r in seen:
            continue
        seen.add(r)
        reps.append(min(x for x in range(lo, hi + 1) if _class_low(x, compared, True) == r))
    return tuple(reps)


def bqst_values(facts: BamFacts, mmq: int) -> list:
    return sorted({min(q, m) for q in facts.compared_quals for m in facts.mapqs if m >= mmq})


def bqst_groups(facts: BamFacts, mmq: int, lo: int, hi: int) -> Tuple[int, ...]:
    capped = bqst_values(facts, mmq)
    lo = max(lo, MIN_USABLE_Q)
    seen = set()
    reps = []
    for v in range(lo, hi + 1):
        r = min(hi, max(_class_low(v, capped, False), lo))
        if r in seen:
            continue
        seen.add(r)
        reps.append(r)
    return tuple(reps)


def snap(facts: BamFacts, name: str, value, *, mmq: int, lo: int, hi: int):
    """The representative of `value` inside [lo, hi] on this BAM."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return value
    v = int(value)
    if name == MBQ:
        v = min(max(v, lo), hi)
        r = _class_low(v, facts.compared_quals, True)
        for rep in mbq_groups(facts.compared_quals, lo, hi):
            if _class_low(rep, facts.compared_quals, True) == r:
                return rep
        return v
    if name == BQST:
        lo = max(lo, MIN_USABLE_Q)
        v = min(max(v, lo), hi)
        capped = bqst_values(facts, mmq)
        r = min(hi, max(_class_low(v, capped, False), lo))
        reps = bqst_groups(facts, mmq, lo, hi)
        if r in reps:
            return r
        for rep in reps:
            if min(hi, max(_class_low(rep, capped, False), lo)) == r:
                return rep
        return v                                          # unknown class: keep the requested value
    if name == MRPAS:
        cap = facts.max_reads_per_start
        if cap <= lo:
            return lo
        if v >= cap:
            return min(cap, hi)
        return min(max(v, lo), hi)
    return value


def bounds_from(facts: BamFacts, ranges: ParamRanges, tool: str = "gatk") -> BamBounds:
    """Groups and caps for one round. A range that does not shrink is left alone."""
    specs = ranges.tools[tool]
    caps: Dict[str, int] = {}
    groups: Dict[str, Tuple[int, ...]] = {}
    full: Dict[str, int] = {}
    count: Dict[str, int] = {}

    span = _span(specs.get(MBQ))
    if span and facts.quals:
        lo, hi, n = span
        top = max(facts.quals)
        if top < hi:
            caps[MBQ] = top
            hi = max(lo, top)
        reps = mbq_groups(facts.compared_quals, lo, hi)
        full[MBQ], count[MBQ] = n, len(reps)
        if len(reps) < n:
            groups[MBQ] = reps

    span = _span(specs.get(BQST))
    mmq_spec = specs.get(MMQ)
    if span and facts.mapqs and mmq_spec is not None and mmq_spec.min is not None:
        lo, hi, _n = span
        usable = hi - max(lo, MIN_USABLE_Q) + 1     # GATK rejects a threshold below 6
        reps = bqst_groups(facts, int(mmq_spec.min), lo, hi)
        full[BQST], count[BQST] = usable, len(reps)
        if 1 < len(reps) < usable:
            groups[BQST] = reps

    span = _span(specs.get(MRPAS))
    if span and facts.max_reads_per_start > 0:
        lo, hi, n = span
        cap = facts.max_reads_per_start
        full[MRPAS] = n
        if cap <= lo:
            groups[MRPAS] = (lo,)
            count[MRPAS] = 1
        elif cap < hi:
            caps[MRPAS] = cap
            count[MRPAS] = cap - lo + 1
        else:
            count[MRPAS] = n

    return BamBounds(caps, groups, full, count, facts)


def bam_bounds(bam: str, ranges: ParamRanges, tool: str = "gatk") -> Optional[BamBounds]:
    """One pass over the BAM. None when it cannot be read: the search keeps the published ranges."""
    import os
    if not bam or not os.path.isfile(bam):
        return None
    try:
        import pysam
        seen_q = bytearray(256)
        chunk = bytearray()
        mapqs = set()
        best = run = 0
        last = None
        with pysam.AlignmentFile(bam) as handle:
            for rec in handle.fetch(until_eof=True):
                q = rec.query_qualities
                if q is not None:
                    chunk.extend(q)
                    if len(chunk) > 1 << 22:
                        for v in set(chunk):
                            seen_q[v] = 1
                        chunk.clear()
                mapqs.add(rec.mapping_quality)
                if rec.reference_id < 0 or rec.reference_start < 0:
                    continue
                pos = (rec.reference_id, rec.reference_start)
                run = run + 1 if pos == last else 1
                last = pos
                best = max(best, run)
        for v in set(chunk):
            seen_q[v] = 1
        quals = frozenset(i for i in range(256) if seen_q[i])
        if not quals:
            return None
        return bounds_from(BamFacts(quals, frozenset(mapqs), best), ranges, tool)
    except Exception:  # noqa: BLE001 - no BAM, no index, no pysam: published ranges
        return None
