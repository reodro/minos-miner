"""Synthetic (injected) variants: seen in the round's reads, absent from the donor's own reads (the parent).

1. candidates  bcftools mpileup+call on the round BAM (chunked, parallel), normalized and split to one allele per
               record; candidates equal to a germline truth allele are dropped
2. parent      alt fraction in the parent at the same site: ~0 -> injected, clearly >0 -> a real donor variant
               the platform left out of the truth (a caller would find it; the validator counts it as a FP)
3. genotype    from the round's allele fraction
Each kept candidate gets p = confidence it is injected (1.0 clear; lower when the parent is thin or ambiguous).
"""
from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Iterable, List, Optional, Set, Tuple

from answer_set.parent import Parent
from common.types import Round, SyntheticCandidate, Variant

MIN_ALT = 3                 # round reads supporting the candidate
MIN_VAF = 0.15
HOM_VAF = 0.85              # round allele fraction at or above this -> 1/1
PARENT_MIN_DEPTH = 8
PARENT_ABSENT = 0.03        # parent alt fraction at or below this -> injected
PARENT_PRESENT = 0.12       # at or above this -> real donor variant (dropped)


@dataclass(frozen=True)
class Candidate:
    pos: int
    ref: str
    alt: str
    alt_reads: int
    depth: int

    @property
    def vaf(self) -> float:
        return self.alt_reads / max(self.depth, 1)


def call_candidates(round: Round, chunk_bp: int = 500_000, workers: int = 14) -> List[Candidate]:
    spans = [(a, min(a + chunk_bp - 1, round.end)) for a in range(round.start, round.end + 1, chunk_bp)]

    def one(span):
        a, b = span
        cmd = (f"set -o pipefail; bcftools mpileup -f {round.reference} -r {round.chrom}:{a}-{b} -q 20 -Q 13 "
               f"-a AD,DP -d 1000 {round.bam} 2>/dev/null | bcftools call -mv 2>/dev/null | "
               f"bcftools norm -f {round.reference} -m -both 2>/dev/null | "
               f"bcftools query -f '%POS\\t%REF\\t%ALT[\\t%AD]\\n'")
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable="/bin/bash")
        if p.returncode != 0:                              # a lost chunk would silently drop its synthetic truth
            raise RuntimeError(f"candidate calling failed on {round.chrom}:{a}-{b} (rc={p.returncode})")
        out = []
        for line in p.stdout.splitlines():
            pos, ref, alt, ad = line.split("\t")
            if alt == "*" or "," in alt:
                continue
            ads = [int(x) for x in ad.split(",") if x != "."]
            if len(ads) < 2:
                continue
            out.append(Candidate(int(pos), ref, alt, ads[1], sum(ads)))
        return out

    with ThreadPoolExecutor(workers) as pool:
        cands = [c for chunk in pool.map(one, spans) for c in chunk]
    seen, uniq = set(), []
    for c in sorted(cands, key=lambda c: (c.pos, c.ref, c.alt)):
        if (c.pos, c.ref, c.alt) not in seen:
            seen.add((c.pos, c.ref, c.alt)); uniq.append(c)
    return uniq


def germline_keys(germline: Iterable[Variant]) -> Set[Tuple[int, str, str]]:
    """(pos, ref, alt) of every germline allele, multi-allelic records expanded and trimmed like bcftools norm."""
    keys = set()
    for v in germline:
        for alt in v.alt.split(","):
            r, a, p = v.ref, alt, v.pos
            while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:           # right-trim shared suffix
                r, a = r[:-1], a[:-1]
            while len(r) > 1 and len(a) > 1 and r[0] == a[0]:             # left-trim shared prefix
                r, a, p = r[1:], a[1:], p + 1
            keys.add((p, r, a))
    return keys


def classify(round: Round, cands: List[Candidate], parent: Optional[Parent], germline: Iterable[Variant]) -> Tuple[List[SyntheticCandidate], dict]:
    gk = germline_keys(germline)
    out, stats = [], {"candidates": len(cands), "germline": 0, "weak": 0, "donor_real": 0, "thin_parent": 0,
                      "ambiguous": 0, "injected": 0}
    for c in cands:
        if (c.pos, c.ref, c.alt) in gk:
            stats["germline"] += 1; continue
        if not round.start <= c.pos <= round.end:
            continue
        if c.alt_reads < MIN_ALT or c.vaf < MIN_VAF:
            stats["weak"] += 1; continue
        if parent is None:                           # no parent reads: every candidate is uncertain
            out.append(SyntheticCandidate(Variant(round.chrom, c.pos, c.ref, c.alt, "1/1" if c.vaf >= HOM_VAF else "0/1"), 0.5))
            stats["no_parent"] = stats.get("no_parent", 0) + 1
            continue
        pa, pn = parent.counts(c.pos, c.ref, c.alt)
        frac = pa / pn if pn else 0.0
        if pn and frac >= PARENT_PRESENT:
            stats["donor_real"] += 1; continue
        if pn < PARENT_MIN_DEPTH:
            p = 0.5; stats["thin_parent"] += 1
        elif frac <= PARENT_ABSENT:
            p = 1.0; stats["injected"] += 1
        else:
            p = 0.5; stats["ambiguous"] += 1
        gt = "1/1" if c.vaf >= HOM_VAF else "0/1"
        out.append(SyntheticCandidate(Variant(round.chrom, c.pos, c.ref, c.alt, gt), p))
    return out, stats
