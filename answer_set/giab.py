"""GIAB regime: which of HG001-HG007 is the donor, and its germline truth.

Germline truth = the donor's NIST v4.2.1 benchmark calls in the round region, every record, no BED filter
(truth minus SYNTHETIC records == the donor's NIST records in the region, both ways: 65 released truths 09-23..27,
and all 24 stored GIAB rounds with BAMs).
Donor: genotype the round BAM at SNPs where the seven NIST call sets disagree and take the best-matching donor.
A donor's genotype at a site is known only inside its confident BED (absent there = hom-ref).
"""
from __future__ import annotations

import bisect
import random
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Tuple

from answer_set.fingerprint import genotype_sites
from answer_set.germline import canonical_gt  # noqa: F401  (re-exported: giab.canonical_gt)
from common.settings import SETTINGS
from common.types import Round, Variant

DONORS = ("HG001", "HG002", "HG003", "HG004", "HG005", "HG006", "HG007")
DECISIVE = 0.95


@dataclass(frozen=True)
class GiabMatch:
    donor: str
    score: float
    runner_up: str
    runner_up_score: float
    sites: int

    @property
    def decisive(self) -> bool:
        return self.score >= DECISIVE


def _sh(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout


def _query(donor: str, round: Round) -> List[Tuple[int, str, str, str]]:
    txt = _sh(f"bcftools query -r {round.region} -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' {SETTINGS.data_dir}/giab/{donor}.vcf.gz")
    return [tuple(line.split("\t")) for line in txt.splitlines()]


def nist_germline(round: Round, donor: str) -> List[Variant]:
    """The donor's NIST records in the region exactly as written (multi-allelic records stay one record)."""
    return [Variant(round.chrom, int(pos), ref, alt, canonical_gt(gt)) for pos, ref, alt, gt in _query(donor, round)]


def _confident(donor: str, round: Round) -> Tuple[List[int], List[int]]:
    starts, ends = [], []
    with open(f"{SETTINGS.data_dir}/giab/{donor}.bed") as f:
        for line in f:
            c, a, b = line.split("\t")[:3]
            if c == round.chrom and int(b) >= round.start and int(a) <= round.end:
                starts.append(int(a)); ends.append(int(b))
    return starts, ends


def discriminating_sites(round: Round, max_sites: int = 1500) -> Dict[int, tuple]:
    """pos -> (ref, alt, [alt-allele count per donor, -1 if unknown]) for biallelic SNPs where donors differ."""
    raw = {d: _query(d, round) for d in DONORS}
    calls = {d: {int(p): (r, a, g) for p, r, a, g in raw[d] if len(r) == 1 and len(a) == 1} for d in DONORS}
    other = {d: {int(p) for p, r, a, g in raw[d] if not (len(r) == 1 and len(a) == 1)} for d in DONORS}
    beds = {d: _confident(d, round) for d in DONORS}

    def count(d, pos, ref, alt):
        if pos in other[d]:                                # an indel / multi-allelic record there: not a clean SNP
            return -1
        c = calls[d].get(pos)
        if c is not None:
            if (c[0], c[1]) != (ref, alt):
                return -1
            g = c[2].replace("|", "/").split("/")
            return sum(x == "1" for x in g)
        s, e = beds[d]
        i = bisect.bisect_right(s, pos - 1) - 1
        return 0 if i >= 0 and s[i] < pos <= e[i] else -1

    sites = {}
    for pos, (ref, alt, _) in {p: v for d in DONORS for p, v in calls[d].items()}.items():
        per = [count(d, pos, ref, alt) for d in DONORS]
        known = {x for x in per if x >= 0}
        if len(known) > 1:
            sites[pos] = (ref, alt, per)
    keys = sorted(sites)
    random.Random(0).shuffle(keys)
    return {p: sites[p] for p in keys[:max_sites]}


def identify(round: Round) -> GiabMatch:
    sites = discriminating_sites(round)
    called = genotype_sites(round.bam, sites)
    scores = []
    for i, d in enumerate(DONORS):
        pairs = [(sites[p][2][i], g) for p, g in called.items() if sites[p][2][i] >= 0]
        scores.append((sum(a == b for a, b in pairs) / max(len(pairs), 1), d))
    scores.sort(reverse=True)
    return GiabMatch(scores[0][1], scores[0][0], scores[1][1], scores[1][0], len(called))
