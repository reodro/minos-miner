"""Identify the round's donor from its BAM: genotype common 1000G panel SNPs by pileup, pick the best-matching sample.

Measured on 36 1000G-regime rounds: correct donor scores >= 0.991, best wrong sample <= 0.776.
GIAB-regime rounds score <= ~0.69 against every 1000G sample (except HG001 = NA12878).
"""
from __future__ import annotations

import bisect
import gzip
import random
import struct
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Tuple

from common.settings import SETTINGS
from common.types import Round

DECISIVE = 0.95        # best-sample score needed to call a 1000G donor
_NT = "=ACMGRSVTWYHKDBN"


@dataclass(frozen=True)
class Fingerprint:
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


def panel_sites(round: Round, max_sites: int = 1500) -> Tuple[List[str], Dict[int, tuple]]:
    """Common biallelic SNPs in the region: pos -> (ref, alt, per-sample alt-allele count or -1)."""
    panel = f"{SETTINGS.kg_dir}/phased/{round.chrom}.vcf.gz"
    samples = _sh(f"bcftools query -l {panel}").split()
    txt = _sh(f"bcftools view -r {round.region} -m2 -M2 -v snps -i 'INFO/AF>=0.1 && INFO/AF<=0.9' {panel} "
              f"| bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' | head -n 20000")
    rows = [line.split("\t") for line in txt.splitlines()]
    random.Random(0).shuffle(rows)
    sites = {int(p[0]): (p[1], p[2], [-1 if g[0] == "." else (g[0] != "0") + (g[2] != "0") for g in p[3:]])
             for p in rows[:max_sites]}
    return samples, sites


def genotype_sites(bam: str, sites: Dict[int, tuple], min_depth: int = 8) -> Dict[int, int]:
    """Alt-allele count (0/1/2) at each site with a clear read pileup; ambiguous sites are skipped."""
    pos_sorted = sorted(sites)
    counts: Dict[int, List[int]] = {}
    f = gzip.open(bam, "rb")
    f.read(4)
    f.read(struct.unpack("<i", f.read(4))[0])
    for _ in range(struct.unpack("<i", f.read(4))[0]):
        f.read(struct.unpack("<i", f.read(4))[0] + 4)
    while True:
        h = f.read(4)
        if len(h) < 4:
            break
        b = f.read(struct.unpack("<i", h)[0])
        _, pos, lrn, mq, _, ncig, flag, lseq = struct.unpack("<iiBBHHHi", b[:20])
        if flag & (4 | 256 | 1024 | 2048) or mq < 20:
            continue
        cig = struct.unpack("<%dI" % ncig, b[32 + lrn:32 + lrn + ncig * 4])
        rlen = sum(c >> 4 for c in cig if (c & 15) in (0, 2, 3, 7, 8))
        i, j = bisect.bisect_left(pos_sorted, pos + 1), bisect.bisect_right(pos_sorted, pos + rlen)
        if i == j:
            continue
        seq = b[32 + lrn + ncig * 4:32 + lrn + ncig * 4 + (lseq + 1) // 2]
        targets, rp, qp, ti = pos_sorted[i:j], pos + 1, 0, 0
        for c in cig:
            op, ln = c & 15, c >> 4
            if op in (0, 7, 8):
                while ti < len(targets) and targets[ti] < rp + ln:
                    t = targets[ti]
                    if t >= rp:
                        k = qp + t - rp
                        base = _NT[(seq[k // 2] >> 4) if k % 2 == 0 else (seq[k // 2] & 15)]
                        ref, alt, _ = sites[t]
                        cnt = counts.setdefault(t, [0, 0])
                        if base == ref:
                            cnt[0] += 1
                        elif base == alt:
                            cnt[1] += 1
                    ti += 1
                rp += ln
                qp += ln
            elif op in (1, 4):
                qp += ln
            elif op in (2, 3):
                while ti < len(targets) and targets[ti] < rp + ln:
                    ti += 1
                rp += ln
    called = {}
    for p, (r, a) in counts.items():
        if r + a >= min_depth:
            fr = a / (r + a)
            g = 0 if fr < 0.1 else 1 if 0.25 <= fr <= 0.75 else 2 if fr > 0.9 else None
            if g is not None:
                called[p] = g
    return called


def best_match(samples: List[str], sites: Dict[int, tuple], called: Dict[int, int]) -> Fingerprint:
    scores = []
    for si, name in enumerate(samples):
        match = total = 0
        for p, g in called.items():
            sg = sites[p][2][si]
            if sg >= 0:
                total += 1
                match += sg == g
        scores.append((match / max(total, 1), name))
    scores.sort(reverse=True)
    return Fingerprint(scores[0][1], scores[0][0], scores[1][1], scores[1][0], len(called))


def identify(round: Round) -> Fingerprint:
    samples, sites = panel_sites(round)
    return best_match(samples, sites, genotype_sites(round.bam, sites))
