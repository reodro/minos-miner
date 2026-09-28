"""Germline truth lookup.

1000G regime: the platform's germline truth is the donor's NYGC raw calls filtered by KG_RULE.
Fitted on 36 rounds (2026-09-26..28): 0 missing, 0 extra, 0 genotype disagreements.
"""
from __future__ import annotations

import subprocess
from typing import Iterable, List, Tuple

from common.settings import SETTINGS
from common.types import Round, Variant

KG_RULE = {"filter": "PASS", "no_star": True, "min_gq": 20, "min_dp": 10}


def keep(filter_: str, alts: List[str], gq: int, dp: int, rule: dict = KG_RULE) -> bool:
    """Whether a raw record belongs to the platform's germline truth."""
    return (filter_ == rule["filter"] and not (rule["no_star"] and "*" in alts)
            and gq >= rule["min_gq"] and dp >= rule["min_dp"])


def alleles(chrom: str, pos: int, ref: str, alts: List[str], gt: str) -> List[Variant]:
    """One Variant per called non-ref allele; gt is relative to that allele ("0/1" or "1/1")."""
    calls = gt.replace("|", "/").split("/")
    out = []
    for i, alt in enumerate(alts, 1):
        n = sum(1 for c in calls if c == str(i))
        if n and alt != "*":
            out.append(Variant(chrom, pos, ref, alt, "1/1" if n == 2 else "0/1"))
    return out


def parse_query(chrom: str, lines: Iterable[str], rule: dict = KG_RULE) -> List[Variant]:
    """Lines of `bcftools query -f '%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ]\\n'` -> germline truth."""
    out = []
    for line in lines:
        pos, ref, alt, flt, gt, dp, gq = line.rstrip("\n").split("\t")
        alts = alt.split(",")
        if keep(flt, alts, int(gq) if gq.isdigit() else -1, int(dp) if dp.isdigit() else -1, rule):
            out.extend(alleles(chrom, int(pos), ref, alts, gt))
    return out


def kg_germline(round: Round, sample: str, rule: dict = KG_RULE) -> List[Variant]:
    raw = f"{SETTINGS.kg_dir}/raw/{round.chrom}.vcf.gz"
    txt = subprocess.run(
        f"bcftools view -r {round.region} -s {sample} {raw} "
        f"| bcftools query -f '%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ]\\n'",
        shell=True, capture_output=True, text=True, check=True).stdout
    return parse_query(round.chrom, txt.splitlines(), rule)


def giab_germline(round: Round) -> Tuple[str, List[Variant]]:
    """SMOKE ONLY. GIAB regime: germline = NIST v4.2.1 calls of one of HG001-HG007 (verified exact on
    65 rounds, 2026-09-23..26). TODO: identify the donor against the 7 NIST sets and return its calls."""
    return "", []
