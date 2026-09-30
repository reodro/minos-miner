"""Germline truth lookup.

1000G regime: the platform's germline truth is the donor's NYGC raw calls filtered by KG_RULE.
Fitted on 36 rounds (2026-09-26..28): 0 missing, 0 extra, 0 genotype disagreements.
"""
from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor
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


def kg_germline(round: Round, sample: str, rule: dict = KG_RULE, chunks: int = 8) -> List[Variant]:
    """The sample's raw calls in the region, filtered by the rule. Queried in parallel chunks (one query of a
    3202-sample file takes ~30 s); records overlapping two chunks are kept once, so the result equals one query."""
    raw = f"{SETTINGS.kg_dir}/raw/{round.chrom}.vcf.gz"
    step = (round.end - round.start) // chunks + 1
    spans = [(a, min(a + step - 1, round.end)) for a in range(round.start, round.end + 1, step)]

    def one(span):
        return subprocess.run(
            f"bcftools query -r {round.chrom}:{span[0]}-{span[1]} -s {sample} "
            f"-f '%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ]\\n' {raw}",
            shell=True, capture_output=True, text=True, check=True).stdout.splitlines()

    with ThreadPoolExecutor(len(spans)) as pool:
        lines = [line for part in pool.map(one, spans) for line in part]
    seen, uniq = set(), []
    for line in lines:
        k = line.split("\t", 3)[:3]
        if tuple(k) not in seen:
            seen.add(tuple(k)); uniq.append(line)
    return parse_query(round.chrom, uniq, rule)
