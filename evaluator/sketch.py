"""A small sketch of a callset, so the search can measure how different two callsets are without the VCFs.

One bottom-k MinHash per score class (common.types.CLASSES): every record is hashed (32 bits) on (chrom, pos, ref,
alt, unphased genotype) and the k smallest hashes of each class are kept, packed as base64 (about 4 KB a callset: it
travels with every fleet result). Two callsets' Jaccard similarity per class is
estimated from their sketches alone (search.distance), exactly when a class holds fewer than k records.

The record's class follows the score: an indel by its largest length change among its alleles, a SNP by its
genotype (homozygous ALT or not). Built on the machine that ran GATK (a fleet node keeps the VCF), so the sketch
travels with the CallSet.
"""
from __future__ import annotations

import base64
import gzip
import hashlib
import heapq
import struct
from typing import Dict, List

from common.types import CLASSES

K = 128


def _class_of(ref: str, alts: List[str], gt: str) -> str:
    change = max((abs(len(ref) - len(a)) for a in alts if a and a != "*"), default=0)
    if change >= 8:
        return "indel_8bp"
    if change >= 4:
        return "indel_4_7bp"
    if change >= 2:
        return "indel_2_3bp"
    if change == 1:
        return "indel_1bp"
    alleles = gt.replace("|", "/").split("/")
    hom = len(alleles) == 2 and alleles[0] == alleles[1] and alleles[0] not in ("0", ".")
    return "snp_hom" if hom else "snp_het"


def _gt(fmt: str, sample: str) -> str:
    keys = fmt.split(":")
    vals = sample.split(":")
    if "GT" not in keys or keys.index("GT") >= len(vals):
        return "."
    alleles = vals[keys.index("GT")].replace("|", "/").split("/")
    return "/".join(sorted(alleles))


def callset_sketch(vcf: str, k: int = K) -> Dict[str, object]:
    """{"k": k, "n": {class: records}, "max_alts": most ALT alleles in one record,
    "h": {class: base64 of the k smallest 32-bit hashes, ascending}}."""
    hashes: Dict[str, set] = {c: set() for c in CLASSES}
    max_alts = 0
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            alts = [a for a in p[4].split(",") if a and a not in ("*", "<NON_REF>")]
            if len(alts) > max_alts:
                max_alts = len(alts)
            gt = _gt(p[8], p[9]) if len(p) > 9 else "."
            cls = _class_of(p[3], p[4].split(","), gt)
            digest = hashlib.blake2b(f"{p[0]}:{p[1]}:{p[3]}:{p[4]}:{gt}".encode(), digest_size=4).digest()
            hashes[cls].add(int.from_bytes(digest, "big"))
    return {"k": k, "n": {c: len(h) for c, h in hashes.items()}, "max_alts": max_alts,
            "h": {c: pack(heapq.nsmallest(k, h)) for c, h in hashes.items()}}


def pack(values: List[int]) -> str:
    return base64.b64encode(struct.pack(f">{len(values)}I", *values)).decode()


def unpack(text: str) -> List[int]:
    raw = base64.b64decode(text)
    return list(struct.unpack(f">{len(raw) // 4}I", raw))
