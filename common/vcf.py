"""Minimal dependency-free VCF reading shared by all parts."""
from __future__ import annotations

import gzip
from typing import Iterator, Optional, Tuple

from common.types import Variant


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def normalize_gt(gt: str) -> str:
    return gt.replace("|", "/")


def records(path: str, region: Optional[Tuple[str, int, int]] = None) -> Iterator[Tuple[Variant, str]]:
    """Yield (Variant, INFO) per record, first sample's GT. Multi-allelic ALT kept as written."""
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            pos = int(p[1])
            if region and (p[0] != region[0] or not region[1] <= pos <= region[2]):
                continue
            gt = normalize_gt(p[9].split(":")[0]) if len(p) > 9 else "./."
            yield Variant(p[0], pos, p[3], p[4], gt), p[7]
