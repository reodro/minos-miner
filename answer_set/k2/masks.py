"""Public 'easy / hard to read' region lists, used as measurements for 1000 Genomes-person rounds.

Minos keeps a DRAGEN+GATK call almost always inside the 1000 Genomes strict mask (SNPs 99.9%,
indels 98.9%; 16 rounds) and nearly all of its unexplained drops are outside it.
Files (in <data_dir>/answer_set/masks/):
  20160622.allChr.mask.bed                         1000G strict accessibility mask (GRCh38)
  20160622.allChr.pilot_mask.bed                   1000G pilot accessibility mask
  GRCh38_notinalldifficultregions.bed.gz           GIAB v3.3 stratification: not in any difficult region
  GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz   GIAB v3.3: tandem repeats + homopolymers
Each is converted once per chromosome into a small .npz (sorted, merged intervals), cached in
<cache_dir>/answer_set/k2/masks/.
"""
import bisect
import gzip

import numpy as np

from answer_set.k2.util import AS_DIR, CACHE

MASK_DIR = AS_DIR / "masks"
NPZ_DIR = CACHE / "masks"
MASKS = {
    "in_strict_mask": "20160622.allChr.mask.bed",
    "in_pilot_mask": "20160622.allChr.pilot_mask.bed",
    "in_giab_easy": "GRCh38_notinalldifficultregions.bed.gz",
    "in_giab_repeat": "GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
}
_cache = {}


def _intervals(name, chrom):
    key = (name, chrom)
    if key in _cache:
        return _cache[key]
    src = MASK_DIR / MASKS[name]
    npz = NPZ_DIR / f"{name}_{chrom}.npz"
    if not npz.exists():
        if not src.exists():
            _cache[key] = None
            return None
        iv = []
        with (gzip.open if src.suffix == ".gz" else open)(src, "rt") as fh:
            for line in fh:
                c, s, e = line.split("\t", 3)[:3]
                if c == chrom:
                    iv.append((int(s), int(e)))
        iv.sort()
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        arr = np.array(merged, dtype=np.int64).reshape(-1, 2)
        NPZ_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(npz, starts=arr[:, 0], ends=arr[:, 1])
    z = np.load(npz)
    _cache[key] = (z["starts"].tolist(), z["ends"].tolist())
    return _cache[key]


def inside(name, chrom, pos, ref):
    """1.0 if the whole variant (1-based pos, REF span) lies in one interval, 0.0 if not, NaN if no file."""
    iv = _intervals(name, chrom)
    if iv is None:
        return float("nan")
    starts, ends = iv
    a, b = pos - 1, pos - 1 + len(ref)            # 0-based half-open span
    i = bisect.bisect_right(starts, a) - 1
    return float(i >= 0 and starts[i] <= a and ends[i] >= b)
