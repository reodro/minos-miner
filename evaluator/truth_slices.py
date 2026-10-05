"""Truth slices built once per (truth file, region) instead of on every scoring. Owner: evaluator.

Upstream's HappyScorer.score_vcf slices the truth VCF to the round's region on every call: two bcftools containers
(view + index) per scoring. `install()` wraps upstream's `slice_truth_vcf` so the first call for a (truth file,
region) runs it into a cache and every later call copies the cached slice to the path upstream asked for. The bytes
are bcftools' own output for the same input and region, so the scoring that follows is unchanged.

Safe across threads and processes: a slice is built under a per-process lock into a unique temporary name and moved
into place atomically; two processes racing build the same file twice at worst.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import sys
import threading
import uuid
from pathlib import Path

from common.settings import SETTINGS

if SETTINGS.minos_subnet_dir not in sys.path:
    sys.path.append(SETTINGS.minos_subnet_dir)

CACHE_DIR = Path(SETTINGS.cache_dir) / "truth_slices"
_lock = threading.Lock()


def install() -> None:
    """Idempotent: route upstream's slice_truth_vcf through the cache (in this process)."""
    import utils.scoring as upstream
    if getattr(upstream.slice_truth_vcf, "_cached", False):
        return
    real = upstream.slice_truth_vcf

    def cached(source_vcf, target_vcf, region):
        source, target = Path(source_vcf).resolve(), Path(target_vcf).resolve()
        if not source.exists():
            return real(source_vcf, target_vcf, region)            # upstream reports the error
        st = source.stat()
        key = hashlib.sha1(f"{source}|{st.st_size}|{st.st_mtime_ns}|{region}".encode()).hexdigest()[:16]
        slice_vcf = CACHE_DIR / f"{key}.vcf.gz"
        with _lock:
            if not slice_vcf.exists():
                CACHE_DIR.mkdir(parents=True, exist_ok=True)
                tmp = CACHE_DIR / f"tmp_{uuid.uuid4().hex[:12]}.vcf.gz"
                if not real(source, tmp, region):
                    return False
                if Path(str(tmp) + ".csi").exists():
                    os.replace(str(tmp) + ".csi", str(slice_vcf) + ".csi")
                os.replace(tmp, slice_vcf)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(slice_vcf, target)
        if Path(str(slice_vcf) + ".csi").exists():
            shutil.copyfile(str(slice_vcf) + ".csi", str(target) + ".csi")
        return True

    cached._cached = True
    upstream.slice_truth_vcf = cached
