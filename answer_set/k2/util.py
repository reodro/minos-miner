"""Paths and shared helpers of the model-based builder. Owner: answer_set.

Inputs, under SETTINGS.data_dir (shareable between clones):
  reference/<chrom>/<chrom>.fa        reference genome (SETTINGS.reference_fasta)
  giab/<donor>.vcf.gz, <donor>.bed    GIAB NIST v4.2.1 calls and confident region (SETTINGS.giab_dir)
  1kg/phased/<chrom>.vcf.gz           1000 Genomes phased panel (SETTINGS.kg_dir; scripts/fetch_reference_data.sh)
  1kg/index/                          panel lookup index, index_<chrom>.npz (python -m answer_set.k2.kg --index all)
  1kg/gatk_sites/                     1000G raw GATK callset, sites + annotations, <chrom>.sites.vcf.gz (EBI
                                      20201028_3202_raw_GT_with_annot)
  1kg/masks/                          region lists used for 1000G rounds (1000G accessibility masks, GIAB v3.3
                                      stratifications)
Models: released ones in git (answer_set/k2/models, see tracks.py); trained ones in the cache.
Caches, under SETTINGS.cache_dir (this clone, rebuildable):
  answer_set/k2/<round tag>/v1/       one build per round
  answer_set/k2/models/, prepared/    trained model versions, prepared training rounds (tracks.py, train.py)
  answer_set/k2/inputs.json           manifest of the inputs above (inputs.py)
  answer_set/k2/giab_origin/, kg_gatk_person/, masks/   downloads and per-chromosome mask intervals
"""
from __future__ import annotations

import bisect
import os
import subprocess
import time
import uuid
from pathlib import Path

from common.settings import SETTINGS

KG_DIR = Path(SETTINGS.kg_dir)
CACHE = Path(SETTINGS.cache_dir) / "answer_set" / "k2"
BCFTOOLS_IMG = "quay.io/biocontainers/bcftools:1.20--h8b25389_0"    # the image the validator's scoring uses
BCFTOOLS_TIMEOUT_S = 900       # one bcftools call; a whole GIAB build takes ~1 min, its longest call less
DONORS = ["HG001", "HG002", "HG003", "HG004", "HG005", "HG006", "HG007"]


def ref_path(chrom):
    return Path(SETTINGS.reference_fasta(chrom))


def donor_vcf(d):
    return Path(SETTINGS.giab_dir) / f"{d}.vcf.gz"


def donor_bed(d):
    return Path(SETTINGS.giab_dir) / f"{d}.bed"


# A bug gives the same answer on every try: not retried. Everything else (network, Docker, partial downloads,
# timeouts) is worth another try.
PERMANENT_ERRORS = (TypeError, AttributeError, NameError, KeyError, IndexError, ZeroDivisionError, AssertionError,
                    NotImplementedError)


def retry(fn, attempts, wait_s, sleep=time.sleep, budget_s=None, clock=time.monotonic):
    """fn() up to `attempts` times, waiting wait_s x attempt number between tries. Returns (result, messages of the
    failed tries); if every try fails, raises RuntimeError listing them all. A permanent error (a bug) fails at
    once; with budget_s, no new try starts when the wait would end past budget_s from the first try."""
    errors, end = [], None if budget_s is None else clock() + budget_s
    for i in range(attempts):
        try:
            return fn(), errors
        except PERMANENT_ERRORS as e:
            errors.append(f"{type(e).__name__}: {str(e)[:300]}")
            raise RuntimeError("failed (a bug, not retried): " + " | ".join(errors)) from e
        except Exception as e:  # noqa: BLE001
            errors.append(f"{type(e).__name__}: {str(e)[:300]}")
            if i + 1 < attempts:
                if end is not None and clock() + wait_s * (i + 1) >= end:
                    raise RuntimeError(f"failed, out of time after {i + 1} attempts: " + " | ".join(errors)) from e
                sleep(wait_s * (i + 1))
    raise RuntimeError(f"failed after {attempts} attempts: " + " | ".join(errors))


def parse_region(region):
    chrom, span = region.split(":")
    a, b = span.replace(",", "").split("-")
    return chrom, int(a), int(b)


def bcftools(cmd, mounts, timeout_s=BCFTOOLS_TIMEOUT_S):
    """Run a bcftools shell pipeline in Docker. mounts: [(host_dir, container_dir), ...]. Past timeout_s the
    container is removed and RuntimeError raised (killing only the docker client would leave it running)."""
    name = f"k2-bcftools-{os.getpid()}-{uuid.uuid4().hex[:12]}"
    args = ["docker", "run", "--rm", "--name", name, "-u", f"{os.getuid()}:{os.getgid()}"]
    for host, cont in mounts:
        args += ["-v", f"{host}:{cont}"]
    args += [BCFTOOLS_IMG, "sh", "-c", cmd]
    try:
        r = subprocess.run(args, capture_output=True, text=True, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        try:
            subprocess.run(["docker", "rm", "-f", name], capture_output=True, timeout=60)
        except subprocess.TimeoutExpired:
            pass
        raise RuntimeError(f"bcftools timed out after {timeout_s} s: {cmd[:200]}") from None
    if r.returncode != 0:
        raise RuntimeError(f"bcftools failed: {r.stderr[-2000:]}")
    return r.stdout


def normalize_vcf(src, dst, ref, region=None):
    """Split multi-allelic records and left-align against the reference."""
    src, dst, ref = Path(src).resolve(), Path(dst).resolve(), Path(ref).resolve()
    reg = f"-r {region}" if region else ""
    idx = "" if Path(str(src) + ".tbi").exists() or Path(str(src) + ".csi").exists() else \
        f"bcftools index -t /i/{src.name} 2>/dev/null; "
    bcftools(f"{idx}bcftools view {reg} /i/{src.name} -Ou | bcftools norm -m -both -f /r/{ref.name} "
             f"-Oz -o /o/{dst.name} 2>/dev/null && bcftools index -f -t /o/{dst.name}",
             [(str(src.parent), "/i"), (str(ref.parent), "/r"), (str(dst.parent), "/o")])
    return dst


def trim(pos, ref, alt):
    """Minimal representation (drop shared suffix, then shared prefix)."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def read_alleles(vcf_path, region=None, keep_info=False):
    """Yield (key, genotype, info) for each ALT allele actually carried by the sample.

    key = (pos, ref, alt) in minimal representation. genotype = 'het' | 'hom'.
    """
    import pysam
    f = pysam.VariantFile(str(vcf_path))
    it = f.fetch(*_region_args(region)) if region else f.fetch()
    for rec in it:
        s = rec.samples[0] if len(rec.samples) else None
        gt = s.get("GT") if s is not None else None
        if not gt or any(g is None for g in gt):
            continue
        carried = {g for g in gt if g and g > 0}
        zyg = "hom" if len(set(gt)) == 1 else "het"
        for i in carried:
            alt = rec.alleles[i]
            if alt.startswith("<") or alt == "*":
                continue
            key = trim(rec.pos, rec.ref, alt)
            yield key, zyg, (dict(rec.info) if keep_info else None)


def _region_args(region):
    c, a, b = parse_region(region)
    return c, a - 1, b


class Bed:
    """Fast point-in-interval lookup for a BED file restricted to one chromosome."""

    def __init__(self, path, chrom):
        iv = []
        with open(path) as fh:
            for line in fh:
                c, a, b = line.split()[:3]
                if c == chrom:
                    iv.append((int(a), int(b)))
        iv.sort()
        self.starts = [a for a, _ in iv]
        self.iv = iv

    def contains(self, pos1):
        i = bisect.bisect_right(self.starts, pos1 - 1) - 1
        return i >= 0 and self.iv[i][0] < pos1 <= self.iv[i][1]

    def covered_bp(self, a, b):
        return sum(max(0, min(y, b) - max(x, a - 1)) for x, y in self.iv if y > a - 1 and x < b)


def donor_records(d, chrom, a, b):
    """{pos: (REF, ALTs, GT)} of the donor's GIAB records in the window, written the way Minos's truth writes them:
    NIST's own REF, ALT order and genotype order, phased '|' as '/'. The validator classes a two-allele record by
    the allele its genotype names first (2/1 vs 1/2), so the order matters for the score (2026-10-04: writing 1/2 for
    NIST's 2/1 moved ~14 records per round from the >= 8 bp class to shorter classes)."""
    import pysam
    out = {}
    f = pysam.VariantFile(str(donor_vcf(d)))
    for rec in f.fetch(chrom, max(0, a - 1), b):
        s = rec.samples[0]
        gt = s.get("GT")
        if not gt or any(g is None for g in gt):
            continue
        out[rec.pos] = (rec.ref, tuple(rec.alts or ()), "/".join(str(g) for g in gt))
    return out


def load_donor_alleles(chrom, a, b):
    """{donor: {key: zygosity}} for each GIAB donor in the window."""
    out = {}
    for d in DONORS:
        out[d] = {k: z for k, z, _ in read_alleles(donor_vcf(d), f"{chrom}:{a}-{b}")}
    return out
