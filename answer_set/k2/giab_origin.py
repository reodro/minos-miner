"""The ORIGINAL reads a GIAB round was built from, fetched for the round's region only.

Every GIAB round BAM (all 7 people, checked 2026-09-30) consists of 128 bp pieces of reads from
Google's public NovaSeq PCR-free WGS BAMs (names hashed, qualities changed, some reads edited).
The 50x files match best. A whole file is 63-69 GB, so only the round region is fetched
(parallel pieces, retries), kept in <cache_dir>/answer_set/k2/giab_origin/ and removable after use.
"""
import concurrent.futures as cf
import time
from pathlib import Path

from answer_set.k2.util import CACHE as K2_CACHE

URL = ("https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/"
       "wgs_pcr_free/50x/{p}.novaseq.pcr-free.50x.dedup.grch38.bam")
CACHE = K2_CACHE / "giab_origin"


def path(person, chrom, a, b):
    return CACHE / f"{person}_{chrom}_{a}_{b}.bam"


def remote_index(person, attempts=6):
    """The remote BAM's index, downloaded once into the cache (htslib would otherwise save it in the current
    working directory on every first open)."""
    import urllib.request
    local = CACHE / f"{person}.novaseq.pcr-free.50x.dedup.grch38.bam.bai"
    if local.exists():
        return local
    CACHE.mkdir(parents=True, exist_ok=True)
    part = Path(str(local) + ".part")
    for t in range(attempts):
        try:
            urllib.request.urlretrieve(URL.format(p=person) + ".bai", part)
            part.rename(local)
            return local
        except OSError:
            part.unlink(missing_ok=True)
            time.sleep(5 * (t + 1))
    raise RuntimeError(f"index of the original reads for {person} failed after {attempts} attempts")


def fetch(person, chrom, a, b, pad=500, chunks=8, attempts=6):
    """Local BAM (indexed) with the person's original reads over [a-pad, b+pad]."""
    import pysam
    out = path(person, chrom, a, b)
    if out.exists() and Path(str(out) + ".bai").exists():
        return out
    bai = remote_index(person)
    lo, hi = max(1, a - pad), b + pad
    step = (hi - lo + chunks) // chunks
    pieces = [(lo + i * step, min(hi, lo + (i + 1) * step - 1)) for i in range(chunks) if lo + i * step <= hi]
    tmp = CACHE / f".tmp_{person}_{chrom}_{a}_{b}"
    tmp.mkdir(exist_ok=True)

    def one(i):
        x, y = pieces[i]
        dest = tmp / f"p{i}.bam"
        for t in range(attempts):
            try:
                src = pysam.AlignmentFile(URL.format(p=person), index_filename=str(bai))
                with pysam.AlignmentFile(str(dest), "wb", template=src) as o:
                    n = 0
                    for r in src.fetch(chrom, x - 1, y):
                        if r.reference_start >= x - 1:        # each read once (by start)
                            o.write(r)
                            n += 1
                return n
            except (OSError, ValueError):
                dest.unlink(missing_ok=True)
                time.sleep(5 * (t + 1))
        raise RuntimeError(f"original reads {person} {chrom}:{x}-{y} failed after {attempts} attempts")

    with cf.ThreadPoolExecutor(len(pieces)) as ex:
        counts = list(ex.map(one, range(len(pieces))))
    parts = [str(tmp / f"p{i}.bam") for i in range(len(pieces))]
    pysam.cat("-o", str(out) + ".part", *parts)
    Path(str(out) + ".part").rename(out)
    pysam.index(str(out))
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    if sum(counts) == 0:
        out.unlink()
        raise RuntimeError(f"no original reads for {person} {chrom}:{a}-{b}")
    return out
