"""The ORIGINAL reads a GIAB round was built from, fetched for the round's region only.

Every GIAB round BAM (all 7 people, checked 2026-09-30) consists of 128 bp pieces of reads from
Google's public NovaSeq PCR-free WGS BAMs (names hashed, qualities changed, some reads edited).
The 50x files match best. A whole file is 63-69 GB, so only the round region is fetched
(parallel pieces, retries), kept in <cache_dir>/answer_set/k2/giab_origin/ and removable after use.

Nothing here may wait forever: main builds the answer set before it can submit anything. htslib's remote reads have
no timeout and a hung thread cannot be stopped, so each piece is fetched in its own process (killed past
PIECE_TIMEOUT_S), and the whole fetch stops at budget_s.
"""
import concurrent.futures as cf
import os
import subprocess
import sys
import time
from pathlib import Path

from answer_set.k2.util import CACHE as K2_CACHE

URL = ("https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/"
       "wgs_pcr_free/50x/{p}.novaseq.pcr-free.50x.dedup.grch38.bam")
CACHE = K2_CACHE / "giab_origin"
INDEX_TIMEOUT_S = 60            # socket timeout of the index download (~9 MB)
PIECE_TIMEOUT_S = 300           # one piece (1/8 of a round region, ~35 MB) takes ~6 s
FETCH_BUDGET_S = 600            # all tries of one fetch together


def path(person, chrom, a, b):
    return CACHE / f"{person}_{chrom}_{a}_{b}.bam"


def remote_index(person, attempts=4, deadline=None):
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
            with urllib.request.urlopen(URL.format(p=person) + ".bai", timeout=INDEX_TIMEOUT_S) as src, \
                    open(part, "wb") as o:
                while True:
                    block = src.read(1 << 20)
                    if not block:
                        break
                    o.write(block)
                    if deadline is not None and time.monotonic() > deadline:
                        raise TimeoutError("fetch budget spent")
            part.rename(local)
            return local
        except OSError:                                   # includes socket and budget timeouts
            part.unlink(missing_ok=True)
            if t + 1 == attempts or (deadline is not None and time.monotonic() + 5 * (t + 1) >= deadline):
                break
            time.sleep(5 * (t + 1))
    raise RuntimeError(f"index of the original reads for {person} failed")


def _piece(person, chrom, x, y, dest, bai):
    """Reads starting in [x, y] into dest (run in its own process by fetch)."""
    import pysam
    src = pysam.AlignmentFile(URL.format(p=person), index_filename=str(bai))
    n = 0
    with pysam.AlignmentFile(str(dest), "wb", template=src) as o:
        for r in src.fetch(chrom, x - 1, y):
            if r.reference_start >= x - 1:                # each read once (by start)
                o.write(r)
                n += 1
    return n


def fetch(person, chrom, a, b, pad=500, chunks=8, attempts=4, budget_s=None):
    """Local BAM (indexed) with the person's original reads over [a-pad, b+pad]. RuntimeError if it cannot be had
    within budget_s (default FETCH_BUDGET_S)."""
    import pysam
    out = path(person, chrom, a, b)
    if out.exists() and Path(str(out) + ".bai").exists():
        return out
    deadline = time.monotonic() + (FETCH_BUDGET_S if budget_s is None else budget_s)
    bai = remote_index(person, deadline=deadline)
    lo, hi = max(1, a - pad), b + pad
    step = (hi - lo + chunks) // chunks
    pieces = [(lo + i * step, min(hi, lo + (i + 1) * step - 1)) for i in range(chunks) if lo + i * step <= hi]
    tmp = CACHE / f".tmp_{person}_{chrom}_{a}_{b}"
    tmp.mkdir(exist_ok=True)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [str(Path(__file__).resolve().parents[2]),
                                                                     os.environ.get("PYTHONPATH")])))

    def one(i):
        x, y = pieces[i]
        dest, last = tmp / f"p{i}.bam", "no try started"
        for t in range(attempts):
            left = deadline - time.monotonic()
            if left <= 0:
                break
            try:
                p = subprocess.run([sys.executable, "-m", "answer_set.k2.giab_origin", person, chrom, str(x), str(y),
                                    str(dest), str(bai)], env=env, capture_output=True, text=True,
                                   timeout=min(PIECE_TIMEOUT_S, left))
                if p.returncode == 0:
                    return int(p.stdout.split()[-1])
                last = (p.stderr.strip().splitlines() or [f"exit {p.returncode}"])[-1][:300]
            except subprocess.TimeoutExpired:                 # subprocess.run kills the child
                last = f"no answer in {min(PIECE_TIMEOUT_S, left):.0f} s"
            dest.unlink(missing_ok=True)
            if t + 1 < attempts and time.monotonic() + 5 * (t + 1) < deadline:
                time.sleep(5 * (t + 1))
        raise RuntimeError(f"original reads {person} {chrom}:{x}-{y} failed within the fetch budget ({last})")

    try:
        with cf.ThreadPoolExecutor(len(pieces)) as ex:
            counts = list(ex.map(one, range(len(pieces))))
        parts = [str(tmp / f"p{i}.bam") for i in range(len(pieces))]
        pysam.cat("-o", str(out) + ".part", *parts)
        Path(str(out) + ".part").rename(out)
        pysam.index(str(out))
    finally:
        for f in tmp.iterdir():
            f.unlink()
        tmp.rmdir()
    if sum(counts) == 0:
        out.unlink()
        raise RuntimeError(f"no original reads for {person} {chrom}:{a}-{b}")
    return out


if __name__ == "__main__":                                # one piece, for fetch: person chrom x y dest bai
    _p, _c, _x, _y, _d, _b = sys.argv[1:7]
    print(_piece(_p, _c, int(_x), int(_y), _d, _b))
