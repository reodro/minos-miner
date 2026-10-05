"""The ORIGINAL reads a 1000 Genomes round was built from, fetched for the round's region only.

A 1000G round BAM consists of the person's public 1000 Genomes high-coverage reads (2x150, ~34x, NYGC CRAMs at EBI):
~82% of read pairs kept, each read cut to its first 128 bases in sequencing order, names hashed, some reads edited
(the inserted mutations, a little noise), realigned (measured 2026-10-04 on HG01519 chr22: 96% of reads with
distinctive qualities have an original with the same qualities and sequence). So readdiff.py works unchanged against
these originals, as for GIAB rounds (giab_origin.py).

The person's CRAM is found in EBI's sequence indexes of the 2,504 + 698 high-coverage sets (cached). A round region
(~10 Mb, ~2.5M reads) takes ~10-20 s in 8 parallel pieces. Same rules as giab_origin: each piece in its own process
(killed past PIECE_TIMEOUT_S; htslib's remote reads have no timeout), the whole fetch stops at budget_s.
Owner: answer_set.
"""
import concurrent.futures as cf
import os
import subprocess
import sys
import time
from pathlib import Path

from answer_set.k2.util import CACHE as K2_CACHE, ref_path

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
INDEXES = ("1000G_2504_high_coverage.sequence.index", "1000G_698_related_high_coverage.sequence.index")
CACHE = K2_CACHE / "kg_origin"
INDEX_TIMEOUT_S = 60
PIECE_TIMEOUT_S = 300
FETCH_BUDGET_S = 600


def path(person, chrom, a, b):
    return CACHE / f"{person}_{chrom}_{a}_{b}.bam"


def _download(url, dest, deadline):
    import urllib.request
    left = deadline - time.monotonic()
    if left <= 0:
        raise RuntimeError(f"no time left to download {url}")
    tmp = Path(str(dest) + ".part")
    with urllib.request.urlopen(url, timeout=min(INDEX_TIMEOUT_S, left)) as r, open(tmp, "wb") as o:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            o.write(chunk)
            if time.monotonic() > deadline:
                raise RuntimeError(f"download of {url} ran past the fetch budget")
    tmp.rename(dest)


def cram_url(person, deadline=None, attempts=4):
    """https URL of the person's CRAM, from EBI's sequence indexes (downloaded once)."""
    deadline = deadline or time.monotonic() + FETCH_BUDGET_S
    CACHE.mkdir(parents=True, exist_ok=True)
    for name in INDEXES:
        local = CACHE / name
        for t in range(attempts):
            if local.exists():
                break
            try:
                _download(BASE + name, local, deadline)
            except Exception as e:                        # network trouble: retry within the budget
                if t + 1 == attempts or time.monotonic() + 5 * (t + 1) > deadline:
                    raise RuntimeError(f"1000G sequence index {name}: {e}") from e
                time.sleep(5 * (t + 1))
        for line in local.read_text().splitlines():
            if line.startswith("#"):
                continue
            u = line.split("\t", 1)[0]
            if u.endswith(f"/{person}.final.cram"):
                return u.replace("ftp://", "https://")
    raise RuntimeError(f"no 1000G high-coverage CRAM listed for {person}")


def remote_index(person, url, deadline, attempts=4):
    local = CACHE / f"{person}.final.cram.crai"
    for t in range(attempts):
        if local.exists():
            return local
        try:
            _download(url + ".crai", local, deadline)
        except Exception as e:
            if t + 1 == attempts or time.monotonic() + 5 * (t + 1) > deadline:
                raise RuntimeError(f"CRAM index of {person}: {e}") from e
            time.sleep(5 * (t + 1))
    return local


def _piece(url, chrom, x, y, dest, crai):
    """Reads starting in [x, y] into dest (run in its own process by fetch)."""
    import pysam
    src = pysam.AlignmentFile(url, "rc", index_filename=str(crai), reference_filename=str(ref_path(chrom)))
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
    url = cram_url(person, deadline)
    crai = remote_index(person, url, deadline)
    lo, hi = max(1, a - pad), b + pad
    step = (hi - lo + chunks) // chunks
    pieces = [(lo + i * step, min(hi, lo + (i + 1) * step - 1)) for i in range(chunks) if lo + i * step <= hi]
    tmp = CACHE / f".tmp_{person}_{chrom}_{a}_{b}"
    tmp.mkdir(parents=True, exist_ok=True)
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
                p = subprocess.run([sys.executable, "-m", "answer_set.k2.kg_origin", url, chrom, str(x), str(y),
                                    str(dest), str(crai)], env=env, cwd=str(tmp), capture_output=True, text=True,
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
        # a remote read can end early without an error and leave a piece nearly empty (seen once in 216 downloads:
        # 31 reads instead of ~195,000, so a twelfth of a round's reads had no original). Such a piece is fetched
        # again within the budget; if it stays nearly empty, it is reported (the build goes on with the rest).
        for _ in range(attempts):
            med = sorted(counts)[len(counts) // 2]
            weak = [i for i, c in enumerate(counts) if med > 1000 and c < 0.2 * med]
            if not weak or time.monotonic() >= deadline:
                break
            for i in weak:
                counts[i] = one(i)
        med = sorted(counts)[len(counts) // 2]
        weak = [i for i, c in enumerate(counts) if med > 1000 and c < 0.2 * med]
        if weak:
            print(f"WARNING original reads {person} {chrom}:{a}-{b}: piece(s) {weak} still nearly empty "
                  f"({[counts[i] for i in weak]} reads vs median {med})", file=sys.stderr)
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


if __name__ == "__main__":                                # one piece, for fetch: url chrom x y dest crai
    _u, _c, _x, _y, _d, _i = sys.argv[1:7]
    print(_piece(_u, _c, int(_x), int(_y), _d, _i))
