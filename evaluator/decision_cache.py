"""Exact hap.py scoring with a per-cluster decision cache. Owner: evaluator.

The region's truth and query variant spans are split into clusters: maximal runs of spans separated by gaps
< `gap` bp. A cluster's hap.py decisions depend only on the truth and query records inside it, so the annotated
hap.py lines around a cluster are cached under (span, query records in it). A callset is scored by assembling the
cached lines of its clusters into an annotated VCF and running the UNCHANGED upstream parsers on it
(parse_happy_vcf_assessed_metrics, compute_synthetic_only_metrics, parse_region_overcall_metrics, parse_happy_vcf,
difficulty_class_counts), in the same order and with the same fallbacks as upstream HappyScorer.score_vcf.
hap.py runs only for clusters never seen before, on a query subset holding just those clusters' records.

Verified (test/decision_cache.py prototype, test/dcache_validate.py): ScoreResults identical (Python ==) to
evaluator.gatk.HappyScorer on 174 callsets x 3 rounds at gap 1000 and 300.

Limitations:
  - Only FULL fidelity.
  - A callset with any new cluster pays one hap.py run (~8 s, the fixed cost dominates), however few clusters are
    new; savings come from callsets whose clusters are all cached.
  - Ratios hap.py takes from summary.csv are not in the annotated VCF. When the assessed-only parser cannot supply
    one (no assessed SNP tv / hom calls, practically never), the callset is scored by the full hap.py path instead.
  - A callset hap.py rejects is only detected when one of its clusters is new (then RuntimeError, as HappyScorer).
"""
from __future__ import annotations

import bisect
import gzip
import hashlib
import os
import shutil
import sqlite3
import sys
import threading
import uuid
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from common.scoring import finalize
from common.settings import SETTINGS
from common.types import CLASSES, FULL, AnswerSet, CallSet, ClassCounts, Round, ScoreResult

if SETTINGS.minos_subnet_dir not in sys.path:
    sys.path.append(SETTINGS.minos_subnet_dir)

GAP = 1000
RATIO_KEYS = ("titv_truth_snp", "titv_query_snp", "hethom_truth_snp", "hethom_query_snp")


class Record(NamedTuple):
    pos: int
    end: int        # pos + len(REF) - 1
    ref: str
    alt: str
    gt: str         # sample GT with phasing dropped ("0|1" -> "0/1"); "" without a sample column
    filter: str


class Cluster(NamedTuple):
    start: int
    end: int
    query: Tuple[Tuple[int, str, str, str, str], ...]   # (pos, ref, alt, gt, filter) of the query records inside

    @property
    def key(self) -> tuple:
        return (self.start, self.end, self.query)


# ---- pure logic ----

def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_records(vcf: str, chrom: str, start: int, end: int) -> List[Record]:
    """Records of `vcf` on `chrom` with start <= POS <= end, sorted by POS (stable)."""
    out = []
    with _open(vcf) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            pos = int(p[1])
            if p[0] == chrom and start <= pos <= end:
                gt = p[9].split(":")[0].replace("|", "/") if len(p) > 9 else ""
                out.append(Record(pos, pos + len(p[3]) - 1, p[3], p[4], gt, p[6]))
    out.sort(key=lambda r: r.pos)
    return out


def make_clusters(truth_spans: Sequence[Tuple[int, int]], query: Sequence[Record], gap: int = GAP) -> List[Cluster]:
    """Maximal runs of truth+query spans whose gaps are < gap, each with the query records starting inside it.

    `query` must be sorted by pos. Clusters come back sorted and disjoint (consecutive ones >= gap apart).
    """
    runs: List[List[int]] = []
    for s, e in sorted(list(truth_spans) + [(q.pos, q.end) for q in query]):
        if runs and s - runs[-1][1] < gap:
            runs[-1][1] = max(runs[-1][1], e)
        else:
            runs.append([s, e])
    qpos = [q.pos for q in query]
    out = []
    for s, e in runs:
        i, j = bisect.bisect_left(qpos, s), bisect.bisect_right(qpos, e)
        out.append(Cluster(s, e, tuple((q.pos, q.ref, q.alt, q.gt, q.filter) for q in query[i:j])))
    return out


def pack_clusters(clusters: Sequence[Cluster], gap: int = GAP) -> List[List[Cluster]]:
    """Bins of clusters that can share one hap.py run: inside a bin, spans are disjoint and >= gap apart, the same
    independence the cache relies on. First-fit over clusters sorted by start; order of the result is stable."""
    bins: List[List[Cluster]] = []
    for c in sorted(clusters, key=lambda c: (c.start, c.end)):
        for b in bins:
            if c.start - b[-1].end >= gap:
                b.append(c)
                break
        else:
            bins.append([c])
    return bins


def window_pad(gap: int = GAP) -> int:
    """Padding around a cluster's span when taking its hap.py lines: windows of distinct clusters never overlap."""
    return gap // 2 - 1


def lines_in_window(lines: Sequence[Tuple[int, str]], start: int, end: int) -> Tuple[str, ...]:
    """The lines with start <= pos <= end from (pos, line) pairs sorted by pos."""
    pos = [p for p, _ in lines]
    return tuple(l for _, l in lines[bisect.bisect_left(pos, start):bisect.bisect_right(pos, end)])


def in_spans(pos: int, spans: Sequence[Tuple[int, int]]) -> bool:
    """Whether pos lies in one of the sorted, disjoint (start, end) spans."""
    i = bisect.bisect_right(spans, (pos, float("inf"))) - 1
    return i >= 0 and spans[i][0] <= pos <= spans[i][1]


def _line_pos(line: str) -> int:
    return int(line.split("\t", 2)[1])


def _bgzip_index(plain: Path) -> Path:
    import pysam
    return Path(pysam.tabix_index(str(plain), preset="vcf", force=True))


# ---- the cache ----

class _Store:
    """One round's cached cluster lines on disk, shared by every process on the machine (SQLite, WAL). A cluster's
    lines depend only on the records inside it, so concurrent writers of one cluster write the same value."""

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self._local = threading.local()

    def _db(self) -> sqlite3.Connection:
        db = getattr(self._local, "db", None)
        if db is None:
            db = sqlite3.connect(str(self.path), timeout=120, isolation_level=None, check_same_thread=False)
            db.execute("PRAGMA journal_mode=WAL")
            db.execute("PRAGMA synchronous=NORMAL")
            db.execute("CREATE TABLE IF NOT EXISTS clusters (k TEXT PRIMARY KEY, v TEXT NOT NULL)")
            db.execute("CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT NOT NULL)")
            self._local.db = db
        return db

    def get_many(self, skeys: Sequence[str]) -> Dict[str, str]:
        out: Dict[str, str] = {}
        db = self._db()
        for i in range(0, len(skeys), 500):
            chunk = list(skeys[i:i + 500])
            q = f"SELECT k, v FROM clusters WHERE k IN ({','.join('?' * len(chunk))})"
            out.update(db.execute(q, chunk).fetchall())
        return out

    def put_many(self, items: Dict[str, str]) -> None:
        if items:
            self._db().executemany("INSERT OR IGNORE INTO clusters (k, v) VALUES (?, ?)", list(items.items()))

    def get_meta(self, k: str) -> Optional[str]:
        row = self._db().execute("SELECT v FROM meta WHERE k = ?", (k,)).fetchone()
        return row[0] if row else None

    def put_meta(self, k: str, v: str) -> None:
        self._db().execute("INSERT OR IGNORE INTO meta (k, v) VALUES (?, ?)", (k, v))


_CONTENT_IDS: Dict[tuple, str] = {}


def _content_id(path) -> str:
    """sha1 of a (small) file, cached by path, size and mtime."""
    try:
        st = os.stat(path)
    except OSError:
        return "missing"
    k = (str(path), st.st_size, st.st_mtime_ns)
    if k not in _CONTENT_IDS:
        h = hashlib.sha1()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        _CONTENT_IDS[k] = h.hexdigest()[:16]
    return _CONTENT_IDS[k]


def _skey(key: tuple) -> str:
    return hashlib.sha1(repr(key).encode()).hexdigest()


class _RoundCache:
    """Cached hap.py lines for one (round, answer-set version, truth): in memory, backed by the machine-wide store."""

    def __init__(self, truth: str, mutations: str, region: str, reference: str, work_dir: Path, gap: int):
        self.truth, self.mutations, self.region, self.reference = truth, mutations, region, reference
        self.chrom, span = region.split(":")
        self.start, self.end = (int(x) for x in span.replace(",", "").split("-"))
        self.sdf = str(Path(reference).parent / f"{self.chrom}.sdf")
        self.work_dir, self.gap = work_dir, gap
        self.truth_spans = [(r.pos, r.end) for r in read_records(truth, self.chrom, self.start, self.end)]
        self.lines: Dict[tuple, Tuple[str, ...]] = {}
        self.header: Optional[str] = None
        self.store = _Store(work_dir / "clusters.sqlite")

    def preload(self, keys: Sequence[tuple]) -> None:
        """Bring clusters other processes on this machine already scored into memory (one query)."""
        want = {_skey(k): k for k in keys if k not in self.lines}
        if want:
            for sk, v in self.store.get_many(list(want)).items():
                self.lines[want[sk]] = tuple(v.splitlines(keepends=True))
        if self.header is None:
            self.header = self.store.get_meta("header")

    def _persist(self, keys: Sequence[tuple]) -> None:
        self.store.put_many({_skey(k): "".join(self.lines[k]) for k in keys if k in self.lines})
        if self.header is not None:
            self.store.put_meta("header", self.header)

    def clusters(self, query: Sequence[Record]) -> List[Cluster]:
        return make_clusters(self.truth_spans, query, self.gap)

    def run_happy(self, query_vcf: Path) -> Optional[str]:
        """Upstream hap.py path on one query VCF; the annotated VCF's path or None if hap.py produced nothing."""
        from evaluator import truth_slices
        from evaluator.sandbox import install_scoring_shim
        from utils.scoring import HappyScorer as UpstreamHappy
        install_scoring_shim()                  # hap.py and bcftools from the image files when the sandbox is on
        truth_slices.install()                 # the truth is sliced once per round
        m = UpstreamHappy().score_vcf(truth_vcf=self.truth, query_vcf=str(query_vcf), reference_fasta=self.reference,
                                      confident_bed=None, region=self.region, reference_sdf=self.sdf,
                                      mutations_vcf=self.mutations)
        path = (m or {}).get("happy_vcf_path")
        return path if path and os.path.exists(path) else None

    def fill(self, callset_vcf: str, missing: List[Cluster]) -> int:
        """Run hap.py on the callset restricted to the missing clusters and cache their lines. Returns hap.py runs.

        If hap.py fails on the subset (e.g. it holds no records) it is retried once on the whole callset; if that
        fails too the callset is unscorable, as in HappyScorer.
        """
        spans = sorted((c.start, c.end) for c in missing)
        tmp = self.work_dir / f"fill_{uuid.uuid4().hex[:12]}"
        tmp.mkdir(parents=True)
        try:
            header, kept, total = [], [], 0
            with _open(callset_vcf) as f:
                for line in f:
                    if line[0] == "#":
                        header.append(line)
                        continue
                    total += 1
                    p = line.split("\t", 2)
                    if p[0] == self.chrom and in_spans(int(p[1]), spans):
                        kept.append(line)
            runs = 0
            happy_vcf = None
            for body in ([kept, None] if len(kept) < total else [None]):
                q = tmp / f"q{runs}.vcf"
                with open(q, "w") as g:
                    g.writelines(header)
                    if body is None:
                        with _open(callset_vcf) as f:
                            g.writelines(l for l in f if l[0] != "#")
                    else:
                        g.writelines(body)
                runs += 1
                happy_vcf = self.run_happy(_bgzip_index(q))
                if happy_vcf:
                    break
            if not happy_vcf:
                raise RuntimeError("hap.py produced no result")
            hdr, lines = [], []
            with gzip.open(happy_vcf, "rt") as f:
                for line in f:
                    if line[0] == "#":
                        hdr.append(line)
                    else:
                        lines.append((_line_pos(line), line))
            lines.sort(key=lambda x: x[0])
            if self.header is None:
                self.header = "".join(hdr)
            pad = window_pad(self.gap)
            for c in missing:
                self.lines[c.key] = lines_in_window(lines, c.start - pad, c.end + pad)
            self._persist([c.key for c in missing])
            return runs
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def fill_packed(self, items: Dict[tuple, Tuple[Cluster, List[str]]], header: str, workers: int = 4) -> int:
        """Cache many clusters (from any callsets) with as few hap.py runs as possible. `items` maps cluster key ->
        (cluster, its query lines). Clusters >= gap apart share a run. Returns hap.py runs; raises if any run fails."""
        bins = pack_clusters([c for c, _ in items.values()], self.gap)
        pad = window_pad(self.gap)

        def run_bin(i_bin):
            i, b = i_bin
            tmp = self.work_dir / f"pack_{uuid.uuid4().hex[:12]}"
            tmp.mkdir(parents=True)
            try:
                q = tmp / "q.vcf"
                body = sorted((l for c in b for l in items[c.key][1]), key=_line_pos)
                q.write_text(header + "".join(body))
                happy_vcf = self.run_happy(_bgzip_index(q))
                if not happy_vcf:
                    raise RuntimeError("hap.py produced no result for a packed run")
                hdr, lines = [], []
                with gzip.open(happy_vcf, "rt") as f:
                    for line in f:
                        (hdr if line[0] == "#" else lines).append(line if line[0] == "#" else (_line_pos(line), line))
                lines.sort(key=lambda x: x[0])
                return "".join(hdr), {c.key: lines_in_window(lines, c.start - pad, c.end + pad) for c in b}
            finally:
                shutil.rmtree(tmp, ignore_errors=True)

        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max(1, min(workers, len(bins)))) as pool:
            results = list(pool.map(run_bin, enumerate(bins)))
        for hdr, got in results:           # only after every run succeeded: all or nothing
            if self.header is None:
                self.header = hdr
            self.lines.update(got)
            self._persist(list(got))
        return len(bins)

    def assemble(self, clusters: Sequence[Cluster], out: Path) -> Path:
        """The annotated VCF for these (all cached) clusters, bgzipped."""
        body = sorted((l for c in clusters for l in self.lines[c.key]), key=_line_pos)
        plain = out.with_suffix("")     # x.vcf.gz -> x.vcf
        with open(plain, "w") as g:
            g.write(self.header or "")
            g.writelines(body)
        return _bgzip_index(plain)


class DecisionCacheScorer:
    """Scorer returning exactly evaluator.gatk.HappyScorer's ScoreResult, running hap.py only for new clusters.

    One cache per (round_id, answer_set.version, truth file): a new answer-set version starts a new cache. It lives
    on disk (work_dir/<round tag>/<cache>/clusters.sqlite) and is shared by every scorer process on the machine that
    uses the same work_dir - a cluster scored by one is never sent to hap.py by another - with each process keeping
    an in-memory copy. Safe from several threads and processes; hap.py fills run concurrently (no lock held while
    hap.py runs): a cluster's lines depend only on the records inside it, so two fills of one cluster are identical.
    """

    def __init__(self, work_dir: Optional[str] = None, gap: int = GAP):
        self.work_dir = Path(work_dir or os.path.join(SETTINGS.cache_dir, "dcache"))
        self.gap = gap
        self._caches: Dict[tuple, _RoundCache] = {}
        self._lock = threading.Lock()
        self._stats = {"callsets": 0, "cached_callsets": 0, "hits": 0, "misses": 0, "happy_runs": 0,
                       "full_fallbacks": 0, "packed_runs": 0, "packed_failures": 0}

    # ---- public ----
    @property
    def stats(self) -> Dict[str, int]:
        """callsets scored, cached_callsets (no hap.py needed), cluster hits/misses, happy_runs, full_fallbacks."""
        with self._lock:
            return dict(self._stats)

    def drop(self, round_id: str) -> None:
        """Forget every cache of a round (e.g. once the round is over)."""
        with self._lock:
            for k in [k for k in self._caches if k[0] == round_id]:
                del self._caches[k]

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        if callset.fidelity != FULL:
            raise NotImplementedError(f"DecisionCacheScorer supports only FULL fidelity, got {callset.fidelity!r}")
        cache = self._cache_for(answer_set)
        query = read_records(callset.vcf, cache.chrom, cache.start, cache.end)
        clusters = cache.clusters(query)
        cache.preload([c.key for c in clusters])
        n_missing, runs = self._ensure_cached(cache, callset.vcf, clusters)
        self._count(callsets=1, cached_callsets=int(n_missing == 0), hits=len(clusters) - n_missing,
                    misses=n_missing, happy_runs=runs)

        tmp = self.work_dir / f"asm_{uuid.uuid4().hex[:12]}"
        tmp.mkdir(parents=True)
        try:
            happy_vcf = str(cache.assemble(clusters, tmp / "happy.vcf.gz"))
            metrics = self._metrics(happy_vcf, cache.mutations)
            if metrics is None:      # a value only hap.py's summary.csv has: take the full path
                self._count(full_fallbacks=1)
                from evaluator.gatk import HappyScorer
                return HappyScorer().score(callset, answer_set)
            from utils.scoring import difficulty_class_counts, parse_happy_vcf
            cc = difficulty_class_counts(parse_happy_vcf(happy_vcf, truth_vcf_path=cache.truth))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        ratios = {k: metrics[k] for k in RATIO_KEYS if k in metrics}
        result = ScoreResult(callset.round_id, callset.config_key, answer_set.version,
                             {k: ClassCounts(**cc[k]) for k in CLASSES}, int(metrics.get("region_fp_total", 0)),
                             answer_set.n_targets or 1, ratios, fidelity=callset.fidelity, job_id=callset.job_id)
        return finalize(result)

    def score_many(self, callsets: Sequence[CallSet], answer_set: AnswerSet, workers: int = 4) -> List[ScoreResult]:
        """Score a batch: new clusters of all callsets are deduplicated and packed into as few hap.py runs as
        possible (run in parallel), then each callset is scored from the cache. Results equal score() per callset.
        If a packed run fails, falls back to per-callset fills so one bad callset cannot fail the others."""
        cache = self._cache_for(answer_set)
        items: Dict[tuple, Tuple[Cluster, List[str]]] = {}
        header = None
        for cs in callsets:
            if cs.fidelity != FULL or not cs.vcf:
                continue
            query = read_records(cs.vcf, cache.chrom, cache.start, cache.end)
            clusters = cache.clusters(query)
            cache.preload([c.key for c in clusters])
            missing = [c for c in clusters if c.key not in cache.lines and c.key not in items]
            if not missing:
                continue
            spans = sorted((c.start, c.end) for c in missing)
            by_span: Dict[Tuple[int, int], List[str]] = {sp: [] for sp in spans}
            with _open(cs.vcf) as f:
                for line in f:
                    if line[0] == "#":
                        continue
                    p = line.split("\t", 2)
                    if p[0] != cache.chrom:
                        continue
                    pos = int(p[1])
                    i = bisect.bisect_right(spans, (pos, float("inf"))) - 1
                    if i >= 0 and spans[i][0] <= pos <= spans[i][1]:
                        by_span[spans[i]].append(line)
            if header is None:
                with _open(cs.vcf) as f:
                    header = "".join(l for l in f if l[0] == "#")
            for c in missing:
                items[c.key] = (c, by_span[(c.start, c.end)])
        items = {k: v for k, v in items.items() if k not in cache.lines}   # another thread may have filled some
        if items:
            try:
                runs = cache.fill_packed(items, header, workers)
                self._count(happy_runs=runs, packed_runs=runs)
            except RuntimeError:
                self._count(packed_failures=1)             # per-callset fills below isolate the bad callset
        out = []
        for cs in callsets:
            out.append(self.score(cs, answer_set))
        return out

    # ---- internals ----
    def _count(self, **inc: int) -> None:
        with self._lock:
            for k, v in inc.items():
                self._stats[k] += v

    def _cache_for(self, answer_set: AnswerSet) -> _RoundCache:
        """The cache of this answer set; keyed by the truth files' CONTENT as well as their paths, so two answer sets
        written to the same path (a rebuilt truth, another builder on the same round) never share hap.py lines."""
        n = answer_set.notes
        truth, mutations, region, reference = (n.get(k) for k in ("truth_vcf", "mutations_vcf", "region", "reference"))
        if not (truth and mutations and region and reference):
            raise NotImplementedError("DecisionCacheScorer needs truth/mutations VCFs, region and reference in "
                                      "answer_set.notes (writing them from reconstructed records is not implemented)")
        key = (answer_set.round_id, answer_set.version, str(truth), _content_id(truth), str(mutations),
               _content_id(mutations), region, str(reference))
        with self._lock:
            cache = self._caches.get(key)
            if cache is None:
                name = hashlib.sha1(repr(key).encode()).hexdigest()[:12]
                d = self.work_dir / Round.tag_of(answer_set.round_id) / name   # machine-wide: every scorer shares it
                cache = _RoundCache(str(truth), str(mutations), region, str(reference), d, self.gap)
                self._caches[key] = cache
            return cache

    @staticmethod
    def _ensure_cached(cache: _RoundCache, vcf: str, clusters: List[Cluster]) -> Tuple[int, int]:
        """(clusters that were missing, hap.py runs) after making sure every cluster is cached."""
        missing = [c for c in clusters if c.key not in cache.lines]
        if not missing:
            return 0, 0
        return len(missing), cache.fill(vcf, missing)     # no lock: concurrent fills of a cluster are identical

    @staticmethod
    def _metrics(happy_vcf: str, mutations: str) -> Optional[dict]:
        """The metrics upstream score_vcf derives from its annotated VCF, same calls and fallbacks, same order.

        None when the result would still hold a ratio taken from hap.py's summary.csv (not reproducible here).
        """
        from utils.scoring import (compute_synthetic_only_metrics, parse_happy_vcf_assessed_metrics,
                                   parse_region_overcall_metrics)
        metrics = {}
        assessed = parse_happy_vcf_assessed_metrics(happy_vcf)
        if not assessed or any(k not in assessed for k in RATIO_KEYS):
            return None
        metrics.update(assessed)
        syn = compute_synthetic_only_metrics(happy_vcf, mutations)
        if syn is None:
            raise RuntimeError("hap.py produced no result")   # upstream returns no scores here
        metrics.update(syn)
        target_total = syn.get("target_total_snp", 0.0) + syn.get("target_total_indel", 0.0)
        target_snp = syn.get("target_total_snp", 0.0)
        if target_total <= 0:
            target_total = syn.get("truth_total_snp", 0.0) + syn.get("truth_total_indel", 0.0)
            target_snp = syn.get("truth_total_snp", 0.0)
        over = parse_region_overcall_metrics(
            happy_vcf, target_total, target_snp,
            v1_truth_total=syn.get("truth_total_snp", 0.0) + syn.get("truth_total_indel", 0.0),
            v1_snp_truth_total=syn.get("truth_total_snp", 0.0))
        if over:
            metrics.update(over)
        return metrics
