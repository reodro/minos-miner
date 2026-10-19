"""N for a round the platform gives no N for, from the round's own pool of inserted mutations and past rounds.

N (the mutations.vcf entries, the v2 score's false-call divisor) is a random draw, and the scored mutations are a
sample of the ones Minos inserts (the pool), never more than 1/8 of them. 488 rounds (2026-09-23..10-19) show how:
  - a region starts with a base pool of ~600-615 inserted per Mb;
  - when 8 x N does not fit in it, the pool is raised until it does: above 660 per Mb, N = q x pool / 8 with
    q ~ 0.97 (the median over past raised rounds);
  - otherwise N is a draw that fitted the base pool: the median of past draws of the same round type (binned
    quality <= 9 distinct base qualities, or full) that fit this pool.
What looked like batches of the task window with their own N range (binned 26-28 Sep median 578 against 432 and
460 later) were raised pools. Live, the pool is the inserted mutations our answer set flags SYNTHETIC (median 0.2%
under the truth's, at most 5% on healthy builds).
Walk-forward over 465 rounds (each sees only rounds revealed before it started), error in N_true / N_est:
team median of the last 10 rounds 17.5%, type median capped at pool/8 15.4%, this 12.6% (test/answer_set/n_study).
"""
from __future__ import annotations

import gzip
import json
import os
import statistics
from typing import Optional, Tuple

import pysam

from answer_set.k2.util import CACHE

ENABLED = True            # False: the caller's N (main.n_estimate) as before
MIN_ROUNDS = 10           # fewer revealed rounds of this type: no estimate (the caller keeps its N)
MIN_RAISED = 5            # fewer revealed raised-pool rounds: no pool rule
CAP_SHARE = 0.125         # N <= pool / 8
RAISED_PER_MB = 660       # pool density above which the pool was raised to fit N
POOL_SCALE = 0.9976       # our pool count / the truth's (median, healthy builds)
POOL_UNDER = 0.95         # ... and at most 5% under
POOL_MIN_PER_MB = 540     # every round had >= 598 inserted per Mb; fewer flagged means the build lost some
ROUND_CACHE = CACHE / "n_rounds.json"


def bam_type(bam: str, reads: int = 5000) -> Optional[str]:
    """'binned' (<= 9 distinct base qualities in the first reads) or 'full'; None when unreadable."""
    try:
        quals, n = set(), 0
        with pysam.AlignmentFile(str(bam)) as f:
            for r in f.fetch(until_eof=True):
                if r.query_qualities is not None:
                    quals.update(r.query_qualities[:50])
                n += 1
                if n >= reads:
                    break
        return ("binned" if len(quals) <= 9 else "full") if quals else None
    except (OSError, ValueError):
        return None


def _span(region: str) -> Optional[int]:
    try:
        a, b = map(int, region.split(":")[1].split("-"))
        return b - a + 1
    except (AttributeError, IndexError, ValueError):
        return None


def _round_facts(d: str) -> Optional[list]:
    """[type, N, pool (truth's SYNTHETIC records), span] of a revealed round directory; None when incomplete."""
    try:
        with open(os.path.join(d, "round.json")) as f:
            region = json.load(f)["region"]
        with gzip.open(os.path.join(d, "mutations.vcf.gz"), "rt") as f:
            n = sum(1 for line in f if not line.startswith("#"))
        with gzip.open(os.path.join(d, "truth.vcf.gz"), "rt") as f:
            pool = sum(1 for line in f if not line.startswith("#") and "SYNTHETIC" in line.split("\t", 8)[7])
    except (OSError, EOFError, ValueError, KeyError, IndexError):
        return None
    ty, span = bam_type(os.path.join(d, "input.bam")), _span(region)
    return [ty, n, pool, span] if ty and n and pool and span else None


def past_rounds(tag: str, rounds_dir: Optional[str] = None) -> list:
    """Facts of every round whose files were out when round `tag` started (cached: a round's files never change)."""
    from main.n_estimate import REVEAL_S, _start
    from common.settings import SETTINGS
    start, root = _start(tag), rounds_dir or SETTINGS.rounds_dir
    if start is None:
        return []
    try:
        tags = [t for t in os.listdir(root) if _start(t) is not None and (start - _start(t)).total_seconds() >= REVEAL_S]
    except OSError:
        return []
    try:
        cache = json.loads(ROUND_CACHE.read_text())
    except (OSError, ValueError):
        cache = {}
    new = False
    for t in tags:
        d = os.path.join(root, t)
        if t not in cache and all(os.path.exists(os.path.join(d, f)) for f in
                                  ("input.bam", "mutations.vcf.gz", "truth.vcf.gz", "round.json")):
            facts = _round_facts(d)
            if facts:
                cache[t] = facts
                new = True
    if new:
        try:
            ROUND_CACHE.parent.mkdir(parents=True, exist_ok=True)
            tmp = ROUND_CACHE.with_suffix(".part")
            tmp.write_text(json.dumps(cache))
            os.replace(tmp, ROUND_CACHE)
        except OSError:
            pass
    return [cache[t] for t in tags if t in cache]


def estimate(tag: str, bam: str, region: str, pool: int, rounds_dir: Optional[str] = None
             ) -> Tuple[Optional[int], dict]:
    """(N, how) - N None when there are too few past rounds of this type."""
    rtype, span = bam_type(bam), _span(region)
    past = past_rounds(tag, rounds_dir) if rtype else []
    ok = bool(span) and pool >= POOL_MIN_PER_MB * span / 1e6        # a short count (lost flags) says nothing
    pool_t = pool / POOL_SCALE if ok else None                       # as the truth would count it
    how = {"type": rtype, "pool": pool, "per_mb": round(pool_t / span * 1e6) if ok else None}
    if ok and pool_t / span * 1e6 >= RAISED_PER_MB:
        q = [n / (p / 8) for _, n, p, s in past if p / s * 1e6 >= RAISED_PER_MB]
        if len(q) >= MIN_RAISED:
            how.update(rule="raised pool", q=round(statistics.median(q), 4), rounds=len(q))
            return int(round(statistics.median(q) * CAP_SHARE * pool_t)), how
    same = [n for ty, n, _, _ in past if ty == rtype]
    cap = CAP_SHARE * pool_t / POOL_UNDER if ok else None
    fit = [n for n in same if cap is None or n <= cap]
    if len(fit) >= MIN_ROUNDS:
        how.update(rule="draws that fit the pool", cap=round(cap) if cap else None, rounds=len(fit))
        return int(round(statistics.median(fit))), how
    if len(same) >= MIN_ROUNDS:
        how.update(rule="type median, capped", cap=round(cap) if cap else None, rounds=len(same))
        return int(round(min(statistics.median(same), cap) if cap else statistics.median(same))), how
    return None, how
