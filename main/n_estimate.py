"""N for a round the platform gives no N for: the median target count of recent released rounds.

The v2 score divides false calls by N (the number of injected targets). The platform does not publish N live, and the
fixed answer_set.n (548) matched rounds before 2026-10-17 (median 567) but not the rounds since (10-18: median 403,
284-535). On 20261018T130600Z (true N 352) the fixed 548 made false calls look cheap and the pick took a call set with
~90 more false calls: 0.175 below the best reviewed config; with the true N the review's Spearman went from 0.56 to
0.996. Over 47 rounds (10-10..10-18) the median of the last 10 revealed rounds lost least: mean 0.0060 against 0.0089
for the fixed 548, worst round 0.077 against 0.113.

A round counts once its truth was out when the round being estimated started (REVEAL_S after its own start), so a
replay sees what a live round would have seen.
"""
from __future__ import annotations

import gzip
import os
import statistics
from datetime import datetime, timezone
from typing import Optional

from common.settings import SETTINGS
from common.tuning import tune

REVEAL_S = 3.5 * 3600          # a round's truth and mutations come out about this long after it starts
MIN_ROUNDS = 3                 # fewer usable rounds: no estimate (the fixed answer_set.n applies)


def _start(tag: str) -> Optional[datetime]:
    try:
        return datetime.strptime(tag, "%Y%m%dT%H%M%SZ").replace(tzinfo=timezone.utc)
    except ValueError:
        return None


def _count(path: str) -> Optional[int]:
    try:
        with gzip.open(path, "rt") as f:
            return sum(1 for line in f if not line.startswith("#"))
    except (OSError, EOFError):
        return None


def recent_n(tag: str, rounds_dir: Optional[str] = None, k: Optional[int] = None) -> Optional[int]:
    """Median N of the last `k` (main.n_recent_rounds) rounds whose mutations were out when round `tag` started;
    None when the setting is 0 or fewer than MIN_ROUNDS such rounds are on disk."""
    k = int(tune("main.n_recent_rounds")) if k is None else int(k)
    start = _start(tag)
    if k <= 0 or start is None:
        return None
    root = rounds_dir or SETTINGS.rounds_dir
    try:
        tags = sorted(t for t in os.listdir(root) if _start(t) is not None)
    except OSError:
        return None
    counts = []
    for t in reversed(tags):
        s = _start(t)
        if (start - s).total_seconds() < REVEAL_S:
            continue                                     # not revealed yet when this round started
        n = _count(os.path.join(root, t, "mutations.vcf.gz"))
        if n:
            counts.append(n)
        if len(counts) >= k:
            break
    if len(counts) < MIN_ROUNDS:
        return None
    return int(round(statistics.median(counts)))
