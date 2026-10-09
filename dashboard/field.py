"""Where a replay's current best would sit on the published field.

The dashboard computes this. Main keeps mining; a slow or failed field read never reaches the search.
"""
from __future__ import annotations

import statistics
import threading
import time
from typing import Callable, List, Optional

from backtest.kpi import rank_in_field

# A finished round on this subnet is about 70 scored submissions. Well under that, the board is still filling.
FILLING_BELOW = 50


def published_scores(round_id: str) -> List[float]:
    """One 0-100 score per submission, from the public scores endpoint."""
    from common.platform import Platform
    payload = Platform().scores(round_id)
    by = {}
    for row in payload.get("scores") or []:
        sid, final = row.get("submission_id"), row.get("combined_final")
        if sid and final is not None:
            by[sid] = 100.0 * float(final)
    return list(by.values())


def place(score: float, field: List[float], *, replay: bool) -> dict:
    """Rank of `score` among the published scores. A replay was never submitted, so it joins the field."""
    scores = [float(s) for s in field]
    if not scores:
        return {"score": score, "rank": None, "field": 0, "field_best": None, "field_median": None,
                "gap_to_best": None, "filling": True}
    best, med = max(scores), statistics.median(scores)
    return {
        "score": score,
        "rank": rank_in_field(score, scores),
        "field": len(scores) + (1 if replay else 0),
        "field_best": round(best, 4),
        "field_median": round(med, 4),
        "gap_to_best": round(best - score, 4),
        "filling": len(scores) < FILLING_BELOW,
    }


class FieldCache:
    """Published scores for a round, reused for `ttl_s`. A failed refresh keeps the previous list."""

    def __init__(self, fetch: Callable[[str], List[float]], ttl_s: float = 60.0):
        self.fetch, self.ttl_s = fetch, ttl_s
        self._lock = threading.Lock()
        self._rows: dict = {}

    def scores(self, round_id: str) -> Optional[List[float]]:
        now = time.monotonic()
        with self._lock:
            hit = self._rows.get(round_id)
            if hit is not None and hit[1] is not None and now - hit[0] < self.ttl_s:
                return list(hit[1])
        try:
            got = self.fetch(round_id)
        except Exception:  # noqa: BLE001 - the board keeps the last good field
            with self._lock:
                old = self._rows.get(round_id)
                return None if old is None or old[1] is None else list(old[1])
        with self._lock:
            self._rows[round_id] = (now, list(got))
            return list(got)


def live_rank(live: Optional[dict], cache: FieldCache) -> Optional[dict]:
    """The replay's best local score placed on its round's published field. Other runs have no live rank."""
    if not isinstance(live, dict):
        return None
    rnd = live.get("round") or {}
    if rnd.get("kind") != "replay":
        return None
    score = (live.get("search") or {}).get("best")
    rid = rnd.get("id")
    if isinstance(score, bool) or not isinstance(score, (int, float)) or not rid:
        return None
    scores = cache.scores(str(rid))
    if not scores:
        return {"score": float(score), "rank": None, "field": 0, "field_best": None, "field_median": None,
                "gap_to_best": None, "filling": True}
    return place(float(score), scores, replay=True)
