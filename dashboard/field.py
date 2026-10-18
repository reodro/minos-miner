"""Where a replay's current best would sit on the published field.

The dashboard computes this. Main keeps mining; a slow or failed field read never reaches the search.
"""
from __future__ import annotations

import json
import os
import statistics
import threading
import time
from pathlib import Path
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


def our_hotkey() -> Optional[str]:
    """The live miner's hotkey (ss58): $MINOS_HOTKEY_SS58, else the public file of wallet $MINOS_WALLET (sn107) /
    hotkey $MINOS_HOTKEY (miner). None when neither exists."""
    if os.environ.get("MINOS_HOTKEY_SS58"):
        return os.environ["MINOS_HOTKEY_SS58"]
    pub = (Path.home() / ".bittensor" / "wallets" / os.environ.get("MINOS_WALLET", "sn107") / "hotkeys" /
           f"{os.environ.get('MINOS_HOTKEY', 'miner')}pub.txt")
    try:
        return json.loads(pub.read_text())["ss58Address"]
    except (OSError, ValueError, KeyError):
        return None


def official_result(round_id: str, hotkey: Optional[str]) -> Optional[dict]:
    """Our submission's score as the validators computed it, and its rank among miners (one row per miner: the mean
    over its validators), once the platform has revealed who submitted what. None before that.
    2026-10-18: the validators' calls for our configs differed from ours (rounds 23:54 and 01:06: 0.31 and 0.89
    lower), so a live round's place comes from here, not from our own TRUE score."""
    if not hotkey:
        return None
    from common.platform import Platform
    payload = Platform().scores(round_id)
    if not payload.get("identity_revealed"):
        return None
    by: dict = {}
    for row in payload.get("scores") or []:
        hk, final = row.get("miner_hotkey"), row.get("combined_final")
        if hk and final is not None:
            by.setdefault(hk, []).append(100.0 * float(final))
    if hotkey not in by:
        return None
    means = {k: statistics.mean(v) for k, v in by.items()}
    ours, scores = means[hotkey], list(means.values())
    best, med = max(scores), statistics.median(scores)
    return {"official": True, "official_score": round(ours, 4), "rank": 1 + sum(1 for x in scores if x > ours),
            "field": len(scores), "field_best": round(best, 4), "field_median": round(med, 4),
            "gap_to_best": round(best - ours, 4), "filling": False}


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


def refresh_rank(row: Optional[dict], cache: "FieldCache", *, replay: bool) -> Optional[dict]:
    """A review's rank row placed again on the field as published now. A live round's truth is out when it closes,
    an hour before its field: main ranked it then against nothing ("1 of 0") and never looked again (2026-10-18).
    No published field yet: no rank, the board waits."""
    if not isinstance(row, dict):
        return row
    true, rid = row.get("true_score"), row.get("round_id")
    if isinstance(true, bool) or not isinstance(true, (int, float)) or not rid:
        return row
    if not replay:                                 # a live round: the validators' own score once it is revealed
        try:
            off = _OFFICIAL.get(str(rid))
        except Exception:  # noqa: BLE001 - the platform being down keeps our own estimate on screen
            off = None
        if off:                                    # shown score = the validators'; our own estimate kept beside it
            return {**row, **off, "our_true": row.get("true_score"), "true_score": off["official_score"]}
    scores = cache.scores(str(rid))
    if not scores:
        return {**row, "rank": None, "field": 0, "field_best": None, "field_median": None, "gap_to_best": None,
                "filling": True}
    p = place(float(true), scores, replay=replay)
    return {**row, **{k: p[k] for k in ("rank", "field", "field_best", "field_median", "gap_to_best", "filling")}}


class _OfficialCache:
    def __init__(self, ttl_s: float = 120.0):
        self.ttl_s, self._rows, self._lock = ttl_s, {}, threading.Lock()

    def get(self, round_id: str) -> Optional[dict]:
        now = time.time()
        with self._lock:
            hit = self._rows.get(round_id)
            if hit and now - hit[0] < self.ttl_s:
                return hit[1]
        val = official_result(round_id, our_hotkey())
        with self._lock:
            self._rows[round_id] = (now, val)
        return val


_OFFICIAL = _OfficialCache()


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
