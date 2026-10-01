"""What the dashboard shows: main writes a live snapshot and a per-round history; the dashboard only reads them.

    <data_dir>/status/live.json      the current round, rewritten every few seconds (atomic replace)
    <data_dir>/status/rounds.jsonl   one line per finished round (+ "true" lines when a true score is known)

Nothing here can slow down or break mining: every method swallows its own errors. Node names are anonymized
("node 1", ...) because the dashboard is public.
"""
from __future__ import annotations

import collections
import json
import logging
import os
import threading
import time
from pathlib import Path
from typing import Optional

from common.settings import SETTINGS

log = logging.getLogger(__name__)
STATUS_DIR = Path(SETTINGS.data_dir) / "status"


class _Recent(logging.Handler):
    """Keeps the last warnings/errors for the dashboard's log panel."""

    def __init__(self, n: int = 40):
        super().__init__(logging.WARNING)
        self.lines = collections.deque(maxlen=n)

    def emit(self, record):
        try:
            self.lines.append([record.created, record.levelname, record.name, record.getMessage()[:300]])
        except Exception:  # noqa: BLE001
            pass


class StatusBoard:
    def __init__(self, directory: Path = STATUS_DIR, every_s: float = 5.0):
        self.dir, self.every_s = directory, every_s
        self.dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._last = 0.0
        self.state: dict = {}
        self._recent = _Recent()
        logging.getLogger().addHandler(self._recent)

    # ---- round lifecycle (called by mine_round) -------------------------------------------------------------------
    def round_started(self, round, deadline_ts: Optional[float]) -> None:
        self._safe(self._started, round, deadline_ts)

    def _started(self, round, deadline_ts):
        self.state = {"round": {"id": round.round_id, "tag": round.tag, "region": round.region,
                                "deadline_ts": deadline_ts, "started_ts": time.time(), "mode": "building answer set"},
                      "answer_set": None, "search": {"evaluated": 0, "failed": 0, "in_flight": 0, "best": None,
                                                     "best_timeline": [], "eval_timeline": []},
                      "submissions": [], "executor": None}
        self._write(force=True)

    def answer_set(self, a) -> None:
        self._safe(self._answer_set, a)

    def _answer_set(self, a):
        n = a.notes or {}
        self.state["answer_set"] = {
            "regime": a.regime.value, "donor": a.donor, "germline": len(a.germline), "synthetic": len(a.synthetic),
            "n_targets": a.n_targets, "n_targets_source": n.get("n_targets_source"), "seconds": n.get("seconds"),
            "degraded": n.get("degraded")}
        self.state["round"]["mode"] = "searching" if a.germline else "static fallback"
        self._write(force=True)

    def tick(self, client, out, executor) -> None:
        """Called from the round loop; writes at most every `every_s` seconds."""
        now = time.time()
        if now - self._last < self.every_s:
            return
        self._safe(self._tick, client, out, executor, now)

    def _tick(self, client, out, executor, now):
        s = self.state.setdefault("search", {})
        best = client.best()
        s.update(evaluated=out.evaluated, failed=out.failed, in_flight=client.in_flight,
                 best=best.score if best else None, best_config=best.config_key if best else None)
        if best:
            s.setdefault("best_timeline", []).append([now, best.score])
        s.setdefault("eval_timeline", []).append([now, out.evaluated])
        for k in ("best_timeline", "eval_timeline"):
            s[k] = s[k][-2000:]
        self.state["submissions"] = [[t, key, score] for t, key, score in out.submissions]
        stats = getattr(executor, "stats", None)
        self.state["executor"] = stats() if stats else None
        self._write()

    def round_ended(self, mined) -> None:
        self._safe(self._ended, mined)

    def _ended(self, mined):
        r = self.state.get("round", {})
        out = mined.outcome
        r["mode"] = f"finished ({mined.mode})"
        row = {"type": "round", "round_id": mined.round_id, "tag": r.get("tag"), "region": r.get("region"),
               "mode": mined.mode, "note": mined.note, "started": r.get("started_ts"), "ended": time.time(),
               "evaluated": out.evaluated if out else 0, "failed": out.failed if out else 0,
               "best_local": out.best.score if out and out.best else None,
               "final_config": out.submissions[-1][1] if out and out.submissions else None,
               "submissions": len(out.submissions) if out else 0,
               "regime": (self.state.get("answer_set") or {}).get("regime"),
               "donor": (self.state.get("answer_set") or {}).get("donor")}
        self._append(row)
        self._write(force=True)

    def true_score(self, round_id: str, score: float, rank: int, field_size: int) -> None:
        """When the released truth is known (replay, or later backtest): the submitted config's true score."""
        self._safe(self._append, {"type": "true", "round_id": round_id, "true_score": score, "rank": rank,
                                  "field": field_size})

    # ---- io -----------------------------------------------------------------------------------------------------------
    def _write(self, force: bool = False):
        with self._lock:
            self.state["updated"] = time.time()
            self.state["log"] = list(self._recent.lines)
            tmp = self.dir / "live.json.tmp"
            tmp.write_text(json.dumps(self.state))
            os.replace(tmp, self.dir / "live.json")
            self._last = time.time()

    def _append(self, row: dict):
        with self._lock, open(self.dir / "rounds.jsonl", "a") as f:
            f.write(json.dumps(row) + "\n")

    @staticmethod
    def _safe(fn, *args):
        try:
            fn(*args)
        except Exception as e:  # noqa: BLE001 - the dashboard must never break mining
            log.debug("status board: %s", e)
