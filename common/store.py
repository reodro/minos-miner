"""Append-only result store: one row per scored callset. The single source of truth for search history,
re-scoring, backtests and portfolio building. Search may read it; only the round loop writes."""
from __future__ import annotations

import glob
import json
import os
import sqlite3
import time
from typing import List, Optional

from common.types import Config, ScoreResult

_SCHEMA = """
CREATE TABLE IF NOT EXISTS results (
    round_id TEXT, config_key TEXT, config TEXT, fidelity TEXT, answer_set_version INTEGER,
    score REAL, score_sd REAL, gate_passed INTEGER, result TEXT, runtime_s REAL, tag TEXT, created REAL
);
CREATE INDEX IF NOT EXISTS results_round ON results (round_id, answer_set_version);
"""


class ResultStore:
    def __init__(self, path: str = ":memory:"):
        self.path = path
        self.db = sqlite3.connect(path, check_same_thread=False)
        self.db.executescript(_SCHEMA)

    def add(self, result: ScoreResult, config: Config, runtime_s: float = 0.0, tag: str = "") -> None:
        self.db.execute(
            "INSERT INTO results VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (result.round_id, result.config_key, json.dumps(config.to_dict()), result.fidelity,
             result.answer_set_version, result.score, result.score_sd,
             None if result.gate_passed is None else int(result.gate_passed),
             json.dumps(result.to_dict()), runtime_s, tag, time.time()))
        self.db.commit()

    def results(self, round_id: str, answer_set_version: Optional[int] = None) -> List[ScoreResult]:
        q, args = "SELECT result FROM results WHERE round_id=?", [round_id]
        if answer_set_version is not None:
            q += " AND answer_set_version=?"
            args.append(answer_set_version)
        return [ScoreResult.from_dict(json.loads(r[0])) for r in self.db.execute(q, args)]

    def portfolio(self, limit: int = 4, tool: str = "gatk") -> List[Config]:
        """Best config of each other round's database in this file's directory, newest file first.

        Scores are not returned: a new round has a new donor and region, so only the configs are reused.
        The in-memory store and this round's own file contribute nothing.
        """
        if limit <= 0 or not self.path or self.path == ":memory:":
            return []
        directory = os.path.dirname(os.path.abspath(self.path))
        me = os.path.abspath(self.path)
        paths = [p for p in glob.glob(os.path.join(directory, "*.sqlite")) if os.path.abspath(p) != me]
        paths.sort(key=os.path.getmtime, reverse=True)
        out: List[Config] = []
        seen = set()
        for path in paths:
            config = _best_config(path, tool)
            if config is not None and config.key not in seen:
                seen.add(config.key)
                out.append(config)
            if len(out) >= limit:
                break
        return out


def _best_config(path: str, tool: str) -> Optional[Config]:
    try:
        db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    except sqlite3.Error:
        return None
    try:
        row = db.execute(
            """SELECT config FROM results
               WHERE score IS NOT NULL AND (gate_passed IS NULL OR gate_passed != 0)
               ORDER BY score DESC LIMIT 1"""
        ).fetchone()
    except sqlite3.Error:
        return None
    finally:
        db.close()
    if row is None:
        return None
    config = Config.from_dict(json.loads(row[0]))
    return config if config.tool == tool else None
