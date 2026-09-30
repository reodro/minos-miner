"""Append-only result store: one row per scored callset. The single source of truth for search history,
re-scoring, backtests and portfolio building. Search may read it; only the round loop writes."""
from __future__ import annotations

import json
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
