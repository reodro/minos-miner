"""Append-only result store: one row per scored callset. The single source of truth for search history,
re-scoring, backtests and portfolio building. Search may read it; only the round loop writes."""
from __future__ import annotations

import glob
import json
import os
import re
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
CREATE TABLE IF NOT EXISTS sketches (
    round_id TEXT, config_key TEXT, fidelity TEXT, sketch TEXT, created REAL
);
CREATE INDEX IF NOT EXISTS sketches_key ON sketches (round_id, config_key, fidelity);
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

    def add_sketch(self, round_id: str, config_key: str, fidelity: str, sketch: Optional[dict]) -> None:
        """The callset sketch beside the score. Older rows have none. Equality of sketches is observed
        sketch equivalence, not a guarantee that the full callsets are identical."""
        if not sketch:
            return
        self.db.execute(
            "INSERT INTO sketches VALUES (?,?,?,?,?)",
            (round_id, config_key, fidelity, json.dumps(sketch), time.time()))
        self.db.commit()

    def results(self, round_id: str, answer_set_version: Optional[int] = None) -> List[ScoreResult]:
        q, args = "SELECT result FROM results WHERE round_id=?", [round_id]
        if answer_set_version is not None:
            q += " AND answer_set_version=?"
            args.append(answer_set_version)
        return [ScoreResult.from_dict(json.loads(r[0])) for r in self.db.execute(q, args)]

    def portfolio(self, limit: int = 4, tool: str = "gatk") -> List[Config]:
        """Best config of each EARLIER round's database in this file's directory, latest round first.

        Files are named <round tag>.sqlite; only rounds before this one count, so a replay of an old round never
        sees configs found in later rounds (a live round only has earlier ones anyway). Scores are not returned:
        a new round has a new donor and region, so only the configs are reused. The in-memory store contributes
        nothing.
        """
        if limit <= 0 or not self.path or self.path == ":memory:":
            return []
        me = os.path.abspath(self.path)
        tag = _round_tag(me)
        if tag is None:
            return []                                   # not a round file: we cannot tell earlier from later
        directory = os.path.dirname(me)
        earlier = [(t, p) for p in glob.glob(os.path.join(directory, "*.sqlite"))
                   for t in [_round_tag(p)] if t is not None and t < tag]
        earlier.sort(reverse=True)
        out: List[Config] = []
        seen = set()
        for _, path in earlier:
            config = _best_config(path, tool)
            if config is not None and config.key not in seen:
                seen.add(config.key)
                out.append(config)
            if len(out) >= limit:
                break
        return out


def clear_round(directory: str, tag: str) -> int:
    """Remove one round tag's saved results, including sqlite sidecars. A replay starts empty: the dashboard
    reads this file, and an earlier replay of the same tag would otherwise appear as this run's top score.
    Other rounds in the directory stay."""
    if not re.fullmatch(r"\d{8}T\d{6}Z", tag):
        raise ValueError(f"not a round tag: {tag}")
    removed = 0
    for path in glob.glob(os.path.join(directory, tag + ".sqlite*")):
        base = os.path.basename(path)
        if base != tag + ".sqlite" and not base.startswith(tag + ".sqlite-"):
            continue
        os.remove(path)
        removed += 1
    return removed


def _round_tag(path: str) -> Optional[str]:
    """The round tag a results file is named after (YYYYMMDDTHHMMSSZ), or None."""
    stem = os.path.basename(path)[:-len(".sqlite")] if path.endswith(".sqlite") else ""
    return stem if re.fullmatch(r"\d{8}T\d{6}Z", stem) else None


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
