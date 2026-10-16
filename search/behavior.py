"""Separate the questions a scored configuration can answer.

Parameter novelty, behavioral novelty, and score novelty are different. A sketch match is observed
equivalence of the stored bottom-k signatures, not a proof that two full callsets are identical.
"""
from __future__ import annotations

import json
import sqlite3
from typing import Optional


def sketch_key(sketch: Optional[dict]) -> Optional[str]:
    if not sketch:
        return None
    return json.dumps(sketch.get("h") or sketch, sort_keys=True)


def summarize(path: str) -> dict:
    """Counts for one results file. Rows without a sketch are omitted from the behavioral counts."""
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        rows = con.execute(
            "SELECT s.config_key, s.sketch, r.score, r.runtime_s "
            "FROM sketches s LEFT JOIN results r ON r.config_key = s.config_key AND r.fidelity = s.fidelity "
            "WHERE s.sketch IS NOT NULL"
        ).fetchall()
    except sqlite3.Error:
        rows = []
    finally:
        con.close()
    seen = {}
    new_behavior = 0
    same_behavior = 0
    for key, raw, score, runtime in rows:
        token = sketch_key(json.loads(raw))
        if token in seen:
            same_behavior += 1
        else:
            new_behavior += 1
            seen[token] = (score, runtime, key)
    return {
        "sketched": len(rows),
        "distinct_sketches": len(seen),
        "repeat_sketches": same_behavior,
        "new_sketches": new_behavior,
    }
