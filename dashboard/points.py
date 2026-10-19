"""Every config scored this round, projected on its three main axes for the dashboard's 3D island view
(GET /api/points). It reads the round's result store read-only (cache/results/<tag>.sqlite), so it never touches main.

The projection (mean and basis) goes with the points, so the page projects the live climbers and island bests the
same way. Encodings are cached per config and the whole answer for CACHE_S, so a page polling it costs little.
"""
from __future__ import annotations

import json
import sqlite3
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from common.settings import SETTINGS

CACHE_S = 20.0
SOURCES = ["climber", "ga", "rf", "step_a", "released", "other"]     # the page's colour per source

_lock = threading.Lock()
_space = None
_enc: Dict[str, List[float]] = {}
_enc_round: Optional[str] = None
_cache: Dict[str, tuple] = {}


def _space_once():
    global _space
    if _space is None:
        from common.params import ParamRanges
        from search.space import Space
        try:
            ranges = ParamRanges.from_platform()
        except Exception:  # noqa: BLE001 - offline: the fixture's ranges still give a stable projection
            ranges = ParamRanges.from_file(str(Path(__file__).resolve().parents[1] / "tests" / "fixtures" /
                                               "parameter_ranges.json"))
        _space = Space(ranges)
    return _space


def source_of(tag: str) -> tuple:
    """(source index, island, climber) of a store tag: B<id> is climber id on island id // 100."""
    tag = tag or ""
    if tag.startswith("B") and tag[1:].isdigit():
        cid = int(tag[1:])
        return 0, cid // 100, cid
    return {"G": 1, "R": 2, "A": 3, "D": 4}.get(tag, 5), -1, -1


def points(tag: str, results_dir: Optional[Path] = None) -> dict:
    """The round's scored configs: names, mean, basis (3 x settings), xyz (flat), score, src, island, climber."""
    path = Path(results_dir or Path(SETTINGS.cache_dir) / "results") / f"{tag}.sqlite"
    with _lock:
        hit = _cache.get(tag)
        if hit is not None and time.time() - hit[0] < CACHE_S:
            return hit[1]
        out = _compute(tag, path)
        _cache.clear()
        _cache[tag] = (time.time(), out)
        return out


def _compute(tag: str, path: Path) -> dict:
    global _enc_round
    empty = {"round": tag, "n": 0, "sources": SOURCES}
    if not path.exists():
        return empty
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        rows = db.execute("SELECT config_key, config, MAX(score), tag FROM results "
                          "WHERE score IS NOT NULL AND fidelity = 'full' GROUP BY config_key").fetchall()
    finally:
        db.close()
    if len(rows) < 4:
        return empty
    from common.types import Config
    space = _space_once()
    if _enc_round != tag:
        _enc.clear()
        _enc_round = tag
    vecs = []
    for key, cfg, _score, _tag in rows:
        v = _enc.get(key)
        if v is None:
            try:
                v = _enc[key] = [float(x) for x in space.encode(Config.from_dict(json.loads(cfg)))]
            except (ValueError, KeyError, TypeError):
                v = _enc[key] = [0.0] * len(space.names)
        vecs.append(v)
    X = np.asarray(vecs, dtype=float)
    mean = X.mean(axis=0)
    _u, _s, vt = np.linalg.svd(X - mean, full_matrices=False)
    basis = vt[:3]
    for k in range(basis.shape[0]):                       # a fixed sign, so the view does not flip between refreshes
        if basis[k, np.argmax(np.abs(basis[k]))] < 0:
            basis[k] = -basis[k]
    xyz = (X - mean) @ basis.T
    src = [source_of(t) for _k, _c, _s, t in rows]
    return {"round": tag, "n": len(rows), "names": list(space.names), "sources": SOURCES,
            "mean": [round(float(x), 5) for x in mean], "basis": [[round(float(x), 5) for x in b] for b in basis],
            "xyz": [round(float(x), 4) for x in xyz.ravel()], "score": [round(float(r[2]), 4) for r in rows],
            "src": [s[0] for s in src], "island": [s[1] for s in src], "climber": [s[2] for s in src]}
