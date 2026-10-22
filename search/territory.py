"""Territories from a scan (generalized Step A, 2026-10-22): which scanned configs climb to the same top.

User's definition: a territory is the set of points that Step B's climb takes to the same vertex. A sparse scan cannot
run that climb from every point, so the territories are read off the scan by nearest-better clustering (Preuss,
"Niching the CMA-ES via nearest-better clustering", 2010): every config links to its nearest config with a higher
score; links much longer than usual (> PHI x the mean link) cross a valley and are cut; each remaining tree is one
territory, its root the best point. Distances are on rung positions (each setting scaled to [0, 1]; a category or a
bool counts 0 or 1), so they follow the steps a climber takes.
"""
from __future__ import annotations

import copy
import dataclasses
from typing import Dict, List, Sequence, Tuple

import numpy as np

from common.types import Config

PHI = 2.0
# The settings whose distance decides a territory: the ones that shape hills (today's key settings, plus the region
# padding). Over all 19 settings a sparse even scan has links of nearly one length and nothing is cut (simulated
# 2026-10-22: 2057 configs -> 1 territory at phi 2; over these: 25, top roots 92.19 / 92.00 / 91.96 / 91.80 / 91.74).
NAMES = ("pcr_indel_model", "min_pruning", "base_quality_score_threshold", "active_probability_threshold",
         "min_base_quality_score", "phred_scaled_global_read_mismapping_rate", "assembly_region_padding")


def vectors(points: Sequence[Tuple[Config, float]], values: Dict[str, list]) -> Tuple[np.ndarray, np.ndarray]:
    """(numeric part, categorical part) of every config: rung positions in [0, 1], and category indices."""
    num_names = [n for n, v in values.items() if len(v) > 2 and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                                                                      for x in v)]
    cat_names = [n for n in values if n not in num_names]
    num = np.zeros((len(points), len(num_names)))
    cat = np.zeros((len(points), len(cat_names)), dtype=int)
    for i, (c, _) in enumerate(points):
        for j, n in enumerate(num_names):
            vals, v = values[n], c.params.get(n)
            k = vals.index(v) if v in vals else int(np.argmin([abs(float(x) - float(v)) for x in vals]))
            num[i, j] = k / max(1, len(vals) - 1)
        for j, n in enumerate(cat_names):
            vals, v = values[n], c.params.get(n)
            cat[i, j] = vals.index(v) if v in vals else -1
    return num, cat


def territories(points: Sequence[Tuple[Config, float]], values: Dict[str, list], phi: float = PHI) -> List[dict]:
    """Territories, best first: {"root": Config, "score": float, "members": [(Config, score), ...], "size": int}."""
    pts = sorted(points, key=lambda t: -t[1])
    n = len(pts)
    if n == 0:
        return []
    num, cat = vectors(pts, values)
    parent = [-1] * n
    dist = np.zeros(n)
    for i in range(1, n):                                    # nearest config with a higher (or equal, earlier) score
        d = np.sqrt(((num[:i] - num[i]) ** 2).sum(axis=1) + (cat[:i] != cat[i]).sum(axis=1))
        j = int(np.argmin(d))
        parent[i], dist[i] = j, float(d[j])
    links = dist[1:]
    cut = phi * float(links.mean()) if len(links) else 0.0
    root = list(range(n))
    for i in range(1, n):                                    # parents come first (higher score), so one pass is enough
        root[i] = i if dist[i] > cut else root[parent[i]]
    groups: Dict[int, List[int]] = {}
    for i, r in enumerate(root):
        groups.setdefault(r, []).append(i)
    out = [{"root": pts[r][0], "score": pts[r][1], "members": [pts[i] for i in ms], "size": len(ms)}
           for r, ms in groups.items()]
    return sorted(out, key=lambda t: -t["score"])


def band(values: Dict[str, list], members: Sequence[Tuple[Config, float]], names: Sequence[str],
         widen: int = 2) -> Dict[str, Tuple[object, object]]:
    """Each named setting's range over the territory's members, widened by `widen` rungs either side: (lo, hi)."""
    out = {}
    for n in names:
        vals = values.get(n)
        if not vals:
            continue
        idx = [vals.index(c.params.get(n)) for c, _ in members if c.params.get(n) in vals]
        if not idx:
            continue
        lo, hi = max(0, min(idx) - widen), min(len(vals) - 1, max(idx) + widen)
        out[n] = (vals[lo], vals[hi])
    return out


def restrict(space, values: Dict[str, list], members: Sequence[Tuple[Config, float]], names: Sequence[str],
             widen: int = 2):
    """A copy of `space` in which each named setting keeps only the territory's values (a category: the ones its
    members use; a number: its members' range widened by `widen` rungs). The island climbs inside it; the other
    settings stay free."""
    new = copy.copy(space)
    new.specs = dict(space.specs)
    rng = band(values, members, [n for n in names if n in values])
    for n in names:
        s, vals = new.specs.get(n), values.get(n)
        if s is None or not vals:
            continue
        if s.allowed and not all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in vals):
            used = tuple(v for v in s.allowed if any(c.params.get(n) == v for c, _ in members))
            new.specs[n] = dataclasses.replace(s, allowed=used or s.allowed)
        elif n in rng:
            lo, hi = rng[n]
            if s.allowed:
                new.specs[n] = dataclasses.replace(s, allowed=tuple(v for v in s.allowed if lo <= v <= hi) or s.allowed)
            else:
                new.specs[n] = dataclasses.replace(s, min=lo, max=hi)
    new.api = {n: {"min": s.min, "max": s.max, "allowed": s.allowed, "type": s.type} for n, s in new.specs.items()}
    new.names = [n for n, s in sorted(new.specs.items()) if not new._fixed(s) and n not in new.pinned and n in space.names]
    return new
