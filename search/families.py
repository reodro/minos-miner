"""Genotyping riders: cheap variants of a config that share its GATK traversal. Owner: search.

A genotyping rider costs ~1% of a run, so a global point (wave 0, a generation, a scout point) carries a few of them:
its calling threshold (the FP-vs-recall knob) and the other genotyping settings, so a new region is judged at its own
best threshold rather than a default one.
"""
from __future__ import annotations

from typing import Callable, List, Optional

from common.types import Config
from search.space import CALL_CONF_GRID, CONTAMINATION_GRID, Space

CONF = "standard_min_confidence_threshold_for_calling"
CONF_STEPS = (3.0, 6.0, 10.0, 15.0, 20.0)


def _num(space: Space, name: str, value: float):
    s = space.specs[name]
    v = min(float(s.max), max(float(s.min), float(value)))
    return int(round(v)) if s.type == "int" else round(v, 6)


def genotyping_riders(space: Space, base: Config, n: int, fresh: Optional[Callable[[Config], bool]] = None) -> List[Config]:
    """Up to `n` genotyping variants of `base`, nearest calling thresholds first, then the global grid, then the
    indel prior, contamination and max alleles. Each differs from `base` in one genotyping axis."""
    if n <= 0:
        return []
    axes = []
    if CONF in space.names:
        s = space.specs[CONF]
        cur = float(base.params.get(CONF, s.default))
        vals: List[float] = []
        for d in CONF_STEPS:                           # alternate above / below, growing
            for v in (cur + d, cur - d):
                v = _num(space, CONF, v)
                if s.min <= v <= s.max and v not in vals and abs(v - cur) > 1e-6:
                    vals.append(v)
        for v in CALL_CONF_GRID:
            if s.min <= v <= s.max and v not in vals and abs(v - cur) > 1e-6:
                vals.append(v)
        axes.append((CONF, vals))
    het = "indel_heterozygosity"
    if het in space.names:
        cur = float(base.params.get(het, space.specs[het].default))
        axes.append((het, [_num(space, het, cur * 2), _num(space, het, cur / 2)]))
    cont = "contamination_fraction_to_filter"
    if cont in space.names:
        axes.append((cont, [v for v in CONTAMINATION_GRID]))
    alt = "max_alternate_alleles"
    if alt in space.names:
        axes.append((alt, [2, 6]))
    # interleave the axes so a short family still covers more than one of them
    rounds = max(len(v) for _, v in axes) if axes else 0
    ordered = []
    for i in range(rounds):
        for name, vals in axes:
            if i < len(vals):
                ordered.append((name, vals[i]))
    out: List[Config] = []
    seen = {base.key}
    for name, v in ordered:
        if name not in space.names or space._same(space.specs[name], base.params.get(name, space.specs[name].default), v):
            continue
        c = space._assemble({**base.params, name: v})
        if c is None or c.key in seen or (fresh is not None and not fresh(c)):
            continue
        seen.add(c.key)
        out.append(c)
        if len(out) >= n:
            break
    return out
