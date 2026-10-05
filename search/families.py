"""Rider families: the cheap variants that share a traversal with a skeleton. Owner: search.

A genotyping rider costs ~1% of a run, so every traversal carries a designed sweep of the calling threshold (the
FP-vs-recall knob) and the other genotyping axes; the skeleton's score is the best of its family. PairHMM riders
(~12%) are added in the finish phase, when latency no longer matters.
"""
from __future__ import annotations

from typing import Callable, List, Optional, Sequence, Tuple

import random

from common.types import Config
from search.space import ASSEMBLY_PARAMS, CALL_CONF_GRID, CONTAMINATION_GRID, PAIRHMM_PARAMS, Space

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
        for c in space._axis_variants(base, [(name, (v,))]):
            if c.key in seen or (fresh is not None and not fresh(c)):
                continue
            seen.add(c.key)
            out.append(c)
            break
        if len(out) >= n:
            break
    return out


STRUCTURAL = tuple(ASSEMBLY_PARAMS) + tuple(PAIRHMM_PARAMS)    # branchable, yet they change the calls' structure


def structural_riders(space: Space, base: Config, n: int, rng: random.Random,
                      fresh: Optional[Callable[[Config], bool]] = None) -> List[Config]:
    """Up to `n` variants of `base` on the branchable structural axes (assembly: pruning, dangling branches,
    haplotypes; PairHMM: base-quality threshold, PCR model, gap penalty, mismapping), one or two axes each: they
    share `base`'s traversal at 0.12-0.55 of a run and are structurally different configs - width for the price
    of riders."""
    pool = [a for a in STRUCTURAL if a in space.names]
    out: List[Config] = []
    seen = {base.key}
    tries = 0
    while len(out) < n and tries < 10 * max(1, n):
        tries += 1
        c = space.neighbour(base, rng, n_params=rng.choice((1, 1, 2)), scale=0.3, pool=pool)
        if c is None or c.key in seen or (fresh is not None and not fresh(c)):
            continue
        seen.add(c.key)
        out.append(c)
    return out


def screen_moves(space: Space, centre: Config, rng: random.Random, n: int, names: Sequence[str],
                 fresh: Optional[Callable[[Config], bool]] = None, k: Tuple[int, int] = (3, 8),
                 scale: float = 0.2) -> List[Config]:
    """Up to `n` configs that move 3-8 of `names` at once from `centre` by moderate steps: the legacy project's
    multi-axis screen, which found optima that no single step reaches (a padding step alone costs FP; with an
    FP-lean move elsewhere it wins)."""
    out: List[Config] = []
    seen = {centre.key}
    live = [a for a in names if a in space.names]
    need = min(k[0], len(live))
    tries = 0
    while len(out) < n and tries < 20 * max(1, n) and live:
        tries += 1
        c = space.neighbour(centre, rng, n_params=rng.randint(need, min(k[1], len(live))), scale=scale, pool=live)
        if c is None or c.key in seen or (fresh is not None and not fresh(c)):
            continue
        if sum(c.params.get(a) != centre.params.get(a) for a in live) < need:   # a step may land on its old value
            continue
        seen.add(c.key)
        out.append(c)
    return out


def pairhmm_riders(space: Space, base: Config, n: int, fresh: Optional[Callable[[Config], bool]] = None) -> List[Config]:
    """Up to `n` PairHMM variants of `base`: the other PCR models, then small steps of the three numeric axes."""
    axes: List[tuple] = []
    pcr = "pcr_indel_model"
    if pcr in space.names:
        axes.append((pcr, tuple(space.specs[pcr].allowed)))
    for name, steps in (("base_quality_score_threshold", (3, -3)), ("pair_hmm_gap_continuation_penalty", (3, -3)),
                        ("phred_scaled_global_read_mismapping_rate", (5, -5))):
        if name in space.names:
            cur = float(base.params.get(name, space.specs[name].default))
            axes.append((name, tuple(_num(space, name, cur + d) for d in steps)))
    out: List[Config] = []
    seen = {base.key}
    for c in space._axis_variants(base, axes):
        if c.key in seen or (fresh is not None and not fresh(c)):
            continue
        seen.add(c.key)
        out.append(c)
        if len(out) >= n:
            break
    return out
