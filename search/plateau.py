"""A second score that only orders the plateau probe (2026-10-14). Owner: search.

The blended score stays the only number that accepts a move; this one is consulted only when the score is exactly
equal (a flat), to choose which tied config the probe is centred on and which rung of a locked setting goes first.
Background: Brockhoff, Friedrich, Hebbinghaus, Klein, Neumann and Zitzler (2009), "On the Effects of Adding Objectives
to Plateau Functions": a second measure removes a flat, and helps only when it points at the exit; added into the
score or allowed to accept a move, it can reverse a real gain. So it never does either here.

The second score is a distance to an interval, one interval per setting, estimated only from runs that already moved
the score:
  slice     one setting varying with every other setting fixed (scored configs grouped by the other settings).
            Kept only when its highest and lowest scores differ by more than SLICE_GAP (0.10, the gap that already
            closes a useless probe): a slice inside that gap is flat and says nothing about where a maximum sits.
  gain      per rung of a kept slice, its score minus the slice's own centre, so a strong centre does not make every
            rung look good. The stored runs do not say which config a line was sent from; the centre here is the
            slice's median score.
  interval  per setting, every rung whose average gain over the kept slices is within NEAR_BEST (0.05, the margin
            the compass treats as no real face) of the best average.
  score     minus the sum over settings of how many rungs a config sits outside that setting's interval, each divided
            by the length of the ladder. 0 inside every interval. A setting with no kept slice is left out: no
            evidence of an interval, so no direction is invented.
Use (search/climb.py, plateau probe): on a flat, the probe is centred on the tied config with the highest second
score; on a locked setting the first rung sent is the one toward the setting's interval, and the opposite rung is
still sent in the same wave when a worker is free. It never closes a direction and never accepts a move.
Limit: it points at where one-setting moves have already raised the score. A higher region no one-setting move has
reached is invisible to it.
"""
from __future__ import annotations

import statistics
import time
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from common.tuning import tune
from common.types import Config

SLICE_GAP = 0.10
NEAR_BEST = 0.05


def _canon(v):
    return round(v, 9) if isinstance(v, float) else v


def intervals(ctx) -> Dict[str, Set[int]]:
    """Per setting, the rung indices of its interval (only settings with at least one kept slice)."""
    scored = list(ctx.scored())
    cache = getattr(ctx, "_plateau_iv", None)
    if cache is not None and cache[0] == len(scored) and time.time() - cache[1] < 60:
        return cache[2]
    names = [n for n in ctx.space.names if ctx.ladders.movable(n)]
    out: Dict[str, Set[int]] = {}
    for n in names:
        others = [m for m in ctx.space.names if m != n]
        groups: Dict[tuple, Dict[int, float]] = defaultdict(dict)
        for c, s in scored:
            if n not in c.params:
                continue
            key = tuple(_canon(c.params.get(m)) for m in others)
            i = int(round(ctx.ladders.index(n, c.params[n])))
            g = groups[key]
            if i not in g or s > g[i]:
                g[i] = float(s)
        gains: Dict[int, List[float]] = defaultdict(list)
        for g in groups.values():
            if len(g) < 2 or max(g.values()) - min(g.values()) <= SLICE_GAP:
                continue
            centre = statistics.median(g.values())
            for i, s in g.items():
                gains[i].append(s - centre)
        if not gains:
            continue
        avg = {i: statistics.mean(v) for i, v in gains.items()}
        top = max(avg.values())
        out[n] = {i for i, a in avg.items() if a >= top - NEAR_BEST}
    ctx._plateau_iv = (len(scored), time.time(), out)
    return out


def _outside(i: int, interval: Set[int]) -> int:
    return min(abs(i - j) for j in interval) if interval else 0


def second_score(ctx, config: Config, iv: Optional[Dict[str, Set[int]]] = None) -> float:
    iv = intervals(ctx) if iv is None else iv
    total = 0.0
    for n, interval in iv.items():
        if n not in config.params:
            continue
        rungs = ctx.ladders.rungs.get(n) or []
        i = int(round(ctx.ladders.index(n, config.params[n])))
        total += _outside(i, interval) / max(1, len(rungs))
    return -total


def direction(ctx, centre: Config, name: str, iv: Dict[str, Set[int]]) -> int:
    """+1 / -1: the side of `name`'s interval from the centre; 0 inside it or with no interval."""
    interval = iv.get(name)
    if not interval or name not in centre.params:
        return 0
    i = int(round(ctx.ladders.index(name, centre.params[name])))
    if i in interval:
        return 0
    nearest = min(interval, key=lambda j: (abs(j - i), j))
    return 1 if nearest > i else -1


def enabled() -> bool:
    return bool(int(tune("climb.plateau_interval")))


def centre_of_flat(ctx, base: float, default: Config) -> Config:
    """On a flat (configs scoring exactly `base`), the tied config with the highest second score; `default` on a tie."""
    if not enabled():
        return default
    tied = [c for c, s in ctx.scored() if s == base]
    if len(tied) < 2:
        return default
    iv = intervals(ctx)
    if not iv:
        return default
    best_c, best_s = default, second_score(ctx, default, iv)
    for c in tied:
        s = second_score(ctx, c, iv)
        if s > best_s + 1e-12:
            best_c, best_s = c, s
    if best_c.key != default.key:
        ctx.say("plateau: %d configs tie at %.4f; the probe centres on %s (second score %.3f)",
                len(tied), base, best_c.key[:8], best_s)
    return best_c
