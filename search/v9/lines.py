"""Single-setting lines through a config: for each setting, the values already run with every other setting at that
config's values, and the values still to try (search/v9/DESIGN.md, B1 lines and Step C).

    gaps     every untried rung between the config and the nearest tried value on each side (what step growth jumped
             over: 3 -> 4 -> 6 -> 9 lost, so 7 and 8)
    beyond   past the nearest tried value, to the end of the ladder: `coarse` rungs spread evenly on each side
    cheap    the genotyping and PairHMM settings take every rung (they share the config's GATK run)

Step B sends the lines of each grid's settings with the grid; Step C sends every setting's line from each peak.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from common.tuning import tune
from common.types import Config
from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS

CHEAP = set(GENOTYPING_PARAMS) | set(PAIRHMM_PARAMS)
DENSER = (1, 2, 4, 8, 16, 1000)          # far.coarse times this, the largest that still fits the budget


def tried(ctx, peak: Config) -> Dict[str, Dict[int, float]]:
    """setting -> {rung: score} of the scored configs that differ from `peak` in that setting only."""
    names = [n for n in ctx.space.names if ctx.ladders.movable(n)]
    pp = peak.params
    out: Dict[str, Dict[int, float]] = {n: {} for n in names}
    for c, s in ctx.scored():
        diff = None
        for n in names:
            if c.params.get(n) != pp.get(n):
                if diff is not None:
                    diff = "two"
                    break
                diff = n
        if diff is None or diff == "two":
            continue
        out[diff][int(round(ctx.ladders.index(diff, c.params[diff])))] = s
    return out


def spread(lo: int, hi: int, n: int) -> List[int]:
    """Up to n rungs from lo to hi inclusive, evenly spread, both ends included."""
    if hi < lo:
        return []
    if hi - lo + 1 <= n:
        return list(range(lo, hi + 1))
    if n <= 1:
        return [lo]
    return sorted({int(round(lo + k * (hi - lo) / (n - 1))) for k in range(n)})


def plan(ctx, peak: Config, known: Dict[str, Dict[int, float]], coarse: int) -> Dict[str, List[int]]:
    """setting -> the rungs to try: the gaps on both sides, and `coarse` rungs beyond each side's nearest tried one."""
    out: Dict[str, List[int]] = {}
    for n, have in known.items():
        rungs = ctx.ladders.rungs[n]
        p = int(round(ctx.ladders.index(n, peak.params.get(n))))
        last = len(rungs) - 1
        if n in CHEAP:
            want = [i for i in range(len(rungs)) if i != p]
        else:
            want = []
            up = min((i for i in have if i > p), default=None)
            down = max((i for i in have if i < p), default=None)
            if up is None:
                want += [p + 1] + spread(p + 2, last, coarse) if p < last else []
            else:
                want += list(range(p + 1, up)) + spread(up + 1, last, coarse)
            if down is None:
                want += [p - 1] + spread(0, p - 2, coarse) if p > 0 else []
            else:
                want += list(range(down + 1, p)) + spread(0, down - 1, coarse)
        out[n] = sorted({i for i in want if 0 <= i <= last and i != p and i not in have})
    return out


def gap(ctx, peak: Config, known: Dict[str, Dict[int, float]], name: str, rung: int) -> bool:
    """Case 1: `rung` lies between the peak and the nearest tried value on its side."""
    p = int(round(ctx.ladders.index(name, peak.params.get(name))))
    have = known.get(name, {})
    if rung > p:
        up = min((i for i in have if i > p), default=None)
        return up is not None and rung < up
    down = max((i for i in have if i < p), default=None)
    return down is not None and rung > down


def single(ctx, peak: Config, name: str, rung: int) -> Optional[Config]:
    return ctx.space._assemble({**peak.params, name: ctx.ladders.rungs[name][rung]})




def densest(ctx, centre: Config, known: Dict[str, Dict[int, float]], budget: int) -> Dict[str, List[int]]:
    """plan() with the densest coarse spacing whose expensive points (not the cheap riders) fit in `budget`; never
    sparser than far.coarse."""
    base = int(tune("far.coarse"))
    chosen = plan(ctx, centre, known, base)
    for f in DENSER:
        p = plan(ctx, centre, known, base * f)
        if sum(len(v) for n, v in p.items() if n not in CHEAP) > budget:
            break
        chosen = p
    return chosen
