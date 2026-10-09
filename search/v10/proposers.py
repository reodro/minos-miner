"""The proposers whose candidates the forest ranks: QGA, GA-crossover, InjectDrain and local moves.

All of them work on the settings' ladders (search/v9/steps.py), so a move can be one rung: the legacy proposers
(8 bins per setting, Gaussian mutation of 10% of the range) could not express the 1-4 unit moves that won.
Settings QLock has frozen keep the base config's value.

    QGA     a probability vector per setting over min(rf.qga_bins, rungs) bins of its ladder; a candidate draws a bin
            per setting, then a rung inside it. Rotation moves rf.qga_theta of the mass (decaying to rf.qga_theta_end)
            toward the bins of every candidate (engine.py), not the single best.
    GA      tournament selection (3) from the elites, uniform crossover, then each setting mutated with probability
            rf.ga_mut_rate by 1..rf.ga_mut_rungs rungs.
    Inject  fresh Latin-hypercube configs over the live settings (InjectDrain's inflow; the drain is the elites' cap).
    Local   around each candidate: every setting's line (search/v9/lines.py: the gaps
            and coarse rungs beyond the nearest tried value) and random 2-3 setting moves of 1-3 rungs. The legacy
            lesson: the model finds the heads, sweeping them finds the gains.
"""
from __future__ import annotations

import random
from typing import Dict, List, Optional, Sequence, Tuple

from common.types import Config
from search.axes import latin_hypercube
from search.space import PAIRHMM_PARAMS
from search.v9.lines import CHEAP, plan, single, tried  # noqa: F401 (CHEAP: the engine's riders)

PAIRHMM = set(PAIRHMM_PARAMS)

Scored = Tuple[Config, float]


class QGA:
    def __init__(self, space, ladders, bins: int, rng: random.Random):
        self.space, self.ladders, self.rng = space, ladders, rng
        self.names = [n for n in space.names if ladders.movable(n)]
        self.nb = {n: max(1, min(bins, len(ladders.rungs[n]))) for n in self.names}
        self.amp: Dict[str, List[float]] = {}
        self.reset()

    def reset(self) -> None:
        self.amp = {n: [1.0 / self.nb[n]] * self.nb[n] for n in self.names}

    def bin_of(self, name: str, value) -> int:
        """The bin whose rungs (rungs_of) contain this value: the exact inverse of rungs_of."""
        k, nb = len(self.ladders.rungs[name]), self.nb[name]
        i = int(round(self.ladders.index(name, value)))
        return min(nb - 1, ((i + 1) * nb - 1) // k)

    def rungs_of(self, name: str, b: int) -> range:
        k, nb = len(self.ladders.rungs[name]), self.nb[name]
        lo, hi = (b * k) // nb, ((b + 1) * k) // nb
        return range(lo, max(lo + 1, hi))

    def sample(self, n: int, base: Config, live: Sequence[str]) -> List[Config]:
        out = []
        for _ in range(n):
            p = dict(base.params)
            for name in live:
                b = self.rng.choices(range(self.nb[name]), weights=self.amp[name], k=1)[0]
                p[name] = self.ladders.rungs[name][self.rng.choice(self.rungs_of(name, b))]
            c = self.space._assemble(p)
            if c is not None:
                out.append(c)
        return out

    def rotate(self, elites: Sequence[Config], theta: float, mut: float, floor: float) -> None:
        if not elites:
            return
        share = theta / len(elites)
        for name in self.names:
            w = list(self.amp[name])
            for e in elites:
                j = self.bin_of(name, e.params.get(name))
                w = [(1 - share) * x + (share if i == j else 0.0) for i, x in enumerate(w)]
            if len(w) > 1 and self.rng.random() < mut:
                i, j = self.rng.sample(range(len(w)), 2)
                w[i], w[j] = w[j], w[i]
            f = floor * 8.0 / len(w)                         # the same uniform mass at any bin count
            w = [max(f, x) for x in w]
            t = sum(w)
            self.amp[name] = [x / t for x in w]


def shift(ladders, name: str, value, k: int, space):
    """`value` moved k rungs, clamped to the ladder; a bool flips."""
    if space.specs[name].type == "bool":
        return not value
    v = ladders.shift(name, value, k)
    if v is None:
        rungs = ladders.rungs[name]
        v = rungs[-1] if k > 0 else rungs[0]
    return v


def breed(space, ladders, rng: random.Random, elites: Sequence[Scored], n: int, base: Config, live: Sequence[str],
          mut_rate: float, mut_rungs: int) -> List[Config]:
    if len(elites) < 2:
        return []
    out = []
    for _ in range(n):
        pa = max(rng.sample(list(elites), min(3, len(elites))), key=lambda t: t[1])[0]
        pb = max(rng.sample(list(elites), min(3, len(elites))), key=lambda t: t[1])[0]
        p = dict(base.params)
        for name in live:
            p[name] = (pa if rng.random() < 0.5 else pb).params.get(name, p.get(name))
            if rng.random() < mut_rate:
                k = rng.randint(1, max(1, mut_rungs)) * rng.choice((1, -1))
                p[name] = shift(ladders, name, p[name], k, space)
        c = space._assemble(p)
        if c is not None:
            out.append(c)
    return out


def inject(space, rng: random.Random, n: int, base: Config, live: Sequence[str]) -> List[Config]:
    out = []
    for c in latin_hypercube(space, rng, n, names=list(live)):
        full = space._assemble({**base.params, **{k: c.params[k] for k in live}})
        if full is not None:
            out.append(full)
    return out


def local(ctx, rng: random.Random, centres: Sequence[Config], coarse: int, moves: int) -> List[Config]:
    """Every setting's line through each centre, and `moves` random 2-3 setting moves of 1-3 rungs around it."""
    names = [n for n in ctx.space.names if ctx.ladders.movable(n)]
    out: List[Config] = []
    for centre in centres:
        known = tried(ctx, centre)
        for name, rungs in plan(ctx, centre, known, coarse).items():
            out += [c for c in (single(ctx, centre, name, i) for i in rungs) if c is not None]
        for _ in range(moves):
            p = dict(centre.params)
            for name in rng.sample(names, min(len(names), rng.randint(2, 3))):
                p[name] = shift(ctx.ladders, name, p[name], rng.randint(1, 3) * rng.choice((1, -1)), ctx.space)
            c = ctx.space._assemble(p)
            if c is not None:
                out.append(c)
    return out


def pareto(points: Sequence[Tuple[Config, float, float, float]], window: float, cap: int) -> List[Scored]:
    """(core up, false positives down) front among points within `window` of the best score: (config, score)."""
    if not points:
        return []
    top = max(p[1] for p in points)
    front, best_fp = [], float("inf")
    for c, s, core, fp in sorted((p for p in points if p[1] >= top - window), key=lambda p: (-p[2], p[3])):
        if fp < best_fp:
            front.append((c, s))
            best_fp = fp
    return sorted(front, key=lambda t: -t[1])[:cap]


def trade_off_ends(points: Sequence[Tuple[Config, float, float, float]], window: float) -> List[Tuple[Config, str]]:
    """The best score, and the two ends of the core/false-positive trade-off among configs within `window` of it: the
    highest core and the fewest false positives. Each round's winners may sit at either end, so both are kept, not
    only the one this round's king happened to use."""
    if not points:
        return []
    top = max(points, key=lambda p: p[1])
    near = [p for p in points if p[1] >= top[1] - window]
    core = max(near, key=lambda p: (p[2], p[1]))
    fp = min(near, key=lambda p: (p[3], -p[1]))
    return [(top[0], "best score"), (core[0], "best core"), (fp[0], "best FP")]


def distinct(ctx, ranked: Sequence[Scored], n: int, min_callset: float, min_settings: float) -> List[Config]:
    """Up to n configs from `ranked` (best first), each far from those taken (callset distance, or settings)."""
    out: List[Config] = []
    for c, _ in ranked:
        if len(out) >= n:
            break
        if all(ctx.dist.far(c, o, min_callset, min_settings) for o in out):
            out.append(c)
    return out
