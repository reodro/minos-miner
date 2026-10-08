"""Step C: what Step B's grid did not look at, one setting at a time, from each peak (search/DESIGN.md, Step C).

For every setting at a peak, the values already run with every other setting at the peak's values are looked up in
the results. Then, on each side of the peak:

    gaps     (case 1) every untried rung between the peak and the nearest tried value: what step growth jumped over
             (3 -> 4 -> 6 -> 9 lost: 7 and 8)
    beyond   (case 2) past the nearest tried value, to the end of the ladder: far.coarse rungs spread evenly (more when
             the fleet has room). The cheap settings take every rung (they share the peak's GATK run)

All of it goes out as one wave of single-setting changes. Then:

    a point beats the peak by climb.accept_min   it becomes the centre and Step B climbs from there; then another pass
    a second high point on a setting's profile   (higher than its profile neighbours, a lower point between it and
    within far.second_margin of the peak         the peak) every rung around it is tried, the best far.landings of
                                                 those get a short climb (far.landing_iterations grids); one that ends
                                                 above the peak becomes the centre, Step B climbs on, another pass
    neither                                      the peak is final
"""
from __future__ import annotations

import asyncio
from typing import Dict, List, Optional, Sequence, Tuple

from common.tuning import tune
from common.types import Config
from search.climb import Climber
from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS

CHEAP = set(GENOTYPING_PARAMS) | set(PAIRHMM_PARAMS)
DENSER = (1, 2, 4, 8, 16, 1000)          # far.coarse times this, the largest that still fits the wave's budget


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


def second_highs(ctx, peak: Config, score: float, known: Dict[str, Dict[int, float]]) -> List[Tuple[float, str, int, int, int]]:
    """(score, setting, rung, low, high) of each second high point within far.second_margin of the peak: higher than its
    profile neighbours, with a lower point between it and the peak. low..high: the untried rungs around it."""
    margin = float(tune("far.second_margin"))
    out = []
    for n, have in known.items():
        p = int(round(ctx.ladders.index(n, peak.params.get(n))))
        prof = sorted({**have, p: score}.items())
        idx = [i for i, _ in prof]
        for j, (i, s) in enumerate(prof):
            if i == p or s < score - margin:
                continue
            left = prof[j - 1][1] if j > 0 else float("-inf")
            right = prof[j + 1][1] if j + 1 < len(prof) else float("-inf")
            if s < left or s < right:
                continue
            between = [x for k, x in prof if min(i, p) < k < max(i, p)]
            if not between or min(between) >= s:
                continue
            low = idx[j - 1] + 1 if j > 0 else 0
            high = idx[j + 1] - 1 if j + 1 < len(prof) else len(ctx.ladders.rungs[n]) - 1
            out.append((s, n, i, low, high))
    return sorted(out, reverse=True)


async def send(ctx, configs: Sequence[Config], tag: str, until: float) -> None:
    futs = []
    heads = [c for c in configs if c is not None and ctx.fresh(c)]
    for c in heads:
        futs += ctx.propose([c], tag, ["small"])
    await ctx.wait_share(futs, float(tune("climb.quorum")), until)


async def one_pass(ctx, cl: Climber, budget: int, until: float) -> bool:
    """One Step C pass on climber `cl`'s peak. True when it moved the peak (and climbed on), False when final."""
    peak, score = cl.centre, cl.score
    known = tried(ctx, peak)
    base = int(tune("far.coarse"))
    chosen = plan(ctx, peak, known, base)
    for f in DENSER:                                   # the densest coarse pass that fits this peak's share of the wave
        p = plan(ctx, peak, known, base * f)
        if sum(len(v) for n, v in p.items() if n not in CHEAP) > budget:
            break
        chosen = p
    configs = [single(ctx, peak, n, i) for n, rungs in chosen.items() for i in rungs]
    n_gap = sum(len(v) for v in chosen.values())
    entry = {"peak": cl.id, "pass": len([e for e in ctx.c_log if e["peak"] == cl.id]) + 1,
             "started_s": round(ctx.elapsed(), 1), "from": round(score, 4), "points": len(configs),
             "gaps": sum(1 for n, rungs in chosen.items() for i in rungs if n not in CHEAP and gap(ctx, peak, known, n, i)),
             "outcome": "running", "detail": ""}
    ctx.c_log.append(entry)
    ctx.say("step C c%d: %d single-setting points around the peak %.4f", cl.id, n_gap, score)
    await send(ctx, configs, f"C{cl.id}", until)
    cl.refresh()
    score = cl.score
    known = tried(ctx, peak)
    best = max(((s, n, i) for n, have in known.items() for i, s in have.items()), default=None)
    accept = float(tune("climb.accept_min"))
    if best is not None and best[0] >= score + accept:
        s, n, i = best
        c = single(ctx, peak, n, i)
        ctx.say("step C c%d: %s -> %s gains %+.4f; Step B climbs from there", cl.id, n, ctx.ladders.rungs[n][i], s - score)
        entry.update(outcome="moved", detail=f"{n} -> {ctx.ladders.rungs[n][i]} {s - score:+.4f}")
        cl.resume(c, s, f"C {n}")
        await cl.run(until)
        return True

    highs = second_highs(ctx, peak, score, known)[:2 * int(tune("far.landings"))]
    if not highs or ctx.clock() >= until:
        ctx.say("step C c%d: nothing beyond the gaps; the peak %.4f is final", cl.id, score)
        entry.update(outcome="final", detail="no better point, no second high point")
        return False
    refine = []
    for _, n, i, low, high in highs:
        refine += [single(ctx, peak, n, r) for r in range(low, high + 1) if r != i]
    ctx.say("step C c%d: %d second high points (%s); %d points around them", cl.id, len(highs),
            ", ".join(f"{n}={ctx.ladders.rungs[n][i]} {s:.4f}" for s, n, i, _, _ in highs), len(refine))
    entry.update(outcome="second highs", detail=", ".join(f"{n}={ctx.ladders.rungs[n][i]} {s:.4f}" for s, n, i, _, _ in highs))
    await send(ctx, refine, f"C{cl.id}", until)
    known = tried(ctx, peak)
    landings: List[Tuple[float, Config]] = []
    for _, n, i, low, high in highs:
        area = [(s, r) for r, s in known[n].items() if low <= r <= high]
        if area:
            s, r = max(area)
            landings.append((s, single(ctx, peak, n, r)))
    landings = sorted({c.key: (s, c) for s, c in landings if c is not None}.values(), key=lambda t: -t[0])
    landings = landings[:int(tune("far.landings"))]
    if not landings or ctx.clock() >= until:
        entry.update(outcome="final")
        return False
    short = []
    for s, c in landings:
        d = Climber(ctx, len(ctx.climbers), c, s, origin=f"step C landing of c{cl.id}")
        d.parent = cl.id
        ctx.climbers.append(d)
        short.append(d)
    n_it = int(tune("far.landing_iterations"))
    await asyncio.gather(*(d.run(until, iterations=n_it, jumps=False) for d in short))
    for d in short:
        if d.state == "climbing":
            d.state = "landed"
    top = max(short, key=lambda d: d.score)
    ctx.say("step C c%d: short climbs %s; the peak is %.4f", cl.id,
            ", ".join(f"c{d.id} {d.opened:.4f} -> {d.score:.4f}" for d in short), cl.score)
    entry["landings"] = [d.id for d in short]
    if top.score >= cl.score + accept:
        entry.update(outcome="landing won", detail=f"c{top.id} {top.opened:.4f} -> {top.score:.4f}")
        cl.resume(top.centre, top.score, f"C landing c{top.id}")
        await cl.run(until)
        return True
    entry.update(outcome="final", detail="short climbs " + ", ".join(f"c{d.id} {d.opened:.4f} -> {d.score:.4f}" for d in short))
    return False


async def run(ctx, climbers: Sequence[Climber], until: float) -> None:
    """Step C on every peak at once, each pass repeated while it moves its peak."""
    peaks: List[Climber] = []
    for c in climbers:
        if c.state in ("peak", "stopped") and all(c.centre.key != o.centre.key for o in peaks):
            peaks.append(c)
    if not peaks:
        ctx.say("step C: no peak")
        return
    budget = max(1, ctx.cap // len(peaks))

    async def improve(cl: Climber) -> None:
        passes = 0
        while ctx.clock() < until and not ctx.over():
            passes += 1
            if not await one_pass(ctx, cl, budget, until):
                break
        if cl.state == "climbing":
            cl.state = "stopped"
        ctx.say("step C c%d: %d pass(es), final %.4f", cl.id, passes, cl.score)

    ctx.say("step C on %d peak(s): %s", len(peaks), ", ".join(f"c{c.id} {c.score:.4f}" for c in peaks))
    await asyncio.gather(*(improve(c) for c in peaks))
