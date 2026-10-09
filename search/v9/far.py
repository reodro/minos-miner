"""Step C: what Step B's grid did not look at, one setting at a time, from each peak (search/v9/DESIGN.md, Step C).
It starts on a climber as soon as Step B leaves that climber at its peak (or at climb.end_min), while the other
climbers may still be climbing.

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
from search.v9.climb import Climber
from search.v9.lines import CHEAP, densest, gap, single, tried

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
    chosen = densest(ctx, peak, known, budget)
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


def share(ctx) -> int:
    """A peak's share of the workers: the capacity over the climbers still at work (climbing or in Step C)."""
    active = [c for c in ctx.climbers if c.parent is None and c.state not in ("joined", "landed")]
    return max(1, ctx.cap // max(1, len(active)))


async def improve(ctx, cl: Climber, until: float) -> None:
    """Step C on one climber's peak, as soon as Step B has left it there; each pass repeated while it moves the peak."""
    passes = 0
    ctx.say("step C c%d: from the peak %.4f", cl.id, cl.score)
    while ctx.clock() < until and not ctx.over():
        passes += 1
        if not await one_pass(ctx, cl, share(ctx), until):
            break
    if cl.state == "climbing":
        cl.state = "stopped"
    ctx.say("step C c%d: %d pass(es), final %.4f", cl.id, passes, cl.score)


async def run(ctx, climbers: Sequence[Climber], until: float) -> None:
    """Step C on every distinct peak at once."""
    peaks: List[Climber] = []
    for c in climbers:
        if c.state in ("peak", "stopped") and all(c.centre.key != o.centre.key for o in peaks):
            peaks.append(c)
    if not peaks:
        ctx.say("step C: no peak")
        return
    await asyncio.gather(*(improve(ctx, c, until) for c in peaks))
