"""Certify that the round's pick is a 1-rung local optimum of the blended GATK score.

A stall used to fire random 2-4 setting jumps. Those do not prove a local peak, and they often
ran out the clock on the same shelf. After Step C, every untried one-rung neighbor of the current
best is scored. A neighbor that gains climb.accept_min becomes the centre and the climb resumes.
When every neighbor is scored and none gains, the pick is a 1-rung local optimum and the round can stop.
A neighbor that crashed, failed the gate, or is refused as a known tie never gets a score; it is settled, not
untried. Counting it as untried looped here until the pick (044200Z, 2026-10-16: 9,600 empty climbers while the
fleet sat idle for 38 minutes).
"""
from __future__ import annotations

import asyncio

from common.tuning import tune
from search.climb import TICK, Climber, held_by_climber, is_one_rung_local_opt, one_rung_neighbors


def open_neighbors(ctx, centre):
    """One-rung neighbors that can still change the certificate: (the ones that can be sent now, how many run)."""
    client = getattr(ctx, "client", None)
    running = client.running() if callable(getattr(client, "running", None)) else set()
    send, waiting = [], 0
    for c in one_rung_neighbors(ctx, centre):
        if ctx.score(c) is not None:
            continue
        if c.key in running:
            waiting += 1
        elif ctx.fresh(c) and held_by_climber(ctx, c) is None:
            send.append(c)
    return send, waiting


def higher_neighbor(ctx, centre, score):
    """The best scored one-rung neighbor that gains climb.accept_min, or None."""
    accept = float(tune("climb.accept_min"))
    best = None
    for c in one_rung_neighbors(ctx, centre):
        s = ctx.score(c)
        if s is not None and s - score >= accept and (best is None or s > best[1]):
            best = (c, s)
    return best


async def certify_round(ctx) -> None:
    ctx.certifying = True                              # island mode: the neighbours may leave the band by one rung
    try:
        await _certify_round(ctx)
    finally:
        ctx.certifying = False


async def _certify_round(ctx) -> None:
    """Until the pick, keep the round's best at a 1-rung local optimum."""
    pick = getattr(ctx, "pick_at", None)
    until = float(pick()) if callable(pick) else ctx.clock()
    if ctx.clock() >= until or ctx.over():
        return
    top = ctx.best() if callable(getattr(ctx, "best", None)) else None
    if top is None:
        ctx.say("certify: nothing scored")
        return
    ctx.phase = "C: certify"
    centre, score = top
    from search.candidates import in_indel_family
    from search.climb import family_in_reach
    if not in_indel_family(centre):
        fam = family_in_reach(ctx, score)
        if fam is not None:
            ctx.say("certify: this shelf is a local peak; remaining time climbs the edge family %s %.4f",
                    fam[0].key[:8], fam[1])
            centre, score = fam
    n = len(one_rung_neighbors(ctx, centre))
    ctx.say("certify: 1-rung neighborhood of %s %.4f (%d neighbors)", centre.key[:8], score, n)
    while ctx.clock() < until and not ctx.over():
        top = ctx.best()
        if top is not None and in_indel_family(top[0]):
            centre, score = top
        else:
            s = ctx.score(centre) if callable(getattr(ctx, "score", None)) else None
            if s is not None:
                score = s
        if is_one_rung_local_opt(ctx, centre, score):
            ctx.say("certify: 1-rung local optimum %.4f %s (%d neighbors)",
                    score, centre.key[:8], len(one_rung_neighbors(ctx, centre)))
            return
        up = higher_neighbor(ctx, centre, score)
        if up is not None:
            ctx.say("certify: scored 1-rung neighbor %s gains %+.4f; it becomes the centre", up[0].key[:8], up[1] - score)
            centre, score = up
            continue
        send, waiting = open_neighbors(ctx, centre)
        if not send and not waiting:
            neighbors = one_rung_neighbors(ctx, centre)
            ctx.say("certify: no 1-rung neighbor of %s %.4f gains; %d of its %d neighbors have no score "
                    "(crashed, gate-failed or a known tie)", centre.key[:8], score,
                    sum(1 for c in neighbors if ctx.score(c) is None), len(neighbors))
            return
        if not send:
            executor = getattr(getattr(ctx, "client", None), "executor", None)
            await asyncio.sleep(0 if hasattr(executor, "now") else TICK)
            continue
        climber = Climber(ctx, len(getattr(ctx, "climbers", []) or []), centre, score, origin="certify")
        ctx.climbers.append(climber)
        await climber.jumps(until)
        if climber.state == "climbing":
            await climber.run(until)
        centre, score = climber.centre, climber.score
        if is_one_rung_local_opt(ctx, centre, score):
            ctx.say("certify: 1-rung local optimum %.4f %s (%d neighbors)",
                    score, centre.key[:8], len(one_rung_neighbors(ctx, centre)))
            return
    missing = sum(1 for c in one_rung_neighbors(ctx, centre) if ctx.score(c) is None)
    ctx.say("certify: clock met the pick with %d one-rung neighbors still untried on %s %.4f",
            missing, centre.key[:8], score)
