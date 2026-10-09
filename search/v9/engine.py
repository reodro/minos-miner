"""Search engine v9 (search/v9/DESIGN.md): A candidates -> B climb each to its peak -> C what B's grid did not look at.

    Step A  search/candidates.py  wave 0 (seeds, portfolio, hypercube), 3 GA generations over the both / core / fp
                                  groups, then one candidate per candidates.workers_per_candidate workers, far apart
                                  by callset distance (search/distance.py)
    Step B  search/climb.py       every candidate climbs at once: 5-setting grids at the settings' resolution
                                  (search/steps.py), cheap riders, step growth, 20 random jumps at a
                                  stall. Ends when every climber is at its peak, or at climb.end_min
    Step C  search/far.py         from every peak, one setting at a time: the rungs step growth jumped over and the
                                  rungs beyond the nearest tried value; a better point or a second high point (after
                                  a short climb) moves the peak and Step B climbs on
    after   the best config off every hill climbed so far (candidates.min_distance from each centre) is climbed with
            Step B's rules, again and again until the pick; the pick
            (far.pick_s before the end) names the highest local score in client.pick_key

The schedule is written for a 72-minute round and scales with the real one. Notes for main and the review:
client.notes["engine"] (version, phase, basins: one row per climber).
"""
from __future__ import annotations

import asyncio
import json
import logging
import math
import os
import random
import time
from typing import Dict, List, Optional, Sequence, Set, Tuple

from common.params import ParamRanges
from common.tuning import tune
from common.types import FULL, Config
from search import portfolio
from search.v9 import far
from search.bamgroups import bam_bounds
from search.base import Strategy
from search.v9.candidates import StepA
from search.client import EvalClient
from search.v9.climb import Climber
from search.v9.distance import Distance
from search.space import Space, bam_caps
from search.v9.steps import Ladders

log = logging.getLogger(__name__)

VERSION = "v9"
ROUND_S = 72 * 60                    # the schedule is written for this much search time; it scales with it
TICK_S = 2.0


async def settled_capacity(client: EvalClient, workers: int, wait_s: float = 60.0, settle_s: float = 20.0) -> int:
    """The workers the round has, once a fleet's nodes have joined (local and virtual executors answer at once)."""
    cap = client.capacity()
    nodes = getattr(client.executor, "nodes", None)
    if nodes is None:
        return cap or workers or 8
    t0 = stable = time.time()
    while time.time() - t0 < wait_s:
        await asyncio.sleep(5.0)
        c = client.capacity() or 0
        if c != cap:
            cap, stable = c, time.time()
            continue
        joining = [n for n in nodes() if not n["ready"] and not n["rejected"]]
        if cap and not joining and time.time() - stable >= settle_s:
            break
    log.info("engine: fleet capacity settled at %d workers after %.0f s", cap or 0, time.time() - t0)
    return cap or workers or 8


class Context:
    """What the steps share: the client, the round's space and ladders, distances, the clock, proposing and waiting."""

    def __init__(self, client: EvalClient, space: Space, rng: random.Random, profile: Optional[str]):
        self.client, self.space, self.rng, self.profile = client, space, rng, profile
        self.ladders = Ladders(space)
        self.dist = Distance(client, space)
        self.clock = client.clock
        self.t0 = self.clock()
        tl = client.time_left()
        self.length = tl if tl is not None and tl > 0 else ROUND_S
        self.scale = self.length / ROUND_S
        self.cap = 8
        self.phase = "starting"
        self.proposed: Set[str] = set()
        self.climbers: List[Climber] = []
        self.events: List[list] = []
        self.step_a: Optional[StepA] = None
        self.c_log: List[dict] = []                  # Step C: one row per pass on a peak
        self.marks: Dict[str, float] = {}            # seconds into the search each step ended
        self.picked: Optional[dict] = None
        self.started_ts = time.time()
        self.last_publish = -1e18

    # ---- time ------------------------------------------------------------------------------------------------------
    def T(self, seconds: float) -> float:
        return seconds * self.scale

    def elapsed(self) -> float:
        return self.clock() - self.t0

    def pick_at(self) -> float:
        return self.t0 + self.length - float(tune("far.pick_s"))

    def over(self) -> bool:
        t = self.client.time_left()
        return t is not None and t <= 0

    # ---- configs ---------------------------------------------------------------------------------------------------
    def fresh(self, c: Optional[Config]) -> bool:
        return c is not None and c.key not in self.proposed and (c.key, FULL) not in self.client.latest

    def score(self, c: Config) -> Optional[float]:
        r = self.client.latest.get((c.key, FULL))
        if r is None or r.score is None or r.gate_passed is False:
            return None
        return r.score

    def scored(self):
        """(config, score) of every scored config."""
        for (key, fid), r in self.client.latest.items():
            if fid == FULL and r.score is not None and r.gate_passed is not False:
                c = self.client.configs.get(key)
                if c is not None:
                    yield c, r.score

    def best(self) -> Optional[Tuple[Config, float]]:
        top = None
        for (key, fid), r in self.client.latest.items():
            if fid != FULL or r.score is None or r.gate_passed is False:
                continue
            if top is None or r.score > top[1]:
                c = self.client.configs.get(key)
                if c is not None:
                    top = (c, r.score)
        return top

    def propose(self, batch: Sequence[Config], tag: str, batches: Sequence[str]) -> List[asyncio.Future]:
        """One GATK run's worth: a head and its riders."""
        for c in batch:
            self.proposed.add(c.key)
        return self.client.propose_many(list(batch), tag=tag, batches=list(batches),
                                        priority=getattr(self, "priority", 0))

    async def wait_share(self, futs: Sequence[asyncio.Future], share: float, until: float) -> None:
        """Until `share` of `futs` are done, `until` passes, or the search's time is over."""
        if not futs:
            return
        need = math.ceil(share * len(futs))
        pending = {f for f in futs if not f.done()}
        while len(futs) - len(pending) < need and pending:
            if self.clock() >= until or self.over():
                break
            _, pending = await asyncio.wait(pending, timeout=TICK_S, return_when=asyncio.FIRST_COMPLETED)
            self.publish()

    def owner(self, c: Config, me: Climber) -> Optional[Climber]:
        """Another live climber already centred on `c`."""
        for o in self.climbers:
            if o is not me and o.state in ("climbing", "peak") and o.centre.key == c.key:
                return o
        return None

    # ---- telling -----------------------------------------------------------------------------------------------------
    def say(self, msg: str, *args) -> None:
        log.info("engine: " + msg, *args)
        try:
            text = msg % args if args else msg
        except (TypeError, ValueError):
            text = msg
        self.events.append([round(time.time(), 1), text])
        del self.events[:-200]

    def name_pick(self) -> None:
        """The highest local score goes in client.pick_key, from far.pick_s before the end on. Not earlier: main submits a
        named pick at once, past its resubmission interval and margin (main/loop.py maybe_submit)."""
        if self.clock() < self.pick_at():
            return
        top = self.best()
        if top is not None and top[0].key != self.client.pick_key:
            self.client.pick_key = top[0].key
            self.picked = {"key": top[0].key, "score": round(top[1], 4), "at_s": round(self.elapsed(), 1)}
            self.say("pick %s (%.4f) at t+%.1fm", top[0].key, top[1], self.elapsed() / 60)

    def publish(self, force: bool = False) -> None:
        """client.notes["engine"]: what main and the review read (the basins ledger)."""
        now = self.clock()
        if now >= self.pick_at():                    # a climb still waiting for its wave must not miss the pick
            self.name_pick()
        if not force and now - self.last_publish < TICK_S:
            return
        self.last_publish = now
        try:
            top = self.best()
            self.client.notes["engine"] = {
                "version": VERSION, "phase": self.phase, "elapsed_s": round(self.elapsed(), 1),
                "started_ts": round(self.started_ts, 1), "capacity": self.cap,
                "best": None if top is None else round(top[1], 4),
                "schedule": {"length_s": round(self.length, 1), "scale": round(self.scale, 4),
                             "a_max_s": round(self.T(60 * float(tune("candidates.max_min"))), 1),
                             "b_end_s": round(self.T(60 * float(tune("climb.end_min"))), 1),
                             "pick_s": round(self.length - float(tune("far.pick_s")), 1),
                             **{f"{k}_done_s": round(v, 1) for k, v in self.marks.items()}},
                "waves": self.step_a.wave_rows() if self.step_a is not None else [],
                "candidates": list(self.step_a.picked) if self.step_a is not None else [],
                "basins": [c.row() for c in self.climbers],
                "step_c": list(self.c_log),
                "pick": self.picked,
                "events": list(self.events[-200:]),
            }
        except Exception:  # noqa: BLE001 - the notes must never end the round
            log.exception("engine: publishing the notes failed")


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: int = 0) -> Strategy:
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        r = client.round
        bounds = bam_bounds(r.bam, ranges, tool) if r.bam else None
        caps = {} if bounds is not None else (bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {})
        profile = portfolio.bam_profile(r.bam, r.chrom, r.start, r.end) if r.bam else None
        space = Space(ranges, tool, caps=caps, bounds=bounds)
        if bounds is not None:
            log.info("engine: bam groups: %s", bounds.report)
        ctx = Context(client, space, rng, profile)
        ctx.cap = await settled_capacity(client, workers)
        ctx.say("%s: %d workers, %d settings searched, a %d-minute round (schedule x%.2f)", VERSION, ctx.cap,
                len(space.names), round(ctx.length / 60), ctx.scale)

        # ---- Step A ----------------------------------------------------------------------------------------------------
        # SEARCH_START is a one-replay override: climb that frozen config and skip the GA. The normal path is unchanged.
        ctx.phase = "A: candidates"
        ctx.step_a = StepA(ctx)
        forced = os.environ.get("SEARCH_START")
        known_score: Dict[str, float] = {}
        if forced:
            raw = json.loads(open(forced).read())
            seat = space._assemble(raw["params"])
            starts = [seat] if seat is not None else []
            if seat is not None and seat.key == raw.get("key"):
                known_score[seat.key] = float(raw["score"])
            ctx.say("step A skipped: climbing frozen seat %s opened at %.4f (%s)",
                    seat.key[:8] if seat is not None else "?", float(raw.get("score") or 0), raw.get("group") or "")
        else:
            try:
                starts = await ctx.step_a.run()
            except Exception:  # noqa: BLE001 - a bug in a step must not end the round
                log.exception("engine: step A failed")
                top = ctx.best()
                starts = [top[0]] if top else []
        ctx.marks["a"] = ctx.elapsed()
        ctx.publish(force=True)

        # ---- Step B ----------------------------------------------------------------------------------------------------
        ctx.phase = "B: climbing"
        end_b = min(ctx.t0 + ctx.T(60 * float(tune("climb.end_min"))), ctx.pick_at())
        for c in starts:
            s = ctx.score(c)
            if s is None:
                s = known_score.get(c.key)
            if s is not None:
                origin = "frozen seat" if c.key in known_score else "step A"
                ctx.climbers.append(Climber(ctx, len(ctx.climbers), c, s, origin=origin))
        n_main = len(ctx.climbers)
        in_b = [n_main]

        async def lifecycle(c: Climber) -> None:
            """Step B to the peak (or climb.end_min), then Step C on that peak at once: no climber waits for another."""
            try:
                await c.run(end_b)
            except Exception:  # noqa: BLE001 - one climber's bug ends that climber, not the round
                log.exception("engine: climber %d failed", c.id)
                c.state = "stopped"
            in_b[0] -= 1
            if in_b[0] == 0:
                ctx.marks["b"] = ctx.elapsed()
                ctx.say("step B done at t+%.1fm: %s", ctx.elapsed() / 60,
                        " | ".join(f"c{d.id} {d.opened:.4f} -> {d.score:.4f} ({d.state}, {d.iteration} it)"
                                   for d in ctx.climbers if d.parent is None))
            ctx.publish(force=True)
            if c.state in ("peak", "stopped") and not ctx.over() and all(
                    o.centre.key != c.centre.key or o.id >= c.id for o in ctx.climbers if o.parent is None):
                if ctx.phase.startswith("B"):
                    ctx.phase = "C: around the peaks" if in_b[0] == 0 else "B+C: climbing, and around the first peaks"
                try:
                    await far.improve(ctx, c, ctx.pick_at())
                except Exception:  # noqa: BLE001
                    log.exception("engine: step C failed on climber %d", c.id)
            if in_b[0] == 0:
                ctx.phase = "C: around the peaks"
            ctx.publish(force=True)

        await asyncio.gather(*(lifecycle(c) for c in list(ctx.climbers)))
        if not ctx.climbers:
            ctx.marks["b"] = ctx.elapsed()
        ctx.marks["c"] = ctx.elapsed()
        ctx.publish(force=True)

        # ---- the rest of the time: climb again from the best start not yet climbed from -------------------------------
        ctx.phase = "B: climbing again"
        while ctx.clock() < ctx.pick_at() and not ctx.over():
            ranked = sorted(ranked_configs(ctx), key=lambda t: -t[1])
            mc, ms = float(tune("candidates.min_distance")), float(tune("candidates.min_setting_distance"))
            c = s = None
            for x, sx in ranked:                         # the best config off every hill climbed so far
                if all(x.key not in (d.start.key, d.centre.key) and ctx.dist.far(x, d.centre, mc, ms)
                       for d in ctx.climbers):
                    c, s = x, sx
                    break
            if c is None:
                break
            d = Climber(ctx, len(ctx.climbers), c, s, origin="climbing again after step C")
            ctx.climbers.append(d)
            ctx.say("climbing again from %s %.4f", c.key[:8], s)
            await climb_all(ctx, [d], ctx.pick_at())

        # ---- the pick --------------------------------------------------------------------------------------------------
        ctx.phase = "picked"
        ctx.name_pick()
        ctx.publish(force=True)
        while not ctx.over() and client.in_flight:      # what is still running may still be submitted by main
            await asyncio.sleep(0 if hasattr(client.executor, "now") else TICK_S)
            ctx.name_pick()
            ctx.publish()

    return run


async def climb_all(ctx: Context, climbers: Sequence[Climber], until: float) -> None:
    async def one(c: Climber) -> None:
        try:
            await c.run(until)
        except Exception:  # noqa: BLE001 - one climber's bug ends that climber, not the round
            log.exception("engine: climber %d failed", c.id)
            c.state = "stopped"
        ctx.publish(force=True)
    await asyncio.gather(*(one(c) for c in climbers))


def ranked_configs(ctx: Context) -> List[Tuple[Config, float]]:
    out = []
    for (key, fid), r in ctx.client.latest.items():
        if fid == FULL and r.score is not None and r.gate_passed is not False and key in ctx.client.configs:
            out.append((ctx.client.configs[key], r.score))
    return out
