"""Step B: climb one candidate to the peak of its hill (search/DESIGN.md, Step B).

One iteration, from the centre:

    grid        v9.17: compass then a cube on the live axes. Each of climb.group_size settings is tried at -d and
                +d from the centre (~10 points). Axes whose face beats the centre are live. Zero or one live axis
                ends the iteration (move to the best point). Two or more: the 3^m cube on those axes only, skipping
                the centre and the faces already scored. d is the setting's resolution: one rung of its
                ladder (search/steps.py), whatever its range.
                    odd iterations   the read-pass settings no batch can share (FULL_RUN_PARAMS).
                                     Once those settings are already dead, that grid is skipped.
                    even iterations  5 of the 9 assembly and activity settings, in turn. A skipped read
                                     grid takes the next assembly group instead.
    riders      v9.2: every grid point carries the genotyping settings at -d and +d, inside its own GATK run
                (about 1% of a run each, so they fit). PairHMM, the cheap pairs, and every other cheap rung
                stay on the centre only, and only go out once every main point of the grid is already running.
                A cheap change must not start the centre's heavy work over on every grid point.
    lines       each grid setting alone over the rest of its ladder (search/lines.py: the gaps step growth jumped over,
                and rungs beyond the nearest tried value), as many as the climber's share of the fleet leaves room for.
                A setting whose score zig-zags (both neighbours lose, a value further out wins: min_mapping_quality_score
                on 20261008T175300Z) shows its best value in the same wave
    move        the best point of the iteration becomes the centre. Nothing is combined.
    growth      a grid setting whose move won is offered a larger step next time in that direction (and d back).
                A long ladder (30 rungs or more) grows 1d, 2d, 3d, 5d and then stays at 5d, so a wide axis
                reaches its hill in fewer grids. A short ladder grows 1d, 2d, 3d, 4d... When the grown step
                stops winning, the setting goes back to d. Example: A+1 wins the grid over (A, B, C); the next
                grid is around (A+1, B, C) with A at -1, 0, +2 and B, C at -1, 0, +1

An iteration that gains less than climb.accept_min has failed. That failure is remembered for those settings only:
a later gain clears the settings in that grid, and the others stay failed. When every expensive setting (read-pass
and all 7 assembly/activity ones, which rotate 5 at a time) has failed since its own last gain, the climb stalls:
every untried one-rung neighbor is scored. A neighbor that gains climb.accept_min becomes the centre
and the grid resumes; otherwise the centre is a 1-rung local optimum of the blended score.
Far points are Step C's (search/far.py). After the fine-tune, search/local_opt.py certifies the round's best.
"""
from __future__ import annotations

import asyncio
import itertools
import logging
import math
import random
import time
from typing import Dict, List, Optional, Sequence, Tuple

log = logging.getLogger(__name__)
TICK = 2.0
CLIMB_PRIORITY = 1          # a candidate update is taken before a GA child (priority 0) waiting in the queue

from common.tuning import tune
from search.duplicates import PARKED
from common.types import Config
from search.candidates import (EDGE_SCOUT_AXES, EDGE_SCOUT_CHOICES, FAMILY_SEAT_FLOOR, SEAT_GAP,
                               in_indel_family, left_the_rim, n_seats)
from search.engine_tag import ELITE_LEFTOVER, VERSION
from search.piles import CANDIDATES, GAP, MAX_PILES, RADIUS, already_peaked, best_member, centre_allowed, claim_peak, enabled as piles_enabled, follow, mating_pool, own_pile, refresh as refresh_piles, seal_peak, unseated
from search.lines import densest, plan, single, tried
from search.space import ACTIVITY_PARAMS, ASSEMBLY_PARAMS, FULL_RUN_PARAMS, GENOTYPING_PARAMS, PAIRHMM_PARAMS

# A ladder this long is walked 1, 2, 3, then 5 rungs at a time. Shorter ladders still grow by one rung per win:
# a 5-rung jump would leave them.
FILL_CUTOFF_S = 150.0        # Step C sends no fill this close to the pick (a GATK run takes about 100 s)
LEAP_SHARE_MAX = 0.25        # centre jumps never take more of a climber's workers: at 0.5 the waves were the slowest
LARGE_RUNGS = 30
LARGE_STEPS = (1, 2, 3, 5)
LINE_COARSE = 1               # one rung past the compass, once GA candidate updates have gone quiet


LATE_AXIS_MIN = 12          # the compass reaches one rung further for this many minutes before the pick


def late_axis(ctx) -> bool:
    """True in the final 12 minutes, when the climb owns the fleet and the grid can look one rung further."""
    pick = getattr(ctx, "pick_at", None)
    if not callable(pick):
        return False
    scale = getattr(ctx, "T", None)
    window = scale(LATE_AXIS_MIN * 60) if callable(scale) else LATE_AXIS_MIN * 60.0
    try:
        return ctx.clock() >= float(pick()) - float(window)
    except Exception:  # noqa: BLE001 - a board without a schedule keeps the one-rung compass
        return False


def _rung_gap(ctx, centre: Config, config: Config, name: str):
    ladders = getattr(ctx, "ladders", None)
    if ladders is None:
        return None
    try:
        a = ladders.index(name, centre.params.get(name))
        b = ladders.index(name, config.params.get(name))
    except Exception:  # noqa: BLE001
        return None
    if a != int(a) or b != int(b):
        return None
    return abs(int(b) - int(a))


def ga_hot(ctx) -> bool:
    """A GA child or a candidate has just moved a climber. The grid stays a compass until that stops."""
    moved = getattr(ctx, "ga_moved", None)
    if moved is None:
        return False
    minutes = float(tune("climb.ga_quiet_min"))
    scale = getattr(ctx, "T", None)
    quiet = scale(60.0 * minutes) if callable(scale) else 60.0 * minutes
    return ctx.clock() - float(moved) < quiet


def leftover_breeds() -> bool:
    """v9.18: the compass and its centre packs go out, and every core they leave keeps breeding."""
    return VERSION.startswith("v9.18")


def line_coarse(ctx) -> int:
    """No line while GA is still moving candidates. After that, one rung, so a wave can finish.

    v9.18 runs the compass on the same pass as the leftover GA, so the quiet window does not drop the line.
    """
    if ga_hot(ctx) and not leftover_breeds():
        return 0
    return LINE_COARSE


def one_rung_neighbors(ctx, centre: Config) -> List[Config]:
    """Every distinct config one ladder rung from `centre` on a single movable setting."""
    out, keys = [], {centre.key}
    space, ladders = ctx.space, ctx.ladders
    for name in space.names:
        if not ladders.movable(name):
            continue
        for step in (-1, 1):
            value = ladders.shift(name, centre.params.get(name), step)
            if value is None:
                continue
            config = space._assemble({**centre.params, name: value})
            if config is None or config.key in keys:
                continue
            keys.add(config.key)
            out.append(config)
    return out


def is_one_rung_local_opt(ctx, centre: Config, score: float) -> bool:
    """True when every one-rung neighbor is scored and none gains climb.accept_min."""
    accept = float(tune("climb.accept_min"))
    neighbors = one_rung_neighbors(ctx, centre)
    if not neighbors:
        return True
    for config in neighbors:
        known = ctx.score(config)
        if known is None:
            return False
        if known - score >= accept:
            return False
    return True


SMALL_LADDER, LARGE_LADDER = 40, 150   # rungs: the three step plans of climb.big_steps


def step_plan(n_rungs: int) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    """(widening steps, narrowing steps) of a setting with this many rungs (2026-10-18, the user's schedule)."""
    if n_rungs <= SMALL_LADDER:
        return (1, 2, 3, 5), (1,)
    if n_rungs <= LARGE_LADDER:
        return (3, 5, 7, 9), (3, 1)
    top = max(19, int(tune("climb.big_step_large_max")))
    return (3, 9, 18, top), (9, 3, 1)


class Climber:
    def __init__(self, ctx, cid: int, start: Config, score: float, origin: str = "step A"):
        # An island's climbers are numbered from 100 x its index, so every B{id} tag in the store is one climber.
        self.ctx, self.id, self.origin = ctx, cid + int(getattr(ctx, "climber_base", 0) or 0), origin
        self.start, self.centre, self.score = start, start, score
        self.found = score                                 # highest score this climber's own runs have returned
        self.opened = score
        self.first: Optional[Tuple[float, str]] = None        # (score, key) after the first iteration
        self.iteration = 0
        self.turn = 0                                          # odd turns read-pass, unless that grid is already dead
        self.growth: Dict[str, Tuple[int, int]] = {}           # setting -> (direction, rungs) offered next time
        self.reach: Dict[str, Tuple[str, int]] = {}            # setting -> (phase, index) of its grid step (step_plan)
        self.flat: set = set()                                 # settings gridded without a gain since their own last one
        self.recent: List[float] = []                          # gain of each finished grid, oldest first
        self.rotation = 0
        self.state = "climbing"                                # climbing | peak | joined | stopped | landed
        self.evaluated = 0
        self.waiting: Dict[str, Config] = {}                   # sent, not yet seen scored by a decision
        self.owed: Dict[str, Config] = {}                      # sent by this climber; a higher return still takes the centre
        self.waves: List[Tuple[str, int, int, float]] = []    # (label, heads, configs, minutes) per wave
        self.history: List[Tuple[float, float, str]] = []     # (seconds into the search, score, why)
        self.born_s = round(ctx.elapsed(), 1)
        self.ended_s: Optional[float] = None
        self.entry_gap: Optional[float] = None
        self.entry_peers: List[dict] = []
        self.led_round = False
        self.parent: Optional[int] = None                      # a Step C landing: the climber whose peak it came from
        self.kind = ""                                         # the grid in flight: read | assembly | jumps
        self.group: List[str] = []
        self.inflight: list = []                               # the heads of the wave in flight
        self.pending: List[Config] = []                        # this wave's heads not yet sent
        self.extra: List[Config] = []                          # line points added to fill free workers; they do not gate
        self.later: List[Config] = []                          # centre-only cheap work; sent after every head is out
        self.leaps: List[Config] = []                          # centre jumps (search/leaps.py); after the packs
        self.leap_futs: list = []                              # leap runs still out, across waves
        self.leap_sent: set = set()                            # keys of every centre jump this climber sent
        self.leap_meta: Dict[str, tuple] = {}                  # jump key -> ((setting, class), centre score when sent)
        self.cover_done: set = set()                           # centres whose one-per-class certificate went out
        self.pair_done: set = set()                            # centres whose two-setting jumps went out
        self.leap_share = float(tune("climb.leap_share"))      # this climber's share; TuRBO-style counters move it
        self.momentum: Optional[Tuple[str, int, int]] = None    # (setting, direction, rungs) of the last winning jump
        self._sent: Dict[object, Tuple[str, str]] = {}          # future -> (config key, head | rider | line | pack | jump)
        self.sent_wave: Dict[str, int] = {}                     # config key -> the wave that sent it
        self._wave_id = 0
        self._wave_t0, self._wave_n = 0.0, 0
        self.wave_results: list = []                           # a wave may close early (early_gain) before it filled
        self.phase = ""
        self.crash_edge: Dict[Tuple[str, int], Tuple[float, float]] = {}  # (setting, direction) -> (centre rung, crash rung)
        self.leap_succ = self.leap_fail = 0
        self._wave_leaps = self._wave_leap_wins = 0
        self.gate_futs: list = []
        self.extra_futs: list = []
        self.batch_live = False
        self.jumps_enabled = True
        self.wave_before = score
        self.wave_centre: Optional[Config] = None
        self.line_keys: set = set()
        self.climb_until = float("inf")
        self._wave_n = 0

        names = [n for n in ctx.space.names if ctx.ladders.movable(n)]
        self.read_names = [n for n in FULL_RUN_PARAMS if n in names]
        self.assembly_names = [n for n in ASSEMBLY_PARAMS + ACTIVITY_PARAMS if n in names]
        self.geno_names = [n for n in GENOTYPING_PARAMS if n in names]
        self.pair_names = [n for n in PAIRHMM_PARAMS if n in names]
        self.cheap_names = self.geno_names + self.pair_names
        self.all_names = names

    # ---- the climb ------------------------------------------------------------------------------------------------
    async def run(self, until: float, iterations: Optional[int] = None, jumps: bool = True) -> None:
        """Climb until the peak, `until`, or `iterations` more grids (a short climb: no random jumps)."""
        ctx = self.ctx
        if free_workers(ctx) is not None:
            await self.run_filled(until, iterations, jumps)
            return
        stop_at = None if iterations is None else self.iteration + iterations
        while ctx.clock() < until and self.state == "climbing":
            if stop_at is not None and self.iteration >= stop_at:
                break
            kind, group = self.next_group()
            if not group:
                if self.stalled():
                    if not jumps:
                        break
                    await self.jumps(until)
                continue
            self.iteration += 1
            before = self.score
            await self.grid(kind, group, until)
            if self.first is None:
                self.first = (self.score, self.centre.key)
            gain = self.score - before
            self.note_gain(group, gain)
            ctx.say("climber %d it %d (%s: %s): %+.4f -> %.4f; waves %s", self.id, self.iteration, kind,
                    ",".join(group), gain, self.score,
                    ", ".join(f"{w[1]}/{w[2]} in {w[3]:.1f}m" for w in self.waves[-2:]))
            if self.state == "climbing" and self.trailing():
                self.state = "stopped"
                ctx.say("climber %d is %.4f behind the best seat after two flat grids; the seat moves on",
                        self.id, self.gap_to_best())
                break
            if self.state == "climbing" and self.stalled():
                if not jumps:
                    break
                await self.jumps(until)
        if self.state == "climbing" and stop_at is None:
            self.state = "stopped"
        self.ended_s = round(ctx.elapsed(), 1)
        self.kind, self.group, self.inflight = "", [], []

    def stalled(self) -> bool:
        return set(self.read_names + self.assembly_names) <= self.flat

    def flat_waves(self) -> bool:
        """The last climb.stall_waves finished grids each gained less than accept_min. The certificate (begin_jumps)
        then scores every untried one-rung neighbour at once instead of one group per wave (the v4 fleet run: last
        gain at t+48.2, certificate only at t+60.1, after every group had gone flat one wave at a time)."""
        n = int(tune("climb.stall_waves"))
        if n <= 0 or len(self.recent) < n:
            return False
        accept = float(tune("climb.accept_min"))
        return all(g < accept for g in self.recent[-n:])

    def _peaked(self) -> None:
        """This fingerprint's peak is in the history, so a similar pile is not opened."""
        self.state = "peak"
        seal_peak(self)

    def note_gain(self, group: Sequence[str], gain: float) -> None:
        """A failed grid marks its own settings. A gain clears those settings and leaves every other failure in place."""
        if gain < float(tune("climb.accept_min")):
            self.flat.update(group)
        else:
            self.flat.difference_update(group)
        self.recent.append(gain)

    def gap_to_best(self) -> float:
        others = [cl.score for cl in self.ctx.climbers if cl is not self]
        return 0.0 if not others else max(others) - self.score

    def trailing(self) -> bool:
        """The last two grids gained nothing and this seat is SEAT_GAP or more behind the best other seat.

        A pile stays seated when it is behind. Falling behind is not a reason to leave its fingerprint.
        """
        if piles_enabled() and (self.origin or "").startswith("pile"):
            return False
        accept = float(tune("climb.accept_min"))
        if len(self.recent) < 2 or self.recent[-1] >= accept or self.recent[-2] >= accept:
            return False
        return self.gap_to_best() >= SEAT_GAP

    def next_group(self) -> Tuple[str, List[str]]:
        """The next grid. The read grid is skipped once every read setting in it is already dead."""
        size = int(tune("climb.group_size"))
        self.turn += 1
        read = self.read_names[:size]
        if self.turn % 2 == 1 and any(name not in self.flat for name in read):
            return "read", read
        return "assembly", self._assembly_group(size)

    def _assembly_group(self, size: int) -> List[str]:
        n = len(self.assembly_names)
        if n <= size:
            return list(self.assembly_names)
        group = [self.assembly_names[(self.rotation + i) % n] for i in range(size)]
        self.rotation = (self.rotation + size) % n
        return group

    def offsets(self, name: str) -> List[int]:
        d, k = self.growth.get(name, (1, 1))
        s = self.step_size(name)
        return [0, d * max(k, s), -d * s]

    def step_size(self, name: str) -> int:
        """The rungs this setting's grid step spans now (climb.big_steps; 1 when off, for a bool or a list)."""
        if not int(tune("climb.big_steps")):
            return 1
        spec = self.ctx.space.specs.get(name)
        if spec is not None and (spec.type == "bool" or spec.allowed):
            return 1
        widen, narrow = step_plan(len(self.ctx.ladders.rungs.get(name) or ()))
        phase, i = self.reach.get(name, ("w", 0))
        if phase == "w":
            return widen[min(i, len(widen) - 1)]
        if phase == "n":
            return narrow[min(i, len(narrow) - 1)]
        return 1

    def adapt_reach(self, group: Sequence[str], old: Optional[Config], new: Optional[Config]) -> None:
        """After a grid iteration: a setting that moved keeps its step; one that did not widens it, past the widest it
        narrows, and after one rung it stays at one rung."""
        if not int(tune("climb.big_steps")):
            return
        for n in group:
            if new is not None and old is not None and new.params.get(n) != old.params.get(n):
                continue
            widen, narrow = step_plan(len(self.ctx.ladders.rungs.get(n) or ()))
            phase, i = self.reach.get(n, ("w", 0))
            if phase == "w":
                self.reach[n] = ("w", i + 1) if i + 1 < len(widen) else ("n", 0)
            elif phase == "n":
                self.reach[n] = ("n", i + 1) if i + 1 < len(narrow) else ("f", 0)

    def compass_offsets(self, name: str) -> List[int]:
        """The compass rungs. In the last 12 minutes each side also takes the next rung out.

        The joint cube stays on the inner offsets, so the extra rung does not grow the cube.
        """
        offs = [off for off in self.offsets(name) if off != 0]
        if not (getattr(self, "compass_only", False) and late_axis(self.ctx)):
            return offs
        further = []
        for off in offs:
            step = 1 if off > 0 else -1
            nxt = off + step
            if nxt != 0 and nxt not in offs and nxt not in further:
                further.append(nxt)
        return offs + further

    def values(self, name: str, base) -> list:
        out = []
        for off in self.offsets(name):
            v = base if off == 0 else self.ctx.ladders.shift(name, base, off)
            if v is not None and v not in out:
                out.append(v)
        return out

    def with_(self, base: Config, changes: Dict[str, object]) -> Optional[Config]:
        return self.ctx.space._assemble({**base.params, **changes})

    def plan_leaps(self, centre: Config, n: Optional[int], taken=()) -> List[Config]:
        """Up to n centre jumps from `centre`, each from the island's jump arm with the highest Thompson sample
        (search/leaps.py). n None is the certificate: one jump per (setting, distance class, direction) the first
        time, two-setting jumps the second time, nothing after that.

        A jump already planned in this wave, already proposed, a tied callset, or on another climber's hill is
        left out."""
        from search import leaps
        if (n is not None and n <= 0) or float(tune("climb.leap_share")) <= 0:
            return []
        ctx = self.ctx
        taken = set(taken)

        def skip(c: Config) -> bool:
            return (c.key in taken or not ctx.fresh(c) or held_by_climber(ctx, c, me=self) is not None
                    or self.beyond_crash(centre, c))

        rng = getattr(ctx, "rng", None) or random.Random(self.id)
        names = self.read_names + self.assembly_names
        try:
            if n is not None:
                got = []
                if self.momentum is not None:                  # the winning move again, then twice as far
                    name, sign, k = self.momentum
                    self.momentum = None
                    got = leaps.momentum(ctx, centre, name, sign, k, skip, taken)
                    taken |= {g[3].key for g in got}
                got += leaps.pick(ctx, centre, names, skip, max(0, n - len(got)), rng, taken)
            elif centre.key not in self.cover_done:
                # One pass: the one-per-class jumps and the two-setting jumps go out together. Sent apart, the second
                # pass was a whole wave for two configs (044200Z islands, 2026-10-16: 14 minutes to name one peak).
                self.cover_done.add(centre.key)
                self.pair_done.add(centre.key)
                got = leaps.coverage(ctx, centre, names, skip, rng, taken)
                got += leaps.pairs(ctx, centre, names, skip, rng, taken=taken | {g[3].key for g in got})
                island = getattr(ctx, "island", None)
                elsewhere = island.elsewhere() if island is not None and hasattr(island, "elsewhere") else ()
                got += leaps.relink(ctx, centre, self.score, skip, rng, taken=taken | {g[3].key for g in got},
                                    elsewhere=elsewhere)
            else:
                got = []
        except Exception:  # noqa: BLE001 - a jump that cannot be planned must not stop the climb
            log.exception("climb: centre jumps for climber %d failed", self.id)
            return []
        for name, ci, k, config in got:
            self.leap_meta[config.key] = ((name, ci, 1 if k > 0 else -1), self.score, centre, config)
        return [config for _n, _c, _k, config in got]

    def fill_step_c(self, n: int) -> int:
        """Step C: workers the fine-tune leaves idle take centre jumps from its centre, then two- and three-setting moves
        of the leading arms (an elite-anchored GA child of the centre), until FILL_CUTOFF_S before the pick: a run sent
        later cannot score before it. 044200Z (2026-10-16): Step C kept 6-34% of Step B's cores busy. Returns how many
        configs were queued."""
        from search import leaps
        ctx = self.ctx
        pick = getattr(ctx, "pick_at", None)
        if n <= 0 or not callable(pick) or ctx.clock() > float(pick()) - FILL_CUTOFF_S:
            return 0
        taken = ({c.key for c in self.pending} | {c.key for c in self.later} | {c.key for c in self.leaps}
                 | set(self.waiting))
        got = self.plan_leaps(self.centre, n, taken)
        if len(got) < n:
            taken |= {c.key for c in got}

            def skip(c: Config) -> bool:
                return (c.key in taken or not ctx.fresh(c) or held_by_climber(ctx, c, me=self) is not None
                        or self.beyond_crash(self.centre, c))

            rng = getattr(ctx, "rng", None) or random.Random(self.id)
            try:
                more = leaps.mutations(ctx, self.centre, self.read_names + self.assembly_names, skip, n - len(got), rng)
            except Exception:  # noqa: BLE001 - a fill that cannot be planned must not stop the fine-tune
                log.exception("climb: step C fill for climber %d failed", self.id)
                more = []
            for name, ci, k, config in more:
                self.leap_meta[config.key] = ((name, ci, 1 if k > 0 else -1), self.score, self.centre, config)
            got += [config for _n, _c, _k, config in more]
        self.leaps.extend(got)
        if got:
            self._filled = getattr(self, "_filled", 0) + len(got)
        return len(got)

    def beyond_crash(self, centre: Config, c: Config, margin: int = 0) -> bool:
        """A jump at or past a value that crashed GATK in the same direction, from a centre within two rungs of the one
        it crashed from. 044200Z (2026-10-16): crashes doubled to 3.6-3.9% with jumps, mostly assembly_region_padding
        jumped down to 5-12 and max_assembly_region_size jumped up, because a crash never reached the arms."""
        if not self.crash_edge:
            return False
        lad = self.ctx.ladders
        for n in self.all_names:
            a, b = centre.params.get(n), c.params.get(n)
            if a == b:
                continue
            try:
                ia, ib = lad.index(n, a), lad.index(n, b)
            except Exception:  # noqa: BLE001
                continue
            sign = 1 if ib > ia else -1
            edge = self.crash_edge.get((n, sign))
            if edge is not None and abs(ia - edge[0]) <= 2 and (ib - edge[1]) * sign >= -margin:
                return True
        return False

    def drop_failed(self) -> int:
        """Runs this climber sent that came back without a score: a crash, or a gate failure. They leave the waiting
        list; a centre jump among them is a loss for its arm, and a crash marks that setting's edge (beyond_crash).
        The context counts them by kind (Context.failures). Returns how many crashed."""
        from search import leaps
        crashed = 0
        for fut, (key, kind) in list(self._sent.items()):
            if not fut.done():
                continue
            del self._sent[fut]
            try:
                out = fut.result()
            except Exception:  # noqa: BLE001 - a cancelled run is not a score
                out = None
            if out is None or not hasattr(out, "result"):
                continue
            r = out.result
            if r is not None and getattr(r, "score", None) is not None and getattr(r, "gate_passed", True) is not False:
                continue
            if r is None and str(getattr(out, "error", "") or "").startswith(PARKED):
                self.waiting.pop(key, None)             # a predicted repeat (search/duplicates.py): not a run, not a loss
                self.owed.pop(key, None)
                self.leap_meta.pop(key, None)
                continue
            crash = r is None
            self.waiting.pop(key, None)
            self.owed.pop(key, None)
            note = getattr(self.ctx, "note_failure", None)
            if callable(note) and crash:
                note(kind, getattr(out, "error", "") or "")
            meta = self.leap_meta.pop(key, None)
            if meta is None:
                continue
            arm, base, centre, config = meta
            if "+" not in arm[0]:
                leaps.arms_of(self.ctx).learn(arm, 0.0, float(tune("climb.accept_min")))
            if crash:
                crashed += 1
                for n in self.all_names:
                    a, b = centre.params.get(n), config.params.get(n)
                    if a == b:
                        continue
                    try:
                        ia, ib = self.ctx.ladders.index(n, a), self.ctx.ladders.index(n, b)
                    except Exception:  # noqa: BLE001
                        continue
                    sign = 1 if ib > ia else -1
                    old = self.crash_edge.get((n, sign))
                    if old is None or (old[1] - ib) * sign > 0 or abs(old[0] - ia) > 2:
                        self.crash_edge[(n, sign)] = (ia, ib)
        return crashed

    def note_momentum(self, old: Optional[Config], new: Config) -> None:
        """A centre jump just won: remember its setting, direction and size for the next wave (search/leaps.py)."""
        if old is None:
            return
        moved = [n for n in self.all_names if new.params.get(n) != old.params.get(n)]
        if len(moved) != 1:
            return
        k = self.ctx.ladders.rungs_moved(moved[0], old.params.get(moved[0]), new.params.get(moved[0]))
        if abs(k) >= 2 and k == int(k):
            self.momentum = (moved[0], 1 if k > 0 else -1, int(abs(k)))

    def learn_leaps(self, results) -> int:
        """Credit every returned centre jump to its arm (a win beats its centre by climb.accept_min); the wins."""
        from search import leaps
        if not self.leap_meta:
            return 0
        arms = leaps.arms_of(self.ctx)
        accept = float(tune("climb.accept_min"))
        won = 0
        for c, s in results:
            meta = self.leap_meta.pop(getattr(c, "key", None), None)
            if meta is None or s is None:
                continue
            arm, base = meta[0], meta[1]
            if "+" in arm[0]:                                  # a two-setting jump is not an arm
                ok = (s - base) >= accept
                won += ok
                two = getattr(self.ctx, "_leap2", None)
                if two is None:
                    two = [0, 0, 0.0]
                    try:
                        self.ctx._leap2 = two
                    except AttributeError:
                        pass
                two[0] += 1
                two[1] += int(ok)
                two[2] += max(0.0, s - base) if ok else 0.0
                continue
            won += arms.learn(arm, max(0.0, s - base), accept)
        self._wave_leap_wins += won
        return won

    def adapt_leap_share(self) -> None:
        """After a grid wave that sent jumps: three such waves with a jump win double the share (at most LEAP_SHARE_MAX),
        two without one halve it (at least 0.1). TuRBO's success and failure counters (Eriksson et al. 2019)."""
        if self._wave_leaps <= 0:
            return
        if self._wave_leap_wins > 0:
            self.leap_succ, self.leap_fail = self.leap_succ + 1, 0
        else:
            self.leap_succ, self.leap_fail = 0, self.leap_fail + 1
        if self.leap_succ >= 3:
            self.leap_share, self.leap_succ = min(LEAP_SHARE_MAX, 2.0 * self.leap_share), 0
        elif self.leap_fail >= 2:
            self.leap_share, self.leap_fail = max(0.1, self.leap_share / 2.0), 0
        self._wave_leaps = self._wave_leap_wins = 0
    def _known_score(self, config: Config):
        """The score this config already has, or None when it has not been scored."""
        scored = getattr(self.ctx, "scored", None)
        if not callable(scored):
            return None
        for other, score in scored():
            if other is not None and other.key == config.key and score is not None:
                return score
        return None

    def land(self, centre: Config, name: str, offset: int):
        """The value `offset` rungs from `centre` along `name`, stepped past a flat.

        A rung that assembles into the same GATK run, or that has already tied the centre's
        score, is one operation with the centre. The grid's slot keeps walking until it finds
        a rung that has not been scored. A rung that has already scored differently stops the
        walk: that direction is measured, and the line, not this slot, looks further out.
        An untried value of the setting is still returned. The ladder itself is unchanged.
        """
        if offset == 0:
            return None
        direction = 1 if offset > 0 else -1
        k = offset
        limit = len(self.ctx.ladders.rungs.get(name, ())) + 2
        tied = getattr(self.ctx, "tied", None)
        for _ in range(limit):
            value = self.ctx.ladders.shift(name, centre.params.get(name), k)
            if value is None:
                return None
            config = self.with_(centre, {name: value})
            same_run = (config is None or config.key == centre.key
                        or config.params.get(name) == centre.params.get(name))
            if same_run:
                k += direction
                continue
            known = self._known_score(config)
            if known is not None and known == self.score:
                k += direction
                continue
            if tied is not None and tied.holds(config.params):
                k += direction
                continue
            if known is not None:
                return None
            return value
        return None

    def compass_heads(self, centre: Config, group: Sequence[str]) -> List[Config]:
        """Each grid setting at −d and +d, the others held. A flat step lands on the next untried rung."""
        heads, keys = [], {centre.key}
        for n in group:
            for off in self.compass_offsets(n):
                if off == 0:
                    continue
                v = self.land(centre, n, off)
                if v is None:
                    continue
                c = self.with_(centre, {n: v})
                if c is not None and c.key not in keys:
                    keys.add(c.key)
                    heads.append(c)
        return heads

    def live_axes(self, centre: Config, group: Sequence[str],
                  results: Sequence[Tuple[Config, float]], base: float) -> List[str]:
        """Grid settings whose compass face beats the centre. The full system ignores a face inside 0.05."""
        from search.separate import FACE_MARGIN
        margin = FACE_MARGIN if getattr(self, "compass_only", False) or (self.origin or "") == "step C" else 0.0
        live = []
        for n in group:
            best_n = None
            for item in results:
                if not isinstance(item, tuple) or len(item) != 2:
                    continue
                c, s = item
                if c is None or not hasattr(c, "params") or s is None:
                    continue
                if c.key in self.line_keys:
                    continue
                moved = [x for x in group if c.params.get(x) != centre.params.get(x)]
                if moved != [n]:
                    continue
                inner = min((abs(off) for off in self.offsets(n) if off), default=1)
                apart = _rung_gap(self.ctx, centre, c, n)
                if apart is not None and apart > inner:
                    continue
                if best_n is None or s > best_n:
                    best_n = s
            if best_n is not None and best_n > base + margin:
                live.append(n)
        return live

    def confirm_heads(self, centre: Config, live: Sequence[str], seen: set) -> List[Config]:
        """The 3^m cube on the live axes, minus the centre and the one-setting faces already scored."""
        if len(live) < 2:
            return []
        axes = [self.values(n, centre.params.get(n)) for n in live]
        heads, keys = [], set(seen) | {centre.key}
        for combo in itertools.product(*axes):
            changes = dict(zip(live, combo))
            if sum(centre.params.get(n) != changes[n] for n in live) < 2:
                continue
            c = self.with_(centre, changes)
            if c is not None and c.key not in keys:
                keys.add(c.key)
                heads.append(c)
        return heads


    def geno_riders(self, head: Config) -> List[Config]:
        """Genotyping at one rung either side. A rung that has already tied steps out to the next untried value.

        About 1% of a run each, so they share the head's GATK run. The head is the config the rung moves from;
        a tie is still judged against this climber's centre score, which is the plateau the grid is on.
        """
        out = []
        for n in self.geno_names:
            for off in (1, -1):
                v = self.land(head, n, off)
                if v is None:
                    continue
                c = self.with_(head, {n: v})
                if c is not None and c.key != head.key:
                    out.append(c)
        return out

    def extra_cost(self, main: Config, variant: Config) -> float:
        """Added share of a solo run when `variant` rides in `main`'s batch. PairHMM is 0.12, genotyping 0.011."""
        changed = [n for n in self.cheap_names if main.params.get(n) != variant.params.get(n)]
        if any(n in self.pair_names for n in changed):
            return 0.12
        return 0.011 if changed else 0.0

    def packs(self, configs: Sequence[Config]) -> List[List[Config]]:
        """Groups that fit in one big batch, so one heavy run is not started per cheap config."""
        limit = float(tune("evaluator.big_batch"))
        out: List[List[Config]] = []
        pack: List[Config] = []
        cost = 0.0
        for c in configs:
            add = 1.0 if not pack else self.extra_cost(pack[0], c)
            if pack and cost + add > limit:
                out.append(pack)
                pack, cost = [], 0.0
                add = 1.0
            pack.append(c)
            cost += add
        if pack:
            out.append(pack)
        return out

    def centre_later(self, centre: Config, known: Dict[str, Dict[int, float]]) -> List[Config]:
        """PairHMM, cheap pairs, and the rest of each cheap ladder. Not the genotyping steps the head already carries."""
        ctx = self.ctx
        seen = {centre.key, *(c.key for c in self.geno_riders(centre))}
        out: List[Config] = []

        def add(c: Optional[Config]) -> None:
            if c is None or c.key in seen or not ctx.fresh(c) or held_by_climber(ctx, c, me=self) is not None:
                return
            seen.add(c.key)
            out.append(c)

        lad = ctx.ladders
        singles: List[Tuple[str, object]] = []
        for n in self.cheap_names:
            for off in (1, -1):
                v = lad.shift(n, centre.params.get(n), off)
                if v is None:
                    continue
                singles.append((n, v))
                if n in self.pair_names:
                    add(self.with_(centre, {n: v}))
        for (n1, v1), (n2, v2) in itertools.combinations(singles, 2):
            if n1 != n2:
                add(self.with_(centre, {n1: v1, n2: v2}))
        cheap_plan = plan(ctx, centre, {n: known[n] for n in self.cheap_names if n in known}, 0)
        for n, rungs in cheap_plan.items():
            for i in rungs:
                add(single(ctx, centre, n, i))
        return out

    async def send(self, heads: Sequence[Config], riders_of, label: str, until: float) -> List[Tuple[Config, float]]:
        """Each head in one GATK run with its riders. Waits for climb.quorum of the heads, and of the riders too when
        climb.wait_riders is 1. Returns (config, score) of this wave's scored configs plus any of this climber's
        earlier configs that came back since the last decision (a late rider)."""
        ctx = self.ctx
        head_futs, all_futs, seen = [], [], set()
        t0 = ctx.clock()
        for h in heads:
            if h.key in seen or held_by_climber(ctx, h, me=self) is not None:
                continue
            seen.add(h.key)
            batch = [h] + [r for r in riders_of(h) if r.key not in seen]
            seen.update(r.key for r in batch)
            fs = ctx.propose(batch, f"B{self.id}", ["small"] + ["big"] * (len(batch) - 1))
            head_futs.append(fs[0])
            all_futs += fs
            for c in batch:
                self._watch(c)
        self.evaluated += len(all_futs)
        self.inflight = head_futs
        wait_all = int(tune("climb.wait_riders")) == 1
        await ctx.wait_share(all_futs if wait_all else head_futs, float(tune("climb.quorum")), until)
        out = []
        for k, c in list(self.waiting.items()):
            s = ctx.score(c)
            if s is not None:
                out.append((c, s))
                del self.waiting[k]
        self.waves.append((label, len(head_futs), len(all_futs), round((ctx.clock() - t0) / 60, 2)))
        return out

    async def grid(self, kind: str, group: List[str], until: float) -> None:
        ctx = self.ctx
        self.kind, self.group = kind, list(group)
        centre = self.centre
        axes = [self.values(n, centre.params.get(n)) for n in group]
        heads, keys = [], {centre.key}
        self.wave_centre = centre
        heads = self.compass_heads(centre, group)
        keys = {centre.key, *(c.key for c in heads)}
        # The lines: every grid setting alone over the rest of its ladder (the gaps and beyond the nearest tried value),
        # with the fleet's room left after the compass; every rung of the cheap settings rides on the centre's run. A
        # setting whose score zig-zags (both neighbours lose, a value further out wins) shows its best value here.
        known = tried(ctx, centre)
        room = max(0, self.budget() - max(len(heads), 1))
        line_plan = densest(ctx, centre, {n: known[n] for n in group if n in known}, room)
        line = [c for c in (single(ctx, centre, n, i) for n, rungs in line_plan.items() for i in rungs)
                if c is not None and c.key not in keys]
        keys.update(c.key for c in line)
        self.line_keys = {c.key for c in line}
        def riders_of(h: Config) -> List[Config]:
            return [] if h.key in self.line_keys else self.geno_riders(h)

        results = await self.send(heads + line, riders_of, kind, until)
        live = self.live_axes(centre, group, results, self.score)
        cube = self.confirm_heads(centre, live, keys)
        if cube:
            keys.update(c.key for c in cube)
            results.extend(await self.send(cube, riders_of, kind, until))
        for pack in self.packs(self.centre_later(centre, known)):
            more = await self.send(pack[:1], lambda h, pack=pack: pack[1:] if h.key == pack[0].key else [], kind, until)
            results.extend(more)
        self.refresh()
        best = max(results, key=lambda t: t[1], default=None)
        old = self.centre
        if best is not None and best[1] > self.score and centre_allowed(self, best[0]):
            on_line = best[0].key in self.line_keys
            if on_line:
                for n in group:
                    self.growth.pop(n, None)
            else:
                self.grow(group, old, best[0])
            self.move(best[0], best[1], f"{kind} {'line' if on_line else 'grid'}")
        else:
            for n in group:
                self.growth.pop(n, None)

    def budget(self) -> int:
        """This climber's share of the workers: the capacity over the climbers still climbing."""
        climbing = sum(1 for c in getattr(self.ctx, "climbers", []) if c.state == "climbing") or 1
        return max(1, int(getattr(self.ctx, "cap", 0)) // climbing)

    def queued(self) -> bool:
        return bool(self.pending or self.later or self.extra or self.leaps)

    def submit_one(self) -> None:
        """One unsent head, then one more line point, then one packed centre batch, then one centre jump."""
        if self.pending:
            self._submit_head()
        elif self.extra:
            self._submit_extra()
        elif self.later:
            self._submit_pack()
        elif self.leaps:
            self._submit_leap()

    def _submit_leap(self) -> None:
        """One centre jump, for a worker the grid left free. It does not hold the wave open; its score is read when
        the wave closes, or at the next one if it is still running."""
        ctx = self.ctx
        while self.leaps and (not ctx.fresh(self.leaps[0]) or held_by_climber(ctx, self.leaps[0], me=self) is not None):
            self.leaps.pop(0)
        if not self.leaps:
            return
        h = self.leaps.pop(0)
        fs = ctx.propose([h], f"B{self.id}", ["small"], priority=CLIMB_PRIORITY)
        self._note_sent(fs, [h], "jump")
        self.leap_futs = [f for f in self.leap_futs if not f.done()] + list(fs)
        self.leap_sent.add(h.key)
        self._wave_leaps += 1
        self._watch(h)
        self.evaluated += len(fs)
        self._wave_n += len(fs)

    def _submit_extra(self) -> None:
        """One more line point, for a worker the planned grid did not fill. It does not hold the wave open."""
        ctx = self.ctx
        while self.extra and (not ctx.fresh(self.extra[0]) or held_by_climber(ctx, self.extra[0], me=self) is not None):
            self.extra.pop(0)
        if not self.extra:
            return
        h = self.extra.pop(0)
        fs = ctx.propose([h], f"B{self.id}", ["small"], priority=CLIMB_PRIORITY)
        self._note_sent(fs, [h], "line")
        self.extra_futs.extend(fs)
        self._watch(h)
        self.evaluated += len(fs)
        self._wave_n += len(fs)

    def _submit_head(self) -> None:
        """The head and its genotyping riders: one GATK run. PairHMM is not on this run."""
        ctx = self.ctx
        probes = getattr(ctx, "_probe_keys", ())
        while self.pending and (not ctx.fresh(self.pending[0]) or (
                self.pending[0].key not in probes and held_by_climber(ctx, self.pending[0], me=self) is not None)):
            self.pending.pop(0)
        if not self.pending:
            return
        h = self.pending.pop(0)
        guard = getattr(ctx, "guard", None)
        risky = self.beyond_crash(self.wave_centre or self.centre, h, margin=1) or (
            guard is not None and guard.risky(h))        # a crash would take its riders down too
        riders = [] if h.key in self.line_keys or h.key in probes or risky else self.geno_riders(h)
        batch = [h] + [r for r in riders if r is not None and r.key != h.key and ctx.fresh(r)]
        fs = ctx.propose(batch, f"B{self.id}", ["small"] + ["big"] * (len(batch) - 1), priority=CLIMB_PRIORITY)
        self._note_sent(fs[:1], batch[:1], "jump" if h.key in self.leap_meta else ("line" if h.key in self.line_keys else "head"))
        self._note_sent(fs[1:], batch[1:], "rider")
        self.gate_futs.append(fs[0])
        self.inflight = self.gate_futs
        for c in batch:
            self._watch(c)
        self.evaluated += len(fs)
        self._wave_n += len(fs)

    def _submit_pack(self) -> None:
        """The next centre pack that fits in one big batch. Same heavy work, cheap changes inside it."""
        ctx = self.ctx
        limit = float(tune("evaluator.big_batch"))
        pack: List[Config] = []
        cost = 0.0
        while self.later:
            c = self.later[0]
            if not ctx.fresh(c) or held_by_climber(ctx, c, me=self) is not None:
                self.later.pop(0)
                continue
            add = 1.0 if not pack else self.extra_cost(pack[0], c)
            if pack and cost + add > limit:
                break
            self.later.pop(0)
            pack.append(c)
            cost += add
        if not pack:
            return
        fs = ctx.propose(pack, f"B{self.id}", ["big"] * len(pack), priority=CLIMB_PRIORITY)
        self._note_sent(fs, pack, "pack")
        self.extra_futs.extend(fs)
        for c in pack:
            self._watch(c)
        self.evaluated += len(fs)
        self._wave_n += len(fs)

    def early_gain(self) -> bool:
        """Move at once: a run of the live wave already beats the centre by climb.accept_min, so the wave closes now
        and the climber moves there (close_wave). The rest of the wave keeps running; a run that later beats the new
        centre is taken the same way (claim_higher). In Step C a higher round best does the same.

        Four latest rounds (2026-10-16): a better result sat 35-41 minutes per run in total before its climber moved
        (13 on the slowest round): about half waiting for the rest of its own wave, a third for a cube around the old
        centre. No side effects: close_wave makes the move."""
        if not self.batch_live or not int(tune("climb.move_at_once")):
            return False
        accept = float(tune("climb.accept_min"))
        wave = getattr(self, "_wave_id", 0)
        for k, c in self.waiting.items():
            if c is None or not hasattr(c, "params") or self.sent_wave.get(k, wave) != wave:
                continue
            s = self.ctx.score(c)
            if s is not None and s - self.score >= accept and family_allows(self, c, s) and centre_allowed(self, c):
                return True
        if (self.origin or "") == "step C":
            best = getattr(self.ctx, "best", None)
            top = best() if callable(best) else None
            if top is not None and top[1] - self.score >= accept and family_allows(self, top[0], top[1]):
                return True
        return False

    def gate_done(self) -> bool:
        """The planned heads have all been sent and climb.quorum of them are back."""
        if not self.batch_live or self.pending:
            return False
        if not self.gate_futs:
            return True
        share = 1.0 if self.kind == "jumps" else float(tune("climb.quorum"))   # a peak is named on every neighbour
        need = math.ceil(share * len(self.gate_futs))
        return sum(1 for f in self.gate_futs if f.done()) >= need

    def open_futs(self) -> list:
        return [f for f in self.gate_futs + self.extra_futs + self.leap_futs if not f.done()]

    def begin_grid(self, until: float) -> bool:
        """The next compass (and its lines), held until a worker is free. False when there is nothing left to climb."""
        ctx = self.ctx
        if ctx.clock() >= until or ctx.over() or self.state != "climbing":
            return False
        if ((self.origin or "") != "step C" and not getattr(self, "compass_only", False)
                and ga_hot(ctx) and not leftover_breeds()):
            return False
        if getattr(self, "compass_only", False):
            from search.separate import compass_axes
            group = compass_axes(ctx, self)
            kind = "read" if group and group[0] in self.read_names else "assembly"
        else:
            kind, group = self.next_group()
        if not group:
            if not getattr(self, "compass_only", False):
                self.state = "stopped"
            return False
        self.iteration += 1
        self.wave_before = self.score
        self.kind, self.group = kind, list(group)
        seen_axes = set(getattr(self, "tried_axes", ()) or ())
        seen_axes.update(group)
        self.tried_axes = seen_axes
        self.phase = "compass"
        self.wave_results = []
        centre = self.centre
        self.wave_centre = centre
        heads = self.compass_heads(centre, group)
        if getattr(self, "compass_only", False) and late_axis(ctx) and heads:
            ctx.say("climber %d: compass reaches one rung further (%d points)", self.id, len(heads))
        if getattr(self, "compass_only", False) and not heads:
            self.flat.update(group)
            group = compass_axes(ctx, self)
            if not group:
                return False
            kind = "read" if group[0] in self.read_names else "assembly"
            self.kind, self.group = kind, list(group)
            seen_axes = set(self.tried_axes)
            seen_axes.update(group)
            self.tried_axes = seen_axes
            heads = self.compass_heads(centre, group)
            if not heads:
                self.flat.update(group)
                return False
        keys = {centre.key, *(c.key for c in heads)}
        known = tried(ctx, centre)
        # A hot GA sends nothing from this climber: begin_grid already returned. Once the updates
        # go quiet, one rung past the compass walks to the peak before the clock.
        coarse = 0 if (self.origin or "") == "step C" or getattr(self, "compass_only", False) else line_coarse(ctx)
        if coarse <= 0:
            line = []
            self.later = []
        else:
            line_plan = plan(ctx, centre, {n: known[n] for n in group if n in known}, coarse)
            line = [c for c in (single(ctx, centre, n, i) for n, rungs in line_plan.items() for i in rungs)
                    if c is not None and c.key not in keys]
            self.later = self.centre_later(centre, known)
        self.line_keys = {c.key for c in line}
        self.pending, self.extra = heads + line, []
        # Centre jumps take climb.leap_share of this climber's workers, after its heads and packs (search/leaps.py).
        if float(tune("climb.leap_share")) > 0 and not getattr(self, "compass_only", False):
            taken = {c.key for c in self.pending} | {c.key for c in self.later}
            self.leaps = self.plan_leaps(centre, int(self.leap_share * self.budget()), taken)
        else:
            self.leaps = []
        self.gate_futs, self.extra_futs = [], []
        self.batch_live = True
        self._wave_t0 = ctx.clock()
        self._wave_n = 0
        self._wave_id += 1
        return True

    def begin_confirm(self, live: Sequence[str]) -> bool:
        """Queue the joint cube on the live axes. False when every combo is already scored."""
        centre = self.wave_centre or self.centre
        seen = {centre.key, *(c.key for c, _ in self.wave_results)}
        cube = self.confirm_heads(centre, live, seen)
        if not cube:
            return False
        self.phase = "confirm"
        self.pending = cube
        self.gate_futs = []
        self.inflight = []
        self._wave_t0 = self.ctx.clock()
        self._wave_n = 0
        self.batch_live = True
        self.ctx.say("climber %d it %d (%s): cube on %d axes", self.id, self.iteration, self.kind, len(live))
        return True
    def begin_jumps(self, until: float) -> bool:
        """Score every untried one-rung neighbor and every untried centre jump. That wave is the peak certificate.

        Centre jumps still running from earlier waves are waited for before the peak is named."""
        ctx = self.ctx
        if not self.jumps_enabled or ctx.clock() >= until or ctx.over():
            return False
        heads = []
        for config in one_rung_neighbors(ctx, self.centre):
            if ctx.score(config) is not None:
                continue
            if not ctx.fresh(config) or held_by_climber(ctx, config, me=self) is not None:
                continue
            heads.append(config)
        n_rung = len(heads)
        self.drop_failed()
        heads += self.plan_leaps(self.centre, None, {c.key for c in heads})
        running = [f for f in self.leap_futs if not f.done()]
        if not heads and not running:
            n = len(one_rung_neighbors(ctx, self.centre))
            self._peaked()
            ctx.say("climber %d: 1-rung local optimum %.4f %s (%d neighbors)",
                    self.id, self.score, self.centre.key[:8], n)
            return False
        self.kind, self.group = "jumps", []
        self.wave_before = self.score
        self.wave_centre = self.centre
        self.phase = "jumps"
        self.line_keys = set()
        self.pending, self.extra, self.leaps = heads, [], []
        self.gate_futs, self.extra_futs = ([] if heads else running), []
        self.batch_live = True
        self._wave_t0 = ctx.clock()
        self._wave_n = 0
        self._wave_id += 1
        if heads:
            ctx.say("climber %d: certifying 1-rung neighborhood, %d untried, %d centre jumps",
                    self.id, n_rung, len(heads) - n_rung)
        else:
            ctx.say("climber %d: waiting for %d centre jumps before naming a peak", self.id, len(running))
        return True

    def extend(self, n: int) -> int:
        """Up to n more single-setting values of this wave, for workers the planned grid did not fill."""
        if n <= 0 or not self.batch_live or self.kind == "jumps" or self.wave_centre is None or self.pending:
            return 0
        ctx, centre = self.ctx, self.wave_centre
        known = tried(ctx, centre)
        have = len(self.line_keys) + len(self.extra) + len(self.extra_futs)
        line_plan = densest(ctx, centre, {n: known[n] for n in self.group if n in known}, have + n)
        seen = set(self.line_keys) | {c.key for c in self.extra} | set(self.waiting)
        added = []
        for name, rungs in line_plan.items():
            for i in rungs:
                c = single(ctx, centre, name, i)
                if c is not None and c.key not in seen:
                    seen.add(c.key)
                    added.append(c)
                if len(added) >= n:
                    break
            if len(added) >= n:
                break
        self.line_keys.update(c.key for c in added)
        self.extra.extend(added)
        return len(added)

    def close_wave(self) -> None:
        """The planned heads are back: after a compass, maybe a cube; then move, then jumps if this grid failed."""
        ctx = self.ctx
        self.drop_failed()
        results = []
        for k, c in list(self.waiting.items()):
            if c is None or not hasattr(c, "params"):
                del self.waiting[k]
                continue
            s = ctx.score(c)
            if s is not None:
                results.append((c, s))
                del self.waiting[k]
        self.learn_leaps(results)
        self.waves.append((self.kind, len(self.gate_futs), self._wave_n, round((ctx.clock() - self._wave_t0) / 60, 2)))
        self.refresh()
        if self.kind == "jumps":
            self.batch_live = False
            self.pending, self.extra, self.leaps = [], [], []
            best = max(results, key=lambda t: t[1], default=None)
            if (best is not None and best[1] - self.score >= float(tune("climb.accept_min"))
                    and centre_allowed(self, best[0])):
                moved = [n for n in self.all_names if best[0].params.get(n) != self.centre.params.get(n)]
                one_rung = len(moved) == 1 and _rung_gap(ctx, self.centre, best[0], moved[0]) in (None, 1)
                ctx.say("climber %d: a %s gains %+.4f; the grid resumes", self.id,
                        "1-rung neighbor" if one_rung else "centre jump (" + ", ".join(moved) + ")",
                        best[1] - self.score)
                if not one_rung:
                    self.note_momentum(self.centre, best[0])
                self.growth.clear()
                self.reach.clear()
                self.move(best[0], best[1], "jump")
                self.flat.clear()
            elif self.jumps_enabled and float(tune("climb.leap_share")) > 0:
                self.begin_jumps(self.climb_until)            # the rest of the certificate, or the peak
            elif is_one_rung_local_opt(ctx, self.centre, self.score):
                n = len(one_rung_neighbors(ctx, self.centre))
                self._peaked()
                ctx.say("climber %d: 1-rung local optimum %.4f %s (%d neighbors)",
                        self.id, self.score, self.centre.key[:8], n)
            elif self.jumps_enabled:
                self.begin_jumps(self.climb_until)
            else:
                self.state = "stopped"
            return
        self.wave_results.extend(results)
        if self.phase == "compass":
            centre = self.wave_centre or self.centre
            live = self.live_axes(centre, self.group, self.wave_results, self.wave_before)
            accept = float(tune("climb.accept_min"))
            cleared = int(tune("climb.move_at_once")) and any(
                s - self.wave_before >= accept and family_allows(self, c, s) for c, s in self.wave_results)
            # A compass point that clears the bar is moved to now; the next compass, from there, tries the other live
            # settings. The cube on the old centre is for faces that only beat it by less than the bar.
            if len(live) >= 2 and not cleared and self.begin_confirm(live):
                return
        self._finish_grid_wave()

    def _finish_grid_wave(self) -> None:
        ctx = self.ctx
        kind = self.kind
        old = self.wave_centre
        allowed = [r for r in self.wave_results if family_allows(self, r[0], r[1])]
        refused = max((r for r in self.wave_results if r not in allowed), key=lambda t: t[1], default=None)
        if refused is not None and refused[1] > self.score:
            self.move(refused[0], refused[1], "refused")               # says why, once per config
        best = max(allowed, key=lambda t: t[1], default=None)
        moved_to = None
        if best is not None and best[1] > self.score and follow(ctx, self, best[0], best[1]):
            moved_to = best[0]
            on_line = best[0].key in self.line_keys
            if on_line:
                for name in self.group:
                    self.growth.pop(name, None)
            else:
                self.grow(self.group, old, best[0])
            if old is not None:
                leaps = [(n, ctx.ladders.rungs_moved(n, old.params.get(n), best[0].params.get(n)))
                         for n in self.group]
                leaps = [(n, r) for n, r in leaps if abs(r) > 1]
                if leaps:
                    best_of = getattr(ctx, "best", None)
                    top = best_of() if callable(best_of) else None
                    holds = top is not None and best[0].key == top[0].key
                    ctx.say("long step c%d %s %+.4f -> %.4f %s", self.id,
                            ", ".join(f"{n} {r:+.0f}" for n, r in leaps),
                            best[1] - self.score, best[1], "round best" if holds else "seat only")
            if best[0].key in self.leap_sent and old is not None:
                jumped = [(n, ctx.ladders.rungs_moved(n, old.params.get(n), best[0].params.get(n)))
                          for n in self.all_names if best[0].params.get(n) != old.params.get(n)]
                ctx.say("centre jump c%d %s %+.4f -> %.4f", self.id,
                        ", ".join(f"{n} {r:+.0f}" for n, r in jumped), best[1] - self.score, best[1])
                self.note_momentum(old, best[0])
            self.move(best[0], best[1], f"{kind} {'line' if on_line else 'grid'}")
        else:
            for name in self.group:
                self.growth.pop(name, None)
        self.adapt_reach(self.group, old, moved_to)
        if self.first is None:
            self.first = (self.score, self.centre.key)
        gain = self.score - self.wave_before
        self.note_gain(self.group, gain)
        self.adapt_leap_share()
        ctx.say("climber %d it %d (%s: %s): %+.4f -> %.4f; waves %s", self.id, self.iteration, kind,
                ",".join(self.group), gain, self.score,
                ", ".join(f"{w[1]}/{w[2]} in {w[3]:.1f}m" for w in self.waves[-2:]))
        self.phase = ""
        self.wave_results = []
        self.batch_live = False
        self.pending, self.extra, self.later, self.leaps = [], [], [], []
        self.gate_futs, self.extra_futs = [], []
        if self.state == "climbing" and self.trailing():
            self.state = "stopped"
            ctx.say("climber %d is %.4f behind the best seat after two flat grids; the seat moves on",
                    self.id, self.gap_to_best())
        elif self.state == "climbing" and (self.stalled() or self.flat_waves()):
            if getattr(self, "compass_only", False):
                pass
            elif self.jumps_enabled:
                self.begin_jumps(self.climb_until)
            else:
                self.state = "stopped"

    async def run_filled(self, until: float, iterations: Optional[int], jumps: bool) -> None:
        """One climber (the Step C fine-tune). A free worker takes its next unsent value. Idle workers do not breed."""
        ctx = self.ctx
        self.jumps_enabled = jumps
        self.climb_until = until
        self._ga: list = []
        stop_at = None if iterations is None else self.iteration + iterations
        while ctx.clock() < until and self.state == "climbing" and not ctx.over():
            if stop_at is not None and self.iteration >= stop_at:
                break
            guided = getattr(self, "compass_only", False)
            if (self.origin or "") == "step C" and not guided and not self.batch_live and self.follow_best():
                continue
            if not self.batch_live and not self.queued():
                if guided:
                    from search.separate import aim_compass
                    await aim_compass(ctx, self, valley=(self.origin or "") != "step C")
                if not self.begin_grid(until):
                    if not (guided and self.state == "climbing" and ctx.clock() < until):
                        break
            if self.gate_done() or self.early_gain():
                self.close_wave()
                if (self.origin or "") == "step C" and self.state != "climbing" and self.follow_best():
                    continue
                continue
            idle = free_workers(ctx)
            guard = 0
            while idle is not None and idle > 0:
                if not self.queued():
                    break
                self.submit_one()
                idle = free_workers(ctx)
                guard += 1
                if guard > int(getattr(ctx, "cap", 1)) + 5:
                    break
            if (self.origin or "") == "step C" and not self.queued():
                idle = free_workers(ctx)
                if idle and self.fill_step_c(idle):
                    continue
            from search.separate import enabled as separate_enabled, fill_free
            if separate_enabled() and (getattr(self, "compass_only", False) or (self.origin or "") == "step C"):
                idle = free_workers(ctx)
                if idle:
                    await fill_free(ctx, self, idle)
            elif separate_enabled():
                idle = free_workers(ctx)
                if idle:
                    await submit_ga(ctx, idle)
            if ctx.clock() >= until or ctx.over():
                break
            if (not self.batch_live and not self.queued() and not self.open_futs()
                    and guided and self.state == "climbing"
                    and not [f for _, f in (getattr(ctx, "ga_live", None) or [])
                             if f is not None and not f.done()]
                    and not [f for f in getattr(self, "_ga", []) if not f.done()]):
                await asyncio.sleep(TICK)
                continue
            await self._wait_one(until)
        if self.batch_live and ctx.clock() >= until:
            self.close_wave()
        if self.state == "climbing" and stop_at is None and not self.batch_live:
            self.state = "stopped"
        self.ended_s = round(ctx.elapsed(), 1)
        self.kind, self.group, self.inflight = "", [], []
        self.batch_live = False

    async def _wait_one(self, until: float) -> None:
        ctx = self.ctx
        pending = self.open_futs() + [f for f in getattr(self, "_ga", []) if not f.done()]
        pending += [f for _, f in (getattr(ctx, "ga_live", None) or []) if f is not None and not f.done()]
        if not pending:
            await asyncio.sleep(0)
            return
        await asyncio.wait(pending, timeout=TICK, return_when=asyncio.FIRST_COMPLETED)
        ctx.publish()

    def next_step(self, name: str, used: int) -> int:
        """The forward step after a win that used `used` rungs. Long ladders: 1, 2, 3, 5, then 5."""
        if len(self.ctx.ladders.rungs.get(name) or []) >= LARGE_RUNGS:
            if used in LARGE_STEPS:
                i = LARGE_STEPS.index(used)
                return LARGE_STEPS[min(i + 1, len(LARGE_STEPS) - 1)]
            return LARGE_STEPS[-1]
        return used + 1

    def grow(self, group: Sequence[str], old: Config, new: Config) -> None:
        lad = self.ctx.ladders
        for n in group:
            r = lad.rungs_moved(n, old.params.get(n), new.params.get(n))
            if r == 0 or self.ctx.space.specs[n].type == "bool":
                self.growth.pop(n, None)
                continue
            d = 1 if r > 0 else -1
            pd, pk = self.growth.get(n, (d, 1))
            if pd != d:
                nxt = 2                                    # the first win in this direction was one rung
            elif len(lad.rungs.get(n) or []) >= LARGE_RUNGS:
                nxt = self.next_step(n, pk)
            else:
                nxt = max(pk + 1, int(round(abs(r))) + 1)
            self.growth[n] = (d, nxt)

    async def jumps(self, until: float) -> None:
        """The stall: score every untried one-rung neighbor. A gain resumes the grid; none is a local optimum."""
        ctx = self.ctx
        accept = float(tune("climb.accept_min"))
        while ctx.clock() < until and self.state == "climbing":
            self.refresh()
            if is_one_rung_local_opt(ctx, self.centre, self.score):
                n = len(one_rung_neighbors(ctx, self.centre))
                self._peaked()
                ctx.say("climber %d: 1-rung local optimum %.4f %s (%d neighbors)",
                        self.id, self.score, self.centre.key[:8], n)
                return
            heads = [c for c in one_rung_neighbors(ctx, self.centre)
                     if ctx.score(c) is None and ctx.fresh(c) and held_by_climber(ctx, c, me=self) is None]
            if not heads:
                n = len(one_rung_neighbors(ctx, self.centre))
                self._peaked()
                ctx.say("climber %d: 1-rung local optimum %.4f %s (%d neighbors, none left to send)",
                        self.id, self.score, self.centre.key[:8], n)
                return
            self.kind, self.group = "jumps", []
            ctx.say("climber %d: certifying 1-rung neighborhood, %d untried", self.id, len(heads))
            results = await self.send(heads, lambda h: [], "jumps", until)
            self.refresh()
            best = max(results, key=lambda t: t[1], default=None)
            if (best is not None and best[1] - self.score >= accept and centre_allowed(self, best[0])):
                ctx.say("climber %d: a 1-rung neighbor gains %+.4f; the grid resumes", self.id, best[1] - self.score)
                self.growth.clear()
                self.reach.clear()
                self.move(best[0], best[1], "jump")
                self.flat.clear()
                return
        if self.state == "climbing" and is_one_rung_local_opt(ctx, self.centre, self.score):
            n = len(one_rung_neighbors(ctx, self.centre))
            self._peaked()
            ctx.say("climber %d: 1-rung local optimum %.4f %s (%d neighbors)",
                    self.id, self.score, self.centre.key[:8], n)

    # ---- bookkeeping ----------------------------------------------------------------------------------------------
    def refresh(self) -> None:
        """The centre's score as the client has it now (an answer-set update rescores everything)."""
        s = self.ctx.score(self.centre)
        if s is not None:
            self.score = s
            if s > self.found:
                self.found = s

    def best_found(self) -> float:
        """The best score this climber has found: its centre, or a higher run of its own that is already back."""
        best = self.found
        if self.score > best:
            best = self.score
        for c in self.waiting.values():
            s = self.ctx.score(c)
            if s is not None and s > best:
                best = s
        if best > self.found:
            self.found = best
        return best

    def move(self, c: Config, s: float, why: str) -> None:
        if not family_allows(self, c, s):
            refused = getattr(self, "_refused", None)
            if refused is None:
                refused = self._refused = set()
            if c.key not in refused:
                refused.add(c.key)
                others = [o.best_found() for o in getattr(self.ctx, "climbers", ())
                          if o is not self and getattr(o, "state", "") == "climbing"]
                self.ctx.say("climber %d keeps its edge-family box: %s %.4f is outside it and not above the other "
                             "seats (best %.4f)", self.id, c.key[:8], s, max(others) if others else float("nan"))
            return
        if family_locked(self) and not in_indel_family(c):
            self.released = True
            self.ctx.say("climber %d leaves its edge-family box for %s %.4f, above every other seat; the lock is "
                         "released", self.id, c.key[:8], s)
        self.centre, self.score = c, s
        if getattr(self, "compass_only", False):
            self.tried_axes = set()
            self.flat.clear()
        if s > self.found:
            self.found = s
        self.history.append((round(self.ctx.elapsed(), 1), round(s, 4), why))
        top = self.ctx.best() if callable(getattr(self.ctx, "best", None)) else None
        if top is not None and c.key == top[0].key:
            self.led_round = True
        other = self.ctx.owner(c, self)
        if other is not None:
            self.state = "joined"
            self.ctx.say("climber %d reached climber %d's centre %s; it stops", self.id, other.id, c.key[:8])

    def _watch(self, c: Config) -> None:
        self.waiting[c.key] = c
        self.owed[c.key] = c
        self.sent_wave[c.key] = getattr(self, "_wave_id", 0)

    def _note_sent(self, futs, configs, kind: str) -> None:
        for f, c in zip(futs, configs):
            if f is not None and hasattr(f, "done"):
                self._sent[f] = (c.key, kind)

    def claim_higher(self) -> bool:
        """A run this climber sent that has come back above the centre takes the centre.

        A live grid keeps the centre until the wave finishes, except for a run sent in an earlier wave: that one
        would have moved the centre had it come back a moment sooner (044200Z islands, 2026-10-16: climber 101's
        82.5627 came back a second after its wave closed and waited nine minutes for the next). A run of the live wave
        moves the centre through early_gain (move at once).
        """
        wave = getattr(self, "_wave_id", 0)
        best = None
        for c in list(self.owed.values()):
            if self.batch_live and self.sent_wave.get(c.key, wave) >= wave:
                continue
            s = self.ctx.score(c)
            if s is None:
                continue
            if s - self.score < float(tune("climb.accept_min")):
                self.owed.pop(c.key, None)
                continue
            if best is None or s > best[1]:
                best = (c, s)
        if best is None or not family_allows(self, best[0], best[1]) or not follow(self.ctx, self, best[0], best[1]):
            return False
        self.learn_leaps([best])
        self.ctx.say("climber %d takes its returned run %s %.4f (centre was %.4f)",
                     self.id, best[0].key[:8], best[1], self.score)
        old = self.centre
        if best[0].key in self.leap_sent and old is not None:
            ctx = self.ctx
            jumped = [(n, ctx.ladders.rungs_moved(n, old.params.get(n), best[0].params.get(n)))
                      for n in self.all_names if best[0].params.get(n) != old.params.get(n)]
            ctx.say("centre jump c%d (returned) %s %+.4f -> %.4f", self.id,
                    ", ".join(f"{n} {r:+.0f}" for n, r in jumped), best[1] - self.score, best[1])
            self.note_momentum(old, best[0])
        self.owed.pop(best[0].key, None)
        self.waiting.pop(best[0].key, None)
        self.pending, self.extra, self.later, self.leaps = [], [], [], []
        self.gate_futs, self.extra_futs, self.inflight = [], [], []
        self.batch_live = False
        self.resume(best[0], best[1], "returned")
        if self.state == "climbing":
            self.begin_grid(getattr(self, "climb_until", 1e9))
        return True

    def follow_best(self) -> bool:
        """The fine-tune climbs the round's best. A higher score restarts the climb from that config."""
        best = getattr(self.ctx, "best", None)
        if not callable(best):
            return False
        top = best()
        if top is None:
            return False
        config, score = top
        if not isinstance(score, (int, float)) or isinstance(score, bool):
            return False
        score = float(score)
        if config.key == self.centre.key or score <= self.score:
            return False
        if not family_allows(self, config, score):
            return False
        self.ctx.say("climber %d restarts from the best %s %.4f (centre was %.4f)",
                     self.id, config.key[:8], score, self.score)
        self.pending, self.extra, self.later, self.leaps = [], [], [], []
        self.waiting.clear()
        self.gate_futs, self.extra_futs, self.inflight = [], [], []
        self.batch_live = False
        self.growth.clear()
        self.reach.clear()
        self.flat.clear()
        self.recent.clear()
        self.kind, self.group = "", []
        self.ended_s = None
        self.state = "climbing"
        self.centre, self.score = config, score
        if score > self.found:
            self.found = score
        self.history.append((round(self.ctx.elapsed(), 1), round(score, 4), "best"))
        self.begin_grid(getattr(self, "climb_until", 1e9))
        return True

    def adopt(self, c: Config, s: float, why: str = "the GA best") -> None:
        """Leave this hill for a higher config. Runs already out still score, and a higher one takes the centre back."""
        if not family_allows(self, c, s):
            return
        was = self.best_found()
        self.ctx.say("climber %d takes %s %s %.4f (its best was %.4f, centre %s)",
                     self.id, why, c.key[:8], s, was, self.centre.key[:8])
        self.pending, self.extra, self.later, self.leaps = [], [], [], []
        self.waiting.clear()
        self.gate_futs, self.extra_futs, self.inflight = [], [], []
        self.batch_live = False
        self.growth.clear()
        self.reach.clear()
        self.flat.clear()
        self.kind, self.group = "", []
        self.resume(c, s, "ga")
        if self.state == "climbing" and not getattr(self.ctx, "ga_opening", False):
            self.begin_grid(self.climb_until)

    def resume(self, c: Config, s: float, why: str) -> None:
        """Climb again from `c` (Step C found it): fresh grids, no growth carried over."""
        self.ended_s = None
        self.growth.clear()
        self.reach.clear()
        self.flat.clear()
        self.recent.clear()
        self.state = "climbing"
        self.move(c, s, why)

    def row(self) -> dict:
        """The basin row main's review reads (first_key, best_key) and what the dashboard draws."""
        back = sum(1 for f in self.inflight if f.done())
        return {"id": self.id, "origin": self.origin, "parent": self.parent, "opened": round(self.opened, 4),
                "first": None if self.first is None else round(self.first[0], 4),
                "first_key": None if self.first is None else self.first[1],
                "best": round(self.score, 4), "best_key": self.centre.key, "status": self.state,
                "born_s": self.born_s, "ended_s": self.ended_s,
                "best_s": self.history[-1][0] if self.history else self.born_s,
                "iterations": self.iteration, "evaluated": self.evaluated,
                "wave": {"kind": self.kind, "group": self.group, "sent": len(self.inflight), "back": back}
                if self.inflight and back < len(self.inflight) else None,
                "growth": {n: d * k for n, (d, k) in self.growth.items() if k > 1},
                "flat": len(self.flat), "expensive": len(self.read_names) + len(self.assembly_names),
                "waves": self.waves[-10:], "history": self.history[-20:],
                "entry_gap": self.entry_gap, "entry_peers": list(self.entry_peers),
                "led_round": self.led_round}


SIMILAR_SETTINGS = 4          # fewer than this many searched-setting differences is the same hill
PROBE_LOSS = 0.10          # a one-rung probe this far below the plateau closes that direction
REPLACEMENT_FLOOR = 1.0    # a scored config further than this below the best is not seated
SCOUT_AXES = 6             # settings a scout moves, so the mean distance clears the seat gap
SCOUT_WAVES = 2            # far-scout waves for one empty seat; then breeding continues alone


def same_hill(ctx, config: Config, other: Config) -> bool:
    """True when two configs are the same key, differ in fewer than four searched settings, or share a callset."""
    if config.key == other.key:
        return True
    diffs = 0
    for n in ctx.space.names:
        if config.params.get(n) != other.params.get(n):
            diffs += 1
            if diffs >= SIMILAR_SETTINGS:          # the rest of the settings cannot change this answer
                break
    if diffs < SIMILAR_SETTINGS:
        return True
    dist = getattr(ctx, "dist", None)
    if dist is None:
        return False
    d, kind = dist.between(config, other)
    return kind == "callset" and d < float(tune("candidates.min_distance"))


def prune_ga_live(ctx) -> None:
    """Keep only GA children whose GATK run is still in flight. A finished child is not a seat."""
    live = []
    for item in getattr(ctx, "ga_live", None) or []:
        c, f = item if isinstance(item, tuple) and len(item) == 2 else (item, None)
        if f is not None and f.done():
            continue
        live.append((c, f))
    ctx.ga_live = live
    ctx.ga_hills = [c for c, _ in live]


def held_by_ga(ctx, config: Config) -> bool:
    """A GA child still in flight that is the same hill as `config` (v9.18). Exact keys stay in `fresh`."""
    for other in getattr(ctx, "ga_hills", ()) or ():
        if same_hill(ctx, config, other):
            return True
    return False


async def held_by_ga_yielding(ctx, config: Config) -> bool:
    """held_by_ga, yielding every 50 ms so a long hill list does not freeze the dashboard."""
    mark = time.monotonic()
    for other in getattr(ctx, "ga_hills", ()) or ():
        if same_hill(ctx, config, other):
            return True
        now = time.monotonic()
        if now - mark > 0.05:
            mark = now
            await asyncio.sleep(0)
    return False


def held_by_climber(ctx, config: Config, me=None):
    """Another climber already has this config, one that differs in fewer than four settings, or the same callset.

    A climber's own grid is not held: `me` is that climber. On chr12 the three centres differed in two settings
    and called the same variants. On chr22 climber 1 and climber 2 differed in three settings and still called
    the same variants, so three is not a new hill either.
    """
    for cl in getattr(ctx, "climbers", ()):
        if me is not None and cl is me:
            continue
        # A peaked or landed climber is finished. Keeping its ball closed stops a Step C landing
        # one setting away, and it makes certify spawn a new climber for every neighbor it cannot send.
        if getattr(cl, "state", "") != "climbing":
            continue
        centre = getattr(cl, "centre", None)
        if centre is None:
            continue
        if same_hill(ctx, config, centre):
            return cl
    return None


def _differs(ctx, centre: Config, other: Config) -> Optional[str]:
    """The one setting where `other` differs from `centre`, or None when it is the centre or differs in two."""
    name = None
    for n in ctx.space.names:
        if other.params.get(n) != centre.params.get(n):
            if name is not None:
                return None
            name = n
    return name


def _probe_value(ctx, centre: Config, name: str, base: float, known: Dict[int, float]):
    """The nearest untried rung of `name`. A side that already lost PROBE_LOSS is closed; a tie or a smaller loss
    steps one rung further. The step is ladders.shift, so a value between rungs moves to the next rung."""
    current = centre.params.get(name)
    best = None
    for direction in (1, -1):
        k = direction
        while True:
            v = ctx.ladders.shift(name, current, k)
            if v is None:
                break
            i = int(round(ctx.ladders.index(name, v)))
            if i not in known:
                if best is None or abs(k) < best[0]:
                    best = (abs(k), v)
                break
            if base - known[i] >= PROBE_LOSS:
                break
            k += direction
    return None if best is None else best[1]


def next_probes(ctx, centre: Config, base: float) -> List[Config]:
    """One untried rung of each setting the plateau holds at a single value.

    The plateau is every config with score `base`. A setting that already varies inside it is the same callset
    and is not probed. A probe changes one locked setting and leaves the rest of the centre alone.
    """
    members = [(c, s) for c, s in ctx.scored() if s == base]
    if len(members) < 2:
        return []
    # A one-setting neighbour that tied is a rung to step past, not a second variety.
    varieties = [c for c, _ in members
                 if sum(c.params.get(n) != centre.params.get(n) for n in ctx.space.names) != 1]
    names = [n for n in ctx.space.names if ctx.ladders.movable(n)]
    locked = []
    for n in names:
        if len({repr(c.params.get(n)) for c in varieties}) == 1:
            locked.append(n)
    known: Dict[str, Dict[int, float]] = {n: {} for n in locked}
    for c, s in ctx.scored():
        n = _differs(ctx, centre, c)
        if n in known:
            known[n][int(round(ctx.ladders.index(n, c.params.get(n))))] = s
    seen = {centre.key}
    for cl in getattr(ctx, "climbers", ()):
        seen.update(c.key for c in getattr(cl, "pending", ()))
        seen.update(getattr(cl, "waiting", {}))
        seen.update(c.key for c in getattr(cl, "later", ()))
    out = []
    for n in locked:
        value = _probe_value(ctx, centre, n, base, known[n])
        if value is None:
            continue
        c = ctx.space._assemble({**centre.params, n: value})
        if c is None or c.key in seen or not ctx.fresh(c):
            continue
        seen.add(c.key)
        out.append(c)
    return out


def _queue_probes(ctx, n: int) -> Optional[list]:
    """On a plateau, queue one-rung probes once, then breed on the workers those probes do not use.

    None means the caller breeds: the best score is still one config, the probes for this score are
    already queued or finished, or there is no climbing climber on that score. [] means send nothing.
    """
    if n <= 0:
        return []
    top = ctx.best()
    if top is None:
        return None
    centre, base = top
    owners = [cl for cl in getattr(ctx, "climbers", ())
              if cl.state == "climbing" and cl.score == base]
    if not owners:
        return None
    owner = min(owners, key=lambda cl: cl.id)
    if getattr(ctx, "_probe_base", None) == base and getattr(ctx, "_probe_done", False):
        return None
    keys = getattr(ctx, "_probe_keys", ())
    if getattr(ctx, "_probe_base", None) == base and (
            any(c.key in keys for c in owner.pending) or any(k in keys for k in owner.waiting)):
        return None
    probes = next_probes(ctx, owner.centre, base)
    ctx._probe_base = base
    if not probes:
        if len([1 for _, s in ctx.scored() if s == base]) >= 2:
            ctx._probe_done = True
        return None
    take = probes[:n]
    owner.pending.extend(take)
    ctx._probe_keys = {c.key for c in take}
    ctx._probe_done = False
    ctx.say("plateau: %d one-rung probes on climber %d; the other free workers breed", len(take), owner.id)
    return None


def far_from_climbers(ctx, config: Config) -> bool:
    """A bred child has to be a different hill from every seated climber.

    Callset distance when both configs have a sketch, otherwise the mean setting difference.
    The thresholds are the same ones that keep two seats apart. With no distance object the
    two-setting rule in `held_by_climber` is the check.
    """
    dist = getattr(ctx, "dist", None)
    if dist is None:
        return True
    mc = float(tune("candidates.min_distance"))
    ms = float(tune("candidates.min_setting_distance"))
    for cl in getattr(ctx, "climbers", ()):
        centre = getattr(cl, "centre", None)
        if centre is None:
            continue
        if config.key == centre.key or not dist.far(config, centre, mc, ms):
            return False
    return True


def ga_keep(ctx, config: Config) -> bool:
    """A GA child is sent when it is not a seated climber's neighbourhood and not a tied callset."""
    if held_by_climber(ctx, config) is not None:
        return False
    if not far_from_climbers(ctx, config):
        return False
    tied = getattr(ctx, "tied", None)
    return tied is None or not tied.holds(config.params)


async def submit_ga(ctx, n: int, gen: int = 1) -> list:
    """Children for workers no candidate update can use. Priority 0, behind every candidate update.

    Step B breeds within the both / core / fp groups (gen 1). On a plateau it also sends one-rung probes,
    once, and the workers those probes do not use keep breeding. The fine-tune does not call this.
    A child that is the same or similar to a climber's centre, or to another GA child still in flight, or closer
    than the seat distance, is not sent.
    """
    if n <= 0 or ctx.over():
        return []
    if gen < 3:
        held = _queue_probes(ctx, n)
        if held is not None:
            return held
    step = getattr(ctx, "step_a", None)
    if step is None:
        return []
    skipped = repeats = close = similar = 0
    prune_ga_live(ctx)
    hills = getattr(ctx, "ga_hills", None)
    if hills is None:
        hills = []
        ctx.ga_hills = hills

    async def accept(c: Config) -> bool:
        nonlocal skipped, repeats, close, similar
        if piles_enabled():
            if not ctx.fresh(c):
                skipped += 1
                return False
            tied = getattr(ctx, "tied", None)
            if tied is not None and tied.holds(c.params):
                repeats += 1
                return False
            names = list(ctx.space.names)
            if getattr(ctx, "ga_opening", False) or ga_hot(ctx):
                top = ctx.best() if callable(getattr(ctx, "best", None)) else None
                if top is not None:
                    d = sum(c.params.get(n) != top[0].params.get(n) for n in names)
                    if 1 <= d <= 3:
                        hills.append(c)
                        return True
            if await held_by_ga_yielding(ctx, c):
                similar += 1
                return False
            for cl in getattr(ctx, "climbers", ()):
                centre = getattr(cl, "centre", None)
                if centre is None:
                    continue
                if sum(c.params.get(n) != centre.params.get(n) for n in names) < 4:
                    skipped += 1
                    return False
            hills.append(c)
            return True
        if held_by_climber(ctx, c) is not None:
            skipped += 1
            return False
        if not far_from_climbers(ctx, c):
            close += 1
            return False
        tied = getattr(ctx, "tied", None)
        if tied is not None and tied.holds(c.params):
            repeats += 1
            return False
        if await held_by_ga_yielding(ctx, c):
            similar += 1
            return False
        hills.append(c)
        return True

    # The peak interval is Generation 1's parent pool. Step B breeds from every scored hill, so a
    # compass that leaves most of the fleet free still has parents outside that interval.
    prev_bound = getattr(step, "_bound_pool", True)
    step._bound_pool = False
    k = max(1, int(tune("climb.model_filter")))      # breed k times as many; the model keeps the best 1/k
    try:
        children = step.breed(gen, n * k, why="GA" if gen < 3 else "GA cross", accept=accept,
                              elite=bool(ELITE_LEFTOVER) and gen < 3)
        if hasattr(children, "__await__"):
            children = await children
    finally:
        step._bound_pool = prev_bound
    if k > 1:
        from search.model_filter import keep_best
        children = keep_best(ctx, list(children), n)
    if skipped:
        ctx.say("GA: %d children left out, same or similar to a climber", skipped)
    if close:
        ctx.say("GA: %d children left out, too close to a climber", close)
    if repeats:
        ctx.say("GA: %d children left out, a callset already scored", repeats)
    if similar:
        ctx.say("GA: %d children left out, same or similar to another GA child", similar)
    keys = getattr(ctx, "ga_keys", None)
    if keys is None:
        keys = set()
        ctx.ga_keys = keys
    futs = []
    live = list(getattr(ctx, "ga_live", None) or [])
    for c in children:
        keys.add(c.key)
        fs = ctx.propose([c], "G", ["big"], priority=0)
        futs += fs
        live.append((c, fs[0] if fs else None))
    ctx.ga_live = live
    ctx.ga_hills = [c for c, _ in live]
    held = getattr(ctx, "ga_futs", None)
    if held is not None:
        held.extend(futs)
    return futs


def free_workers(ctx) -> Optional[int]:
    """Workers with no run queued or in flight, or None when the executor cannot say (the unit tests)."""
    client = getattr(ctx, "client", None)
    fn = getattr(client, "idle_slots", None) if client is not None else None
    if fn is None:
        return None
    try:
        n = fn()
    except Exception:  # noqa: BLE001
        return None
    return None if n is None else int(n)


def ranked_ga(ctx) -> List[Tuple[Config, float]]:
    """Scored GA configs, highest first. A config only the `best()` helper knows still counts."""
    keys = getattr(ctx, "ga_keys", None) or ()
    if not keys:
        return []
    found: Dict[str, Tuple[Config, float]] = {}
    best = getattr(ctx, "best", None)
    if callable(best):
        hit = best()
        if hit is not None and hit[0].key in keys:
            found[hit[0].key] = hit
    scored = getattr(ctx, "scored", None)
    if scored is not None:
        for c, s in scored():
            if c.key in keys and s is not None and (c.key not in found or s > found[c.key][1]):
                found[c.key] = (c, s)
    return sorted(found.values(), key=lambda t: (-t[1], t[0].key))


def highest_ga(ctx) -> Optional[Tuple[Config, float]]:
    """The highest scored config the GA proposed, even when a climber's own run is higher."""
    rows = ranked_ga(ctx)
    return rows[0] if rows else None


def family_centre(ctx) -> Optional[Tuple[Config, float]]:
    """The best edge-family config already seated, else the best scored family config."""
    hit: Optional[Tuple[Config, float]] = None
    for cl in getattr(ctx, "climbers", ()):
        c = getattr(cl, "centre", None)
        if c is None or not in_indel_family(c):
            continue
        s = float(getattr(cl, "score", 0.0) or 0.0)
        if hit is None or s > hit[1]:
            hit = (c, s)
    if hit is not None:
        return hit
    nxt = far_replacement(ctx)
    if nxt is None:
        return None
    return nxt[0], nxt[1]


def is_far_origin(cl: Climber) -> bool:
    """True when this seat was opened as a far hill (B27) and must stay off the best hill (B30)."""
    return (cl.origin or "").lower().startswith("far")


def family_locked(cl: Climber) -> bool:
    """Opened as the reserved edge-family hill. Leftover GA and follow_best may not walk it back to this shelf.

    A seat whose lock was released (see `family_allows`) is not locked any more."""
    if getattr(cl, "released", False):
        return False
    if is_far_origin(cl):
        return True
    origin = (cl.origin or "").lower()
    start = getattr(cl, "start", None)
    if start is None or not in_indel_family(start):
        return False
    return origin.startswith("step c") or origin.startswith("certify")


FAMILY_REACH = 0.05          # Step C and certify climb the edge family only when it is this close to the round's best


def family_in_reach(ctx, best_score: float):
    """The best edge-family config when it is within FAMILY_REACH of `best_score`, else None. Further below, a
    fine-tune of eight minutes does not close the gap (Step C gained 0.02-0.03 in the 2026-10-16 runs)."""
    fam = family_centre(ctx)
    if fam is None or fam[1] < best_score - FAMILY_REACH:
        return None
    return fam


def family_allows(cl: Climber, c: Config, s: Optional[float] = None) -> bool:
    """A family seat may only take a config that is still in that box, or one scoring above every other climbing
    seat of its search (the lock is then released when it moves: `Climber.move`).

    044200Z islands (2026-10-16): island 1's far seat found 82.5627 with min_pruning 2, inside the API box, while the
    other seat was at 82.4066; the lock kept it on 82.3514 and logged +0.0000."""
    if not family_locked(cl) or in_indel_family(c):
        return True
    return s is not None and above_other_seats(cl, s)


def above_other_seats(cl: Climber, s: float) -> bool:
    others = [o.best_found() for o in getattr(cl.ctx, "climbers", ())
              if o is not cl and getattr(o, "state", "") == "climbing"]
    return not others or s > max(others)


def near_hill_ga(ctx, config: Config, near: Sequence[Climber]) -> bool:
    """True when `config` sits on a near climber's hill (four settings, or callset / setting distance)."""
    if not near:
        return False
    dist = getattr(ctx, "dist", None)
    mc = float(tune("candidates.min_distance"))
    ms = float(tune("candidates.min_setting_distance"))
    for cl in near:
        centre = getattr(cl, "centre", None)
        if centre is None:
            continue
        if same_hill(ctx, config, centre):
            return True
        if dist is not None and (config.key == centre.key or not dist.far(config, centre, mc, ms)):
            return True
    return False


def adopt_reference(ctx, live: Sequence[Climber]) -> bool:
    """While the GA owns the fleet, the highest score within 0.10 of a pile becomes that climber's centre.

    Once the grid is walking, the pile does not move the centre. The grid writes the setting it found.
    The highest score becomes the fingerprint, and a pile it overlaps is merged into this one.
    """
    if not getattr(ctx, "_piles", None) or not ga_hot(ctx):
        return False
    climbing = [cl for cl in live if cl.state == "climbing"]
    if not climbing:
        return False
    refresh_piles(ctx)
    for cl in climbing:
        own = own_pile(ctx, cl)
        if own is None:
            continue
        config, score = best_member(own)
        if config.key == cl.centre.key or score <= cl.score:
            continue
        if not centre_allowed(cl, config):
            continue
        piles = getattr(ctx, "_piles", None) or []
        if own in piles:
            claim_peak(ctx, piles, piles.index(own), config, score)
        cl.adopt(config, score, "the highest within 0.10")
        return True
    return False


def _free_member(pile, held: set):
    """The highest config in this pile that no climbing or peaked climber already holds."""
    best = None
    for config in pile.members:
        if config.key in held:
            continue
        score = float(pile.scores.get(config.key, pile.score))
        if best is None or score > best[0]:
            best = (score, config)
    return best


def switch_to_candidate(ctx, live: Sequence[Climber]) -> bool:
    """The climber on the lowest pile moves when a candidate's score is higher.

    A config another climber already holds, including one who has peaked, is not that candidate's
    score. Moving onto it would join that climber and free the seat, and the next seat would do the same.
    """
    if not getattr(ctx, "_piles", None):
        return False
    climbing = [cl for cl in live if cl.state == "climbing"]
    if not climbing:
        return False
    piles = list(refresh_piles(ctx))
    if not piles:
        return False
    held = {cl.centre.key for cl in getattr(ctx, "climbers", ())
            if getattr(cl, "state", "") in ("climbing", "peak") and getattr(cl, "centre", None) is not None}
    lowest = None
    taken = set()
    for cl in climbing:
        own = own_pile(ctx, cl)
        if own is None or own not in piles:
            continue
        i = piles.index(own)
        if i in taken:
            continue
        taken.add(i)
        score = best_member(own)[1]
        if lowest is None or score < lowest[0] or (score == lowest[0] and cl.id < lowest[1]):
            lowest = (score, cl.id, cl)
    if lowest is None:
        return False
    low_score, _, cl = lowest
    best = None
    for i, pile in enumerate(piles):
        if i in taken or already_peaked(ctx, pile):
            continue
        free = _free_member(pile, held)
        if free is None:
            continue
        score, config = free
        if score > low_score and (best is None or score > best[0]):
            best = (score, i, config)
    if best is None:
        return False
    score, i, config = best
    cl.origin = piles[i].label(i)
    cl.adopt(config, score, "the higher candidate")
    return True


def adopt_ga_best(ctx, live: Sequence[Climber]) -> bool:
    """A near seat takes a near-hill GA; a far seat takes only a GA that is far from every near centre.

    v9.21: the highest score within 0.10 of a pile is that pile's reference, and its climber develops from it.
    A child outside every ball does not move an existing climber; an empty seat takes it when it is farther than 0.10.
    """
    if piles_enabled() and getattr(ctx, "_piles", None):
        return _adopt_pile(ctx, live)
    return _adopt_ga_best_legacy(ctx, live)


def _adopt_pile(ctx, live: Sequence[Climber]) -> bool:
    from search.piles import home
    climbing = [cl for cl in live if cl.state == "climbing"]
    if not climbing:
        return False
    seated = {cl.centre.key for cl in climbing}
    piles = refresh_piles(ctx)
    ignored = getattr(ctx, "ga_ignored", None)
    if ignored is None:
        ignored = set()
        ctx.ga_ignored = ignored
    for c, s in ranked_ga(ctx):
        if c.key in seated:
            continue
        i = home(ctx, c)
        if i is None:
            continue
        pool = [cl for cl in climbing if (cl.origin or "") == piles[i].label(i) and cl.best_found() < s]
        if not pool:
            pool = [cl for cl in climbing if own_pile(ctx, cl) is piles[i] and cl.best_found() < s]
        if not pool:
            continue
        taker = min(pool, key=lambda cl: (cl.best_found(), cl.id))
        if not family_allows(taker, c) or not centre_allowed(taker, c):
            continue
        other = held_by_climber(ctx, c, me=taker)
        if other is not None and (other.origin or "") != piles[i].label(i):
            if c.key not in ignored and other.best_found() < s:     # silent when that climber is already as high
                ignored.add(c.key)
                ctx.say("GA best %s is the same or similar to climber %d; it stays there", c.key[:8], other.id)
            continue
        taker.adopt(c, s)
        return True
    return False


def _adopt_ga_best_legacy(ctx, live: Sequence[Climber]) -> bool:
    """A near seat takes a near-hill GA; a far seat takes only a child still in the edge family.

    Best means the highest score that climber's own runs have returned. Step B checks this on every pass,
    so the handoff does not wait for a grid to end. A GA child already sitting on a climber is skipped.
    On v9.17 the far seat was the lowest climber, so every free GA — including the best hill — landed on it.
    A four-setting AGGRESSIVE child of this shelf looks 'far' by setting count; without the family lock
    it walks the reserved seat back in one handoff (replay 20261010T044200Z, t+11m).
    """
    climbing = [cl for cl in live if cl.state == "climbing"]
    if not climbing:
        return False
    seated = {cl.centre.key for cl in climbing}
    near = [cl for cl in climbing if not is_far_origin(cl)]
    far = [cl for cl in climbing if is_far_origin(cl)]
    if not near:
        near, far = climbing, []
    ignored = getattr(ctx, "ga_ignored", None)
    if ignored is None:
        ignored = set()
        ctx.ga_ignored = ignored
    # Highest GA score first: once a score is no higher than every climber's best, no later one can be taken.
    # Walking the whole list on every pass grew with every GA result.
    floor = min(cl.best_found() for cl in climbing)
    for c, s in ranked_ga(ctx):
        if s <= floor:
            break
        if c.key in seated:
            continue
        if near_hill_ga(ctx, c, near):
            pool = [cl for cl in near if cl.best_found() < s]
        else:
            pool = [cl for cl in far if cl.best_found() < s and family_allows(cl, c, s)]
            if not pool:
                pool = [cl for cl in near if cl.best_found() < s]
        if not pool:
            continue
        taker = min(pool, key=lambda cl: (cl.best_found(), cl.id))
        if not family_allows(taker, c, s):
            continue
        other = held_by_climber(ctx, c, me=taker)
        if other is not None:
            # A higher GA config on another climber's own hill becomes that climber's centre. 044200Z islands
            # (2026-10-16): island 1's GA scored 82.5737 on climber 101's hill (callsets under 0.05 apart) while
            # that climber climbed from 82.3514; "it stays there" left it for good, 22 times in that run.
            if int(tune("climb.own_hill_adopt")) and other.best_found() < s and family_allows(other, c, s):
                other.adopt(c, s, why="the GA best on its own hill")
                return True
            if c.key not in ignored and other.best_found() < s:     # silent when that climber is already as high
                ignored.add(c.key)
                ctx.say("GA best %s is the same or similar to climber %d; it stays there", c.key[:8], other.id)
            continue
        taker.adopt(c, s)
        return True
    return False


def far_replacement(ctx) -> Optional[Tuple[Config, float, str]]:
    """The best scored config that has left the rim of the best and is a different hill from every seat.

    The Step A backup list is fixed at the start of the climb. A GA result, or a scout, can be a better far
    seat and is not on that list. Callset distance when both configs have a sketch, otherwise the mean
    setting difference. A config that still shares the best's PCR model, or that differs
    in fewer than four settings from a centre, or more than REPLACEMENT_FLOOR below the best, is left out.
    The scan is reused until a new score arrives or a centre moves.
    """
    top = ctx.best()
    dist = getattr(ctx, "dist", None)
    if top is None or dist is None:
        return None
    _, best = top
    centres = [cl.centre for cl in getattr(ctx, "climbers", ()) if getattr(cl, "centre", None) is not None]
    if not centres:
        return None
    latest = getattr(getattr(ctx, "client", None), "latest", None)
    stamp = (round(best, 6), tuple(c.key for c in centres), None if latest is None else len(latest))
    old = getattr(ctx, "_far_stamp", None)
    if stamp[2] is not None and old is not None and old[0] == stamp[0] and old[1] == stamp[1]:
        hit = getattr(ctx, "_far_hit", None)
        # A full pass over the scored pool is the slow part. Reuse a hit until the best or a centre
        # changes, and reuse a miss until dozens of new scores exist.
        if hit is not None or stamp[2] - old[2] < 50:
            return hit
    seated = set()
    for cl in ctx.climbers:
        for c in (getattr(cl, "start", None), cl.centre):
            if c is not None:
                seated.add(c.key)
    names = list(ctx.space.names)
    mc = float(tune("candidates.min_distance"))
    ms = float(tune("candidates.min_setting_distance"))
    tied = getattr(ctx, "tied", None)
    pool = [(s, c) for c, s in ctx.scored()
            if c.key not in seated and s >= best - FAMILY_SEAT_FLOOR]
    pool.sort(key=lambda t: -t[0])
    peak = top[0]
    hit = None
    for s, c in pool:
        if not left_the_rim(c, peak):
            continue
        if any(sum(c.params.get(n) != centre.params.get(n) for n in names) < SIMILAR_SETTINGS
               for centre in centres):
            continue
        if tied is not None and tied.holds(c.params):
            continue
        if any(c.key == centre.key or not dist.far(c, centre, mc, ms) for centre in centres):
            continue
        kind = "GA" if c.key in getattr(ctx, "ga_keys", ()) else "scored"
        hit = (c, s, f"far ({kind})")
        break
    if stamp[2] is not None:
        ctx._far_stamp = stamp
        ctx._far_hit = hit
    return hit


def scout_bases(ctx) -> List[Config]:
    """Centres a far scout starts from: the core and fewest-FP parents, else the best seated centre."""
    step = getattr(ctx, "step_a", None)
    bases: List[Config] = []
    groups = getattr(step, "groups", None)
    if callable(groups):
        g = groups()
        if asyncio.iscoroutine(g):
            g.close()
            g = None
        if g:
            for name in ("core", "fp"):
                pts = g.get(name) or []
                if pts:
                    c = pts[0].config
                    if all(c.key != b.key for b in bases):
                        bases.append(c)
    if bases:
        return bases
    centres = [(cl.centre, cl.score) for cl in getattr(ctx, "climbers", ())
               if getattr(cl, "centre", None) is not None]
    if not centres:
        return []
    return [max(centres, key=lambda t: t[1])[0]]


def far_scouts(ctx, n: int) -> List[Config]:
    """Configs built far from every climber, starting from the far parent ends.

    A GA child changes two settings by one rung and stays inside the seat distance. A scout sets
    EDGE_SCOUT_AXES validator-edge values, which is enough for the mean setting difference to
    clear the seat distance and to sit outside the published API box.
    """
    if n <= 0:
        return []
    centres = [(cl.centre, cl.score) for cl in getattr(ctx, "climbers", ())
               if getattr(cl, "centre", None) is not None]
    starts = scout_bases(ctx)
    if not centres or not starts:
        return []
    names = [name for name in ctx.space.names if ctx.ladders.movable(name)]
    if len(names) < 2:
        return []
    dist = getattr(ctx, "dist", None)
    mc = float(tune("candidates.min_distance"))
    ms = float(tune("candidates.min_setting_distance"))
    axes = [name for name in EDGE_SCOUT_CHOICES if name in ctx.space.names]
    out, seen = [], set()
    tries = 0
    while len(out) < n and tries < 20 * n:
        tries += 1
        base = starts[tries % len(starts)]
        params = dict(base.params)
        k = min(EDGE_SCOUT_AXES, len(axes))
        chosen = ctx.rng.sample(axes, k) if k else []
        for name in chosen:
            params[name] = ctx.rng.choice(EDGE_SCOUT_CHOICES[name])
        c = ctx.space._assemble(params)
        if c is None or c.key in seen or not ctx.fresh(c):
            continue
        if held_by_climber(ctx, c) is not None:
            continue
        if dist is not None and any(not dist.far(c, centre, mc, ms) for centre, _ in centres):
            continue
        peak = max(centres, key=lambda t: t[1])[0]
        if not left_the_rim(c, peak):
            continue
        seen.add(c.key)
        out.append(c)
    return out


def ensure_far_seat(ctx, live: List[Climber]) -> bool:
    """If no climber is still in the edge family, the lowest non-best seat takes the best family config.

    Two AGGRESSIVE centres that differ in four genotyping settings look like two hills, so the old
    same-hill check left the stolen far seat there. Family-or-not is the test.
    """
    climbing = [cl for cl in live if cl.state == "climbing" and getattr(cl, "centre", None) is not None]
    if len(climbing) < 2:
        return False
    if any(in_indel_family(cl.centre) for cl in climbing):
        return False
    # A released far seat left the box because something outside it scored higher: the box has had its seat, and no
    # other seat is moved into it (044200Z islands, 2026-10-16: climber 0 went from 82.5557 to 82.5274 that way).
    if any(getattr(cl, "released", False) for cl in getattr(ctx, "climbers", ())):
        return False
    nxt = far_replacement(ctx)
    if nxt is None:
        return False
    best = max(climbing, key=lambda cl: (cl.score, -cl.id))
    taker = min((cl for cl in climbing if cl.id != best.id), key=lambda cl: (cl.best_found(), cl.id))
    ctx.say("step B: far seat %s %.4f (%s) on climber %d", nxt[0].key[:8], nxt[1], nxt[2], taker.id)
    taker.origin = f"far seat ({nxt[2]})"
    taker.adopt(nxt[0], nxt[1], why="the edge-family seat")
    return True


def _next_far(ctx, queue: List[Tuple[Config, float, str]]) -> Optional[Tuple[Config, float, str]]:
    mc, ms = float(tune("candidates.min_distance")), float(tune("candidates.min_setting_distance"))
    while queue:
        c, s, origin = queue.pop(0)
        if all(c.key not in (d.start.key, d.centre.key) and ctx.dist.far(c, d.centre, mc, ms) for d in ctx.climbers):
            return c, s, origin
    return None


LATE_SEAT_MIN = 20.0         # no new far seat this close to the end of Step B (scaled minutes)


def late_for_seat(ctx, until: float) -> bool:
    """Too late in Step B to open a far seat: 044200Z islands (2026-10-16) seated two far configs 0.25-0.33 below the best
    with 10 minutes left; climbs needed 30 or more minutes to gain that much."""
    scale = getattr(ctx, "T", None)
    window = scale(60.0 * LATE_SEAT_MIN) if callable(scale) else 60.0 * LATE_SEAT_MIN
    return until - ctx.clock() < window


def head_reserve(live: Sequence[Climber]) -> int:
    """Workers the GA leaves free for the climbers' next heads: a climber whose wave is at least half back keeps as
    many as that wave had heads, up to climb.head_reserve. 044200Z islands (2026-10-16): a new wave's heads queued
    behind GA runs, and the median grid wave took 4.0 minutes against 2.8 before islands."""
    cap = int(tune("climb.head_reserve"))
    if cap <= 0:
        return 0
    keep = 0
    for cl in live:
        if cl.state != "climbing" or not cl.batch_live or cl.pending or not cl.gate_futs:
            continue
        back = sum(1 for f in cl.gate_futs if f.done())
        if back * 2 >= len(cl.gate_futs):
            keep += min(cap, len(cl.gate_futs))
    return keep


def note_reserve(ctx, held: int) -> None:
    """Worker-seconds the head reserve kept idle (free workers the GA would have taken), for the notes and the log."""
    now = ctx.clock()
    last = getattr(ctx, "_reserve_last", None)
    if last is not None and now > last[0]:
        ctx.reserve_ws = getattr(ctx, "reserve_ws", 0.0) + last[1] * (now - last[0])
    ctx._reserve_last = (now, held)


def ga_beat_lowest(ctx, live: Sequence[Climber]) -> bool:
    """A free GA score above the lowest climbing seat, whether or not a seat may take it (an edge-family lock)."""
    climbing = [cl for cl in live if cl.state == "climbing"]
    top = highest_ga(ctx)
    return bool(climbing) and top is not None and top[1] > min(cl.best_found() for cl in climbing)


async def open_with_ga(ctx, live: List[Climber], until: float) -> None:
    """The start of Step B: the whole fleet breeds until a free GA beats the lowest climber, or the cap."""
    if not live or ctx.over() or ctx.clock() >= until:
        return
    if free_workers(ctx) is None:
        return
    if float(tune("climb.ga_open_min")) <= 0:          # off (the default): Step B starts with the grids
        return
    cap_s = 60 * float(tune("climb.ga_open_min"))
    T = getattr(ctx, "T", None)
    if callable(T):
        cap_s = T(cap_s)
    stop = min(ctx.clock() + cap_s, until)
    ctx.ga_opening = True
    ctx.say("step B: GA on the whole fleet for up to %.0f min, then the grids", cap_s / 60)
    ga: list = []
    reason = "cap"
    try:
        while live and ctx.clock() < stop and not ctx.over():
            if adopt_ga_best(ctx, live):
                reason = "seat"
                break
            if ga_beat_lowest(ctx, live):              # 044200Z islands: a lock kept island 0 waiting 4.5 minutes
                reason = "beat"
                break
            idle = free_workers(ctx)
            if idle is not None and idle > 0:
                ga = [f for f in ga if not f.done()] + await submit_ga(ctx, idle)
            pending = [f for f in ga if not f.done()]
            if not pending:
                if idle is not None and idle > 0:
                    reason = "idle"
                    break
                await asyncio.sleep(TICK)
                continue
            await asyncio.wait(pending, timeout=TICK, return_when=asyncio.FIRST_COMPLETED)
            pub = getattr(ctx, "publish", None)
            if callable(pub):
                pub()
    finally:
        ctx.ga_opening = False
    marks = getattr(ctx, "marks", None)
    if isinstance(marks, dict):
        marks["ga_open"] = ctx.elapsed()
    ctx.ga_open_reason = reason
    ctx.say("step B: GA opening ended by %s at t+%.1fm", reason, ctx.elapsed() / 60)
    if reason == "seat":
        ctx.ga_moved = ctx.clock()
        ctx.say("step B: a free GA beat the lowest climber; the grids stay a compass until the GA goes quiet")
    else:
        ctx.say("step B: the GA opening ended; the short grids start")


def note_seat(ctx, cl: Climber) -> None:
    """Record why this seat was opened: its gap to the best score, and its distance from every other seat."""
    top = ctx.best() if callable(getattr(ctx, "best", None)) else None
    best = top[1] if top is not None else cl.score
    peers = []
    for other in ctx.climbers:
        if other is cl:
            continue
        distance, kind = ctx.dist.between(cl.centre, other.centre)
        settings = sum(cl.centre.params.get(n) != other.centre.params.get(n) for n in ctx.space.names)
        peers.append({"seat": other.id, "distance": round(float(distance), 4), "kind": kind, "settings": settings})
    cl.entry_gap = round(best - cl.score, 4)
    cl.entry_peers = peers
    ctx.say("seat c%d %s %.4f gap %.4f from %s", cl.id, cl.origin, cl.score, cl.entry_gap,
            ", ".join(f"c{p['seat']} {p['kind']} {p['distance']:.3f} / {p['settings']} settings" for p in peers) or "no other seat")


def best_unclimbed(ctx) -> Optional[Tuple[Config, float, str]]:
    """The round's best, when no climber has stood on it and it beats every climb so far by climb.accept_min.

    A peaked config is not climbed again (B33); a GA score above every climber's peak is a new hill to climb."""
    top = ctx.best()
    if top is None:
        return None
    c, s = top
    accept = float(tune("climb.accept_min"))
    for cl in getattr(ctx, "climbers", ()):
        if any(d is not None and d.key == c.key for d in (getattr(cl, "start", None), getattr(cl, "centre", None))):
            return None
        if s < float(cl.score) + accept:
            return None
    return c, s, "round best"


async def occupy(ctx, ready: Sequence[Tuple[Config, float, str]], until: float) -> None:
    """Step B. A free worker takes the next unsent value of a climbing candidate.

    The first workers // candidates.workers_per_candidate candidates are loaded before anything is sent, so
    one of them cannot fill the fleet before the others exist. Step B opens with the whole fleet on the GA.
    Grids start once a free GA score beats the lowest climber, or climb.ga_open_min has passed. After that a
    free worker takes the next unsent value of a climbing candidate, and a hill hands its workers to the next
    Step A result. Workers still free after every ready candidate update has been sent run the GA.
    """
    queue = list(ready)
    from search.separate import enabled as separate_enabled
    if separate_enabled():
        from search.separate import survey
        await survey(ctx, queue, until)
        return
    n = n_seats(ctx)
    if piles_enabled():
        n = min(n, MAX_PILES - CANDIDATES)
    live: List[Climber] = []

    def spawn(item: Tuple[Config, float, str]) -> None:
        c, s, origin = item
        cl = Climber(ctx, len(ctx.climbers), c, s, origin=origin)
        cl.climb_until = until
        ctx.climbers.append(cl)
        note_seat(ctx, cl)
        if ga_hot(ctx) or cl.begin_grid(until):
            live.append(cl)
        else:
            cl.ended_s = round(ctx.elapsed(), 1)

    # Every opening candidate exists before any grid is built. budget() is the fleet divided by who is
    # climbing; building the first grid first would hand it the whole fleet.
    opening = queue[:n]
    del queue[:n]
    opened: List[Climber] = []
    for item in opening:
        c, s, origin = item
        cl = Climber(ctx, len(ctx.climbers), c, s, origin=origin)
        cl.climb_until = until
        ctx.climbers.append(cl)
        opened.append(cl)
    for cl in opened:
        note_seat(ctx, cl)
        live.append(cl)
    if live:
        ctx.say("step B: %d candidate(s); GA first, then grids; leftover workers keep breeding",
                len(live))
        await open_with_ga(ctx, live, until)
        if not piles_enabled():
            ensure_far_seat(ctx, live)
        else:
            for item in unseated(ctx, refresh_piles(ctx)):
                if len([cl for cl in live if cl.state == "climbing"]) >= n:
                    break
                spawn(item)
        started: List[Climber] = []
        for cl in live:
            if cl.state != "climbing":
                cl.ended_s = round(ctx.elapsed(), 1)
                continue
            if ga_hot(ctx) or cl.batch_live or cl.begin_grid(until):
                started.append(cl)
            else:
                cl.ended_s = round(ctx.elapsed(), 1)
        live = started
        if live:
            if leftover_breeds():
                ctx.say("step B: %d candidate(s); compass then a cube, and leftover cores keep breeding", len(live))
            else:
                ctx.say("step B: %d candidate(s); the GA keeps the fleet while it is still moving them, then one rung",
                        len(live))
    rr = 0
    ga: list = []
    # No climber is not the end of Step B: the GA and the scouts keep the fleet, and the round's best takes a seat
    # when it is a new hill (044200Z, 2026-10-16: the one seat peaked at t+29.6m and Step B returned with half
    # the round left).
    while ctx.clock() < until and not ctx.over():
        if any(cl.state == "climbing" and cl.claim_higher() for cl in live):
            continue
        if adopt_reference(ctx, live):
            continue
        if switch_to_candidate(ctx, live):
            continue
        if adopt_ga_best(ctx, live):
            continue
        if ensure_far_seat(ctx, live):
            continue
        hold = ga_hot(ctx) and not leftover_breeds()
        if hold:
            for cl in live:
                if cl.state != "climbing" or (cl.origin or "") == "step C":
                    continue
                cl.pending, cl.extra, cl.later, cl.leaps = [], [], [], []
                if not cl.gate_futs and not cl.inflight:
                    cl.batch_live = False
        else:
            for cl in live:
                if cl.state == "climbing" and not cl.batch_live and not cl.queued():
                    cl.begin_grid(until)
        idle = free_workers(ctx)
        guard = 0
        while idle is not None and idle > 0 and not hold:
            # Heads, then the planned line, then a centre pack. v9.18 stops there and breeds on the
            # cores that remain. Later tags walk further up the line until the fleet is full.
            heads = [cl for cl in live if cl.state == "climbing" and cl.pending]
            extras = [cl for cl in live if cl.state == "climbing" and cl.extra]
            choices = heads or extras
            if not choices:
                if leftover_breeds():
                    choices = [cl for cl in live if cl.state == "climbing" and cl.later]
                else:
                    climbing = [cl for cl in live if cl.state == "climbing" and cl.batch_live and cl.kind != "jumps"]
                    share = max(1, idle // len(climbing)) if climbing else 0
                    if any(cl.extend(share) for cl in climbing):
                        choices = [cl for cl in live if cl.state == "climbing" and cl.extra]
                    else:
                        choices = [cl for cl in live if cl.state == "climbing" and cl.later]
            if not choices:
                # Centre jumps (search/leaps.py) take their share before the GA breeds on what is left.
                choices = [cl for cl in live if cl.state == "climbing" and cl.leaps]
            if not choices:
                break
            rr += 1
            choices[rr % len(choices)].submit_one()
            idle = free_workers(ctx)
            guard += 1
            if guard > int(ctx.cap) + 5:
                break
        acted = False
        for cl in list(live):
            if cl.state != "climbing" or not (cl.gate_done() or cl.early_gain()):
                continue
            try:
                cl.close_wave()
            except Exception:  # noqa: BLE001 - one climber's bug ends that climber
                log.exception("climb: climber %d failed", cl.id)
                cl.state = "stopped"
                cl.batch_live = False
            if cl.state == "climbing" and not cl.batch_live and not hold:
                cl.begin_grid(until)
            acted = True
        for cl in list(live):
            if cl.state != "climbing":
                cl.ended_s = round(ctx.elapsed(), 1)
                cl.kind, cl.group, cl.inflight = "", [], []
        finished = [cl for cl in live if cl.state != "climbing"]
        live = [cl for cl in live if cl.state == "climbing"]
        late = late_for_seat(ctx, until)
        if late and (finished or len(live) < n) and not getattr(ctx, "_late_seat_said", False):
            ctx._late_seat_said = True
            ctx.say("step B: no new far seat with %.0f min left; the climbers still seated keep the workers",
                    (until - ctx.clock()) / 60)
        if finished and ctx.clock() < until and not ctx.over() and not late:
            for _ in finished:
                # A peaked shelf seat goes to the edge family first. The backup queue is this table.
                nxt = None if piles_enabled() else far_replacement(ctx)
                if nxt is None and queue and not piles_enabled():
                    nxt = _next_far(ctx, queue)
                if nxt is None and piles_enabled():
                    seats = unseated(ctx, refresh_piles(ctx))
                    nxt = seats[0] if seats else None
                if nxt is None:
                    break
                ctx.say("step B: seat %s %.4f (%s)", nxt[0].key[:8], nxt[1], nxt[2])
                spawn(nxt)
                acted = True
        while len(live) < n and ctx.clock() < until and not ctx.over() and not late:
            if piles_enabled():
                seats = unseated(ctx, refresh_piles(ctx))
                nxt = seats[0] if seats else None
            else:
                nxt = far_replacement(ctx)
            if nxt is None:
                break
            ctx.say("step B: seat %s %.4f (%s)", nxt[0].key[:8], nxt[1], nxt[2])
            spawn(nxt)
            acted = True
        if not live and ctx.clock() < until and not ctx.over():
            nxt = best_unclimbed(ctx)
            if nxt is not None:
                ctx.say("step B: seat %s %.4f (%s)", nxt[0].key[:8], nxt[1], nxt[2])
                spawn(nxt)
                acted = True
        if (not piles_enabled()) and len(live) < n and ctx.clock() < until and not ctx.over():
            pending_scouts = [f for f in getattr(ctx, "_scout_futs", []) if not f.done()]
            ctx._scout_futs = pending_scouts
            idle = free_workers(ctx)
            if (not pending_scouts and getattr(ctx, "_scout_waves", 0) < SCOUT_WAVES
                    and idle is not None and idle > 0):
                scouts = far_scouts(ctx, min(idle, 4 * (n - len(live))))
                ctx._scout_waves = getattr(ctx, "_scout_waves", 0) + (1 if scouts else SCOUT_WAVES)
                if scouts:
                    ctx._scout_futs = ctx.propose(scouts, "F", ["big"] * len(scouts), priority=1)
                    ctx.say("step B: %d scouts; nothing scored is far from the climbers", len(scouts))
        idle = free_workers(ctx)
        held = 0
        if (hold or leftover_breeds()) and idle is not None and idle > 0:
            heads_left = any(cl.state == "climbing" and cl.pending for cl in live)
            packs_left = any(cl.state == "climbing" and (cl.later or cl.leaps) for cl in live)
            if not heads_left and not packs_left:
                keep = head_reserve(live)
                spare = idle - keep
                held = min(idle, keep)
                if spare > 0:
                    ga = [f for f in ga if not f.done()] + await submit_ga(ctx, spare)
        note_reserve(ctx, held)
        if not live and not any(not f.done() for f in ga) and \
                not any(not f.done() for f in getattr(ctx, "_scout_futs", [])):
            break       # nothing climbs and nothing runs: a virtual clock would never move again
        pending = ([f for cl in live for f in cl.open_futs()] + [f for f in ga if not f.done()]
                   + [f for f in getattr(ctx, "_scout_futs", []) if not f.done()])
        if not pending:
            await asyncio.sleep(0)
            continue
        await asyncio.wait(pending, timeout=TICK, return_when=asyncio.FIRST_COMPLETED)
        ctx.publish()
    for cl in live:
        if cl.batch_live and ctx.clock() >= until:
            try:
                cl.close_wave()
            except Exception:  # noqa: BLE001
                log.exception("climb: climber %d failed", cl.id)
                cl.state = "stopped"
        if cl.state == "climbing":
            cl.state = "stopped"
        cl.ended_s = round(ctx.elapsed(), 1)
        cl.kind, cl.group, cl.inflight = "", [], []
        cl.batch_live = False
