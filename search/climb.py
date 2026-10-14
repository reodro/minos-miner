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
climb.random_jumps random local jumps (2-4 settings, 2-5 rungs
each). A jump that gains climb.accept_min becomes the centre and the grid resumes; otherwise the centre is the peak.
Far points are Step C's (search/far.py).
"""
from __future__ import annotations

import asyncio
import itertools
import logging
import math
import time
from typing import Dict, List, Optional, Sequence, Tuple

log = logging.getLogger(__name__)
TICK = 2.0
CLIMB_PRIORITY = 1          # a candidate update is taken before a GA child (priority 0) waiting in the queue

from common.tuning import tune
from common.types import Config
from search.candidates import SEAT_GAP
from search.engine_tag import ELITE_LEFTOVER, VERSION
from search.piles import CANDIDATES, GAP, MAX_PILES, RADIUS, already_peaked, best_member, centre_allowed, claim_peak, enabled as piles_enabled, follow, mating_pool, own_pile, refresh as refresh_piles, seal_peak, unseated
from search.lines import densest, plan, single, tried
from search.space import ACTIVITY_PARAMS, ASSEMBLY_PARAMS, FULL_RUN_PARAMS, GENOTYPING_PARAMS, PAIRHMM_PARAMS

# A ladder this long is walked 1, 2, 3, then 5 rungs at a time. Shorter ladders still grow by one rung per win:
# a 5-rung jump would leave them.
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


class Climber:
    def __init__(self, ctx, cid: int, start: Config, score: float, origin: str = "step A"):
        self.ctx, self.id, self.origin = ctx, cid, origin
        self.start, self.centre, self.score = start, start, score
        self.found = score                                 # highest score this climber's own runs have returned
        self.opened = score
        self.first: Optional[Tuple[float, str]] = None        # (score, key) after the first iteration
        self.iteration = 0
        self.turn = 0                                          # odd turns read-pass, unless that grid is already dead
        self.growth: Dict[str, Tuple[int, int]] = {}           # setting -> (direction, rungs) offered next time
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
        return [0, d * k, -d]

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
        return bool(self.pending or self.later or self.extra)

    def submit_one(self) -> None:
        """One unsent head, then one more line point, then one packed centre batch."""
        if self.pending:
            self._submit_head()
        elif self.extra:
            self._submit_extra()
        elif self.later:
            self._submit_pack()

    def _submit_extra(self) -> None:
        """One more line point, for a worker the planned grid did not fill. It does not hold the wave open."""
        ctx = self.ctx
        while self.extra and (not ctx.fresh(self.extra[0]) or held_by_climber(ctx, self.extra[0], me=self) is not None):
            self.extra.pop(0)
        if not self.extra:
            return
        h = self.extra.pop(0)
        fs = ctx.propose([h], f"B{self.id}", ["small"], priority=CLIMB_PRIORITY)
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
        riders = [] if h.key in self.line_keys or h.key in probes else self.geno_riders(h)
        batch = [h] + [r for r in riders if r is not None and r.key != h.key and ctx.fresh(r)]
        fs = ctx.propose(batch, f"B{self.id}", ["small"] + ["big"] * (len(batch) - 1), priority=CLIMB_PRIORITY)
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
        self.extra_futs.extend(fs)
        for c in pack:
            self._watch(c)
        self.evaluated += len(fs)
        self._wave_n += len(fs)

    def gate_done(self) -> bool:
        """The planned heads have all been sent and climb.quorum of them are back."""
        if not self.batch_live or self.pending:
            return False
        if not self.gate_futs:
            return True
        need = math.ceil(float(tune("climb.quorum")) * len(self.gate_futs))
        return sum(1 for f in self.gate_futs if f.done()) >= need

    def open_futs(self) -> list:
        return [f for f in self.gate_futs + self.extra_futs if not f.done()]

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
        self.gate_futs, self.extra_futs = [], []
        self.batch_live = True
        self._wave_t0 = ctx.clock()
        self._wave_n = 0
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
        ctx = self.ctx
        if not self.jumps_enabled or ctx.clock() >= until or ctx.over():
            return False
        n = int(tune("climb.random_jumps"))
        heads, keys = [], {self.centre.key}
        for _ in range(20 * n):
            if len(heads) >= n or len(self.all_names) < 2:
                break
            names = ctx.rng.sample(self.all_names, min(len(self.all_names), ctx.rng.randint(2, 4)))
            changes = {}
            for name in names:
                if ctx.space.specs[name].type == "bool":
                    changes[name] = not self.centre.params.get(name)
                    continue
                k = ctx.rng.randint(2, 5) * ctx.rng.choice((1, -1))
                v = ctx.ladders.shift(name, self.centre.params.get(name), k)
                if v is None:
                    v = ctx.ladders.shift(name, self.centre.params.get(name), -k)
                if v is not None:
                    changes[name] = v
            c = self.with_(self.centre, changes) if changes else None
            if c is not None and c.key not in keys and ctx.fresh(c) and held_by_climber(ctx, c, me=self) is None:
                keys.add(c.key)
                heads.append(c)
        if not heads:
            self._peaked()
            return False
        self.kind, self.group = "jumps", []
        self.wave_before = self.score
        self.wave_centre = self.centre
        self.phase = "jumps"
        self.wave_before = self.score
        self.wave_centre = self.centre
        self.line_keys = set()
        self.pending, self.extra = heads, []
        self.gate_futs, self.extra_futs = [], []
        self.batch_live = True
        self._wave_t0 = ctx.clock()
        self._wave_n = 0
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
        results = []
        for k, c in list(self.waiting.items()):
            if c is None or not hasattr(c, "params"):
                del self.waiting[k]
                continue
            s = ctx.score(c)
            if s is not None:
                results.append((c, s))
                del self.waiting[k]
        self.waves.append((self.kind, len(self.gate_futs), self._wave_n, round((ctx.clock() - self._wave_t0) / 60, 2)))
        self.refresh()
        if self.kind == "jumps":
            self.batch_live = False
            self.pending, self.extra = [], []
            best = max(results, key=lambda t: t[1], default=None)
            if (best is not None and best[1] - self.score >= float(tune("climb.accept_min"))
                    and centre_allowed(self, best[0])):
                ctx.say("climber %d: a jump gains %+.4f; the grid resumes", self.id, best[1] - self.score)
                self.growth.clear()
                self.move(best[0], best[1], "jump")
                self.flat.clear()
            else:
                self._peaked()
                ctx.say("climber %d: peak %.4f %s (%d jumps, none gains %.3f)", self.id, self.score,
                        self.centre.key[:8], len(self.gate_futs), float(tune("climb.accept_min")))
            return
        self.wave_results.extend(results)
        if self.phase == "compass":
            centre = self.wave_centre or self.centre
            live = self.live_axes(centre, self.group, self.wave_results, self.wave_before)
            if len(live) >= 2 and self.begin_confirm(live):
                return
        self._finish_grid_wave()

    def _finish_grid_wave(self) -> None:
        ctx = self.ctx
        kind = self.kind
        old = self.wave_centre
        best = max(self.wave_results, key=lambda t: t[1], default=None)
        if best is not None and best[1] > self.score and follow(ctx, self, best[0], best[1]):
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
            self.move(best[0], best[1], f"{kind} {'line' if on_line else 'grid'}")
        else:
            for name in self.group:
                self.growth.pop(name, None)
        if self.first is None:
            self.first = (self.score, self.centre.key)
        gain = self.score - self.wave_before
        self.note_gain(self.group, gain)
        ctx.say("climber %d it %d (%s: %s): %+.4f -> %.4f; waves %s", self.id, self.iteration, kind,
                ",".join(self.group), gain, self.score,
                ", ".join(f"{w[1]}/{w[2]} in {w[3]:.1f}m" for w in self.waves[-2:]))
        self.phase = ""
        self.wave_results = []
        self.batch_live = False
        self.pending, self.extra, self.later = [], [], []
        self.gate_futs, self.extra_futs = [], []
        if self.state == "climbing" and self.trailing():
            self.state = "stopped"
            ctx.say("climber %d is %.4f behind the best seat after two flat grids; the seat moves on",
                    self.id, self.gap_to_best())
        elif self.state == "climbing" and self.stalled():
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
            if self.gate_done():
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
        """The stall: random local jumps around the centre. One gaining climb.accept_min resumes the grid."""
        ctx = self.ctx
        n = int(tune("climb.random_jumps"))
        heads, keys = [], {self.centre.key}
        for _ in range(20 * n):
            if len(heads) >= n or len(self.all_names) < 2:
                break
            names = ctx.rng.sample(self.all_names, min(len(self.all_names), ctx.rng.randint(2, 4)))
            changes = {}
            for name in names:
                if ctx.space.specs[name].type == "bool":
                    changes[name] = not self.centre.params.get(name)
                    continue
                k = ctx.rng.randint(2, 5) * ctx.rng.choice((1, -1))
                v = ctx.ladders.shift(name, self.centre.params.get(name), k)
                if v is None:
                    v = ctx.ladders.shift(name, self.centre.params.get(name), -k)
                if v is not None:
                    changes[name] = v
            c = self.with_(self.centre, changes) if changes else None
            if c is not None and c.key not in keys and ctx.fresh(c) and held_by_climber(ctx, c, me=self) is None:
                keys.add(c.key)
                heads.append(c)
        if not heads or ctx.clock() >= until:
            self._peaked()
            return
        self.kind, self.group = "jumps", []
        results = await self.send(heads, lambda h: [], "jumps", until)
        self.refresh()
        best = max(results, key=lambda t: t[1], default=None)
        if (best is not None and best[1] - self.score >= float(tune("climb.accept_min"))
                and centre_allowed(self, best[0])):
            ctx.say("climber %d: a jump gains %+.4f; the grid resumes", self.id, best[1] - self.score)
            self.growth.clear()
            self.move(best[0], best[1], "jump")
            self.flat.clear()
        else:
            self._peaked()
            ctx.say("climber %d: peak %.4f %s (%d jumps, none gains %.3f)", self.id, self.score, self.centre.key[:8],
                    len(heads), float(tune("climb.accept_min")))

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

    def claim_higher(self) -> bool:
        """A run this climber sent that has come back above the centre takes the centre.

        A live grid keeps the centre until the wave finishes. The wave writes the setting it found.
        """
        if self.batch_live:
            return False
        best = None
        for c in list(self.owed.values()):
            s = self.ctx.score(c)
            if s is None:
                continue
            if s - self.score < float(tune("climb.accept_min")):
                self.owed.pop(c.key, None)
                continue
            if best is None or s > best[1]:
                best = (c, s)
        if best is None or not follow(self.ctx, self, best[0], best[1]):
            return False
        self.ctx.say("climber %d takes its returned run %s %.4f (centre was %.4f)",
                     self.id, best[0].key[:8], best[1], self.score)
        self.owed.pop(best[0].key, None)
        self.waiting.pop(best[0].key, None)
        self.pending, self.extra, self.later = [], [], []
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
        self.ctx.say("climber %d restarts from the best %s %.4f (centre was %.4f)",
                     self.id, config.key[:8], score, self.score)
        self.pending, self.extra, self.later = [], [], []
        self.waiting.clear()
        self.gate_futs, self.extra_futs, self.inflight = [], [], []
        self.batch_live = False
        self.growth.clear()
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
        was = self.best_found()
        self.ctx.say("climber %d takes %s %s %.4f (its best was %.4f, centre %s)",
                     self.id, why, c.key[:8], s, was, self.centre.key[:8])
        self.pending, self.extra, self.later = [], [], []
        self.waiting.clear()
        self.gate_futs, self.extra_futs, self.inflight = [], [], []
        self.batch_live = False
        self.growth.clear()
        self.flat.clear()
        self.kind, self.group = "", []
        self.resume(c, s, "ga")
        if self.state == "climbing" and not getattr(self.ctx, "ga_opening", False):
            self.begin_grid(self.climb_until)

    def resume(self, c: Config, s: float, why: str) -> None:
        """Climb again from `c` (Step C found it): fresh grids, no growth carried over."""
        self.ended_s = None
        self.growth.clear()
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


def _probe_side(ctx, centre: Config, name: str, base: float, known: Dict[int, float], direction: int):
    """_probe_value on one side only: the nearest untried rung in `direction`, None when that side is closed."""
    current = centre.params.get(name)
    k = direction
    while True:
        v = ctx.ladders.shift(name, current, k)
        if v is None:
            return None
        i = int(round(ctx.ladders.index(name, v)))
        if i not in known:
            return v
        if base - known[i] >= PROBE_LOSS:
            return None
        k += direction


def next_probes(ctx, centre: Config, base: float) -> List[Config]:
    """One untried rung of each setting the plateau holds at a single value (two when search/plateau.py knows the
    setting's interval: the rung toward it first, the opposite rung after every first rung).

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
    # search/plateau.py: a locked setting whose interval lies to one side sends the rung toward it first; the
    # opposite rung goes after every first rung (sent when a worker is free). No interval: the nearest untried rung.
    from search import plateau
    iv = plateau.intervals(ctx) if plateau.enabled() else {}
    first, second = [], []
    for n in locked:
        d = plateau.direction(ctx, centre, n, iv) if iv else 0
        if d == 0:
            picks = [(first, _probe_value(ctx, centre, n, base, known[n]))]
        else:
            toward = _probe_side(ctx, centre, n, base, known[n], d)
            away = _probe_side(ctx, centre, n, base, known[n], -d)
            picks = [(first, toward), (second, away)] if toward is not None else [(first, away)]
        for tier, value in picks:
            if value is None:
                continue
            c = ctx.space._assemble({**centre.params, n: value})
            if c is None or c.key in seen or not ctx.fresh(c):
                continue
            seen.add(c.key)
            tier.append(c)
    return first + second


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
    from search.plateau import centre_of_flat
    probes = next_probes(ctx, centre_of_flat(ctx, base, owner.centre), base)   # the tied config nearest the intervals
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
    try:
        children = step.breed(gen, n, why="GA" if gen < 3 else "GA cross", accept=accept,
                              elite=bool(ELITE_LEFTOVER) and gen < 3)
        if hasattr(children, "__await__"):
            children = await children
    finally:
        step._bound_pool = prev_bound
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


def is_far_origin(cl: Climber) -> bool:
    """True when this seat was opened as a far hill (B27) and must stay off the best hill (B30)."""
    return (cl.origin or "").lower().startswith("far")


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
    if not piles_enabled() or not getattr(ctx, "_piles", None) or not ga_hot(ctx):
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
    if not piles_enabled() or not getattr(ctx, "_piles", None):
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
        if not centre_allowed(taker, c):
            continue
        other = held_by_climber(ctx, c, me=taker)
        if other is not None and (other.origin or "") != piles[i].label(i):
            if c.key not in ignored:
                ignored.add(c.key)
                ctx.say("GA best %s is the same or similar to climber %d; it stays there", c.key[:8], other.id)
            continue
        taker.adopt(c, s)
        return True
    return False


def _adopt_ga_best_legacy(ctx, live: Sequence[Climber]) -> bool:
    """A near seat takes a near-hill GA; a far seat takes only a GA that is far from every near centre.

    Best means the highest score that climber's own runs have returned. Step B checks this on every pass,
    so the handoff does not wait for a grid to end. A GA child already sitting on a climber is skipped.
    On v9.17 the far seat was the lowest climber, so every free GA — including the best hill — landed on it.
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
    for c, s in ranked_ga(ctx):
        if c.key in seated:
            continue
        if near_hill_ga(ctx, c, near):
            pool = [cl for cl in near if cl.best_found() < s]
        else:
            pool = [cl for cl in far if cl.best_found() < s] or [cl for cl in near if cl.best_found() < s]
        if not pool:
            continue
        taker = min(pool, key=lambda cl: (cl.best_found(), cl.id))
        other = held_by_climber(ctx, c, me=taker)
        if other is not None:
            if c.key not in ignored:
                ignored.add(c.key)
                ctx.say("GA best %s is the same or similar to climber %d; it stays there", c.key[:8], other.id)
            continue
        taker.adopt(c, s)
        return True
    return False


def far_replacement(ctx) -> Optional[Tuple[Config, float, str]]:
    """The best scored config that is a different hill from every seated climber.

    The Step A backup list is fixed at the start of the climb. A GA result, or a scout, can be a better far
    seat and is not on that list. Callset distance when both configs have a sketch, otherwise the mean
    setting difference. A config that differs in fewer than four settings from a centre, or more than REPLACEMENT_FLOOR below the
    best, is left out. The scan is reused until a new score arrives or a centre moves.
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
            if c.key not in seated and s >= best - REPLACEMENT_FLOOR]
    pool.sort(key=lambda t: -t[0])
    hit = None
    for s, c in pool:
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

    A GA child changes two settings by one rung and stays inside the seat distance. A scout moves
    SCOUT_AXES settings by about a quarter of each ladder, which is enough for the mean setting
    difference to clear that distance, and keeps the rest of a config that has already scored.
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
    out, seen = [], set()
    tries = 0
    while len(out) < n and tries < 20 * n:
        tries += 1
        base = starts[tries % len(starts)]
        params = dict(base.params)
        for name in ctx.rng.sample(names, min(SCOUT_AXES, len(names))):
            rungs = ctx.ladders.rungs.get(name) or []
            span = max(4, len(rungs) // 4)
            step = ctx.rng.randint(max(2, span // 2), span) * ctx.rng.choice((1, -1))
            v = ctx.ladders.shift(name, params.get(name), step)
            if v is None:
                v = ctx.ladders.shift(name, params.get(name), -step)
            if v is not None:
                params[name] = v
        c = ctx.space._assemble(params)
        if c is None or c.key in seen or not ctx.fresh(c):
            continue
        if held_by_climber(ctx, c) is not None:
            continue
        if dist is not None and any(not dist.far(c, centre, mc, ms) for centre, _ in centres):
            continue
        seen.add(c.key)
        out.append(c)
    return out


def ensure_far_seat(ctx, live: List[Climber]) -> bool:
    """If every climber is the same hill as the best, the lowest takes the best far scored config."""
    climbing = [cl for cl in live if cl.state == "climbing" and getattr(cl, "centre", None) is not None]
    if len(climbing) < 2:
        return False
    best = max(climbing, key=lambda cl: (cl.score, -cl.id))
    if any(not same_hill(ctx, cl.centre, best.centre) for cl in climbing if cl.id != best.id):
        return False
    nxt = far_replacement(ctx)
    if nxt is None:
        return False
    taker = min((cl for cl in climbing if cl.id != best.id), key=lambda cl: (cl.best_found(), cl.id))
    ctx.say("step B: far seat %s %.4f (%s) on climber %d", nxt[0].key[:8], nxt[1], nxt[2], taker.id)
    taker.adopt(nxt[0], nxt[1])
    taker.origin = f"far seat ({nxt[2]})"
    return True


def _next_far(ctx, queue: List[Tuple[Config, float, str]]) -> Optional[Tuple[Config, float, str]]:
    mc, ms = float(tune("candidates.min_distance")), float(tune("candidates.min_setting_distance"))
    while queue:
        c, s, origin = queue.pop(0)
        if all(c.key not in (d.start.key, d.centre.key) and ctx.dist.far(c, d.centre, mc, ms) for d in ctx.climbers):
            return c, s, origin
    return None


async def open_with_ga(ctx, live: List[Climber], until: float) -> None:
    """The start of Step B: the whole fleet breeds until a free GA beats the lowest climber, or the cap."""
    if not live or ctx.over() or ctx.clock() >= until:
        return
    if free_workers(ctx) is None:
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
    n = max(1, int(getattr(ctx, "cap", 1)) // int(tune("candidates.workers_per_candidate")))
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
    while live and ctx.clock() < until and not ctx.over():
        if any(cl.state == "climbing" and cl.claim_higher() for cl in live):
            continue
        if adopt_reference(ctx, live):
            continue
        if switch_to_candidate(ctx, live):
            continue
        if adopt_ga_best(ctx, live):
            continue
        hold = ga_hot(ctx) and not leftover_breeds()
        if hold:
            for cl in live:
                if cl.state != "climbing" or (cl.origin or "") == "step C":
                    continue
                cl.pending, cl.extra, cl.later = [], [], []
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
                break
            rr += 1
            choices[rr % len(choices)].submit_one()
            idle = free_workers(ctx)
            guard += 1
            if guard > int(ctx.cap) + 5:
                break
        acted = False
        for cl in list(live):
            if cl.state != "climbing" or not cl.gate_done():
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
        if finished and ctx.clock() < until and not ctx.over():
            for _ in finished:
                nxt = _next_far(ctx, queue) if queue and not piles_enabled() else None
                if nxt is None and piles_enabled():
                    seats = unseated(ctx, refresh_piles(ctx))
                    nxt = seats[0] if seats else None
                if nxt is None and not piles_enabled():
                    nxt = far_replacement(ctx)
                if nxt is None:
                    break
                ctx.say("step B: seat %s %.4f (%s)", nxt[0].key[:8], nxt[1], nxt[2])
                spawn(nxt)
                acted = True
        while len(live) < n and ctx.clock() < until and not ctx.over():
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
        if (not piles_enabled()) and live and len(live) < n and ctx.clock() < until and not ctx.over():
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
        if (hold or leftover_breeds()) and idle is not None and idle > 0:
            heads_left = any(cl.state == "climbing" and cl.pending for cl in live)
            packs_left = any(cl.state == "climbing" and cl.later for cl in live)
            if not heads_left and not packs_left:
                ga = [f for f in ga if not f.done()] + await submit_ga(ctx, idle)
        if not live:
            break
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
