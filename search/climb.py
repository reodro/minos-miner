"""Step B: climb one candidate to the peak of its hill (search/DESIGN.md, Step B).

One iteration, from the centre:

    grid        every combination of climb.group_size settings at -d, 0, +d (3^5 = 243 points, fewer when a setting
                sits at the end of its ladder or the BAM collapses it). d is the setting's resolution: one rung of its
                ladder (search/steps.py), whatever its range.
                    odd iterations   the read-pass settings no batch can share (FULL_RUN_PARAMS)
                    even iterations  5 of the 7 assembly and activity settings, in turn
    riders      every grid point carries the cheap settings (genotyping, PairHMM) one at a time at -d and +d, in the
                same GATK run; the centre's own run carries every pair of them
    move        the best point of the iteration becomes the centre. Nothing is combined.
    growth      a grid setting whose move won is offered 2d next time in that direction (and d back), then 3d,
                4d... while it keeps winning; when it does not, it goes back to d. Example: A+1 wins the grid over
                (A, B, C); the next grid is around (A+1, B, C) with A at -1, 0, +2 and B, C at -1, 0, +1

An iteration that gains less than climb.accept_min has failed. When every expensive setting (read-pass and all 7
assembly/activity ones, which rotate 5 at a time) has been in a failed grid since the last gain, the climb stalls:
climb.random_jumps random local jumps (2-4 settings, 2-5 rungs
each). A jump that gains climb.accept_min becomes the centre and the grid resumes; otherwise the centre is the peak.
Far points are Step C's (search/far.py).
"""
from __future__ import annotations

import itertools
from typing import Dict, List, Optional, Sequence, Tuple

from common.tuning import tune
from common.types import Config
from search.space import ACTIVITY_PARAMS, ASSEMBLY_PARAMS, FULL_RUN_PARAMS, GENOTYPING_PARAMS, PAIRHMM_PARAMS


class Climber:
    def __init__(self, ctx, cid: int, start: Config, score: float, origin: str = "step A"):
        self.ctx, self.id, self.origin = ctx, cid, origin
        self.start, self.centre, self.score = start, start, score
        self.opened = score
        self.first: Optional[Tuple[float, str]] = None        # (score, key) after the first iteration
        self.iteration = 0
        self.turn = 0                                          # odd turns read-pass, even turns assembly/activity
        self.growth: Dict[str, Tuple[int, int]] = {}           # setting -> (direction, rungs) offered next time
        self.flat: set = set()                                 # settings gridded without a gain since the last one
        self.rotation = 0
        self.state = "climbing"                                # climbing | peak | joined | stopped | landed
        self.evaluated = 0
        self.waiting: Dict[str, Config] = {}                   # sent, not yet seen scored by a decision
        self.waves: List[Tuple[str, int, int, float]] = []    # (label, heads, configs, minutes) per wave
        self.history: List[Tuple[float, float, str]] = []     # (seconds into the search, score, why)
        self.born_s = round(ctx.elapsed(), 1)
        self.ended_s: Optional[float] = None
        self.parent: Optional[int] = None                      # a Step C landing: the climber whose peak it came from
        self.kind = ""                                         # the grid in flight: read | assembly | jumps
        self.group: List[str] = []
        self.inflight: list = []                               # the heads of the wave in flight

        names = [n for n in ctx.space.names if ctx.ladders.movable(n)]
        self.read_names = [n for n in FULL_RUN_PARAMS if n in names]
        self.assembly_names = [n for n in ASSEMBLY_PARAMS + ACTIVITY_PARAMS if n in names]
        self.cheap_names = [n for n in GENOTYPING_PARAMS + PAIRHMM_PARAMS if n in names]
        self.all_names = names

    # ---- the climb ------------------------------------------------------------------------------------------------
    async def run(self, until: float, iterations: Optional[int] = None, jumps: bool = True) -> None:
        """Climb until the peak, `until`, or `iterations` more grids (a short climb: no random jumps)."""
        ctx = self.ctx
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
            if gain < float(tune("climb.accept_min")):
                self.flat.update(group)
            else:
                self.flat.clear()
            ctx.say("climber %d it %d (%s: %s): %+.4f -> %.4f; waves %s", self.id, self.iteration, kind,
                    ",".join(group), gain, self.score,
                    ", ".join(f"{w[1]}/{w[2]} in {w[3]:.1f}m" for w in self.waves[-2:]))
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

    def next_group(self) -> Tuple[str, List[str]]:
        size = int(tune("climb.group_size"))
        self.turn += 1
        if self.turn % 2 == 1:
            return "read", self.read_names[:size]
        n = len(self.assembly_names)
        if n <= size:
            return "assembly", list(self.assembly_names)
        group = [self.assembly_names[(self.rotation + i) % n] for i in range(size)]
        self.rotation = (self.rotation + size) % n
        return "assembly", group

    def offsets(self, name: str) -> List[int]:
        d, k = self.growth.get(name, (1, 1))
        return [0, d * k, -d]

    def values(self, name: str, base) -> list:
        out = []
        for off in self.offsets(name):
            v = base if off == 0 else self.ctx.ladders.shift(name, base, off)
            if v is not None and v not in out:
                out.append(v)
        return out

    def with_(self, base: Config, changes: Dict[str, object]) -> Optional[Config]:
        return self.ctx.space._assemble({**base.params, **changes})

    def riders(self, head: Config, pairs: bool = False) -> List[Config]:
        lad, out = self.ctx.ladders, []
        singles: List[Tuple[str, object]] = []
        for n in self.cheap_names:
            for off in (1, -1):
                v = lad.shift(n, head.params.get(n), off)
                if v is not None:
                    singles.append((n, v))
        for n, v in singles:
            out.append(self.with_(head, {n: v}))
        if pairs:
            for (n1, v1), (n2, v2) in itertools.combinations(singles, 2):
                if n1 != n2:
                    out.append(self.with_(head, {n1: v1, n2: v2}))
        return [c for c in out if c is not None]

    async def send(self, heads: Sequence[Config], riders_of, label: str, until: float) -> List[Tuple[Config, float]]:
        """Each head in one GATK run with its riders. Waits for climb.quorum of the heads, and of the riders too when
        climb.wait_riders is 1. Returns (config, score) of this wave's scored configs plus any of this climber's
        earlier configs that came back since the last decision (a late rider)."""
        ctx = self.ctx
        head_futs, all_futs, seen = [], [], set()
        t0 = ctx.clock()
        for h in heads:
            if h.key in seen:
                continue
            seen.add(h.key)
            batch = [h] + [r for r in riders_of(h) if r.key not in seen]
            seen.update(r.key for r in batch)
            fs = ctx.propose(batch, f"B{self.id}", ["small"] + ["big"] * (len(batch) - 1))
            head_futs.append(fs[0])
            all_futs += fs
            for c in batch:
                self.waiting[c.key] = c
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
        for combo in itertools.product(*axes):
            c = self.with_(centre, dict(zip(group, combo)))
            if c is not None and c.key not in keys:
                keys.add(c.key)
                heads.append(c)
        heads.insert(0, centre)
        results = await self.send(heads, lambda h: self.riders(h, pairs=h.key == centre.key), kind, until)
        self.refresh()
        best = max(results, key=lambda t: t[1], default=None)
        old = self.centre
        if best is not None and best[1] > self.score:
            self.grow(group, old, best[0])
            self.move(best[0], best[1], f"{kind} grid")
        else:
            for n in group:
                self.growth.pop(n, None)

    def grow(self, group: Sequence[str], old: Config, new: Config) -> None:
        lad = self.ctx.ladders
        for n in group:
            r = lad.rungs_moved(n, old.params.get(n), new.params.get(n))
            if r == 0 or self.ctx.space.specs[n].type == "bool":
                self.growth.pop(n, None)
                continue
            d = 1 if r > 0 else -1
            pd, pk = self.growth.get(n, (d, 1))
            self.growth[n] = (d, max(pk + 1 if pd == d else 2, int(round(abs(r))) + 1))

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
            if c is not None and c.key not in keys and ctx.fresh(c):
                keys.add(c.key)
                heads.append(c)
        if not heads or ctx.clock() >= until:
            self.state = "peak"
            return
        self.kind, self.group = "jumps", []
        results = await self.send(heads, lambda h: [], "jumps", until)
        self.refresh()
        best = max(results, key=lambda t: t[1], default=None)
        if best is not None and best[1] - self.score >= float(tune("climb.accept_min")):
            ctx.say("climber %d: a jump gains %+.4f; the grid resumes", self.id, best[1] - self.score)
            self.growth.clear()
            self.move(best[0], best[1], "jump")
            self.flat.clear()
        else:
            self.state = "peak"
            ctx.say("climber %d: peak %.4f %s (%d jumps, none gains %.3f)", self.id, self.score, self.centre.key[:8],
                    len(heads), float(tune("climb.accept_min")))

    # ---- bookkeeping ----------------------------------------------------------------------------------------------
    def refresh(self) -> None:
        """The centre's score as the client has it now (an answer-set update rescores everything)."""
        s = self.ctx.score(self.centre)
        if s is not None:
            self.score = s

    def move(self, c: Config, s: float, why: str) -> None:
        self.centre, self.score = c, s
        self.history.append((round(self.ctx.elapsed(), 1), round(s, 4), why))
        other = self.ctx.owner(c, self)
        if other is not None:
            self.state = "joined"
            self.ctx.say("climber %d reached climber %d's centre %s; it stops", self.id, other.id, c.key[:8])

    def resume(self, c: Config, s: float, why: str) -> None:
        """Climb again from `c` (Step C found it): fresh grids, no growth carried over."""
        self.ended_s = None
        self.growth.clear()
        self.flat.clear()
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
                "waves": self.waves[-10:], "history": self.history[-20:]}
