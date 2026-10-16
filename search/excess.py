"""A second measure beside the score: core above the run's own core-vs-FP curve (2026-10-14). Owner: search.

Why. The published score blends core and false positives, so a better trade-off can sit behind a score dip: on
20261010T044200Z a climb on this measure (test/niche_climb) reached core 0.8723 at FP 1353 (the king: 0.8703 at
1232), +0.0040 above our curve, while its score went 82.7 -> 81.9 -> 82.2 on the way. A climb on the score alone never
crosses such a dip (Brockhoff et al. 2009, "On the Effects of Adding Objectives to Plateau Functions"; Knowles,
Watson and Corne 2001, multiobjectivization: a second measure removes a flat only when it points at the exit).

The measure. J(config) = core - curve(FP clamped to [band_low, band_high] x the anchor's FP)      (0.95 / 1.05)
                         - the score's own price of each % of FP above that band, in core.
  Tight on purpose (2026-10-14, test/conversion_sweep): with a 0.90-1.15 band and a token 0.0003 core per % above it,
  climbs settled at 1.2x our FP, +0.0040 above our curve but where the score charges FP heavily; cutting FP with the
  calling threshold, BQST and contamination took them only from 82.24 to 82.39 (our best 82.73) and halved the edge.
  The field's edge exists at our own FP too (82% of its top 40 within 0.97-1.03x our FP are above our curve, +0.0018),
  so J looks for more core at about the same FP, which is what turns into score.
  curve(x): the best core of every scored config at FP <= x, EXCEPT the lane's own results (they are what is measured
            against the rest of the search; with them in, nothing could ever be above the curve).
  anchor:   the FP of the run's best score.
  WARNING (found the hard way, test/region_hunt): the curve is dense only near the anchor. Comparing a config at its
  own FP far from it let climbs drift to FP 700-780 and look +0.01 to +0.08 "above" a thin curve; ranking every
  out-of-band config below the band (ignoring core) made climbs trade any core for FP. Clamping to the band edge, as
  here, is the version that held.

Where it is used.
  Excess lane (excess.lanes, default 0 = off; v9.18.2 with it off is v9.18.1): during Step B, starting from the best-J scored config at least
     excess.start_far apart in settings from the best config, it climbs on J (every setting at +-1, 2, 4, 8 rungs; the
     best single move or the top-3 / top-6 combination), then from the J peak it climbs on the score (converting the
     better trade-off into points). Everything it scores goes to the store, so the pick (highest score) sees it.
     The score stays the only thing that names the pick.
  Exact score ties are left to the plateau probe (search/plateau.py); J is not used there.
"""
from __future__ import annotations

import asyncio
import bisect
import logging
import math
import time
from typing import Dict, List, Optional, Sequence, Set, Tuple

from common.tuning import tune
from common.types import FULL, Config
from search.distance import setting_distance

log = logging.getLogger(__name__)
LANE_TAG = "X"
OFFSETS = (1, -1, 2, -2, 4, -4, 8, -8)


class Curve:
    """The run's core-vs-FP frontier (lane results excluded) and J against it."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.lane_keys: Set[str] = set()
        self.fps: List[float] = []
        self.cur: List[float] = []
        self.anchor: Optional[float] = None
        self._n, self._t = -1, 0.0

    def refresh(self, force: bool = False) -> None:
        latest = self.ctx.client.latest
        n = len(latest)
        if not force and (n == self._n or (self._n > 0 and n < 1.02 * self._n and time.time() - self._t < 30)):
            return
        pts = []
        for (key, fid), r in latest.items():
            if (fid != FULL or key in self.lane_keys or r.score is None or r.gate_passed is False
                    or r.core is None or r.region_fp is None):
                continue
            pts.append((float(r.region_fp), float(r.core)))
        pts.sort()
        fps, cur, b = [], [], -1.0
        for fp, c in pts:
            b = max(b, c)
            fps.append(fp)
            cur.append(b)
        self.fps, self.cur, self._n, self._t = fps, cur, n, time.time()
        top = self.ctx.client.best(FULL)
        self.anchor = float(top.region_fp) if top is not None and top.region_fp is not None else None

    def n_targets(self) -> float:
        a = getattr(self.ctx.client, "answer_set", None)
        n = getattr(a, "n_targets", None)
        return float(n) if n else 548.0

    def over_cost(self, hi: float) -> float:
        """Core per % of the anchor's FP that the score itself charges at the band's upper edge:
        score = 100 x (0.70 core + 0.30 exp(-FP / 8N)), so one FP costs (30 / 8N) exp(-FP / 8N) points and one point of
        core is worth 70."""
        n = self.n_targets()
        per_fp = 30.0 / (8.0 * n) * math.exp(-hi / (8.0 * n))
        return per_fp * self.anchor / 100.0 / 70.0 * float(tune("excess.over_scale"))

    def j(self, core: float, fp: float) -> Optional[float]:
        if self.anchor is None or not self.fps:
            return None
        lo, hi = float(tune("excess.band_low")) * self.anchor, float(tune("excess.band_high")) * self.anchor
        x = min(max(fp, lo), hi)
        i = bisect.bisect_right(self.fps, x)
        ref = self.cur[i - 1] if i else self.cur[0]
        over = 100.0 * max(0.0, fp - hi) / self.anchor
        return core - ref - self.over_cost(hi) * over

    def of(self, c: Config) -> Optional[float]:
        r = self.ctx.client.latest.get((c.key, FULL))
        if r is None or r.score is None or r.gate_passed is False or r.core is None or r.region_fp is None:
            return None
        return self.j(float(r.core), float(r.region_fp))


def curve(ctx) -> Curve:
    cv = getattr(ctx, "excess_curve", None)
    if cv is None:
        cv = ctx.excess_curve = Curve(ctx)
    cv.refresh()
    return cv


class Lane:
    """Climb on J, then on the score from the J peak."""

    def __init__(self, ctx, lid: int):
        self.ctx, self.id = ctx, lid
        self.centre: Optional[Config] = None
        self.phase, self.passes = "waiting", 0
        self.best_j: Optional[float] = None
        self.score: Optional[float] = None

    def _moves(self, centre: Config) -> Dict[str, Tuple[str, object, Config]]:
        ctx, out = self.ctx, {}
        for n in ctx.space.names:
            if not ctx.ladders.movable(n):
                continue
            spec = ctx.space.specs[n]
            if spec.allowed or spec.type == "bool":
                vals = [v for v in ctx.ladders.rungs.get(n, []) if v != centre.params.get(n)]
            else:
                vals = [ctx.ladders.shift(n, centre.params.get(n), k) for k in OFFSETS]
            for v in vals:
                if v is None or v == centre.params.get(n):
                    continue
                c = ctx.space._assemble({**centre.params, n: v})
                if c is not None and c.key not in out:
                    out[c.key] = (n, v, c)
        return out

    async def _score(self, configs: Sequence[Config], until: float) -> None:
        """Send what is neither scored nor already proposed, and wait for most of it."""
        ctx, send = self.ctx, []
        for c in configs:
            if (c.key, FULL) in ctx.client.latest:
                continue
            if ctx.fresh(c):
                send.append(c)
        if not send:
            return
        cv = curve(self.ctx)
        cv.lane_keys.update(c.key for c in send)
        futs = []
        for c in send:
            futs += ctx.propose([c], f"{LANE_TAG}{self.id}", ["small"], priority=0)
        await ctx.wait_share(futs, 0.95, min(until, ctx.clock() + ctx.T(20 * 60)))

    def _value(self, c: Config, on_j: bool) -> Optional[float]:
        if on_j:
            return curve(self.ctx).of(c)
        return self.ctx.score(c)

    async def climb(self, start: Config, on_j: bool, until: float) -> Config:
        ctx = self.ctx
        self.centre, self.phase = start, ("J" if on_j else "score")
        await self._score([start], until)
        cur = self._value(start, on_j)
        if cur is None:
            return start
        need = 1e-5 if on_j else float(tune("climb.accept_min"))
        while ctx.clock() < until and not ctx.over():
            moves = self._moves(self.centre)
            await self._score([t[2] for t in moves.values()], until)
            gains = []
            for k, (n, v, c) in moves.items():
                val = self._value(c, on_j)
                if val is not None and val > cur + need:
                    gains.append((val, n, v))
            gains.sort(key=lambda t: -t[0])
            picked: List[Tuple[float, Dict[str, object]]] = [(g[0], {g[1]: g[2]}) for g in gains[:1]]
            combos = []
            singles = []
            for g in gains:                                   # one per setting, best first
                if all(g[1] != s[1] for s in singles):
                    singles.append(g)
            for top in (3, 6):
                if len(singles) >= 2 and (top == 3 or len(singles) > 3):
                    step = {s[1]: s[2] for s in singles[:top]}
                    c = ctx.space._assemble({**self.centre.params, **step})
                    if c is not None:
                        combos.append((c, step))
            if combos:
                await self._score([c for c, _ in combos], until)
                for c, step in combos:
                    val = self._value(c, on_j)
                    if val is not None and val > cur + need:
                        picked.append((val, step))
            self.passes += 1
            if not picked:
                break
            val, step = max(picked, key=lambda t: t[0])
            nxt = ctx.space._assemble({**self.centre.params, **step})
            if nxt is None:
                break
            self.centre, cur = nxt, val
            if on_j:
                self.best_j = cur
            self.score = ctx.score(nxt)
            ctx.say("excess lane %d (%s): pass %d -> %s %.4f (score %s) by %s", self.id, self.phase, self.passes,
                    "J" if on_j else "score", cur, f"{self.score:.4f}" if self.score is not None else "?", step)
        return self.centre

    def row(self) -> dict:
        return {"lane": self.id, "phase": self.phase, "passes": self.passes, "best_j": self.best_j,
                "score": self.score, "centre": self.centre.key if self.centre is not None else None}


def _start(ctx, taken: Sequence[Config]) -> Optional[Config]:
    """The best-J scored config at least excess.start_far apart in settings from the best config and from other lanes."""
    top = ctx.best()
    if top is None:
        return None
    best_c = top[0]
    far = float(tune("excess.start_far"))
    cv = curve(ctx)
    ranked = []
    for c, s in ctx.scored():
        if c.key in cv.lane_keys:
            continue
        if setting_distance(ctx.space, c, best_c) < far or any(setting_distance(ctx.space, c, t) < far for t in taken):
            continue
        j = cv.of(c)
        if j is not None:
            ranked.append((j, c))
    if not ranked:
        return None
    return max(ranked, key=lambda t: t[0])[1]


async def run_lanes(ctx, until: float) -> None:
    """Step B's excess lanes: each climbs on J, then converts its J peak into score."""
    n = int(tune("excess.lanes"))
    if n <= 0:
        return
    ctx.excess_lanes = lanes = [Lane(ctx, i) for i in range(n)]
    starts: List[Config] = []
    for lane in lanes:
        s = _start(ctx, starts)
        if s is None:
            ctx.say("excess lane %d: no scored config far enough from the best; not started", lane.id)
            continue
        starts.append(s)
        ctx.say("excess lane %d: starts at %s (J %+.4f, score %.4f)", lane.id, s.key[:8],
                curve(ctx).of(s) or 0.0, ctx.score(s) or 0.0)

    async def one(lane: Lane, s: Config) -> None:
        try:
            peak = await lane.climb(s, on_j=True, until=until)
            ctx.say("excess lane %d: J peak %+.4f at %s after %d passes; climbing its score", lane.id,
                    lane.best_j if lane.best_j is not None else 0.0, peak.key[:8], lane.passes)
            if int(tune("excess.convert")):
                await lane.climb(peak, on_j=False, until=until)
            lane.phase = "done"
        except asyncio.CancelledError:
            raise
        except Exception:  # noqa: BLE001 - a lane failing must not end Step B
            log.exception("engine: excess lane %d failed", lane.id)
            lane.phase = "failed"

    await asyncio.gather(*(one(lane, s) for lane, s in zip(lanes, starts)))
