"""Search engine v12: trust regions searched at once, picked by asynchronous Thompson sampling (search/v12/DESIGN.md).

v10 committed to one hill by minute ~15-20; v11 grouped evaluated configs into hills and its merge rule collapsed them
again. v12 does not group anything: a hill is a search region (TuRBO, Eriksson 2019; CASMOPOLITAN, Wan 2021; Bounce,
Papenmeier 2023), and its size is learned from its own successes and failures.

    opening     v10's (defaults, official seed, portfolio share, Latin hypercube; cheap riders).
    screening   successive halving: tr.sh_regions (16) regions open at the opening's results far apart in calls (going
                well below the top), share the picks in equal turns, and every tr.sh_step_min the better half by best +
                gain since opened is kept, until tr.regions remain. tr.sh_far_starts of the starts are far in settings
                from the opening's top (hills no engine so far converged to); each halving also keeps the
                tr.sh_indel_keep regions with the best indel part at the reference FP (the king's edge) and the
                tr.sh_far_keep best far starts, whatever their score.
    regions     tr.regions regions start at the opening's best configs that are tr.new_far apart in calls. A region
                is a centre (its best config) and a length L (tr.length_init): a candidate changes up to ceil(L x
                settings / 2) settings, each by up to L/2 x its ladder x its box weight (1 / the forest's importance
                of that setting, geometric mean 1, clipped by tr.importance_clip: settings that matter move less).
    steps       every tr.step_results of a region's picks back is one step: a success if its best gained
                tr.success_gain. tr.success_tol successes in a row double L (up to tr.length_max), tr.failure_tol
                failures in a row halve it. L is also capped on the clock (Bounce): the cap falls from tr.length_max at
                the loop's start to tr.length_min at the endgame, so regions end polishing.
    restarts    a region below tr.length_min is finished: its centre goes on a tabu list and a new region starts at
                the best evaluated config tr.new_far (callset distance) from every active and finished centre. Two
                active centres within tr.dup_far: the weaker one is closed. No new region within tr.no_restart_min of
                the endgame (a simulated loss came from restarts at minutes 38-40 that took budget from the leader).
    lines       tr.share_lines of the picks: every setting's untried rungs along its whole ladder through a region's
                centre (cheap settings included), coarse to fine, not ranked and not limited by L (the endgame's
                line scans found the winning region's +17-rung move on 044200Z; its local candidates never did).
    picks       asynchronous Thompson sampling (Kandasamy 2018): every free worker draws one tree of the forest and
                takes the best candidate of all regions under it that is not taken (TuRBO's bandit across regions,
                no fixed shares); candidates the crash model rates more likely bad than not are dropped (SCBO).
                rf.share_inject of the picks are fresh hypercube configs. Every pick carries v10's cheap riders.
    endgame     from rf.endgame_min: the loop's queued picks dropped, v9's grid climb with lines and Step C on the
                centres within tr.endgame_window of the best (at most tr.endgame), then "climbing again".
    pick        the highest local score (unchanged).
"""
from __future__ import annotations

import asyncio
import json
import logging
import math
import random
import time
from collections import deque
from dataclasses import dataclass, field
from typing import Deque, Dict, List, Optional, Sequence, Tuple

import numpy as np

from common.params import ParamRanges
from common.scoring import DIFFICULTY_WEIGHTS
from common.tuning import tune
from common.types import CLASSES, FULL, Config
from search import portfolio
from search.bamgroups import bam_bounds
from search.base import Strategy
from search.client import EvalClient
from search.space import Space, bam_caps
from search.v9 import far as far_step
from search.v9.distance import setting_distance
from search.v9.candidates import OFFICIAL_SEED
from search.v9.climb import Climber
from search.v9.engine import TICK_S, settled_capacity
from search.v10 import proposers
from search.v9.lines import plan, single, tried
from search.v10.engine import Loop, Point, Refit, RFContext, Snapshot, riders

log = logging.getLogger(__name__)

VERSION = "v12"
LOOP_TAGS = {"opening", "tr", "line", "inject"}


@dataclass
class Region:
    id: int
    centre: Config
    best: float
    length: float
    born: float                                    # elapsed seconds
    step_start: float = float("-inf")              # its best when the current step began
    step_n: int = 0
    step_t: float = 0.0                            # elapsed seconds when the current step began
    succ: int = 0
    fail: int = 0
    picks: int = 0
    back: int = 0
    ended: Optional[float] = None
    why: str = ""
    last_gain: float = 0.0                         # elapsed seconds of its last new best
    seeds: Deque[Config] = field(default_factory=deque)   # configs to evaluate first (a fresh region's own centre)
    start_best: float = float("-inf")              # its best when opened (the screen ranks by gain since)
    far: bool = False                              # a screening start far in settings from the opening's top
    indel: Optional[float] = None                  # its best indel part at the reference FP (region_indel)

    def row(self, now: float) -> dict:
        return {"id": self.id, "best": round(self.best, 4), "centre": self.centre.key, "length": round(self.length, 4),
                **({"far": True} if self.far else {}), **({"indel": round(self.indel, 4)} if self.indel is not None else {}),
                **({"start": round(self.start_best, 4)} if self.start_best > float("-inf") else {}),
                "succ": self.succ, "fail": self.fail, "picks": self.picks, "back": self.back,
                "born_m": round(self.born / 60, 1), "age_m": round((now - self.born) / 60, 1),
                **({"ended_m": round(self.ended / 60, 1), "why": self.why} if self.ended is not None else {})}


def box_weights(importances: np.ndarray, clip: float) -> np.ndarray:
    """TuRBO weights a region's box by the model's lengthscales (geometric mean 1); the forest's analogue is 1 /
    importance: a setting that explains much of the score gets a narrow side, one that explains nothing a wide one."""
    imp = np.maximum(np.asarray(importances, dtype=float), 1e-6)
    w = 1.0 / imp
    w = w / math.exp(float(np.mean(np.log(w))))
    return np.clip(w, 1.0 / clip, clip)


def step_update(r: Region, gain: float, succ_tol: int, fail_tol: int, l_max: float) -> str:
    """One step of a region (its tr.step_results picks are back): success or failure, then grow or shrink."""
    r.step_n = 0
    if r.best >= r.step_start + gain:
        r.succ, r.fail = r.succ + 1, 0
        out = "success"
        if r.succ >= succ_tol:
            r.length, r.succ = min(l_max, 2 * r.length), 0
            out = "grow"
    else:
        r.fail, r.succ = r.fail + 1, 0
        out = "failure"
        if r.fail >= fail_tol:
            r.length, r.fail = r.length / 2, 0
            out = "shrink"
    r.step_start = r.best
    return out


def coarse_to_fine(n: int, p: int) -> List[int]:
    """Rung indices 0..n-1 except p, coarse first: both ends and the middle, then halving steps, each level nearest
    to p first. The first few points of a line span the whole ladder."""
    if n <= 1:
        return []
    out, seen = [], {p}
    step = 1
    while step * 2 < n:
        step *= 2
    while step >= 1:
        level = sorted({i for i in range(0, n, step)} | ({n - 1} if step == 1 or not out else set()), key=lambda i: abs(i - p))
        for i in level:
            if i not in seen:
                seen.add(i)
                out.append(i)
        step //= 2
    return out


def region_lines(snap, centre: Config, coarse: int) -> List[Config]:
    """Every setting's untried rungs along its whole ladder through `centre`, the other settings fixed (v9's line
    planner: cheap genotyping/PairHMM settings get every rung, the others their gaps and `coarse` spread rungs),
    ordered coarse to fine and interleaved across settings. Not asked of the forest and not limited by the region's
    length: on 044200Z (2026-10-10) the region whose endgame climb won needed base_quality_score_threshold +17 rungs
    and the gap penalty -7, moves its local candidates never tried (394 configs near its centre, all within 3 rungs)."""
    want = plan(snap, centre, tried(snap, centre), coarse)
    queues = []
    for name, rungs in want.items():
        if not rungs:
            continue
        n = len(snap.ladders.rungs[name])
        p = int(round(snap.ladders.index(name, centre.params.get(name))))
        allowed = set(rungs)
        order = [i for i in coarse_to_fine(n, p) if i in allowed]
        queues.append(deque((name, i) for i in order))
    out = []
    while queues:
        for q in list(queues):
            name, i = q.popleft()
            c = single(snap, centre, name, i)
            if c is not None:
                out.append(c)
            if not q:
                queues.remove(q)
    return out


INDEL_CLASSES = ("indel_1bp", "indel_2_3bp", "indel_4_7bp", "indel_8bp")


def indel_part(counts, weights: Dict[str, float] = DIFFICULTY_WEIGHTS) -> Optional[float]:
    """The four indel terms of core: difficulty-weighted F1 with upstream's class rules (a class without truth drops
    out unless it has calls; weights renormalised over every class kept, SNP classes included), so core = SNP part +
    this. The field's kings lead us here in 18 of 18 rounds (test/field_study/REPORT.md)."""
    total_w, acc = 0.0, 0.0
    for k in CLASSES:
        c = counts.get(k)
        if c is None:
            continue
        tp, fp, fn = c.tp, c.fp, c.fn
        if tp + fn <= 0 and fp <= 0:
            continue
        w = weights[k]
        total_w += w
        if k in INDEL_CLASSES and tp > 0:
            prec, rec = tp / (tp + fp), tp / (tp + fn)
            acc += w * 2 * prec * rec / (prec + rec)
    return acc / total_w if total_w > 0 else None


def screen_keep(regions: Sequence["Region"], keep: int, gain_weight: float, indel_keep: int = 0,
                far_keep: int = 0) -> List["Region"]:
    """Successive halving (Jamieson & Talwalkar 2016; Hyperband): the regions to keep, ranked by best + gain_weight x
    gain since opened - how high it is and how fast it climbs, not where it started (on 044200Z the region that
    started third ended best). Reserved first: the `indel_keep` regions with the best indel part at the reference FP
    (the king's edge; such a hill can trail on score while it is young) and the `far_keep` best far starts (hills
    away from the one every engine so far converged to)."""
    def key(r: "Region") -> float:
        gain = r.best - r.start_best if r.start_best > float("-inf") else 0.0
        return (r.best if r.best > float("-inf") else -1e9) + gain_weight * gain
    ranked = sorted(regions, key=key, reverse=True)
    reserved: List["Region"] = []
    for r in sorted((r for r in regions if r.indel is not None), key=lambda r: -r.indel)[:indel_keep]:
        reserved.append(r)
    for r in [r for r in ranked if r.far and r not in reserved][:far_keep]:
        reserved.append(r)
    reserved = reserved[:keep]
    rest = [r for r in ranked if r not in reserved][:keep - len(reserved)]
    return sorted(reserved + rest, key=key, reverse=True)


def loop_spent(regions: Sequence["Region"], now: float, quiet_s: float) -> bool:
    """Every live region has had no new best for `quiet_s`: the loop has nothing left to gain."""
    return bool(regions) and all(now - r.last_gain >= quiet_s for r in regions)


def clock_cap(frac: float, l_max: float, l_min: float) -> float:
    """Bounce's budget awareness as a cap that only bites late: l_max x (l_min / l_max)^(frac^3) - 0.82 at half the
    loop, 0.105 at 80%, l_min at the endgame. (A plain geometric cap set every region's length in a simulated round:
    the learned sizes never mattered.)"""
    frac = min(1.0, max(0.0, frac))
    return l_max * (l_min / l_max) ** (frac ** 3)


@dataclass
class TRRefit(Refit):
    cands: List[Config] = field(default_factory=list)
    region_of: List[int] = field(default_factory=list)
    order: Optional[np.ndarray] = None            # trees x candidates, each tree's ranking (best first)
    by_region: Dict[int, Deque[Config]] = field(default_factory=dict)   # each region's candidates by mean
    lines: Dict[int, Deque[Config]] = field(default_factory=dict)       # each region's line points


class TRContext(RFContext):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.region_rows: List[dict] = []

    def publish(self, force: bool = False) -> None:
        before = self.last_publish
        super().publish(force)
        eng = self.client.notes.get("engine")
        if eng is not None and (force or self.last_publish != before):
            eng["version"] = VERSION
            eng["regions"] = list(self.region_rows)


class TRLoop(Loop):
    """v10's loop (forest, crash model, riders, calibration log) with trust regions and Thompson-sampled picks."""

    def __init__(self, ctx: TRContext):
        super().__init__(ctx)
        self.regions: Dict[int, Region] = {}
        self.finished: List[Region] = []
        self.next_id = 0
        self.key_region: Dict[str, int] = {}
        self.cands: List[Config] = []
        self.region_of: List[int] = []
        self.order: Optional[np.ndarray] = None
        self.ptr: List[int] = []
        self.by_region: Dict[int, Deque[Config]] = {}
        self.lines: Dict[int, Deque[Config]] = {}
        self.reopened: set = set()
        self.screening = False                    # successive halving over tr.sh_regions regions, until tr.regions
        self.screen_opened = False
        self.last_halving = 0.0
        self.loop_start = 0.0
        self.events: List[str] = []

    def qlock(self, centre: Config, n: int) -> None:   # off: regions set their own size
        return

    def restart_check(self, best_s: float, n_results: int) -> None:   # off: regions restart on their own
        self.best = max(self.best, best_s)

    # ---- distances (event loop: callset sketches read the client) ---------------------------------------------------
    def far_from(self, c: Config, others: Sequence[Config], callset_thr: float, setting_thr: float) -> bool:
        for o in others:
            d, kind = self.ctx.dist.between(c, o)
            if d < (callset_thr if kind == "callset" else setting_thr):
                return False
        return True

    def choose_start(self, pts: Sequence[Point], exclude: Sequence[Config]) -> Optional[Point]:
        """The best evaluated config far (tr.new_far) from every active and finished centre."""
        thr, sthr = float(tune("tr.new_far")), float(tune("tr.setting_new_far"))
        for p in sorted(pts, key=lambda p: -p[1])[:3000]:
            if self.far_from(p[0], exclude, thr, sthr):
                return p
        return None

    def open_region(self, p: Point, why: str, length: Optional[float] = None) -> Region:
        now = self.ctx.elapsed()
        r = Region(self.next_id, p[0], p[1], float(tune("tr.length_init")) if length is None else length, now,
                   step_start=p[1], step_t=now, last_gain=now, start_best=p[1])
        self.next_id += 1
        self.regions[r.id] = r
        self.key_region[p[0].key] = r.id
        self.events.append(f"region r{r.id} opened at {p[0].key[:8]} {p[1]:.4f} ({why})")
        return r

    def fresh_start(self, exclude: Sequence[Config]) -> Optional[Config]:
        """TuRBO restarts from a fresh design: when no evaluated config is far enough from every centre, a new
        hypercube config far (by settings: it has no calls yet) from all of them."""
        ctx = self.ctx
        top = max(((c, s) for c, s, _, _ in self.points()), key=lambda t: t[1], default=None) if hasattr(ctx, "client") else None
        base = top[0] if top else ctx.space._assemble({})
        for c in proposers.inject(ctx.space, ctx.rng, 64, base, list(self.names)):
            if ctx.fresh(c) and self.far_from(c, exclude, float(tune("tr.new_far")), float(tune("tr.setting_new_far"))):
                return c
        return None

    def close_region(self, r: Region, why: str) -> None:
        r.ended, r.why = self.ctx.elapsed(), why
        self.regions.pop(r.id, None)
        self.finished.append(r)
        self.events.append(f"region r{r.id} finished at {r.best:.4f} after {(r.ended - r.born) / 60:.1f}m ({why})")

    def manage(self, pts: Sequence[Point]) -> None:
        """On the event loop before a refit: centres from riders, the clock cap, restarts, duplicates."""
        ctx = self.ctx
        for c, s, _, _ in pts:                       # a rider of a region's pick may beat the pick itself
            rid = self.key_region.get(c.key)
            r = self.regions.get(rid) if rid is not None else None
            if r is not None and s > r.best:
                r.best, r.centre, r.last_gain = s, c, ctx.elapsed()
        l_min, l_max = float(tune("tr.length_min")), float(tune("tr.length_max"))
        end = ctx.T(60 * float(tune("rf.endgame_min")))
        frac = (ctx.elapsed() - self.loop_start) / max(1.0, end - self.loop_start)
        cap = clock_cap(frac, l_max, l_min)
        late = ctx.elapsed() > end - ctx.T(60 * float(tune("tr.no_restart_min")))
        step_s = ctx.T(60 * float(tune("tr.step_max_min")))
        for r in self.regions.values():
            if ctx.elapsed() - r.step_t >= step_s:   # a starved or stuck region is judged on what it has
                what = step_update(r, float(tune("tr.success_gain")), int(tune("tr.success_tol")),
                                   int(tune("tr.failure_tol")), l_max)
                r.step_t = ctx.elapsed()
                if what in ("grow", "shrink"):
                    self.events.append(f"region r{r.id} {what}s to L {r.length:.4f} (best {r.best:.4f}, timed step)")
            r.length = min(r.length, cap)
        if not late:                                 # late, the clock cap shrinks every region: they keep polishing
            for r in sorted(list(self.regions.values()), key=lambda r: r.best):
                if r.length < l_min:
                    self.close_region(r, "collapsed")
        ranked = sorted(self.regions.values(), key=lambda r: -r.best)
        for i, a in enumerate(ranked):
            if a.id not in self.regions:
                continue
            for b in ranked[i + 1:]:
                if b.id in self.regions and not self.far_from(b.centre, [a.centre], float(tune("tr.dup_far")),
                                                               float(tune("tr.setting_dup_far"))):
                    self.close_region(b, f"same hill as r{a.id}")
        if late:                                     # a new region this late cannot catch up: its workers go to the rest
            self.reopen_best(pts)                    # (except a small one on the best)
            return
        if not self.regions and not self.finished and int(tune("tr.sh_regions")) > int(tune("tr.regions")):
            self.screening, self.last_halving = True, ctx.elapsed()
            self.events.append(f"screening: {int(tune('tr.sh_regions'))} regions, halved every "
                               f"{ctx.T(60 * float(tune('tr.sh_step_min'))) / 60:.1f} minutes to {int(tune('tr.regions'))}")
        if self.screening:
            self.screen(pts)
            return
        self.reopen_best(pts)
        while len(self.regions) < int(tune("tr.regions")):
            exclude = [r.centre for r in self.regions.values()] + [r.centre for r in self.finished]
            p = self.choose_start(pts, exclude)
            if p is None:                            # nothing evaluated is far enough: a fresh start (TuRBO)
                c = self.fresh_start(exclude)
                if c is None:
                    break
                r = self.open_region((c, float("-inf"), 0.0, 0.0), "fresh start")
                r.seeds.append(c)
                continue
            self.open_region(p, "restart" if self.finished else "start")

    def screen(self, pts: Sequence[Point]) -> None:
        """Open tr.sh_regions regions once (far apart, going down the opening's results), then halve them every
        tr.sh_step_min by screen_keep until tr.regions remain."""
        ctx, target, keep_n = self.ctx, int(tune("tr.sh_regions")), int(tune("tr.regions"))
        if not self.screen_opened:                   # once: each attempt measures thousands of callset distances
            self.screen_opened = True
            n_far = min(int(tune("tr.sh_far_starts")), max(0, target - 1))
            while len(self.regions) < target - n_far:
                exclude = [r.centre for r in self.regions.values()]
                p = self.choose_start(pts, exclude)
                if p is None:
                    break
                self.open_region(p, "screening start")
            home = [p[0] for p in sorted(pts, key=lambda p: -p[1])[:int(tune("tr.sh_home"))]]
            for _ in range(n_far):
                exclude = [r.centre for r in self.regions.values()]
                p = self.choose_far_start(pts, home, exclude)
                if p is None:
                    break
                self.open_region(p, "far screening start").far = True
            if len(self.regions) <= keep_n:          # too few distinct starts to screen
                self.screening = False
            return
        if ctx.elapsed() - self.last_halving < ctx.T(60 * float(tune("tr.sh_step_min"))):
            return
        self.last_halving = ctx.elapsed()
        live = list(self.regions.values())
        for rid, v in self.region_indel(pts).items():
            self.regions[rid].indel = v
        keep = screen_keep(live, max(keep_n, (len(live) + 1) // 2), float(tune("tr.sh_gain_weight")),
                           int(tune("tr.sh_indel_keep")), int(tune("tr.sh_far_keep")))
        kept = {r.id for r in keep}
        for r in live:
            if r.id not in kept:
                self.close_region(r, "screened out")
        if len(self.regions) <= keep_n:
            self.screening = False
            self.events.append("screening done: " + ", ".join(
                f"r{r.id} {r.best:.4f} (+{r.best - r.start_best:.3f}{', far' if r.far else ''}"
                f"{f', indel {r.indel:.4f}' if r.indel is not None else ''})" for r in sorted(keep, key=lambda r: -r.best)))

    def choose_far_start(self, pts: Sequence[Point], home: Sequence[Config], exclude: Sequence[Config]) -> Optional[Point]:
        """The best evaluated config at least tr.sh_far_setting from every one of the opening's top tr.sh_home in
        settings, and tr.new_far from every opened centre: a start on a hill no engine so far converged to (the
        opening's ranking keeps every other start near our usual hill; the king's is outside it, REPORT.md)."""
        thr = float(tune("tr.sh_far_setting"))
        for p in sorted(pts, key=lambda p: -p[1])[:3000]:
            if all(setting_distance(self.ctx.space, p[0], h) >= thr for h in home) and \
                    self.far_from(p[0], exclude, float(tune("tr.new_far")), float(tune("tr.setting_new_far"))):
                return p
        return None

    def indel_of(self, c: Config) -> Optional[float]:
        client = getattr(self.ctx, "client", None)
        r = client.latest.get((c.key, FULL)) if client is not None else None
        return indel_part(r.counts) if r is not None and r.counts else None

    def region_indel(self, pts: Sequence[Point]) -> Dict[int, float]:
        """Each live region's best indel part among its configs with region FP at most the best config's x (1 +
        tr.indel_fp_tol): indel calling at equal FP, as the field study measured the king's edge."""
        if not pts:
            return {}
        ref = max(pts, key=lambda p: p[1])[3] * (1 + float(tune("tr.indel_fp_tol")))
        out: Dict[int, float] = {}
        for c, _, _, fp in pts:
            rid = self.key_region.get(c.key)
            if rid not in self.regions or fp > ref:
                continue
            v = self.indel_of(c)
            if v is not None and v > out.get(rid, float("-inf")):
                out[rid] = v
        return out

    def reopen_best(self, pts: Sequence[Point]) -> None:
        """The best config in no live region gets a small one (allowed late, one extra slot): on 044200Z the best
        (82.653) sat untouched for 30 minutes after its region collapsed. Once per config."""
        if not pts:
            return
        best = max(pts, key=lambda p: p[1])
        if best[0].key in self.reopened:
            return
        if not self.far_from(best[0], [r.centre for r in self.regions.values()], float(tune("tr.dup_far")),
                             float(tune("tr.setting_dup_far"))):
            return
        if len(self.regions) > int(tune("tr.regions")):
            return
        self.reopened.add(best[0].key)
        self.open_region(best, "the best is in no region", length=float(tune("tr.reopen_length")))

    # ---- the refit (a worker thread on a live fleet) ----------------------------------------------------------------
    def refit(self, pts: Sequence[Point], plan: Sequence[Tuple[int, Config, float]], bad: Sequence[Config] = ()) -> TRRefit:
        ctx, rng = self.ctx, self.ctx.rng
        t0 = time.time()
        out = TRRefit(n=len(pts))
        if len(pts) < 8:
            return out
        have = {p[0].key for p in pts}
        self.forest.fit([p[0] for p in pts], [p[1] for p in pts], float(tune("rf.density_radius")))
        self.forest.fit_bad([p[0] for p in pts], [c for c in bad if c.key not in have])
        out.fit_s = time.time() - t0
        self.refits += 1
        names = self.forest.names
        if int(tune("tr.importance_weights")):
            w = box_weights(self.forest.model.feature_importances_, float(tune("tr.importance_clip")))
        else:                                        # equal: importance is not a lengthscale (044200Z, see region_lines)
            w = np.ones(len(names))
        weight = dict(zip(names, w))
        n_c = int(tune("tr.candidates"))
        cands, region_of, seen = [], [], set()
        for rid, centre, length in plan:
            k_max = max(1, math.ceil(length * len(names) / 2))
            for _ in range(n_c):
                p = dict(centre.params)
                for name in rng.sample(names, rng.randint(1, min(k_max, len(names)))):
                    rungs = ctx.ladders.rungs[name]
                    if len(rungs) <= 4:                    # a switch or a short list: another value
                        p[name] = rng.choice([v for v in rungs if v != p.get(name)] or [p.get(name)])
                        continue
                    reach = max(1, int(round(length / 2 * weight[name] * (len(rungs) - 1))))
                    i = int(round(ctx.ladders.index(name, p.get(name)))) + rng.choice((-1, 1)) * rng.randint(1, reach)
                    p[name] = rungs[min(len(rungs) - 1, max(0, i))]
                c = ctx.space._assemble(p)
                if c is not None and c.key not in seen and ctx.fresh(c):
                    seen.add(c.key)
                    cands.append(c)
                    region_of.append(rid)
        if cands:
            pt = self.forest.per_tree(cands)                       # trees x candidates
            pb = self.forest.p_bad(cands)
            keep = np.where(pb <= 0.5)[0] if pb.any() else np.arange(len(cands))
            cands = [cands[i] for i in keep]
            region_of = [region_of[i] for i in keep]
            pt = pt[:, keep]
            out.order = np.argsort(-pt, axis=1).astype(np.int32)
            mu, sd = pt.mean(axis=0), pt.std(axis=0)
            out.pred = {c.key: (float(mu[i]), float(sd[i]), 0.0) for i, c in enumerate(cands)}
            for i in np.argsort(-mu):
                out.by_region.setdefault(region_of[int(i)], deque()).append(cands[int(i)])
        out.cands, out.region_of = cands, region_of
        snap = Snapshot(ctx, pts)
        for rid, centre, _length in plan:
            out.lines[rid] = deque(c for c in region_lines(snap, centre, int(tune("far.coarse"))) if ctx.fresh(c))
        best_c = max(((c, s) for c, s, _, _ in pts), key=lambda t: t[1])[0]
        out.inject = deque(c for c in proposers.inject(ctx.space, rng, int(int(tune("rf.pool"))
                                                                          * float(tune("rf.share_inject"))),
                                                       best_c, list(self.names)) if ctx.fresh(c))
        out.total_s = time.time() - t0
        return out

    def apply(self, r: Refit) -> None:
        self.inject_q, self.last = r.inject, r
        self.pred = r.pred
        if isinstance(r, TRRefit) and r.order is not None:
            self.cands, self.region_of, self.order = r.cands, r.region_of, r.order
            self.ptr = [0] * r.order.shape[0]
            self.by_region = r.by_region
            self.lines = r.lines
            self.pool = deque((c, "tr") for c in r.cands)          # for "nothing ranked yet" checks only
        now = self.ctx.elapsed()
        self.ctx.region_rows = [x.row(now) for x in sorted(self.regions.values(), key=lambda x: -x.best)]
        self.ctx.rf.update(self.stats())
        for e in self.events:
            self.ctx.say("%s", e)
        self.events = []
        if self.refits % 10 == 1 and self.regions:
            self.ctx.say("regions at t+%.1fm: %s", now / 60, "; ".join(
                f"r{x.id} {x.best:.4f} L {x.length:.3f} ({x.back} back)" for x in
                sorted(self.regions.values(), key=lambda x: -x.best)))

    # ---- picking (event loop) -------------------------------------------------------------------------------------
    def take(self, n: int) -> List[Tuple[Config, str]]:
        """n configs: rf.share_inject a hypercube config, else Thompson sampling: one random tree of the forest, its
        best candidate of all regions not taken yet."""
        ctx = self.ctx
        out: List[Tuple[Config, str]] = []
        taken: set = set()
        share_i = float(tune("rf.share_inject"))
        inject_cap = max(1, int(np.ceil(share_i * n)))
        n_inject, tries = 0, 0
        T = 0 if self.order is None else self.order.shape[0]
        live = [rid for rid in self.regions if self.by_region.get(rid)]
        floor = min(0.5, float(tune("tr.min_share")) * len(live))
        if self.screening:                           # a fair screen: equal turns, not the leader's Thompson draws
            floor = float(tune("tr.sh_floor"))
        for r in self.regions.values():               # a fresh region's own centre first
            while r.seeds and len(out) < n:
                c = r.seeds.popleft()
                if c.key not in taken and ctx.fresh(c):
                    out.append((c, "tr"))
                    taken.add(c.key)
                    self.key_region[c.key] = r.id
                    r.picks += 1
        share_l = float(tune("tr.share_lines"))
        while len(out) < n and tries < 4 * n + 16:
            tries += 1
            with_lines = [rid for rid in self.regions if self.lines.get(rid)]
            if with_lines and ctx.rng.random() < share_l:    # lines through a region's centre (not ranked)
                rid = ctx.rng.choice(with_lines)
                q = self.lines[rid]
                while q:
                    c = q.popleft()
                    if c.key not in taken and ctx.fresh(c):
                        out.append((c, "line"))
                        taken.add(c.key)
                        self.key_region[c.key] = rid
                        self.regions[rid].picks += 1
                        break
                continue
            if live and ctx.rng.random() < floor:      # every region gets enough picks to be judged
                rid = ctx.rng.choice(live)
                q = self.by_region.get(rid) or deque()
                while q:
                    c = q.popleft()
                    if c.key not in taken and ctx.fresh(c):
                        out.append((c, "tr"))
                        taken.add(c.key)
                        self.key_region[c.key] = rid
                        self.regions[rid].picks += 1
                        break
                continue
            if T == 0 or (ctx.rng.random() < share_i and n_inject < inject_cap):
                if n_inject >= inject_cap:
                    break
                c = self.next_inject(n)
                n_inject += 1
                if c is not None and c.key not in taken and ctx.fresh(c):
                    out.append((c, "inject"))
                    taken.add(c.key)
                continue
            t = ctx.rng.randrange(T)
            row, i = self.order[t], self.ptr[t]
            while i < len(row):
                j = int(row[i])
                i += 1
                c = self.cands[j]
                if c.key not in taken and ctx.fresh(c):
                    rid = self.region_of[j]
                    out.append((c, "tr"))
                    taken.add(c.key)
                    self.key_region[c.key] = rid
                    if rid in self.regions:
                        self.regions[rid].picks += 1
                    break
            self.ptr[t] = i
        return out

    def riders(self, c: Config) -> List[Config]:
        out = super().riders(c)
        rid = self.key_region.get(c.key)
        if rid is not None:
            for r in out:
                self.key_region.setdefault(r.key, rid)
        return out

    def note_back(self, c: Config, actual: float) -> None:
        super().note_back(c, actual)
        rid = self.key_region.get(c.key)
        r = self.regions.get(rid) if rid is not None else None
        if r is None:
            return
        r.back += 1
        r.step_n += 1
        if actual > r.best:
            r.best, r.centre, r.last_gain = actual, c, self.ctx.elapsed()
        if r.step_n >= int(tune("tr.step_results")):
            r.step_t = self.ctx.elapsed()
            what = step_update(r, float(tune("tr.success_gain")), int(tune("tr.success_tol")),
                               int(tune("tr.failure_tol")), float(tune("tr.length_max")))
            if what in ("grow", "shrink"):
                self.events.append(f"region r{r.id} {what}s to L {r.length:.4f} (best {r.best:.4f})")

    def stats(self) -> dict:
        out = super().stats()
        now = self.ctx.elapsed()
        out["regions"] = [r.row(now) for r in sorted(self.regions.values(), key=lambda x: -x.best)]
        out["finished"] = [r.row(now) for r in self.finished[-40:]]
        return out


def endgame_centres(regions: Sequence[Region], count: int, window: float) -> List[Region]:
    ranked = sorted(regions, key=lambda r: -r.best)
    if not ranked:
        return []
    return [r for r in ranked if r.best >= ranked[0].best - window][:count]


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: Optional[int] = None) -> Strategy:
    import os
    if seed is None:
        seed = int(os.environ.get("MINOS_SEED") or (time.time_ns() % 1_000_000_007))
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        r = client.round
        bounds = bam_bounds(r.bam, ranges, tool) if r.bam else None
        caps = {} if bounds is not None else (bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {})
        profile = portfolio.bam_profile(r.bam, r.chrom, r.start, r.end) if r.bam else None
        space = Space(ranges, tool, caps=caps, bounds=bounds)
        if bounds is not None:
            log.info("engine: bam groups: %s", bounds.report)
        ctx = TRContext(client, space, rng, profile)
        ctx.cap = await settled_capacity(client, workers)
        virtual = hasattr(client.executor, "now")
        ctx.say("%s trust regions: %d workers, %d settings searched, a %d-minute round (schedule x%.2f), seed %d, "
                "%d regions", VERSION, ctx.cap, len(space.names), round(ctx.length / 60), ctx.scale, seed,
                int(tune("tr.regions")))
        loop = TRLoop(ctx)
        if not virtual:
            from pathlib import Path
            from common.settings import SETTINGS
            d = Path(SETTINGS.status_dir) / "rf_picks"
            try:
                d.mkdir(parents=True, exist_ok=True)
                loop.log_path = d / f"{r.tag}_{int(time.time())}.jsonl"
            except OSError:
                loop.log_path = None
        pending: Dict[asyncio.Future, str] = {}

        def send(c: Config, src: str, extra: Sequence[Config] = ()) -> None:
            batch = [c] + [x for x in extra if x.key != c.key]
            futs = ctx.propose(batch, src, ["small"] + ["big"] * (len(batch) - 1))
            pending[futs[0]] = src
            loop.picks[src] += 1

        async def drain(timeout: float) -> int:
            if not pending:
                await asyncio.sleep(0 if virtual else timeout)
                return 0
            done, _ = await asyncio.wait(list(pending), timeout=timeout, return_when=asyncio.FIRST_COMPLETED)
            t = time.time()
            for f in done:
                src = pending.pop(f)
                if f.cancelled():
                    continue
                o = f.result()
                if o is None or not o.ok or o.result.gate_passed is False or not o.score:
                    if o is not None:
                        loop.sent.pop(o.config.key, None)
                        err = o.error or ""
                        if not err.startswith(("invalid", "budget", "round ended")):
                            loop.bad[o.config.key] = o.config
                    continue
                loop.note_back(o.config, o.score)
                if o.score > loop.best:
                    loop.wins[src] += 1
                    loop.best = o.score
            loop.results += len(done)
            loop.loop_s += time.time() - t
            ctx.publish()
            return len(done)

        # ---- opening (v10's) -------------------------------------------------------------------------------------------
        ctx.phase = "A: opening"
        seeds: List[Config] = []
        for p in ({}, OFFICIAL_SEED):
            c = space._assemble(dict(p))
            if c is not None and all(c.key != x.key for x in seeds):
                seeds.append(c)
        for c in portfolio.pick(space, profile, int(ctx.cap * float(tune("rf.portfolio_share"))), exclude_round=r.tag):
            if all(c.key != x.key for x in seeds):
                seeds.append(c)
        hyper = [c for c in proposers.inject(space, rng, max(0, ctx.cap - len(seeds)), space._assemble({}), loop.names)
                 if all(c.key != x.key for x in seeds)]
        ctx.say("opening: %d seeds (%d from the portfolio) and %d hypercube configs", len(seeds), len(seeds) - 2,
                len(hyper))
        for c in seeds + hyper:
            if ctx.fresh(c):
                send(c, "opening", riders(ctx, c, int(tune("candidates.riders"))))
        n_open = len(pending)

        def endgame() -> float:
            return ctx.t0 + ctx.T(60 * float(tune("rf.endgame_min")))
        while pending and n_open - len(pending) < float(tune("rf.opening_back")) * n_open and ctx.clock() < endgame() \
                and not ctx.over():
            await drain(TICK_S)
        ctx.marks["a"] = ctx.elapsed()
        top = ctx.best()
        ctx.say("opening: %d of %d back at t+%.1fm, best %.4f", n_open - len(pending), n_open, ctx.elapsed() / 60,
                top[1] if top else float("nan"))

        # ---- trust regions -------------------------------------------------------------------------------------------
        ctx.phase = "B: trust regions"
        loop.loop_start = ctx.elapsed()

        def refit_n() -> int:
            return max(8, ctx.cap // int(tune("rf.refit_per_fleet")))
        since = refit_n()
        task: Optional[asyncio.Task] = None
        headroom = max(4, ctx.cap // 20)

        def start_refit():
            pts = loop.points()
            t = time.time()
            if pts:
                loop.restart_check(max(p[1] for p in pts), len(pts))
                loop.manage(pts)
            plan = [(x.id, x.centre, x.length) for x in loop.regions.values()]
            loop.loop_s += time.time() - t
            return pts, plan, list(loop.bad.values())

        try:
            loop.apply(loop.refit(*start_refit()))
            since = 0
        except Exception:  # noqa: BLE001
            log.exception("engine: first refit failed")
        while ctx.clock() < endgame() and not ctx.over():
            try:
                if task is None and (since >= refit_n() or loop.order is None):
                    pts, plan, bad = start_refit()
                    since = 0
                    if virtual:
                        loop.apply(loop.refit(pts, plan, bad))
                    else:
                        task = asyncio.create_task(asyncio.to_thread(loop.refit, pts, plan, bad))
                if task is not None and task.done():
                    t_, task = task, None
                    loop.apply(t_.result())
                idle = client.idle_slots()
                if idle is None:
                    idle = max(0, ctx.cap - len(pending))
                free = idle + headroom if idle > 0 else 0
                if free > 0:
                    t = time.time()
                    batch = loop.take(free)
                    loop.note_sent(batch)
                    for c, src in batch:
                        send(c, src, loop.riders(c))
                    loop.loop_s += time.time() - t
            except Exception:  # noqa: BLE001
                log.exception("engine: trust-region step failed")
                task = None if task is not None and task.done() else task
            got = await drain(TICK_S)
            since += got
            loop.stall += got
            end_s = ctx.T(60 * float(tune("rf.endgame_min")))
            late = ctx.elapsed() > end_s - ctx.T(60 * float(tune("tr.no_restart_min")))
            if late and loop_spent(list(loop.regions.values()), ctx.elapsed(), ctx.T(60 * float(tune("tr.spent_min")))):
                ctx.say("trust regions spent at t+%.1fm: no region gained in %.0f minutes and no new region may open; "
                        "the endgame starts now", ctx.elapsed() / 60, ctx.T(60 * float(tune("tr.spent_min"))) / 60)
                break
        task = None
        ctx.marks["b"] = ctx.elapsed()
        ctx.rf.update(loop.stats())
        ctx.say("trust regions done at t+%.1fm: %s", ctx.elapsed() / 60,
                json.dumps({k: v for k, v in loop.stats().items() if k not in ("centres", "finished", "calibration")}))

        # ---- endgame: v9's grid climb with lines, then Step C, on the best centres ------------------------------------
        ctx.phase = "C: endgame climb"
        ctx.priority = int(tune("rf.endgame_priority"))
        loop_tags = {tag for (_f, _c, tag) in list(client._pending.values()) if tag in LOOP_TAGS}
        dropped = client.drop_queued(loop_tags)
        if dropped:
            ctx.say("endgame: %d of the loop's queued picks dropped; %d jobs still in flight", dropped, client.in_flight)
        chosen = endgame_centres(list(loop.regions.values()), int(tune("tr.endgame")), float(tune("tr.endgame_window")))
        starts = [(x.centre, ctx.score(x.centre) or x.best, f"region {x.id}") for x in chosen]
        top = ctx.best()
        if top is not None and all(top[0].key != c.key for c, _, _ in starts):
            starts.insert(0, (top[0], top[1], "best"))      # the best config may sit in no live region
        for c, s, label in starts[:max(1, int(tune("tr.endgame")))]:
            cl = Climber(ctx, len(ctx.climbers), c, s, origin=f"v12 {label}")
            cl.rider_cap = int(tune("rf.endgame_riders"))
            ctx.climbers.append(cl)

        async def lifecycle(cl: Climber) -> None:
            try:
                await cl.run(ctx.pick_at())
                if cl.state in ("peak", "stopped") and not ctx.over():
                    await far_step.improve(ctx, cl, ctx.pick_at())
            except Exception:  # noqa: BLE001
                log.exception("engine: endgame climber %d failed", cl.id)
            ctx.publish(force=True)

        await asyncio.gather(*(lifecycle(cl) for cl in list(ctx.climbers)))
        while ctx.clock() < ctx.pick_at() and not ctx.over():
            ranked = sorted(((c, s) for c, s, _, _ in loop.points()), key=lambda t: -t[1])
            mc, ms = float(tune("candidates.min_distance")), float(tune("candidates.min_setting_distance"))
            nxt = next(((c, s) for c, s in ranked if all(c.key not in (d.start.key, d.centre.key)
                                                         and ctx.dist.far(c, d.centre, mc, ms) for d in ctx.climbers)),
                       None)
            if nxt is None:
                break
            cl = Climber(ctx, len(ctx.climbers), nxt[0], nxt[1], origin="climbing again")
            cl.rider_cap = int(tune("rf.endgame_riders"))
            ctx.climbers.append(cl)
            ctx.say("endgame: climbing again from %s %.4f", nxt[0].key[:8], nxt[1])
            await lifecycle(cl)
        ctx.marks["c"] = ctx.elapsed()

        # ---- the pick --------------------------------------------------------------------------------------------------
        ctx.phase = "picked"
        ctx.name_pick()
        ctx.publish(force=True)
        while not ctx.over() and client.in_flight:
            await asyncio.sleep(0 if virtual else TICK_S)
            ctx.name_pick()
            ctx.publish()

    return run
