"""Search engine v10: RF_UCB (search/v10/DESIGN.md).

    opening     the defaults, the official seed, rf.portfolio_share of the workers from the portfolio (small, so the
                forest is not fitted around other rounds' configs) and a Latin hypercube for the rest, each with
                calling-threshold riders. The loop starts once rf.opening_back of it is back.
    RF_UCB      one asynchronous loop, no waves: every free worker takes the next config of a pool ranked by the
                forest (forest.py: UCB = mean + kappa x tree spread).
    refits      in the background, throttling themselves: one starts when the previous one has finished and at least
                workers/rf.refit_per_fleet new results are in (the legacy cadence), on every result in by then. In the
                opening's burst they run back to back; the loop never waits for one.
    pool        the proposers (proposers.py) fill one pool, ranked as a whole by UCB:
                    QGA     rf.share_qga of rf.pool   probability vectors over the ladders, rotated toward the centres
                    GA      rf.share_ga of rf.pool    tournament crossover of the elites and the core/FP front
                    Local   every setting's line and rf.local_moves small moves around each centre
                The centres: the best score and the two ends of the core/false-positive trade-off (the highest core
                and the fewest false positives within rf.ends_window of the best), then the best distinct scores, up
                to rf.centres configs far apart by callset distance.
    batch       every free worker takes, by chance: rf.share_inject a fresh hypercube config (unranked); rf.share_mutate
                a mutation (1-3 settings, 1-3 rungs) of one of the rf.mutate_top best-ranked proposals, unranked (a
                forest predicts the same score across a leaf's box, so the top of the ranking is often near-duplicates:
                the mutations spread the batch in and around the best box and hedge the forest's error there);
                rf.qga_picks QGA's best-ranked remaining proposal (in one ranking QGA's many-setting draws never
                reached the top: 0 picks in two simulated rounds); rf.share_sweep the next of every setting's line
                through the actual best, unranked, rebuilt at every refit (legacy's sweep lane); the rest the next
                config of the ranked pool. Every pick carries rf.riders cheap riders.
    QLock       freezes settings the forest finds flat around the best (rf.qlock_min_obs results, rf.qlock_resolution).
    Restart     no gain of rf.restart_min_gain for rf.restart_waves fleets' worth of results: QGA reset, locks lifted;
    Reinit      a second stall in a row: the next fleet's worth of picks are fresh hypercube. The forest keeps every
                result through both.
    endgame     from rf.endgame_min (scaled): v9's grid climb with lines on the best score, the best core and the best FP
                (rf.endgame_hills distinct configs), then
                v9's Step C around each peak, until the pick.
    pick        the highest local score, far.pick_s before the end.
"""
from __future__ import annotations

import asyncio
import json
import logging
import random
import time
from collections import Counter, deque
from dataclasses import dataclass, field
from typing import Deque, Dict, List, Optional, Sequence, Tuple

import numpy as np

from common.params import ParamRanges
from common.tuning import tune
from common.types import FULL, Config
from search import portfolio
from search.bamgroups import bam_bounds
from search.base import Strategy
from search.client import EvalClient
from search.space import CALL_CONF_GRID, Space, bam_caps
from search.v9 import far
from search.v9.candidates import CONF, OFFICIAL_SEED
from search.v9.climb import Climber
from search.v9.engine import TICK_S, Context, settled_capacity
from search.v9.lines import plan, single, tried
from search.v10 import proposers
from search.v10.forest import Forest

log = logging.getLogger(__name__)

VERSION = "v10"
Point = Tuple[Config, float, float, float]            # config, score, core, false positives


class Snapshot:
    """The scored configs as they were when a refit started: what the refit thread reads instead of the live client
    (whose dicts the event loop keeps changing)."""

    def __init__(self, ctx, points: Sequence[Point]):
        self.space, self.ladders, self.dist = ctx.space, ctx.ladders, ctx.dist
        self._scored = [(c, s) for c, s, _, _ in points]

    def scored(self):
        return iter(self._scored)


@dataclass
class Refit:
    """What one refit produced, applied on the event loop."""
    pool: Deque[Tuple[Config, str]] = field(default_factory=deque)      # ranked by UCB, best first
    top: List[Config] = field(default_factory=list)                      # the best-ranked, for mutations
    qga: Deque[Config] = field(default_factory=deque)                   # QGA's own proposals, ranked, for its share
    pred: Dict[str, Tuple[float, float, float]] = field(default_factory=dict)   # key -> (mean, spread, novelty)
    sweep: Deque[Config] = field(default_factory=deque)                 # every setting's line through the actual best
    inject: Deque[Config] = field(default_factory=deque)
    fit_s: float = 0.0
    total_s: float = 0.0
    n: int = 0
    frozen: List[str] = field(default_factory=list)


class RFContext(Context):
    """v9's context (proposing, waiting, scores, distances, the pick) with v10's notes."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.rf: Dict[str, object] = {}

    def publish(self, force: bool = False) -> None:
        now = self.clock()
        if now >= self.pick_at():                    # a climb still waiting for its wave must not miss the pick
            self.name_pick()
        if not force and now - self.last_publish < TICK_S:
            return
        self.last_publish = now
        try:
            top = self.best()
            end = round(self.T(60 * float(tune("rf.endgame_min"))), 1)
            self.client.notes["engine"] = {
                "version": VERSION, "phase": self.phase, "elapsed_s": round(self.elapsed(), 1),
                "started_ts": round(self.started_ts, 1), "capacity": self.cap,
                "best": None if top is None else round(top[1], 4),
                "schedule": {"length_s": round(self.length, 1), "scale": round(self.scale, 4), "a_max_s": end,
                             "b_end_s": end, "pick_s": round(self.length - float(tune("far.pick_s")), 1),
                             **{f"{k}_done_s": round(v, 1) for k, v in self.marks.items()}},
                "rf": dict(self.rf),
                "basins": [c.row() for c in self.climbers],
                "step_c": list(self.c_log),
                "pick": self.picked,
                "events": list(self.events[-200:]),
            }
        except Exception:  # noqa: BLE001 - the notes must never end the round
            log.exception("engine: publishing the notes failed")


class Loop:
    """The RF_UCB loop's state: forest, proposers, the ranked pool, QLock, Restart/Reinitialize."""

    def __init__(self, ctx: RFContext):
        self.ctx = ctx
        self.forest = Forest(ctx.space, ctx.ladders, trees=int(tune("rf.trees")), min_leaf=int(tune("rf.min_leaf")),
                             jobs=int(tune("rf.fit_threads")), sample=float(tune("rf.tree_sample")))
        self.qga = proposers.QGA(ctx.space, ctx.ladders, int(tune("rf.qga_bins")), ctx.rng)
        self.names = list(self.qga.names)
        self.frozen: set = set()
        self.centres: List[Config] = []
        self.pool: Deque[Tuple[Config, str]] = deque()
        self.top: List[Config] = []
        self.qga_q: Deque[Config] = deque()
        self.sweep_q: Deque[Config] = deque()
        self.inject_q: Deque[Config] = deque()
        self.picks = Counter()
        self.wins = Counter()                          # new bests by source
        self.best = float("-inf")
        self.gain_mark = float("-inf")
        self.stall = 0                                 # results since the last gain of rf.restart_min_gain
        self.restarts = self.reinits = 0
        self.last_was_restart = False
        self.reinit_left = 0
        self.qlock_after = 0                           # QLock waits until this many results (after a restart)
        self.refits = 0
        self.last: Optional[Refit] = None
        self.loop_s = 0.0                              # time the event loop spent on the controller's own work
        self.results = 0
        self.pred: Dict[str, Tuple[float, float, float]] = {}     # the current pool's predictions
        self.sent: Dict[str, Tuple[str, float, float, float, int]] = {}   # picked key -> (source, mean, spread, novelty, refit)
        self.rows: List[Tuple[str, float, float, float, float]] = []   # back: (source, mean, spread, novelty, actual)
        self.log_path = None
        self.bad: Dict[str, Config] = {}               # final failures and zero-scored configs: taught as a floor

    # ---- data -----------------------------------------------------------------------------------------------------
    def points(self) -> List[Point]:
        out = []
        client = self.ctx.client
        for (key, fid), r in list(client.latest.items()):
            if fid != FULL or r.score is None or r.gate_passed is False or r.score <= 0:
                continue                               # zero scores reach the forest as bad configs, at the floor
            c = client.configs.get(key)
            if c is not None:
                out.append((c, r.score, r.core if r.core is not None else 0.0, float(r.region_fp or 0)))
        return out

    def elites(self, pts: Sequence[Point]) -> List[Tuple[Config, float]]:
        top = sorted(((c, s) for c, s, _, _ in pts), key=lambda t: -t[1])[:int(tune("rf.elites"))]
        have = {c.key for c, _ in top}
        return top + [t for t in proposers.pareto(pts, float(tune("rf.front_window")), 32) if t[0].key not in have]

    def pick_centres(self, pts: Sequence[Point], n: Optional[int] = None) -> List[Config]:
        """On the event loop (callset distances read the client): the best score and the two ends of the core/false-
        positive trade-off (proposers.trade_off_ends, within rf.ends_window of the best), then the best distinct
        scores, up to n (rf.centres) configs far apart."""
        n = int(tune("rf.centres")) if n is None else n
        mc, ms = float(tune("candidates.min_distance")), float(tune("candidates.min_setting_distance"))
        out: List[Config] = []
        for c, _ in proposers.trade_off_ends(pts, float(tune("rf.ends_window"))):
            # the ends always take their slots, unless the same config: they are different trade-offs even when most
            # of their calls agree (second fleet replay: the best-core config, core 0.8745 like the kings, was
            # dropped as "near" the best score by callset distance)
            if len(out) < n and all(c.key != o.key for o in out):
                out.append(c)
        ranked = sorted(((c, s) for c, s, _, _ in pts), key=lambda t: -t[1])
        for c, _ in ranked:
            if len(out) >= n:
                break
            if all(self.ctx.dist.far(c, o, mc, ms) for o in out):
                out.append(c)
        return out

    # ---- the refit (a worker thread on a live fleet) ----------------------------------------------------------------
    def refit(self, pts: Sequence[Point], centres: Sequence[Config], bad: Sequence[Config] = ()) -> Refit:
        ctx, rng = self.ctx, self.ctx.rng
        t0 = time.time()
        out = Refit(n=len(pts))
        if len(pts) < 8:
            return out
        # Configs that crashed for good (after the evaluator's own solo and cold retries: stock GATK crashes too) or
        # scored nothing are real outcomes on the validator (first fleet replay: 6,246 final failures, 7,790 zero
        # scores). A separate classifier learns them (forest.fit_bad) and UCB sinks likely ones; the score forest is
        # fitted on scored configs only.
        have = {p[0].key for p in pts}
        extra = [c for c in bad if c.key not in have]
        self.forest.fit([p[0] for p in pts], [p[1] for p in pts], float(tune("rf.density_radius")))
        self.forest.fit_bad([p[0] for p in pts], extra)
        out.fit_s = time.time() - t0
        self.refits += 1
        best_c, best_s = max(((c, s) for c, s, _, _ in pts), key=lambda t: t[1])
        self.qlock(best_c, len(pts))
        out.frozen = sorted(self.frozen)
        centres = list(centres) or [best_c]
        progress = min(1.0, ctx.elapsed() / max(1.0, ctx.T(60 * float(tune("rf.endgame_min")))))
        theta = float(tune("rf.qga_theta")) + (float(tune("rf.qga_theta_end")) - float(tune("rf.qga_theta"))) * progress
        self.qga.rotate(centres, theta, float(tune("rf.qga_mut")), float(tune("rf.qga_floor")))
        live = [n for n in self.names if n not in self.frozen]
        size = int(tune("rf.pool"))
        made: List[Tuple[Config, str]] = []
        made += [(c, "qga") for c in self.qga.sample(int(size * float(tune("rf.share_qga"))), best_c, live)]
        made += [(c, "ga") for c in proposers.breed(ctx.space, ctx.ladders, rng, self.elites(pts),
                                                    int(size * float(tune("rf.share_ga"))), best_c, live,
                                                    float(tune("rf.ga_mut_rate")), int(tune("rf.ga_mut_rungs")))]
        made += [(c, "local") for c in proposers.local(Snapshot(ctx, pts), rng, centres, int(tune("far.coarse")),
                                                       int(tune("rf.local_moves")))]
        seen, fresh = set(), []
        for c, src in made:
            if c.key not in seen and ctx.fresh(c):
                seen.add(c.key)
                fresh.append((c, src))
        if fresh:
            u, mu, sd, nov = self.forest.ucb([c for c, _ in fresh], float(tune("rf.kappa")),
                                             float(tune("rf.kappa_dist")))
            order = np.argsort(-u)
            out.pred = {c.key: (float(mu[i]), float(sd[i]), float(nov[i])) for i, (c, _) in enumerate(fresh)}
            out.pool = deque(fresh[int(i)] for i in order)
            out.top = [fresh[int(i)][0] for i in order[:int(tune("rf.mutate_top"))]]
            out.qga = deque(fresh[int(i)][0] for i in order if fresh[int(i)][1] == "qga")
        out.inject = deque(c for c in proposers.inject(ctx.space, rng, int(size * float(tune("rf.share_inject"))),
                                                       best_c, live) if ctx.fresh(c))
        out.sweep = deque(c for c in sweep_lines(Snapshot(ctx, pts), best_c, int(tune("far.coarse"))) if ctx.fresh(c))
        out.total_s = time.time() - t0
        return out

    def restart_check(self, best_s: float, n_results: int) -> None:
        """On the event loop, before a refit starts (no refit running): Restart and Reinitialize. None within
        rf.no_restart_min of the endgame: a Reinitialize there (first fleet replay: minute 49 of 54) spends the loop's
        last minutes on random configs nothing can climb from."""
        ctx = self.ctx
        self.best = max(self.best, best_s)
        if self.best >= self.gain_mark + float(tune("rf.restart_min_gain")):
            self.gain_mark, self.stall, self.last_was_restart = self.best, 0, False
            return
        if self.stall < float(tune("rf.restart_waves")) * ctx.cap:
            return
        endgame = ctx.t0 + ctx.T(60 * float(tune("rf.endgame_min")))
        if ctx.clock() > endgame - ctx.T(60 * float(tune("rf.no_restart_min"))):
            return
        self.stall = 0
        self.qga.reset()
        self.frozen.clear()
        self.qlock_after = n_results + int(tune("rf.qlock_min_obs"))    # or the same forest re-freezes them at once
        if self.last_was_restart:
            self.reinits += 1
            self.reinit_left = ctx.cap
            ctx.say("RF_UCB: reinitialize %d (no gain of %.3f since %.4f): the next %d picks are fresh hypercube",
                    self.reinits, float(tune("rf.restart_min_gain")), self.gain_mark, ctx.cap)
        else:
            self.restarts += 1
            ctx.say("RF_UCB: restart %d (no gain of %.3f since %.4f): QGA reset, locks lifted", self.restarts,
                    float(tune("rf.restart_min_gain")), self.gain_mark)
        self.last_was_restart = not self.last_was_restart

    def qlock(self, centre: Config, n: int) -> None:
        if n < max(int(tune("rf.qlock_min_obs")), self.qlock_after):
            return
        res = float(tune("rf.qlock_resolution"))
        newly = [x for x in self.names if x not in self.frozen and self.forest.flat(centre, x, 12, res)]
        if newly:
            self.frozen.update(newly)
            self.ctx.say("RF_UCB: QLock froze %s (flat within %.3f at the best)", ", ".join(newly), res)

    def apply(self, r: Refit) -> None:
        """On the event loop: the new pool replaces the old one."""
        self.pool, self.top, self.qga_q, self.inject_q, self.last = r.pool, r.top, r.qga, r.inject, r
        self.pred = r.pred
        self.sweep_q = r.sweep
        self.ctx.rf.update(self.stats())
        if self.refits % 10 == 1:
            self.ctx.say("RF_UCB refit %d at t+%.1fm on %d results: best %.4f, fit %.1fs (refit %.1fs), pool %d, "
                         "frozen %d, loop %.2f ms/result", self.refits, self.ctx.elapsed() / 60, r.n, self.best,
                         r.fit_s, r.total_s, len(self.pool), len(self.frozen),
                         1000 * self.loop_s / max(1, self.results))

    # ---- picking (event loop) -------------------------------------------------------------------------------------
    def mutate(self) -> Optional[Config]:
        """A config 1-3 settings and 1-3 rungs away from one of the best-ranked proposals (not asked of the forest)."""
        ctx = self.ctx
        live = [x for x in self.names if x not in self.frozen]
        if not self.top or not live:
            return None
        for _ in range(20):
            p = dict(ctx.rng.choice(self.top).params)
            for name in ctx.rng.sample(live, min(len(live), ctx.rng.randint(1, 3))):
                p[name] = proposers.shift(ctx.ladders, name, p[name], ctx.rng.randint(1, 3) * ctx.rng.choice((1, -1)),
                                          ctx.space)
            c = ctx.space._assemble(p)
            if c is not None and ctx.fresh(c):
                return c
        return None

    def next_inject(self, n: int) -> Optional[Config]:
        ctx = self.ctx
        if not self.inject_q:
            top = ctx.best()
            base = top[0] if top else ctx.space._assemble({})
            self.inject_q.extend(proposers.inject(ctx.space, ctx.rng, max(n, 64), base,
                                                  [x for x in self.names if x not in self.frozen]))
        while self.inject_q:
            c = self.inject_q.popleft()
            if ctx.fresh(c):
                return c
        return None

    def take(self, n: int) -> List[Tuple[Config, str]]:
        """n configs: by chance an Inject (rf.share_inject), a mutation of a top proposal (rf.share_mutate), QGA's best
        remaining proposal (rf.qga_picks), the next of the sweep through the actual best (rf.share_sweep), else the
        next of the ranked pool. Reinitialize: all Inject."""
        ctx = self.ctx
        share_i, share_m, share_q, share_s = (float(tune("rf.share_inject")), float(tune("rf.share_mutate")),
                                              float(tune("rf.qga_picks")), float(tune("rf.share_sweep")))
        out: List[Tuple[Config, str]] = []
        taken: set = set()

        def ok(c: Optional[Config]) -> bool:
            return c is not None and c.key not in taken and ctx.fresh(c)

        def first(q: Deque) -> Optional[Config]:
            while q:
                c = q.popleft()
                c = c[0] if isinstance(c, tuple) else c
                if ok(c):
                    return c
            return None

        tries = 0
        inject_cap = n if self.reinit_left > 0 else max(1, int(np.ceil(share_i * n)))
        n_inject = 0
        while len(out) < n and tries < 4 * n + 16:
            tries += 1
            if not (self.pool or self.qga_q or self.sweep_q or self.top) and self.reinit_left <= 0 \
                    and n_inject >= inject_cap:
                break                                      # nothing ranked yet: Inject keeps its share, the rest waits
            roll = ctx.rng.random()
            if (self.reinit_left > 0 or roll < share_i) and n_inject < inject_cap:
                c = self.next_inject(n)
                if ok(c):
                    out.append((c, "inject"))
                    n_inject += 1
                    self.reinit_left = max(0, self.reinit_left - 1)
                elif c is None:
                    self.reinit_left = 0
                if c is not None:
                    taken.add(c.key)
                continue
            src, c = None, None
            if roll < share_i + share_m:
                src, c = "mutate", self.mutate()
            elif roll < share_i + share_m + share_q:
                src, c = "qga", first(self.qga_q)
            elif roll < share_i + share_m + share_q + share_s:
                src, c = "sweep", first(self.sweep_q)
            if not ok(c):                                  # nothing there: the ranked pool
                while self.pool:
                    c2, src2 = self.pool.popleft()
                    if ok(c2):
                        src, c = src2, c2
                        break
                else:
                    c = None
            if c is not None and ok(c):
                out.append((c, src))
                taken.add(c.key)
        return out

    def riders(self, c: Config) -> List[Config]:
        """rf.riders cheap settings (genotyping, PairHMM) moved one rung each, in the pick's own GATK run."""
        ctx = self.ctx
        cheap = [n for n in self.names if n in proposers.CHEAP]
        ctx.rng.shuffle(cheap)
        chosen, pairhmm = [], 0
        for name in cheap:                         # a small batch's head carries at most 0.3 of a run: two PairHMM
            if len(chosen) >= int(tune("rf.riders")):     # riders (0.12 each) fit, a third would run on its own
                break
            if name in proposers.PAIRHMM:
                if pairhmm >= 2:
                    continue
                pairhmm += 1
            chosen.append(name)
        out = []
        for name in chosen:
            v = ctx.ladders.shift(name, c.params.get(name), ctx.rng.choice((1, -1)))
            r = ctx.space._assemble({**c.params, name: v}) if v is not None else None
            if r is not None and ctx.fresh(r):
                out.append(r)
        return out

    def note_sent(self, batch: Sequence[Tuple[Config, str]]) -> None:
        """Remember what the forest predicted for each pick (unranked picks are predicted here, in one pass)."""
        missing = [c for c, _ in batch if c.key not in self.pred]
        if missing and self.forest.ready:
            _, mu, sd, nov = self.forest.ucb(missing, 0.0, 0.0)
            for i, c in enumerate(missing):
                self.pred[c.key] = (float(mu[i]), float(sd[i]), float(nov[i]))
        for c, src in batch:
            m, sd_, nv = self.pred.get(c.key, (float("nan"),) * 3)
            self.sent[c.key] = (src, m, sd_, nv, self.refits)

    def note_back(self, c: Config, actual: float) -> None:
        got = self.sent.pop(c.key, None)
        if got is None:
            return
        src, m, sd_, nv, refit = got
        self.rows.append((src, m, sd_, nv, actual))
        if self.log_path is not None:
            try:
                with open(self.log_path, "a") as f:
                    f.write(json.dumps({"t": round(self.ctx.elapsed(), 1), "key": c.key, "src": src, "refit": refit,
                                        "mean": round(m, 4), "spread": round(sd_, 4), "novelty": round(nv, 4),
                                        "actual": round(actual, 4)}) + "\n")
            except OSError:
                self.log_path = None

    def calibration(self) -> dict:
        """Per source: picks back, how far the forest's prediction was off on average, how well it ordered them
        (Spearman), and their median novelty: did UCB send workers into new regions, and was it right there?"""
        out = {}
        by: Dict[str, list] = {}
        for src, m, sd_, nv, a in self.rows:
            if m == m:                                     # predicted (not NaN)
                by.setdefault(src, []).append((m, sd_, nv, a))
        for src, rows in by.items():
            arr = np.array(rows)
            rho = None
            if len(arr) >= 10 and np.std(arr[:, 0]) > 0 and np.std(arr[:, 3]) > 0:
                rk = lambda x: np.argsort(np.argsort(x))
                rho = float(np.corrcoef(rk(arr[:, 0]), rk(arr[:, 3]))[0, 1])
            out[src] = {"n": len(arr), "bias": round(float(np.mean(arr[:, 3] - arr[:, 0])), 4),
                        "spearman": None if rho is None else round(rho, 3),
                        "novelty_median": round(float(np.median(arr[:, 2])), 4)}
        return out

    def stats(self) -> dict:
        r = self.last
        return {"refits": self.refits, "fit_s": round(r.fit_s, 1) if r else None,
                "refit_s": round(r.total_s, 1) if r else None, "fit_on": r.n if r else 0,
                "restarts": self.restarts, "reinits": self.reinits, "frozen": sorted(self.frozen),
                "best": round(self.best, 4), "pool": len(self.pool), "picks": dict(self.picks),
                "new_bests": dict(self.wins), "centres": [c.key for c in self.centres], "bad": len(self.bad),
                "calibration": self.calibration(),
                "loop_ms_per_result": round(1000 * self.loop_s / max(1, self.results), 2)}


def sweep_lines(snap, centre: Config, coarse: int) -> List[Config]:
    """The sweep lane (legacy: the coordinate sweeps beside the model produced every winning config after the opening):
    each setting's untried values through `centre`, one setting at a time, the nearest rungs first (search/v9/lines.py:
    the gaps, then rungs beyond the nearest tried value). Not asked of the forest."""
    known = tried(snap, centre)
    out: List[Tuple[float, Config]] = []
    for name, rungs in plan(snap, centre, known, coarse).items():
        if name in proposers.CHEAP:
            continue                                   # the cheap settings ride along anyway (rf.riders)
        p = snap.ladders.index(name, centre.params.get(name))
        for i in rungs:
            c = single(snap, centre, name, i)
            if c is not None:
                out.append((abs(i - p) / max(1, len(snap.ladders.rungs[name]) - 1), c))
    return [c for _, c in sorted(out, key=lambda t: t[0])]


def riders(ctx: RFContext, c: Config, k: int) -> List[Config]:
    out = []
    if CONF not in ctx.space.names:
        return out
    for v in CALL_CONF_GRID:
        if len(out) >= k:
            break
        if v != c.params.get(CONF):
            r = ctx.space._assemble({**c.params, CONF: v})
            if r is not None and ctx.fresh(r):
                out.append(r)
    return out


class Convergence:
    """The loop has converged: no gain of `gain` for `window` seconds, after `after` seconds. The loop has reached its
    summit (second fleet replay: from minute ~30, 82% of picks polished it for +0.009), so the endgame climbs start
    then instead of at rf.endgame_min. Counted in time, not results: riders make a fleet's worth of results about 40 s
    (a first try at that fired at minute 23 of a simulated round whose earlier runs still gained about 5 after it).
    Simulated rounds gain 0.02 at most 12 minutes apart.

    A restart starts a new window: on 175300Z (2026-10-09) a restart at minute 40 was cut off by convergence 1.2 minutes
    later, so it never had the chance to find another hill."""

    def __init__(self):
        self.mark, self.at, self.restarts = float("-inf"), 0.0, 0

    def check(self, best: float, restarts: int, now: float, after: float, window: float, gain: float) -> bool:
        if best >= self.mark + gain:
            self.mark, self.at = best, now
            return False
        if restarts != self.restarts:
            self.restarts, self.at = restarts, now
            return False
        return now >= after and now - self.at >= window


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: Optional[int] = None) -> Strategy:
    # A seed per run (second fleet replay: both runs of one round used the same 972 opening configs, so they were not
    # independent tests); MINOS_SEED=<n> repeats a run.
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
        ctx = RFContext(client, space, rng, profile)
        ctx.cap = await settled_capacity(client, workers)
        virtual = hasattr(client.executor, "now")
        ctx.say("%s RF_UCB: %d workers, %d settings searched, a %d-minute round (schedule x%.2f), seed %d", VERSION,
                ctx.cap, len(space.names), round(ctx.length / 60), ctx.scale, seed)
        loop = Loop(ctx)
        if not virtual:                                # every pick's prediction next to its result, for the review
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
                o = f.result()
                if o is None or not o.ok or o.result.gate_passed is False or not o.score:
                    if o is not None:
                        loop.sent.pop(o.config.key, None)
                        err = o.error or ""
                        if not err.startswith(("invalid", "budget", "round ended")):
                            loop.bad[o.config.key] = o.config  # a final failure or a zero score
                    continue
                loop.note_back(o.config, o.score)
                if o.score > loop.best:
                    loop.wins[src] += 1
                    loop.best = o.score
            loop.results += len(done)
            loop.loop_s += time.time() - t
            ctx.publish()
            return len(done)

        # ---- opening ---------------------------------------------------------------------------------------------------
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
        def endgame() -> float:                        # read at every check: a change on the dashboard takes effect
            return ctx.t0 + ctx.T(60 * float(tune("rf.endgame_min")))
        while pending and n_open - len(pending) < float(tune("rf.opening_back")) * n_open and ctx.clock() < endgame() \
                and not ctx.over():
            await drain(TICK_S)
        ctx.marks["a"] = ctx.elapsed()
        top = ctx.best()
        ctx.say("opening: %d of %d back at t+%.1fm, best %.4f", n_open - len(pending), n_open, ctx.elapsed() / 60,
                top[1] if top else float("nan"))

        # ---- RF_UCB ---------------------------------------------------------------------------------------------------
        ctx.phase = "B: RF_UCB"
        conv = Convergence()

        def converged() -> bool:
            return conv.check(loop.best, loop.restarts, ctx.elapsed(),
                              ctx.T(60 * float(tune("rf.converge_after_min"))), ctx.T(60 * float(tune("rf.converge_min"))),
                              float(tune("rf.converge_gain")))

        def refit_n() -> int:                          # live, as the setting says
            return max(8, ctx.cap // int(tune("rf.refit_per_fleet")))
        since = refit_n()
        task: Optional[asyncio.Task] = None
        headroom = max(4, ctx.cap // 20)

        def start_refit():
            pts = loop.points()
            t = time.time()
            if pts:
                loop.restart_check(max(p[1] for p in pts), len(pts))
            loop.centres = loop.pick_centres(pts)
            loop.loop_s += time.time() - t
            return pts, list(loop.centres), list(loop.bad.values())

        try:                                           # the first refit before any pick: no pool means random picks
            loop.apply(loop.refit(*start_refit()))
            since = 0
        except Exception:  # noqa: BLE001
            log.exception("engine: first refit failed")
        while ctx.clock() < endgame() and not ctx.over():
            try:                                       # a bug in a step must not end the round's proposing
                if task is None and (since >= refit_n() or not loop.pool):
                    pts, centres, bad = start_refit()
                    since = 0
                    if virtual:                        # deterministic: the virtual clock waits for the fit
                        loop.apply(loop.refit(pts, centres, bad))
                    else:                              # the loop keeps filling workers while the forest is fitted
                        task = asyncio.create_task(asyncio.to_thread(loop.refit, pts, centres, bad))
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
                log.exception("engine: RF_UCB step failed")
                task = None if task is not None and task.done() else task
            got = await drain(TICK_S)
            since += got
            loop.stall += got
            if converged():
                ctx.say("RF_UCB converged at t+%.1fm: no gain of %.3f since t+%.1fm (best %.4f); the endgame starts now",
                        ctx.elapsed() / 60, float(tune("rf.converge_gain")), conv.at / 60, loop.best)
                break
        task = None                                    # a refit still running is left to finish unused: no waiting
        ctx.marks["b"] = ctx.elapsed()
        ctx.rf.update(loop.stats())
        ctx.say("RF_UCB done at t+%.1fm: %s", ctx.elapsed() / 60, json.dumps({k: v for k, v in loop.stats().items()
                                                                            if k != "centres"}))

        # ---- endgame: v9's grid climb with lines, then Step C, on the best candidates ---------------------------------
        ctx.phase = "C: endgame climb"
        # The endgame's waves go ahead of the loop's leftovers in every queue (first fleet replay: its first grid wave
        # took 15 minutes behind them, including a Reinitialize's 972 random picks).
        ctx.priority = int(tune("rf.endgame_priority"))
        loop.reinit_left = 0
        # The loop's picks no worker has started are dropped: they would only delay the endgame's waves on a busy
        # fleet. Running ones finish and still count.
        dropped = client.drop_queued({"opening", "qga", "ga", "local", "mutate", "sweep", "inject"})
        if dropped:
            ctx.say("endgame: %d of the loop's queued picks dropped; %d jobs still in flight", dropped, client.in_flight)
        pts = loop.points()
        hills = loop.pick_centres(pts, int(tune("rf.endgame_hills")))      # best score, best core, best FP
        ends = {c.key: label for c, label in proposers.trade_off_ends(pts, float(tune("rf.ends_window")))}
        for c in hills:
            cl = Climber(ctx, len(ctx.climbers), c, ctx.score(c) or 0.0, origin=f"RF_UCB {ends.get(c.key, 'best distinct')}")
            cl.rider_cap = int(tune("rf.endgame_riders"))   # fewer riders per grid point: a wave ends sooner
            ctx.climbers.append(cl)

        async def lifecycle(cl: Climber) -> None:
            try:
                await cl.run(ctx.pick_at())
                if cl.state in ("peak", "stopped") and not ctx.over():
                    await far.improve(ctx, cl, ctx.pick_at())
            except Exception:  # noqa: BLE001
                log.exception("engine: endgame climber %d failed", cl.id)
            ctx.publish(force=True)

        await asyncio.gather(*(lifecycle(cl) for cl in list(ctx.climbers)))
        while ctx.clock() < ctx.pick_at() and not ctx.over():   # hills done early: climb the best config off them
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
