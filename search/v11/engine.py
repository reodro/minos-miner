"""Search engine v11: RF_UCB over several hills at once (search/v11/DESIGN.md).

v10 collapsed onto one hill: on 20261009T153000Z, 48 of its best 50 configs were on the final hill by minute 20 and all
50 by minute 30; every config above 86.3 (24K) was that hill. Every proposer except Inject followed the global best,
the "centres" were the ends of one ridge, a restart kept the same forest and elites, and the endgame climbed that
ridge. Local tests around its best (test/compound_het/README.md) found nothing better within reach, while the field
sat 0.056 higher somewhere else. v11 keeps v10's forest, proposers, riders and endgame climber, and adds the hills:

    hills       at every refit the best hill.points results are grouped by what they call: a config joins the first
                hill whose best (its founder) is within hill.callset_far (0.03) by v10's callset distance (class-weighted
                1 - Jaccard of the sketches), else founds a new hill if within hill.window of the best. Hills keep
                their identity across refits (a new founder near an old one is the same hill) and their own clock.
    valleys     hills with no valley between them are merged (Hills.merge_valleys: a hill-valley test on the configs
                already known between their bests and on the cross midpoints), so a small radius does not fund one hill
                twice; pairs with no verdict yet get their midpoints first in the cross share.
    funding     hill.count hills are funded: the best by score, except that a hill hill.drop below the leader with no
                gain for hill.patience_min gives its place to the next. Each funded hill has its own pool, ranked by
                the forest: GA from its own members, Local around its own trade-off ends, the sweep through its own
                best, mutations of its own top proposals.
    picks       rf.share_inject a hypercube config; hill.share_qga QGA (rotated toward every funded hill: new hills);
                hill.share_cross a child of two funded hills' bests, mixing whole stage blocks (the briefing: the
                largest single source of gains in v4.1); the rest go to the funded hills, the leader hill.lead_share
                and the others equal shares. Within a hill: hill.share_sweep the sweep, hill.share_mutate a mutation,
                else its ranked pool.
    converged   every funded hill has had no gain of rf.converge_gain for rf.converge_min (each on its own clock; a
                new hill starts its clock when founded), after rf.converge_after_min; else the endgame at rf.endgame_min.
    endgame     v9's grid climb with lines and Step C on the best of each funded hill within hill.endgame_window of
                the leader (at most hill.endgame).
    pick        the highest local score (unchanged).

QLock and Reinitialize are off: a setting flat at one hill's best is not flat at another's, and the hills' own
exploration (QGA, Inject, cross) replaces the random restart. A stall still resets QGA.
"""
from __future__ import annotations

import asyncio
import itertools
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
from common.types import Config
from search import portfolio
from search.bamgroups import bam_bounds
from search.base import Strategy
from search.client import EvalClient
from search.space import Space, bam_caps
from search.v9 import far as far_step
from search.v9.candidates import OFFICIAL_SEED
from search.v9.climb import Climber
from search.v9.distance import _unpack, callset_distance, setting_distance
from search.v9.engine import TICK_S, settled_capacity
from search.v10 import proposers
from search.v10.engine import Loop, Point, Refit, RFContext, Snapshot, riders, sweep_lines

log = logging.getLogger(__name__)

VERSION = "v11"
# Whole stage blocks for cross-hill children (the GATK stages, in run order).
BLOCKS = {
    "read": ("min_base_quality_score", "min_mapping_quality_score", "max_reads_per_alignment_start",
             "assembly_region_padding"),
    "activity": ("active_probability_threshold", "min_assembly_region_size", "max_assembly_region_size"),
    "assembly": ("min_pruning", "min_dangling_branch_length", "recover_all_dangling_branches",
                 "max_num_haplotypes_in_population"),
    "pairhmm": ("base_quality_score_threshold", "pcr_indel_model", "pair_hmm_gap_continuation_penalty",
                "phred_scaled_global_read_mismapping_rate"),
    "genotyping": ("standard_min_confidence_threshold_for_calling", "heterozygosity", "indel_heterozygosity",
                   "max_alternate_alleles", "contamination_fraction_to_filter"),
}
LOOP_TAGS = {"opening", "qga", "ga", "local", "mutate", "sweep", "inject", "cross"}


@dataclass
class Hill:
    id: int
    founder: Config                                # its best config
    best: float
    born: float                                    # elapsed seconds
    mark: float                                    # its best at its last gain of rf.converge_gain
    last_gain: float                               # elapsed seconds of that gain
    members: List[Point] = field(default_factory=list)   # its configs among the grouped results, best first

    def row(self, now: float) -> dict:
        return {"id": self.id, "best": round(self.best, 4), "founder": self.founder.key, "members": len(self.members),
                "born_m": round(self.born / 60, 1), "last_gain_m": round(self.last_gain / 60, 1),
                "quiet_m": round((now - self.last_gain) / 60, 1)}


def valley_verdict(s_a: float, s_b: float, between: Sequence[float], drop: float, min_n: int) -> str:
    """Hill-valley test (HillVallEA): two hills are one when some config between them scores within `drop` of the
    lower of the two (no valley on that path); two when `min_n` or more configs between them all score lower than
    that; else unknown."""
    lo = min(s_a, s_b) - drop
    if any(x >= lo for x in between):
        return "same"
    if len(between) >= min_n:
        return "valley"
    return "unknown"


class Hills:
    """The hills among the scored configs, kept across refits. Grouped by what the configs call, not their settings:
    v10's callset distance (search/v9/distance.py: class-weighted 1 - Jaccard from the bottom-k sketches every result
    carries), so settings that change no call on this round (a range edge, a BAM equivalence class) do not split a
    hill and a one-rung move that changes many calls does. A config with no sketch (a simulator) falls back to the mean
    normalised setting difference, against its own threshold. Runs on the refit thread from a snapshot of the
    sketches taken on the event loop; distances are cached per pair.

    The radius (hill.callset_far) is set small on purpose: a hill split in two is still climbed, a real hill merged
    into another is lost. Duplicates are then merged by a hill-valley test on what is already known between two hills'
    founders: configs whose calls lie between them (closer to each founder than the founders are to each other, and
    at least hill.valley_mid of that distance from both), plus the cross children that take 2-3 of the 5 stage blocks
    from each side. No extra evaluations: both are evaluated anyway."""

    def __init__(self, space):
        self.space = space
        self.hills: Dict[int, Hill] = {}
        self.ids = itertools.count()
        self.funded: List[Hill] = []
        self._dec: Dict[str, Optional[dict]] = {}
        self._dist: Dict[Tuple[str, str], Tuple[float, bool]] = {}     # pair -> (distance, by callset)
        self.mid: Dict[str, Tuple[Config, frozenset]] = {}             # cross midpoint key -> (config, its hill pair)
        self.verdicts: Dict[frozenset, str] = {}
        self.undecided: List[frozenset] = []
        self.events: List[str] = []

    def _sketch(self, key: str, raw: Dict[str, Optional[dict]]) -> Optional[dict]:
        if key not in self._dec:
            sk = raw.get(key)
            self._dec[key] = None if not sk else {"k": sk.get("k", 128),
                                                  "h": {cls: _unpack(v) for cls, v in sk["h"].items()}}
        return self._dec[key]

    def dist(self, a: Config, b: Config, raw: Dict[str, Optional[dict]]) -> Tuple[float, bool]:
        k = (a.key, b.key) if a.key < b.key else (b.key, a.key)
        got = self._dist.get(k)
        if got is None:
            sa, sb = self._sketch(a.key, raw), self._sketch(b.key, raw)
            if sa is not None and sb is not None:
                got = (callset_distance(sa, sb), True)
            else:
                got = (setting_distance(self.space, a, b), False)
            if len(self._dist) > 400_000:
                self._dist.clear()
            self._dist[k] = got
        return got

    def near(self, a: Config, b: Config, raw: Dict[str, Optional[dict]]) -> bool:
        d, by_calls = self.dist(a, b, raw)
        return d < float(tune("hill.callset_far" if by_calls else "hill.setting_far"))

    def note_mid(self, child: Config, a: int, b: int) -> None:
        """A cross child taking 2-3 of the 5 stage blocks from each of hills a and b: a midpoint for their valley test."""
        if len(self.mid) > 20_000:
            self.mid.clear()
        self.mid[child.key] = (child, frozenset((a, b)))

    def update(self, pts: Sequence[Point], raw: Dict[str, Optional[dict]], now: float, gain: float) -> None:
        if not pts:
            return
        ranked = sorted(pts, key=lambda p: -p[1])[:int(tune("hill.points"))]
        self._dec = {k: v for k, v in self._dec.items() if k in raw}   # decoded sketches: this refit's configs only
        window, track = float(tune("hill.window")), int(tune("hill.track"))
        top = ranked[0][1]
        founders: List[int] = []                   # indices into ranked
        members: Dict[int, List[int]] = {}
        for i, p in enumerate(ranked):
            home = next((f for f in founders if self.near(ranked[f][0], p[0], raw)), None)
            if home is not None:
                members[home].append(i)
                continue
            if p[1] >= top - window and len(founders) < track:
                founders.append(i)
                members[i] = [i]
        # identity: a founder near an old hill's founder is that hill (the best unmatched one)
        old = sorted(self.hills.values(), key=lambda h: -h.best)
        used: set = set()
        new: Dict[int, Hill] = {}
        for i in founders:
            c, s = ranked[i][0], ranked[i][1]
            match = next((h for h in old if h.id not in used and self.near(h.founder, c, raw)), None)
            if match is None:
                h = Hill(next(self.ids), c, s, now, s, now)
            else:
                used.add(match.id)
                h = match
                h.founder, h.best = c, s
                if s >= h.mark + gain:
                    h.mark, h.last_gain = s, now
            h.members = [ranked[m] for m in members[i]]
            new[h.id] = h
        self.hills = new
        self.merge_valleys(ranked, {p[0].key: p[1] for p in pts}, raw)

    def merge_valleys(self, ranked: Sequence[Point], score: Dict[str, float], raw: Dict[str, Optional[dict]]) -> None:
        """Merge hills with no valley between them (valley_verdict), the lower into the higher; note the pairs still
        unknown, whose midpoints the cross share sends first."""
        drop, min_n, mid = float(tune("hill.valley_drop")), int(tune("hill.valley_min")), float(tune("hill.valley_mid"))
        tops = sorted(self.hills.values(), key=lambda h: -h.best)[:int(tune("hill.valley_hills"))]
        mids_by_pair: Dict[frozenset, List[float]] = {}
        for key, (c, pair) in self.mid.items():
            if key in score:
                mids_by_pair.setdefault(pair, []).append(score[key])
        self.undecided = []
        gone: set = set()
        for a, b in itertools.combinations(tops, 2):
            if a.id in gone or b.id in gone:
                continue
            pair = frozenset((a.id, b.id))
            d_ab, by_calls = self.dist(a.founder, b.founder, raw)
            between: List[float] = list(mids_by_pair.get(pair, []))
            if by_calls and d_ab > 0:
                for p in ranked:
                    c = p[0]
                    if c.key in (a.founder.key, b.founder.key):
                        continue
                    da, ka = self.dist(a.founder, c, raw)
                    db, kb = self.dist(b.founder, c, raw)
                    if ka and kb and da < d_ab and db < d_ab and min(da, db) >= mid * d_ab:
                        between.append(p[1])
            verdict = valley_verdict(a.best, b.best, between, drop, min_n)
            said = frozenset((a.founder.key, b.founder.key))       # logged once per pair of bests, not per refit
            if self.verdicts.get(said) != verdict and verdict != "unknown":
                self.events.append(f"h{a.id} {a.best:.4f} / h{b.id} {b.best:.4f}: "
                                   + ("no valley, merged" if verdict == "same" else "valley, kept apart")
                                   + f" ({len(between)} between, best {max(between, default=float('nan')):.4f}, "
                                     f"distance {d_ab:.3f})")
            self.verdicts[said] = verdict
            if verdict == "same":
                hi, lo = (a, b) if a.best >= b.best else (b, a)
                hi.members = sorted(hi.members + lo.members, key=lambda p: -p[1])
                hi.last_gain = max(hi.last_gain, lo.last_gain)
                gone.add(lo.id)
            elif verdict == "unknown":
                self.undecided.append(pair)
        for i in gone:
            self.hills.pop(i, None)

    def fund(self, now: float, count: int, drop: float, patience: float) -> List[Hill]:
        """The best `count` hills by score, except those `drop` below the leader with no gain for `patience`."""
        ranked = sorted(self.hills.values(), key=lambda h: -h.best)
        if not ranked:
            self.funded = []
            return self.funded
        lead = ranked[0].best
        keep = [h for h in ranked if h is ranked[0] or not (h.best < lead - drop and now - h.last_gain >= patience)]
        self.funded = keep[:count]
        return self.funded


@dataclass
class HillRefit(Refit):
    funded: List[Hill] = field(default_factory=list)
    hill_q: Dict[int, dict] = field(default_factory=dict)      # hill id -> {"pool", "top", "sweep"}
    cross: Deque[Config] = field(default_factory=deque)


class HillContext(RFContext):
    """v10's context; the notes carry v11's version and the hills."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.hill_rows: List[dict] = []

    def publish(self, force: bool = False) -> None:
        before = self.last_publish
        super().publish(force)
        eng = self.client.notes.get("engine")
        if eng is not None and (force or self.last_publish != before):
            eng["version"] = VERSION
            eng["hills"] = list(self.hill_rows)


class HillLoop(Loop):
    """v10's loop with a pool per funded hill, cross-hill children, and picks funded per hill."""

    def __init__(self, ctx: HillContext):
        super().__init__(ctx)
        self.hills = Hills(ctx.space)
        self.funded: List[Hill] = []
        self.hill_q: Dict[int, dict] = {}
        self.cross_q: Deque[Config] = deque()
        self.hill_picks: Counter = Counter()

    def qlock(self, centre: Config, n: int) -> None:   # off: flat at one hill's best is not flat at another's
        return

    def restart_check(self, best_s: float, n_results: int) -> None:
        """A stall resets QGA only (no Reinitialize: Inject, QGA and the cross children keep exploring)."""
        ctx = self.ctx
        self.best = max(self.best, best_s)
        if self.best >= self.gain_mark + float(tune("rf.restart_min_gain")):
            self.gain_mark, self.stall = self.best, 0
            return
        if self.stall < float(tune("rf.restart_waves")) * ctx.cap:
            return
        self.stall = 0
        self.qga.reset()
        self.restarts += 1
        ctx.say("RF_UCB: restart %d (no gain of %.3f since %.4f): QGA reset", self.restarts,
                float(tune("rf.restart_min_gain")), self.gain_mark)

    # ---- the refit (a worker thread on a live fleet) ----------------------------------------------------------------
    def sketches(self, pts: Sequence[Point]) -> Dict[str, Optional[dict]]:
        """On the event loop: the raw sketches of the configs the next refit groups (the client's dicts keep changing;
        the refit thread decodes them)."""
        from common.types import FULL
        ranked = [p[0] for p in sorted(pts, key=lambda p: -p[1])[:int(tune("hill.points"))]]
        ranked += [h.founder for h in self.hills.hills.values()]   # an old hill's founder may have left the top
        ranked += [c for c, _ in list(self.hills.mid.values())]       # midpoints: low scores are the valleys
        out = {}
        for c in ranked:
            cs = self.ctx.client.callsets.get((c.key, FULL))
            out[c.key] = getattr(cs, "sketch", None) if cs is not None else None
        return out

    def refit(self, pts: Sequence[Point], sketches: Dict[str, Optional[dict]], bad: Sequence[Config] = ()) -> HillRefit:
        ctx, rng = self.ctx, self.ctx.rng
        t0 = time.time()
        out = HillRefit(n=len(pts))
        if len(pts) < 8:
            return out
        have = {p[0].key for p in pts}
        self.forest.fit([p[0] for p in pts], [p[1] for p in pts], float(tune("rf.density_radius")))
        self.forest.fit_bad([p[0] for p in pts], [c for c in bad if c.key not in have])
        out.fit_s = time.time() - t0
        self.refits += 1
        now = ctx.elapsed()
        self.hills.update(pts, sketches or {}, now, float(tune("rf.converge_gain")))
        funded = self.hills.fund(now, int(tune("hill.count")), float(tune("hill.drop")),
                                 ctx.T(60 * float(tune("hill.patience_min"))))
        out.funded = list(funded)
        best_c = max(((c, s) for c, s, _, _ in pts), key=lambda t: t[1])[0]
        live = list(self.names)
        size = int(tune("rf.pool"))
        progress = min(1.0, ctx.elapsed() / max(1.0, ctx.T(60 * float(tune("rf.endgame_min")))))
        theta = float(tune("rf.qga_theta")) + (float(tune("rf.qga_theta_end")) - float(tune("rf.qga_theta"))) * progress
        leaders = [h.founder for h in funded] or [best_c]
        self.qga.rotate(leaders, theta, float(tune("rf.qga_mut")), float(tune("rf.qga_floor")))
        made: List[Tuple[Config, str, int]] = []      # config, source, hill id (-1: global)
        made += [(c, "qga", -1) for c in self.qga.sample(int(size * float(tune("rf.share_qga")) / 2), best_c, live)]
        per = max(32, int(size * float(tune("rf.share_ga")) / max(1, len(funded))))
        snap = Snapshot(ctx, pts)
        sweeps: Dict[int, Deque[Config]] = {}
        for h in funded:
            mem = h.members or [(h.founder, h.best, 0.0, 0.0)]
            elites = [(c, s) for c, s, _, _ in mem[:int(tune("rf.elites"))]]
            made += [(c, "ga", h.id) for c in proposers.breed(ctx.space, ctx.ladders, rng, elites, per, h.founder, live,
                                                              float(tune("rf.ga_mut_rate")),
                                                              int(tune("rf.ga_mut_rungs")))]
            ends = [c for c, _ in proposers.trade_off_ends(mem, float(tune("rf.ends_window")))] or [h.founder]
            made += [(c, "local", h.id) for c in proposers.local(snap, rng, ends, int(tune("far.coarse")),
                                                                 int(tune("rf.local_moves")))]
            sweeps[h.id] = deque(c for c in sweep_lines(snap, h.founder, int(tune("far.coarse"))) if ctx.fresh(c))
        mids: set = set()
        for a, b in itertools.permutations(funded, 2):  # a's settings with b's blocks
            for mask in itertools.product((0, 1), repeat=len(BLOCKS)):
                if not any(mask) or all(mask):
                    continue
                p = dict(a.founder.params)
                for blk, m in zip(BLOCKS, mask):
                    if m:
                        for s in BLOCKS[blk]:
                            if s in b.founder.params:
                                p[s] = b.founder.params[s]
                c = ctx.space._assemble(p)
                if c is not None:
                    made.append((c, "cross", -2))
                    if 2 <= sum(mask) <= 3:            # a midpoint: 2-3 of the 5 blocks from each side
                        self.hills.note_mid(c, a.id, b.id)
                        if frozenset((a.id, b.id)) in self.hills.undecided:
                            mids.add(c.key)
        seen, fresh = set(), []
        for c, src, hid in made:
            if c.key not in seen and ctx.fresh(c):
                seen.add(c.key)
                fresh.append((c, src, hid))
        if fresh:
            u, mu, sd, nov = self.forest.ucb([c for c, _, _ in fresh], float(tune("rf.kappa")),
                                             float(tune("rf.kappa_dist")))
            order = np.argsort(-u)
            out.pred = {c.key: (float(mu[i]), float(sd[i]), float(nov[i])) for i, (c, _, _) in enumerate(fresh)}
            out.pool = deque((fresh[int(i)][0], fresh[int(i)][1]) for i in order)    # everything: the fallback
            out.qga = deque(fresh[int(i)][0] for i in order if fresh[int(i)][1] == "qga")
            cross = [fresh[int(i)][0] for i in order if fresh[int(i)][1] == "cross"]
            # midpoints of hill pairs with no verdict yet go first: they decide whether two hills are one
            out.cross = deque([c for c in cross if c.key in mids] + [c for c in cross if c.key not in mids])
            for h in funded:
                ranked = [(fresh[int(i)][0], fresh[int(i)][1]) for i in order if fresh[int(i)][2] == h.id]
                out.hill_q[h.id] = {"pool": deque(ranked), "top": [c for c, _ in ranked[:int(tune("rf.mutate_top"))]],
                                    "sweep": sweeps.get(h.id, deque())}
            out.top = [fresh[int(i)][0] for i in order[:int(tune("rf.mutate_top"))]]
        out.inject = deque(c for c in proposers.inject(ctx.space, rng, int(size * float(tune("rf.share_inject"))),
                                                       best_c, live) if ctx.fresh(c))
        out.sweep = sweeps.get(funded[0].id, deque()) if funded else deque()
        out.total_s = time.time() - t0
        return out

    def apply(self, r: Refit) -> None:
        super().apply(r)
        if isinstance(r, HillRefit):
            self.funded, self.hill_q, self.cross_q = r.funded, r.hill_q, r.cross
            self.centres = [h.founder for h in self.funded]
            now = self.ctx.elapsed()
            self.ctx.hill_rows = [dict(h.row(now), funded=True, picks=self.hill_picks[h.id]) for h in self.funded]
            for e in self.hills.events:
                self.ctx.say("hills: %s", e)
            self.hills.events = []
            if self.refits % 10 == 1 and self.funded:
                self.ctx.say("hills: %s", "; ".join(f"h{h.id} {h.best:.4f} ({len(h.members)} configs, quiet "
                                                     f"{(now - h.last_gain) / 60:.1f}m)" for h in self.funded))

    # ---- picking (event loop) -------------------------------------------------------------------------------------
    def mutate_from(self, top: Sequence[Config]) -> Optional[Config]:
        ctx = self.ctx
        if not top:
            return None
        for _ in range(20):
            p = dict(ctx.rng.choice(top).params)
            for name in ctx.rng.sample(self.names, min(len(self.names), ctx.rng.randint(1, 3))):
                p[name] = proposers.shift(ctx.ladders, name, p[name], ctx.rng.randint(1, 3) * ctx.rng.choice((1, -1)),
                                          ctx.space)
            c = ctx.space._assemble(p)
            if c is not None and ctx.fresh(c):
                return c
        return None

    def take(self, n: int) -> List[Tuple[Config, str]]:
        """n configs: Inject (rf.share_inject), QGA (hill.share_qga), a cross-hill child (hill.share_cross), else a
        funded hill (the leader hill.lead_share, the others equal shares): its sweep (hill.share_sweep), a mutation of
        its top proposals (hill.share_mutate), else its ranked pool. Anything empty falls back to the global ranking."""
        ctx = self.ctx
        if not self.funded:
            return super().take(n)
        s_i, s_q, s_x = float(tune("rf.share_inject")), float(tune("hill.share_qga")), float(tune("hill.share_cross"))
        s_sw, s_mu = float(tune("hill.share_sweep")), float(tune("hill.share_mutate"))
        lead = float(tune("hill.lead_share"))
        k = len(self.funded)
        weights = [lead] + [(1 - lead) / (k - 1)] * (k - 1) if k > 1 else [1.0]
        out: List[Tuple[Config, str]] = []
        taken: set = set()
        inject_cap = max(1, int(np.ceil(s_i * n)))
        n_inject = 0

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
        while len(out) < n and tries < 4 * n + 16:
            tries += 1
            roll = ctx.rng.random()
            src, c, hid = None, None, None
            if roll < s_i:
                if n_inject >= inject_cap:
                    continue
                src, c = "inject", self.next_inject(n)
                n_inject += 1
            elif roll < s_i + s_q:
                src, c = "qga", first(self.qga_q)
            elif roll < s_i + s_q + s_x:
                src, c = "cross", first(self.cross_q)
            else:
                h = ctx.rng.choices(self.funded, weights=weights)[0]
                hid, q = h.id, self.hill_q.get(h.id) or {}
                r2 = ctx.rng.random()
                if r2 < s_sw:
                    src, c = "sweep", first(q.get("sweep", deque()))
                elif r2 < s_sw + s_mu:
                    src, c = "mutate", self.mutate_from(q.get("top", []))
                if not ok(c):
                    pool = q.get("pool", deque())
                    while pool:
                        c2, src2 = pool.popleft()
                        if ok(c2):
                            src, c = src2, c2
                            break
                    else:
                        c = None
            if not ok(c):                                  # nothing there: the global ranking
                hid = None
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
                if hid is not None:
                    self.hill_picks[hid] += 1
        return out

    def stats(self) -> dict:
        out = super().stats()
        now = self.ctx.elapsed()
        out["hills"] = [dict(h.row(now), picks=self.hill_picks[h.id]) for h in self.funded]
        out["hills_tracked"] = len(self.hills.hills)
        return out


def endgame_hills(funded: Sequence[Hill], count: int, window: float) -> List[Hill]:
    """The leader and every funded hill within `window` of it, best first, at most `count`."""
    ranked = sorted(funded, key=lambda h: -h.best)
    if not ranked:
        return []
    return [h for h in ranked if h.best >= ranked[0].best - window][:count]


def converged(funded: Sequence[Hill], now: float, after: float, window: float) -> bool:
    """Every funded hill has had no gain for `window` (each on its own clock), after `after`."""
    return bool(funded) and now >= after and all(now - h.last_gain >= window for h in funded)


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
        ctx = HillContext(client, space, rng, profile)
        ctx.cap = await settled_capacity(client, workers)
        virtual = hasattr(client.executor, "now")
        ctx.say("%s RF_UCB over hills: %d workers, %d settings searched, a %d-minute round (schedule x%.2f), seed %d, "
                "%d hills funded", VERSION, ctx.cap, len(space.names), round(ctx.length / 60), ctx.scale, seed,
                int(tune("hill.count")))
        loop = HillLoop(ctx)
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

        # ---- RF_UCB over hills ---------------------------------------------------------------------------------------
        ctx.phase = "B: RF_UCB hills"

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
            sk = loop.sketches(pts)
            loop.loop_s += time.time() - t
            return pts, sk, list(loop.bad.values())

        try:
            loop.apply(loop.refit(*start_refit()))
            since = 0
        except Exception:  # noqa: BLE001
            log.exception("engine: first refit failed")
        while ctx.clock() < endgame() and not ctx.over():
            try:
                if task is None and (since >= refit_n() or not loop.pool):
                    pts, centres, bad = start_refit()
                    since = 0
                    if virtual:
                        loop.apply(loop.refit(pts, centres, bad))
                    else:
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
            if converged(loop.funded, ctx.elapsed(), ctx.T(60 * float(tune("rf.converge_after_min"))),
                         ctx.T(60 * float(tune("rf.converge_min")))):
                ctx.say("RF_UCB converged at t+%.1fm: no funded hill gained %.3f in %.0f minutes (%s); the endgame "
                        "starts now", ctx.elapsed() / 60, float(tune("rf.converge_gain")),
                        ctx.T(60 * float(tune("rf.converge_min"))) / 60,
                        ", ".join(f"h{h.id} {h.best:.4f}" for h in loop.funded))
                break
        task = None
        ctx.marks["b"] = ctx.elapsed()
        ctx.rf.update(loop.stats())
        ctx.say("RF_UCB done at t+%.1fm: %s", ctx.elapsed() / 60, json.dumps({k: v for k, v in loop.stats().items()
                                                                            if k != "centres"}))

        # ---- endgame: v9's grid climb with lines, then Step C, on each funded hill's best --------------------------------
        ctx.phase = "C: endgame climb"
        ctx.priority = int(tune("rf.endgame_priority"))
        loop_tags = {tag for (_f, _c, tag) in list(client._pending.values()) if tag in LOOP_TAGS}
        dropped = client.drop_queued(loop_tags)
        if dropped:
            ctx.say("endgame: %d of the loop's queued picks dropped; %d jobs still in flight", dropped, client.in_flight)
        hills = endgame_hills(loop.funded, int(tune("hill.endgame")), float(tune("hill.endgame_window")))
        if not hills:
            top = ctx.best()
            starts = [(top[0], top[1], "best")] if top else []
        else:
            starts = [(h.founder, ctx.score(h.founder) or h.best, f"hill {h.id}") for h in hills]
        for c, s, label in starts:
            cl = Climber(ctx, len(ctx.climbers), c, s, origin=f"RF_UCB {label}")
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
