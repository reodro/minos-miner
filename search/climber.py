"""A domain's climb: iterated local search (ParamILS) on one hill. Owner: search.

The domain only proposes moves and takes in results; it never decides its own lifecycle (founding, merging,
dropping are the controller's, engine.py). It reports facts: what is in flight, its progress rate, `state`, and
notes for the dashboard's event feed (`events`).

  scan         a domain founded fresh (not by a merge) first scans every setting over its whole range:
               climb.scan_points points per setting from the centre (the legacy pipeline's "injection", +0.2 on its
               rounds, 7x more top-30 configs per evaluation than the opening wave)
  sweep        the neighbourhood of the centre. COARSE: one setting at a time at its learned step and at half and
               double it, plus pairs of the strongest settings. When a coarse sweep comes back without a gain, a WIDE
               round first: all pairs of the WIDE_AXES strongest settings at two step sizes (WIDE_STEPS of the
               range; at most climb.wide_max configs - the legacy pipeline's all-pairs rung). Then FINE: steps halved
               and never doubled again, settings whose step changed nothing are left out until the centre moves,
               triples join the pairs, and each exhausted fine sweep halves the scale once more
  acceptance   the centre moves to the best result of the sweep once most of it is back (a decisive gain at once); the
               sweep's other winners are re-applied on the new centre first. An exact tie is never a move: identical
               scores mean identical calls (the setting changed nothing here)
  perturbation when a fine sweep is exhausted at every setting's smallest step: shake K settings of the centre (K =
               climb.shake_k in turn), then one local round around the best shaken point; the result replaces the
               centre if it is better OR EQUAL (ParamILS's acceptance). Round after round - it never gives up
  chains       parallel ParamILS walkers on the workers the main line leaves idle, up to climb.max_chains: each
               shakes 2-3 settings of its own incumbent, climbs one local round around the shaken point
               (climb.chain_local moves), keeps the result when it is better or equal, and repeats; sometimes
               (climb.chain_restart) or when it falls behind, it restarts from the domain's centre. A chain result
               better than the domain's best by the resolution moves the domain's centre. In the final the shakes
               are 1-2 settings and smaller
  final        the last minutes on this one point: single steps at halving scales and all pairs of the top settings,
               chains until the end
  siblings     a structural move carries one genotyping and one PairHMM variation of itself (the question "does this
               move work better with another threshold / PCR model"), packed into the move's GATK run
"""
from __future__ import annotations

import itertools
import random
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

from common.types import Config
from common.tuning import tune
from search.axes import AxisBook
from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS, Space

# Tunable (common/tuning.py): climb.accept_min (the answer set's resolution between neighbours: a smaller gain is
# not a gain), climb.decisive (a gain that moves the centre at once), climb.back_frac (a set of moves counts as back
# at this share), climb.equal_accepts, climb.shake_k, climb.shake_points, climb.scan_points, climb.wide_max,
# climb.max_chains, climb.chain_local, climb.chain_restart, riders.siblings.
EQUAL = 1e-9
COARSE_SCALES = (1.0, 0.5, 2.0)
PAIR_AXES, TRIPLE_AXES, TRIPLES_MAX = 8, 5, 40
WIDE_AXES, WIDE_STEPS = 12, (0.04, 0.12)
SHAKE_SCALE = 0.3
CHAIN_KS, CHAIN_SCALE = (2, 3), 0.15
FINAL_CHAIN_KS, FINAL_CHAIN_SCALE = (1, 2), 0.06
STEP_MAX = 0.5
PROGRESS_WINDOW_S = 600.0
CHEAP = tuple(GENOTYPING_PARAMS) + tuple(PAIRHMM_PARAMS)
SWEEP_KINDS = ("scan", "n1", "n2", "n3", "wide", "reapply", "share", "transplant", "final", "pairing")
CHAIN_KINDS = ("chain", "chainl")


@dataclass
class Move:
    kind: str                # scan | n1 | n2 | n3 | wide | reapply | share | transplant | shake | local | final | pairing | chain | chainl
    config: Config
    axes: Tuple[str, ...]
    centre_key: str
    centre_score: float
    cycle: int = -1          # the perturbation step it belongs to
    chain: int = -1          # the chain it belongs to


@dataclass
class Chain:
    """One ParamILS walker of a domain."""
    id: int
    base: Config             # the chain's incumbent
    base_score: float
    stage: str = "shake"     # shake: its shaken point is out | local: the local round around that point is out
    sent: int = 0
    back: int = 0
    point: Optional[Tuple[Config, float]] = None
    best: Optional[Tuple[Config, float]] = None
    axes: Tuple[str, ...] = ()
    iters: int = 0


class Climber:
    def __init__(self, id: int, space: Space, book: AxisBook, centre: Config, score: float, now: float,
                 owner: Optional[Callable[[Config], Optional[int]]] = None, extra_axes: Sequence[str] = (),
                 scan: bool = False):
        self.id, self.space, self.book = id, space, book
        self.centre, self.best = centre, score
        self.birth_score = score
        self.first_scan_score: Optional[float] = None       # best point the first scan measured
        self.first_scan_key: Optional[str] = None
        self.first_scan_held: Optional[bool] = None         # True: that point did not move the centre
        self.born = now
        self.owner = owner                                 # config -> the domain whose hill it is on (None: no hill)
        self.steps: Dict[str, float] = {a: book.axes[a].step for a in book.names}
        for a in extra_axes:                               # the genotyping settings climb too (not in the axis book)
            s = space.specs.get(a)
            if s is not None and a in space.names and a not in self.steps:
                self.steps[a] = max(0.12, 1.0 / max(1, s.max - s.min)) if s.type == "int" else 0.12
        self.tried: Set[str] = {centre.key}
        self.pool: List[Move] = []
        self.inflight = 0
        self.state = "climb"                               # climb | perturb | final | exhausted
        self.mode = "coarse"                               # coarse | fine
        self.scale = 1.0                                   # fine-mode scale, halved per exhausted fine sweep
        self.fresh = True                                  # the centre's sweep is not sent yet
        self.scan_next = scan                              # the first sweep is a full-range scan
        self.sweep_is_scan = False                         # the sweep in flight is that scan
        self.scan_flat = False                             # the scan came back under climb.decisive: the peak is held
        self.wide_next = self.wide_done = False            # the wide round of this centre: queued / done
        self.eager_wide = False                            # set by the engine: few domains, so wide with the first sweep
        self.sweep_sent = self.sweep_back = 0
        self.winners: List[Tuple[Move, Config, float]] = []
        self.flat: Set[str] = set()
        self.recent: List[str] = []                        # settings of accepted moves, most recent last
        self.axis_gain: Dict[str, float] = {}
        self.queued: List[Tuple[str, Dict[str, object]]] = []    # (kind, changes) to try on the next sweep first
        self.cycle = 0                                     # index into climb.shake_k
        self.cycle_stage = "points"                        # points | local
        self.cycle_sent = self.cycle_back = 0
        self.cycle_best: Optional[Tuple[Config, float]] = None
        self.rounds = 0                                    # full perturbation rounds without an acceptance
        self.equal_accepts = 0
        self.chains: Dict[int, Chain] = {}
        self._chain_seq = 0
        self._live: Sequence[str] = ()                     # the live settings and rng of the last refill: results
        self._rng = random.Random(id)                      # start chain moves without waiting for a refill
        self.events: List[str] = []                        # notes for the event feed (the engine drains them)
        self.gains: List[Tuple[float, float]] = []         # (time, gain) of the domain's OWN accepted moves
        self.history: List[Tuple[float, float, str]] = [(now, score, "born")]

    # ---- facts for the controller --------------------------------------------------------------------------------
    def progress(self, now: float) -> float:
        """Own gain per minute over the last PROGRESS_WINDOW_S (or since birth, at least two minutes)."""
        span = max(120.0, min(PROGRESS_WINDOW_S, now - self.born))
        return sum(g for t, g in self.gains if t >= now - span) / (span / 60.0)

    def age(self, now: float) -> float:
        return now - self.born

    def note(self, msg: str) -> None:
        self.events.append(msg)
        del self.events[:-50]

    # ---- values ----------------------------------------------------------------------------------------------------
    def _floor(self, a: str) -> float:
        s = self.space.specs[a]
        return max(0.02, 1.0 / max(1, (s.max - s.min))) if s.type == "int" else 0.02

    def _allowed(self, a: str, v) -> bool:
        return a not in self.book.axes or self.book.allowed(a, v)

    def step_values(self, a: str, base: Config, scale: float = 1.0) -> List[object]:
        s = self.space.specs[a]
        cur = base.params.get(a, s.default)
        if s.allowed or s.type == "bool":
            return [v for v in (s.allowed or (True, False)) if v != cur]
        u = self.space._to_unit(s, cur)
        h = max(self._floor(a), self.steps[a] * scale)
        out = []
        for sgn in (1, -1):
            v = self.space._from_unit(s, min(1.0, max(0.0, u + sgn * h)))
            if v != cur and v not in out and self._allowed(a, v):
                out.append(v)
        return out

    def _variant(self, base: Config, changes: Dict[str, object]) -> Optional[Config]:
        c = self.space._assemble({**base.params, **changes})
        if c is None or c.key in self.tried or c.key == self.centre.key:
            return None
        if self.owner is not None and self.owner(c) not in (self.id, None):
            return None                                    # on another domain's hill
        return c

    def _add(self, kind: str, base: Config, base_score: float, changes: Dict[str, object], out: List[Move],
             cycle: int = -1, chain: int = -1) -> bool:
        c = self._variant(base, changes)
        if c is None:
            return False
        self.tried.add(c.key)
        out.append(Move(kind, c, tuple(sorted(changes)), base.key, base_score, cycle, chain))
        return True

    def strongest(self, live: Sequence[str], n: int) -> List[str]:
        ranked = sorted((a for a in live if a in self.axis_gain), key=lambda a: -self.axis_gain[a])[:n]
        for a in list(reversed(self.recent)) + list(live):
            if len(ranked) >= n:
                break
            if a in live and a not in ranked:
                ranked.append(a)
        return ranked

    # ---- neighbourhoods ----------------------------------------------------------------------------------------------
    def singles(self, base: Config, score: float, live: Sequence[str], kind: str, scale: float, skip=(),
                cycle: int = -1) -> List[Move]:
        out: List[Move] = []
        for a in live:
            if a in skip or a not in self.steps:
                continue
            for v in self.step_values(a, base, scale):
                self._add(kind, base, score, {a: v}, out, cycle)
        return out

    def pairs(self, live: Sequence[str], rng: random.Random, scale: float, kind: str = "n2") -> List[Move]:
        axes = [a for a in self.strongest(live, PAIR_AXES) if a in self.steps and self.step_values(a, self.centre, scale)]
        out: List[Move] = []
        for a, b in itertools.combinations(axes, 2):
            for x in self.step_values(a, self.centre, scale):
                for y in self.step_values(b, self.centre, scale):
                    self._add(kind, self.centre, self.best, {a: x, b: y}, out)
        rng.shuffle(out)
        return out

    def triples(self, live: Sequence[str], rng: random.Random, scale: float) -> List[Move]:
        axes = [a for a in self.strongest(live, TRIPLE_AXES) if a in self.steps and self.step_values(a, self.centre, scale)]
        out: List[Move] = []
        combos = list(itertools.combinations(axes, 3))
        rng.shuffle(combos)
        for t in combos:
            for _ in range(4):
                self._add("n3", self.centre, self.best, {a: rng.choice(self.step_values(a, self.centre, scale)) for a in t}, out)
            if len(out) >= TRIPLES_MAX:
                break
        return out

    def scan_moves(self, live: Sequence[str]) -> List[Move]:
        """Every setting over its whole range: climb.scan_points evenly spaced values (every value of a choice)."""
        n = max(2, tune("climb.scan_points"))
        out: List[Move] = []
        for a in live:
            s = self.space.specs.get(a)
            if s is None:
                continue
            cur = self.centre.params.get(a, s.default)
            if s.allowed or s.type == "bool":
                vals = [v for v in (s.allowed or (True, False)) if v != cur]
            else:
                vals = []
                for i in range(n):
                    v = self.space._from_unit(s, i / (n - 1))
                    if v != cur and v not in vals and self._allowed(a, v):
                        vals.append(v)
            for v in vals:
                self._add("scan", self.centre, self.best, {a: v}, out)
        return out

    def wide_pairs(self, live: Sequence[str], rng: random.Random) -> List[Move]:
        """All pairs of the WIDE_AXES strongest settings, each at WIDE_STEPS of its range either way (at most two
        other values of a choice); a random climb.wide_max of them when there are more."""
        def values(a: str) -> List[object]:
            s = self.space.specs[a]
            cur = self.centre.params.get(a, s.default)
            if s.allowed or s.type == "bool":
                return [v for v in (s.allowed or (True, False)) if v != cur][:2]
            u, out = self.space._to_unit(s, cur), []
            for h in WIDE_STEPS:
                for sgn in (1, -1):
                    v = self.space._from_unit(s, min(1.0, max(0.0, u + sgn * h)))
                    if v != cur and v not in out and self._allowed(a, v):
                        out.append(v)
            return out
        axes = [a for a in self.strongest(live, WIDE_AXES) if a in self.space.specs]
        vals = {a: values(a) for a in axes}
        combos = [(a, b, x, y) for a, b in itertools.combinations(axes, 2) for x in vals[a] for y in vals[b]]
        rng.shuffle(combos)
        out: List[Move] = []
        for a, b, x, y in combos:
            if len(out) >= tune("climb.wide_max"):
                break
            self._add("wide", self.centre, self.best, {a: x, b: y}, out)
        return out

    def _shake(self, base: Config, live: Sequence[str], rng: random.Random, k: int, scale: float) -> Dict[str, object]:
        """Random changes to k settings of `base` (gaussian steps of `scale` of the range; another value of a choice)."""
        axes = [a for a in live if a in self.steps]
        ch: Dict[str, object] = {}
        for a in rng.sample(axes, min(k, len(axes))) if axes else []:
            s = self.space.specs[a]
            cur = base.params.get(a, s.default)
            if s.allowed or s.type == "bool":
                others = [v for v in (s.allowed or (True, False)) if v != cur]
                if others:
                    ch[a] = rng.choice(others)
            else:
                v = self.space._from_unit(s, min(1.0, max(0.0, self.space._to_unit(s, cur) + rng.gauss(0, scale))))
                if v != cur and self._allowed(a, v):
                    ch[a] = v
        return ch

    def shake_points(self, live: Sequence[str], rng: random.Random, k: int, n: int) -> List[Move]:
        out: List[Move] = []
        tries = 0
        while len(out) < n and tries < 10 * n:
            tries += 1
            ch = self._shake(self.centre, live, rng, k, SHAKE_SCALE)
            if len(ch) >= max(1, k // 2):
                self._add("shake", self.centre, self.best, ch, out, self.cycle)
        return out

    def siblings(self, move: Move, rng: random.Random, n: Optional[int] = None) -> List[Config]:
        """One genotyping and one PairHMM variation of a structural move (none for a cheap-tier move)."""
        if n is None:
            n = tune("riders.siblings")
        if n <= 0 or all(a in CHEAP for a in move.axes):
            return []
        out: List[Config] = []
        tiers = ([a for a in GENOTYPING_PARAMS if a in self.steps], [a for a in PAIRHMM_PARAMS if a in self.space.names])
        for i in range(n):
            pool = tiers[i % 2]
            if not pool:
                continue
            for a in sorted(pool, key=lambda x: -self.axis_gain.get(x, 0.0))[:2] + [rng.choice(pool)]:
                vals = self.step_values(a, move.config) if a in self.steps else \
                    [v for v in (self.space.specs[a].allowed or ()) if v != move.config.params.get(a)]
                if not vals:
                    continue
                c = self.space._assemble({**move.config.params, a: rng.choice(vals)})
                if c is not None and c.key not in self.tried and c.key != move.config.key:
                    self.tried.add(c.key)
                    out.append(c)
                    break
        return out

    # ---- the main line ---------------------------------------------------------------------------------------------
    def _sweep(self, live: Sequence[str], rng: random.Random) -> List[Move]:
        out: List[Move] = []
        for kind, ch in self.queued:                       # re-applied winners and settings shared by other domains
            self._add(kind, self.centre, self.best, ch, out)
        self.queued = []
        self.sweep_is_scan = False
        if self.state == "final":
            out += self.singles(self.centre, self.best, live, "final", self.scale)
            out += self.pairs(live, rng, self.scale, kind="final")
        elif self.scan_next:
            self.scan_next = False
            self.sweep_is_scan = True
            out += self.scan_moves(live)
            self.note(f"full-range scan: {len(out)} configs")
        elif self.wide_next:
            self.wide_next = False
            out += self.wide_pairs(live, rng)
            self.note(f"coarse sweep without a gain: wide round of {len(out)} pairs")
        elif self.mode == "coarse":
            for sc in COARSE_SCALES:
                out += self.singles(self.centre, self.best, live, "n1", sc)
            out += self.pairs(live, rng, 1.0)
            if self.eager_wide and not self.wide_done:     # few domains: the workers are there, do not wait for it
                self.wide_done = True
                wide = self.wide_pairs(live, rng)
                out += wide
                self.note(f"wide round of {len(wide)} pairs with the first coarse sweep (few domains)")
        else:
            out += self.singles(self.centre, self.best, live, "n1", self.scale, skip=self.flat)
            out += self.singles(self.centre, self.best, live, "n1", self.scale * 0.5, skip=self.flat)
            out += self.pairs(live, rng, self.scale)
            out += self.triples(live, rng, self.scale)
        return out

    def refill(self, live: Sequence[str], rng: random.Random, room: int) -> List[Move]:
        """More moves, when the pool is empty. `room`: the workers this domain may still fill. The main line
        (sweeps, then perturbation) first; when it is waiting for results or has nothing, chains take the room."""
        self._live, self._rng = live, rng
        out = self._main(live, rng, room)
        if self.scan_flat and not self.scan_next:
            return out                                 # a flat scan is the end of this peak; chains do not follow it
        if not out and room > 0:
            out = self._start_chains(live, rng, room)
        return out

    def _main(self, live: Sequence[str], rng: random.Random, room: int) -> List[Move]:
        out: List[Move] = []
        if self.fresh:
            self.fresh = False
            out = self._sweep(live, rng)
            self.sweep_sent, self.sweep_back = len(out), 0
            self.pool.extend(out)
            if out or self.state == "final":
                return out
        if self.state in ("climb", "final") and self.sweep_sent and \
                self.sweep_back < tune("climb.back_frac") * self.sweep_sent:
            return []                                      # the sweep is still coming back
        if self.state == "climb" and self.scan_flat and not self.scan_next:
            return []                                      # the scan came back under the answer set's resolution
        if self.state == "final":                          # a halved final sweep, down to the floor; then chains
            if any(self.scale * st > self._floor(a) for a, st in self.steps.items()):
                self.scale *= 0.5
                self.fresh = True
                return self._main(live, rng, room)
            return []
        if self.state == "climb":
            if self.mode == "coarse":
                self.wide_done = True                      # no wide pair round: the v4.11 climbs spent it and gained nothing
                self.mode = "fine"
                for a in self.steps:
                    self.steps[a] = max(self._floor(a), self.steps[a] * 0.5)
                self.scale, self.fresh = 1.0, True
                self.note("coarse sweeps exhausted: fine mode")
                return self._main(live, rng, room)
            if any(self.scale * st > self._floor(a) for a, st in self.steps.items()):
                self.scale *= 0.5
                self.fresh = True
                return self._main(live, rng, room)
            self.state, self.cycle, self.cycle_stage = "perturb", 0, "points"
            self.cycle_sent = 0
            self.note("fine sweeps exhausted at the smallest steps: perturbation")
        return self._perturb(live, rng, room)

    def _perturb(self, live: Sequence[str], rng: random.Random, room: int) -> List[Move]:
        out: List[Move] = []
        if self.cycle_sent == 0:                           # start this K: its shaken points
            ks = tune("climb.shake_k")
            for _ in range(len(ks) + 1):
                if self.cycle >= len(ks):                  # a full round of K without acceptance: start over
                    self.cycle, self.rounds = 0, self.rounds + 1
                    self.note(f"perturbation round {self.rounds} without acceptance: starting over")
                out = self.shake_points(live, rng, ks[self.cycle], max(tune("climb.shake_points"), min(room, 60)))
                if out:
                    break
                self.cycle += 1
            if not out:
                return []                                  # nothing can be shaken (every point tried)
            self.cycle_stage, self.cycle_sent, self.cycle_back, self.cycle_best = "points", len(out), 0, None
        elif self.cycle_back < tune("climb.back_frac") * self.cycle_sent:
            return []
        elif self.cycle_stage == "points":                 # the points are back: one local round around the best
            if self.cycle_best is None:
                self.cycle, self.cycle_sent = self.cycle + 1, 0
                return self._perturb(live, rng, room)
            cfg, sc = self.cycle_best
            out = self.singles(cfg, sc, live, "local", 1.0, cycle=self.cycle)
            self.cycle_stage, self.cycle_sent, self.cycle_back = "local", len(out), 0
            if not out:
                return self._judge(live, rng, room)
        else:
            return self._judge(live, rng, room)
        self.pool.extend(out)
        return out

    def _judge(self, live: Sequence[str], rng: random.Random, room: int) -> List[Move]:
        """End of one perturbation step: better or equal replaces the centre (ParamILS), else the next K."""
        cfg, sc = self.cycle_best if self.cycle_best else (None, None)
        if cfg is not None and sc >= self.best - EQUAL and cfg.key != self.centre.key and \
                (sc > self.best + tune("climb.accept_min") or self.equal_accepts < tune("climb.equal_accepts")):
            own = sc > self.best + tune("climb.accept_min")
            if not own:
                self.equal_accepts += 1
            self._recentre(cfg, sc, self.history[-1][0], (), own=own)
            self.state, self.scale, self.fresh = "climb", 1.0, True
            self.note(f"perturbation accepted at {sc:.4f}")
            return self._main(live, rng, room)
        self.cycle, self.cycle_sent = self.cycle + 1, 0
        return self._perturb(live, rng, room)

    # ---- chains --------------------------------------------------------------------------------------------------------
    def _start_chains(self, live: Sequence[str], rng: random.Random, room: int) -> List[Move]:
        out: List[Move] = []
        per = 1 + max(1, tune("climb.chain_local"))
        while room > 0 and len(self.chains) < tune("climb.max_chains"):
            self._chain_seq += 1
            ch = Chain(self._chain_seq, self.centre, self.best)
            mv = self._chain_shake(ch, live, rng)
            if mv is None:
                break
            self.chains[ch.id] = ch
            out.append(mv)
            room -= per
        if out:
            self.pool.extend(out)
        return out

    def _chain_shake(self, ch: Chain, live: Sequence[str], rng: random.Random) -> Optional[Move]:
        final = self.state == "final"
        ks, scale = (FINAL_CHAIN_KS, FINAL_CHAIN_SCALE) if final else (CHAIN_KS, CHAIN_SCALE)
        for _ in range(10):
            changes = self._shake(ch.base, live, rng, rng.choice(ks), scale)
            out: List[Move] = []
            if changes and self._add("chain", ch.base, ch.base_score, changes, out, cycle=ch.iters, chain=ch.id):
                ch.stage, ch.sent, ch.back, ch.point, ch.best = "shake", 1, 0, None, None
                ch.axes = out[0].axes
                return out[0]
        return None

    def _chain_local(self, ch: Chain) -> List[Move]:
        cfg, sc = ch.point
        n = max(1, tune("climb.chain_local"))
        axes = list(ch.axes) + [a for a in self.strongest(self._live, n) if a not in ch.axes]
        scale = 0.5 if self.state == "final" else 1.0
        out: List[Move] = []
        for a in axes:
            if a not in self.steps:
                continue
            for v in self.step_values(a, cfg, scale):
                self._add("chainl", cfg, sc, {a: v}, out, cycle=ch.iters, chain=ch.id)
            if len(out) >= n:
                break
        return out[:n]

    def _chain_result(self, move: Move, config: Config, score: Optional[float]) -> None:
        ch = self.chains.get(move.chain)
        if ch is None or move.cycle != ch.iters:
            return                                         # a chain that ended, or a late result of a past round
        if move.kind == "chain":                           # the shaken point is back: a local round around it
            if ch.stage != "shake":
                return
            if score is None:
                self._chain_next(ch)
                return
            ch.point = ch.best = (config, score)
            moves = self._chain_local(ch)
            if not moves:
                self._chain_judge(ch)
                return
            ch.stage, ch.sent, ch.back = "local", len(moves), 0
            self.pool.extend(moves)
            return
        if ch.stage != "local":
            return
        ch.back += 1
        if score is not None and (ch.best is None or score > ch.best[1]):
            ch.best = (config, score)
        if ch.back >= tune("climb.back_frac") * ch.sent:
            self._chain_judge(ch)

    def _chain_judge(self, ch: Chain) -> None:
        """Better or equal: the chain moves there (ParamILS); then its next shake."""
        if ch.best is not None and ch.best[1] >= ch.base_score - EQUAL:
            ch.base, ch.base_score = ch.best
        self._chain_next(ch)

    def _chain_next(self, ch: Chain) -> None:
        ch.iters += 1                                      # results of the round that ended are ignored from now on
        if ch.base_score < self.best - tune("climb.decisive") or self._rng.random() < tune("climb.chain_restart"):
            ch.base, ch.base_score = self.centre, self.best        # fell behind (or by chance): from the centre again
        mv = self._chain_shake(ch, self._live, self._rng)
        if mv is None:
            self.chains.pop(ch.id, None)
        else:
            self.pool.append(mv)

    # ---- bookkeeping ---------------------------------------------------------------------------------------------------
    def dropped(self, move: Move) -> None:
        """A move that could not be sent counts as back, so a sweep or chain never waits for it."""
        self.inflight = max(0, self.inflight - 1)
        if move.kind in SWEEP_KINDS and move.centre_key == self.centre.key:
            self.sweep_back += 1
        if move.kind in ("shake", "local") and move.cycle == self.cycle and self.state == "perturb":
            self.cycle_back += 1
        if move.kind in CHAIN_KINDS:
            self._chain_result(move, move.config, None)

    @property
    def main_pending(self) -> bool:
        """Main-line moves (sweeps, perturbation, the final) waiting in the pool; chain moves do not count, so many
        chains never starve the main line of a refill."""
        return any(m.kind not in CHAIN_KINDS for m in self.pool)

    def take(self) -> Optional[Move]:
        if not self.pool:
            return None
        self.inflight += 1
        i = next((j for j, m in enumerate(self.pool) if m.kind not in CHAIN_KINDS), 0)   # the main line first
        return self.pool.pop(i)

    # ---- inputs from the controller ------------------------------------------------------------------------------------
    def queue(self, kind: str, changes: Dict[str, object]) -> None:
        """Changes to try on the centre at the next sweep (settings shared by another domain, a final pairing)."""
        self.queued.append((kind, changes))

    def scan_open(self) -> bool:
        """A full-range scan is about to be sent, or is still short of climb.back_frac of its results."""
        if self.scan_next:
            return True
        return bool(self.sweep_is_scan and self.sweep_sent
                    and self.sweep_back < tune("climb.back_frac") * self.sweep_sent)

    def offer(self, config: Config, score: float, now: float) -> bool:
        """A point the registry puts on this hill. A gain of climb.decisive moves the centre at once and the scan
        runs again. A smaller gain does not interrupt a scan, and does not move a peak whose scan already came
        back flat. Not counted as the domain's own progress."""
        if config.key == self.centre.key or score <= self.best + tune("climb.accept_min"):
            return False
        if self.scan_open() and score < self.best + tune("climb.decisive"):
            return False                                       # a small gain does not abort a scan that is still out
        self._recentre(config, score, now, (), own=False)
        return True

    def start_final(self) -> None:
        self.state, self.pool, self.fresh, self.scale = "final", [], True, 0.5
        self.winners = []
        self.scan_next = self.wide_next = False
        self.sweep_is_scan = self.scan_flat = False
        self.chains.clear()                                # queued chain moves were just dropped with the pool

    def transplant(self, other: Config, live: Sequence[str]) -> int:
        """Another centre's values where it differs from this one become moves on this centre: each alone and pairs
        of the first four (v3.2's merge)."""
        diff = [(a, other.params.get(a)) for a in live if a in self.space.specs
                and not self.space._same(self.space.specs[a], self.centre.params.get(a), other.params.get(a))]
        for a, v in diff:
            self.queued.append(("transplant", {a: v}))
        for (a, va), (b, vb) in itertools.combinations(diff[:4], 2):
            self.queued.append(("transplant", {a: va, b: vb}))
        return len(diff) + len(list(itertools.combinations(diff[:4], 2)))

    def inherit(self, tried: Set[str], steps: Dict[str, float], axis_gain: Dict[str, float]) -> None:
        """At a merge: the parents' tested configs, finer steps and per-setting gains carry over."""
        self.tried |= tried
        for a, st in steps.items():
            if a in self.steps:
                self.steps[a] = min(self.steps[a], st)
        for a, g in axis_gain.items():
            self.axis_gain[a] = self.axis_gain.get(a, 0.0) + g

    # ---- results --------------------------------------------------------------------------------------------------------
    def consider(self, move: Move, is_head: bool, config: Config, score: Optional[float], failed: bool,
                 now: float) -> bool:
        """A result for this domain (a move or one of its siblings). True when the centre moved."""
        bad = failed or score is None or score <= 0
        current = move.centre_key == self.centre.key
        if is_head:
            self.inflight = max(0, self.inflight - 1)
            if move.kind in SWEEP_KINDS and current:
                self.sweep_back += 1
            if move.kind in ("shake", "local") and move.cycle == self.cycle and self.state == "perturb":
                self.cycle_back += 1
            if move.kind in ("n1", "local", "final") and len(move.axes) == 1:
                self._learn(move.axes[0], config, move, score, bad, now)
            if move.kind in CHAIN_KINDS:
                self._chain_result(move, config, None if bad else score)
        if bad:
            return self._maybe_recentre(now)
        if move.kind in ("shake", "local") and self.state == "perturb" and move.cycle == self.cycle:
            if self.cycle_best is None or score > self.cycle_best[1]:
                self.cycle_best = (config, score)
        if score > self.best + tune("climb.decisive"):                  # decisive: at once, whatever produced it
            self._credit(move, score)
            self._recentre(config, score, now, move.axes if current else (), own=True)
            return True
        if score > self.best + tune("climb.accept_min") and move.kind in CHAIN_KINDS:
            if self.scan_open() and score < self.best + tune("climb.decisive"):
                return False                                    # a chain under 0.05 does not abort a scan still out
            self._credit(move, score)                                   # a chain beat the domain: its centre moves
            self.note(f"a chain found {score:.4f}")
            self._recentre(config, score, now, (), own=True)
            return True
        if score > self.best + tune("climb.accept_min") and move.kind in SWEEP_KINDS and self.state in ("climb", "final"):
            self._credit(move, score)
            self.winners.append((move, config, score))
        return self._maybe_recentre(now)

    def _record_first_scan(self, score: float, key: str, held: bool) -> None:
        """The first scan's best point, once. Later scans do not overwrite it."""
        if self.first_scan_score is not None:
            return
        self.first_scan_score, self.first_scan_key, self.first_scan_held = score, key, held
        self.note("first scan %.4f (opened %.4f, %s)" % (score, self.birth_score, "peak held" if held else "centre moved"))

    def _credit(self, move: Move, score: float) -> None:
        for a in move.axes:
            if a in self.steps:
                self.axis_gain[a] = self.axis_gain.get(a, 0.0) + (score - move.centre_score)

    def _maybe_recentre(self, now: float) -> bool:
        back = self.sweep_sent and self.sweep_back >= tune("climb.back_frac") * self.sweep_sent
        if self.scan_flat and not self.scan_next:
            self.winners = []                              # the rest of a flat scan does not move the centre
            return False
        if self.sweep_is_scan and back:
            picked = None
            for item in self.winners:
                if picked is None or item[2] > picked[2]:
                    picked = item
            self.winners = []
            if picked is not None and picked[2] > self.best + tune("climb.decisive"):
                self.sweep_is_scan = False
                return False                            # a decisive gain already moved the centre and restarts the scan
            if picked is not None and picked[2] > self.best + tune("climb.accept_min"):
                move, config, score = picked
                gain = score - self.best
                self._recentre(config, score, now, move.axes if move.centre_key == self.centre.key else (), own=True)
                self.note(f"scan back, gain {gain:.4f}: centre moved")
                return True                             # above the answer set's resolution: keep it, do not rescan
            score = picked[2] if picked is not None else self.best
            key = picked[1].key if picked is not None and score > self.best else self.centre.key
            self.sweep_is_scan = False
            self.scan_flat = True
            self._record_first_scan(score if score > self.best else self.best, key, held=True)
            self.note(f"scan back, gain {max(0.0, score - self.best):.4f}: peak held")
            return False
        if not self.winners or self.sweep_back < tune("climb.back_frac") * max(1, self.sweep_sent):
            return False
        move, config, score = max(self.winners, key=lambda t: t[2])
        others = [(m, c, s) for m, c, s in self.winners if c.key != config.key and m.centre_key == move.centre_key]
        self._recentre(config, score, now, move.axes if move.centre_key == self.centre.key else (), own=True)
        changes: List[Dict[str, object]] = []
        for m, c, s in sorted(others, key=lambda t: -t[2])[:8]:        # re-apply the sweep's other winners
            ch = {a: c.params.get(a) for a in m.axes
                  if not self.space._same(self.space.specs[a], self.centre.params.get(a), c.params.get(a))}
            if ch:
                changes.append(ch)
                self.queued.append(("reapply", ch))
        if tune("climb.reapply_full"):                                  # v3.2: added one by one, pairs of the top 4
            acc: Dict[str, object] = dict(changes[0]) if changes else {}
            for ch in changes[1:]:
                acc.update(ch)
                self.queued.append(("reapply", dict(acc)))
            for a, b in itertools.combinations(changes[:4], 2):
                if len(changes) > 2:                                    # with two, the pair is the running sum
                    self.queued.append(("reapply", {**a, **b}))
        return True

    def _learn(self, axis: str, config: Config, move: Move, score: Optional[float], bad: bool, now: float) -> None:
        s = self.space.specs[axis]
        base = self.centre.params.get(axis, s.default) if move.centre_key == self.centre.key else None
        delta = None if bad else score - move.centre_score
        if base is not None and axis in self.book.axes:
            self.book.note_move(axis, config.params.get(axis), base, delta, bad, now)
        if bad or s.allowed or s.type == "bool" or axis not in self.steps:
            return
        if delta is not None and abs(delta) <= EQUAL:     # below the setting's resolution
            if self.mode == "coarse" and self.state == "climb":
                self.steps[axis] = min(STEP_MAX, self.steps[axis] * 2)
            elif base is not None:
                self.flat.add(axis)
        elif delta is not None and delta < 0:
            self.steps[axis] = max(self._floor(axis), self.steps[axis] * 0.5)

    def _recentre(self, config: Config, score: float, now: float, axes: Sequence[str], own: bool) -> None:
        was_scan = self.sweep_is_scan
        gain = score - self.best
        if was_scan:
            self._record_first_scan(score, config.key, held=False)
        if own and gain > 0:
            self.gains.append((now, gain))
            self.equal_accepts = 0
            self.rounds = 0
        self.cycle, self.cycle_sent, self.cycle_best = 0, 0, None   # a new centre: the next perturbation starts at K=2
        self.centre, self.best = config, score
        self.tried.add(config.key)
        self.pool = [m for m in self.pool if m.kind in CHAIN_KINDS]  # the chains walk on; the centre's sweep restarts
        self.winners = []
        self.flat.clear()
        self.fresh = True
        self.wide_done = False
        self.sweep_sent = self.sweep_back = 0
        self.sweep_is_scan = False
        if gain > tune("climb.accept_min") and self.state != "final":
            self.scan_flat = False                     # a resolved gain is a new centre; the neighbourhood climb goes on
            if gain >= tune("climb.decisive") and self.first_scan_score is None:
                self.scan_next = True                  # one full-range scan per basin; a later gain stays on the neighbourhood
            if self.state == "exhausted":
                self.state = "climb"
        if self.state == "perturb":
            self.state = "climb"
        self.scale = 1.0 if self.state != "final" else self.scale
        for a in axes:
            if a in self.recent:
                self.recent.remove(a)
            self.recent.append(a)
        self.recent = self.recent[-PAIR_AXES:]
        self.history.append((now, score, "own" if own else "taken"))
