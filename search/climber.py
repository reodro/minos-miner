"""A domain's climb: iterated local search (ParamILS) on one hill. Owner: search.

The domain only proposes moves and takes in results; it never decides its own lifecycle (founding, merging, freezing,
dropping are the controller's, engine.py). It reports facts: what is in flight, whether it has moves left, its
progress rate, and `state`.

  sweep        the neighbourhood of the centre. COARSE: one setting at a time at its learned step and at half and
               double it, plus pairs of the strongest settings. FINE (entered when a coarse sweep is exhausted without
               a gain): steps halved and never doubled again, settings whose step changed nothing are left out until
               the centre moves, triples join the pairs, and each exhausted fine sweep halves the scale once more
  acceptance   the centre moves to the best result of the sweep once most of it is back (a decisive gain at once); the
               sweep's other winners are re-applied on the new centre first. An equal-score neighbour is a sideways
               move once a sweep is back without a gain (a walk across the plateau, at most climb.sideways_max in a row)
  perturbation when a fine sweep is exhausted at every setting's smallest step: shake K settings of the centre (K = 2,
               3, 5, 8 in turn), then one local round around the best shaken point; the result replaces the centre if
               it is better OR EQUAL (ParamILS's acceptance). A full cycle of K values without that: state "exhausted"
  final        the last minutes on this one point: single steps at halving scales and all pairs of the top settings
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
# at this share), climb.sideways_max, climb.equal_accepts, climb.shake_k, climb.shake_points, riders.siblings.
EQUAL = 1e-9
COARSE_SCALES = (1.0, 0.5, 2.0)
PAIR_AXES, TRIPLE_AXES, TRIPLES_MAX = 8, 5, 40
SHAKE_SCALE = 0.3
STEP_MAX = 0.5
PROGRESS_WINDOW_S = 600.0
CHEAP = tuple(GENOTYPING_PARAMS) + tuple(PAIRHMM_PARAMS)
SWEEP_KINDS = ("n1", "n2", "n3", "reapply", "share", "final", "pairing")


@dataclass
class Move:
    kind: str                # n1 | n2 | n3 | reapply | share | shake | local | final | pairing
    config: Config
    axes: Tuple[str, ...]
    centre_key: str
    centre_score: float
    cycle: int = -1          # the perturbation step it belongs to


class Climber:
    def __init__(self, id: int, space: Space, book: AxisBook, centre: Config, score: float, now: float,
                 owner: Optional[Callable[[Config], Optional[int]]] = None, extra_axes: Sequence[str] = ()):
        self.id, self.space, self.book = id, space, book
        self.centre, self.best = centre, score
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
        self.state = "climb"                               # climb | perturb | exhausted | final
        self.mode = "coarse"                               # coarse | fine
        self.scale = 1.0                                   # fine-mode scale, halved per exhausted fine sweep
        self.fresh = True                                  # the centre's sweep is not sent yet
        self.sweep_sent = self.sweep_back = 0
        self.winners: List[Tuple[Move, Config, float]] = []
        self.equals: List[Config] = []
        self.sideways = 0
        self.flat: Set[str] = set()
        self.recent: List[str] = []                        # settings of accepted moves, most recent last
        self.axis_gain: Dict[str, float] = {}
        self.queued: List[Tuple[str, Dict[str, object]]] = []    # (kind, changes) to try on the next sweep first
        self.cycle = 0                                     # index into climb.shake_k
        self.cycle_stage = "points"                        # points | local
        self.cycle_sent = self.cycle_back = 0
        self.cycle_best: Optional[Tuple[Config, float]] = None
        self.equal_accepts = 0
        self.gains: List[Tuple[float, float]] = []         # (time, gain) of the domain's OWN accepted moves
        self.history: List[Tuple[float, float, str]] = [(now, score, "born")]

    # ---- facts for the controller --------------------------------------------------------------------------------
    def progress(self, now: float) -> float:
        """Own gain per minute over the last PROGRESS_WINDOW_S (or since birth, at least two minutes)."""
        span = max(120.0, min(PROGRESS_WINDOW_S, now - self.born))
        return sum(g for t, g in self.gains if t >= now - span) / (span / 60.0)

    def age(self, now: float) -> float:
        return now - self.born

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
             cycle: int = -1) -> bool:
        c = self._variant(base, changes)
        if c is None:
            return False
        self.tried.add(c.key)
        out.append(Move(kind, c, tuple(sorted(changes)), base.key, base_score, cycle))
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

    def shake_points(self, live: Sequence[str], rng: random.Random, k: int, n: int) -> List[Move]:
        out: List[Move] = []
        axes = [a for a in live if a in self.steps]
        tries = 0
        while len(out) < n and tries < 10 * n and axes:
            tries += 1
            ch: Dict[str, object] = {}
            for a in rng.sample(axes, min(k, len(axes))):
                s = self.space.specs[a]
                cur = self.centre.params.get(a, s.default)
                if s.allowed or s.type == "bool":
                    others = [v for v in (s.allowed or (True, False)) if v != cur]
                    if others:
                        ch[a] = rng.choice(others)
                else:
                    v = self.space._from_unit(s, min(1.0, max(0.0, self.space._to_unit(s, cur) + rng.gauss(0, SHAKE_SCALE))))
                    if v != cur and self._allowed(a, v):
                        ch[a] = v
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

    # ---- the state machine -------------------------------------------------------------------------------------------
    def _sweep(self, live: Sequence[str], rng: random.Random) -> List[Move]:
        out: List[Move] = []
        for kind, ch in self.queued:                       # re-applied winners and settings shared by other domains
            self._add(kind, self.centre, self.best, ch, out)
        self.queued = []
        if self.state == "final":
            out += self.singles(self.centre, self.best, live, "final", self.scale)
            out += self.pairs(live, rng, self.scale, kind="final")
        elif self.mode == "coarse":
            for sc in COARSE_SCALES:
                out += self.singles(self.centre, self.best, live, "n1", sc)
            out += self.pairs(live, rng, 1.0)
        else:
            out += self.singles(self.centre, self.best, live, "n1", self.scale, skip=self.flat)
            out += self.singles(self.centre, self.best, live, "n1", self.scale * 0.5, skip=self.flat)
            out += self.pairs(live, rng, self.scale)
            out += self.triples(live, rng, self.scale)
        return out

    def refill(self, live: Sequence[str], rng: random.Random, room: int) -> List[Move]:
        """More moves, when the pool is empty. `room`: the workers this domain may still fill."""
        if self.state == "exhausted":
            return []
        out: List[Move] = []
        if self.fresh:
            self.fresh = False
            out = self._sweep(live, rng)
            self.sweep_sent, self.sweep_back = len(out), 0
            self.pool.extend(out)
            if out or self.state == "final":
                return out
        if self.state in ("climb", "final") and self.sweep_back < tune("climb.back_frac") * max(1, self.sweep_sent):
            return []                                      # the sweep is still coming back
        if self.state in ("climb", "final") and self.equals and not self.winners and self.sideways < tune("climb.sideways_max"):
            self.sideways += 1
            self._recentre(self.equals[-1], self.best, self.history[-1][0], (), own=False)
            return self.refill(live, rng, room)
        if self.state == "final":                          # a halved final sweep, down to the floor
            if any(self.scale * st > self._floor(a) for a, st in self.steps.items()):
                self.scale *= 0.5
                self.fresh = True
                return self.refill(live, rng, room)
            return []
        if self.state == "climb":
            if self.mode == "coarse":
                self.mode = "fine"
                for a in self.steps:
                    self.steps[a] = max(self._floor(a), self.steps[a] * 0.5)
                self.scale, self.fresh = 1.0, True
                return self.refill(live, rng, room)
            if any(self.scale * st > self._floor(a) for a, st in self.steps.items()):
                self.scale *= 0.5
                self.fresh = True
                return self.refill(live, rng, room)
            self.state, self.cycle, self.cycle_stage = "perturb", 0, "points"
            self.cycle_sent = 0
        return self._perturb(live, rng, room)

    def _perturb(self, live: Sequence[str], rng: random.Random, room: int) -> List[Move]:
        out: List[Move] = []
        if self.cycle_sent == 0:                           # start this K: its shaken points
            if self.cycle >= len(tune("climb.shake_k")):
                self.state = "exhausted"
                return []
            out = self.shake_points(live, rng, tune("climb.shake_k")[self.cycle], max(tune("climb.shake_points"), min(room, 60)))
            self.cycle_stage, self.cycle_sent, self.cycle_back, self.cycle_best = "points", len(out), 0, None
            if not out:
                self.cycle += 1
                return self._perturb(live, rng, room)
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
            return self.refill(live, rng, room)
        self.cycle, self.cycle_sent = self.cycle + 1, 0
        return self._perturb(live, rng, room)

    def dropped(self, move: Move) -> None:
        """A move that could not be sent counts as back, so a sweep never waits for it."""
        self.inflight = max(0, self.inflight - 1)
        if move.kind in SWEEP_KINDS and move.centre_key == self.centre.key:
            self.sweep_back += 1
        if move.kind in ("shake", "local") and move.cycle == self.cycle and self.state == "perturb":
            self.cycle_back += 1

    def take(self) -> Optional[Move]:
        if not self.pool:
            return None
        self.inflight += 1
        return self.pool.pop(0)

    # ---- inputs from the controller ------------------------------------------------------------------------------------
    def queue(self, kind: str, changes: Dict[str, object]) -> None:
        """Changes to try on the centre at the next sweep (settings shared by another domain, a final pairing)."""
        self.queued.append((kind, changes))

    def offer(self, config: Config, score: float, now: float) -> bool:
        """A global point (a wave cell, a scout) the registry puts on this hill: the centre moves there when it is
        better by the resolution. Not counted as the domain's own progress."""
        if config.key == self.centre.key or score <= self.best + tune("climb.accept_min"):
            return False
        self._recentre(config, score, now, (), own=False)
        return True

    def start_final(self) -> None:
        self.state, self.pool, self.fresh, self.scale = "final", [], True, 0.5
        self.winners, self.equals = [], []

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
        if bad:
            return self._maybe_recentre(now)
        if move.kind in ("shake", "local") and self.state == "perturb" and move.cycle == self.cycle:
            if self.cycle_best is None or score > self.cycle_best[1]:
                self.cycle_best = (config, score)
        if score > self.best + tune("climb.decisive"):                  # decisive: at once, whatever produced it
            self._credit(move, score)
            self._recentre(config, score, now, move.axes if current else (), own=True)
            return True
        if score > self.best + tune("climb.accept_min") and move.kind in SWEEP_KINDS and self.state in ("climb", "final"):
            self._credit(move, score)
            self.winners.append((move, config, score))
            return self._maybe_recentre(now)
        if abs(score - self.best) <= EQUAL and is_head and current and move.kind in ("n1", "n2"):
            self.equals.append(config)
        return self._maybe_recentre(now)

    def _credit(self, move: Move, score: float) -> None:
        for a in move.axes:
            if a in self.steps:
                self.axis_gain[a] = self.axis_gain.get(a, 0.0) + (score - move.centre_score)

    def _maybe_recentre(self, now: float) -> bool:
        if not self.winners or self.sweep_back < tune("climb.back_frac") * max(1, self.sweep_sent):
            return False
        move, config, score = max(self.winners, key=lambda t: t[2])
        others = [(m, c, s) for m, c, s in self.winners if c.key != config.key and m.centre_key == move.centre_key]
        self._recentre(config, score, now, move.axes if move.centre_key == self.centre.key else (), own=True)
        for m, c, s in sorted(others, key=lambda t: -t[2])[:8]:        # re-apply the sweep's other winners
            ch = {a: c.params.get(a) for a in m.axes
                  if not self.space._same(self.space.specs[a], self.centre.params.get(a), c.params.get(a))}
            if ch:
                self.queued.append(("reapply", ch))
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
        gain = score - self.best
        if own and gain > 0:
            self.gains.append((now, gain))
            self.sideways = 0
            self.equal_accepts = 0
        self.cycle, self.cycle_sent, self.cycle_best = 0, 0, None   # a new centre: the next perturbation starts at K=2
        self.centre, self.best = config, score
        self.tried.add(config.key)
        self.pool, self.winners, self.equals = [], [], []
        self.flat.clear()
        self.fresh = True
        self.sweep_sent = self.sweep_back = 0
        if self.state == "perturb":
            self.state = "climb"
        self.scale = 1.0 if self.state != "final" else self.scale
        for a in axes:
            if a in self.recent:
                self.recent.remove(a)
            self.recent.append(a)
        self.recent = self.recent[-PAIR_AXES:]
        self.history.append((now, score, "own" if own else "taken"))
