"""A domain: one isolated region of the config space, climbed by iterated local search. Owner: search.

The objective is discrete: the score changes only when the calls change, so steps of one unit often change nothing,
equal scores are common, and categorical axes have no step at all. The domain therefore

  sweep             the neighbourhood of the centre. COARSE mode (a new domain on a steep slope): one axis at a time at
                    the axis's effective step and at half and double it, plus all pairs among the strongest axes.
                    FINE mode (the slope has flattened: a batch gained less than FINE_GAIN): the steps are halved and
                    never doubled, an axis whose step changed nothing is left out until the centre moves, and pairs
                    and triples among the strongest axes join the single steps; every exhausted sweep halves the scale
                    again down to the axis floor
  acceptance        the centre moves to the BEST of the sweep once most of it is back (a decisive gain at once), and
                    the sweep's other winners are merged onto the new centre first (acc2..accN and pairs of the top
                    winners); an equal-score neighbour is a SIDEWAYS move once the sweep is back without a gain (at
                    most SIDEWAYS_MAX in a row: a walk across a plateau); worse only informs the step table
  shaking           only in coarse mode and only once the sweep is exhausted: perturb k axes of the best (k = 2, 3, 5,
                    8 in turn), then one neighbourhood round around the best of the points; the domain moves there
                    only if that beats its best (variable neighbourhood search's escape, started near the good point)
  bracket           at the end of each batch (BATCH tests completed in this domain): 4 values across every live axis
                    around the centre, run in parallel, plus all pairs among the 5 strongest axes - the feasible form
                    of a dense grid; also early, on spare workers, once per centre. In fine mode its offsets shrink
  freeze            the fine sweep exhausted down to every axis's floor, or no gain of RESOLUTION (cumulative) for
                    STALL_S: the domain is frozen - its best stays a candidate, its workers go to the others, its
                    region is not explored again
  siblings          a structural move carries SIBLINGS cheap variations of itself (one genotyping, one PairHMM step):
                    the question "does this move work better with another threshold / PCR model", for 13% of a run,
                    packed into the move's traversal by the evaluator
  inherit           at a merge the survivor takes the loser's centre values on the axes where they differ as
                    transplants, the loser's tested set, and its finer steps; nothing crosses between live domains

The genotyping axes are ordinary axes of the neighbourhoods: a step on one shares the centre's traversal with the
other genotyping steps, so the evaluator packs them into one run.
"""
from __future__ import annotations

import itertools
import random
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

from common.types import Config
from search.axes import AxisBook, far_axes
from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS, Space

ACCEPT_MIN = 0.005                     # the answer set's resolution: smaller gains are not gains
DECISIVE = 0.05                        # a gain this large recentres at once; smaller ones wait for the sweep
EQUAL = 1e-9
SIDEWAYS_MAX = 3
SHAKE_K = (2, 3, 5, 8)
SHAKE_POINTS = 10
SHAKE_SCALE = 0.35                     # unit-space spread of a shaken axis
FINAL_SHAKE_SCALE = 0.12
BRACKET_OFFSETS = (-0.4, -0.15, 0.15, 0.4)
PROBE_POINTS = 7                       # a new domain's first sweep also probes each axis across its whole range
BRACKET_PAIR_AXES = 5
PAIR_AXES = 8
TRIPLE_AXES = 5
TRIPLES_MAX = 40
COARSE_SCALES = (1.0, 0.5, 2.0)
FINE_GAIN = 0.05                       # a batch below this gain: the slope has flattened - fine mode
FREEZE_BATCHES = 2                     # (kept for the record) batches without a gain; freezing is by exhaustion now
RESOLUTION = 0.05                      # the local score resolves the truth only above this: smaller gains are noise
STALL_S = 12 * 60                      # no gain >= RESOLUTION (cumulative) for this long: frozen
BACK_FRAC = 0.8                        # a set of moves counts as back once this share of it is
STEP_MAX = 0.5
SIBLINGS = 2
CHEAP = tuple(GENOTYPING_PARAMS) + tuple(PAIRHMM_PARAMS)


@dataclass
class Move:
    kind: str                            # n1 | n2 | n3 | shake | local | bracket | merge | transplant | final
    config: Config
    axes: Tuple[str, ...]
    centre_key: str
    centre_score: float
    priority: int
    shake_id: int = -1


@dataclass
class Shake:
    id: int
    k: int
    sent: int = 0                        # points sent
    back: int = 0
    best: Optional[Tuple[Config, float]] = None
    local_sent: int = -1                 # -1: the local round not started yet
    local_back: int = 0
    done: bool = False


class Domain:
    def __init__(self, id: int, space: Space, book: AxisBook, centre: Config, score: float, now: float,
                 owner: Optional[Callable[[Config], int]] = None):
        self.id, self.space, self.book = id, space, book
        self.centre, self.best = centre, score
        self.born = self.improved_at = now
        self.owner = owner                               # config -> domain id that owns it (isolation)
        self.steps: Dict[str, float] = {a: book.axes[a].step for a in book.names}
        for a in GENOTYPING_PARAMS:                      # the genotyping axes climb too (the book tracks structural ones)
            s = space.specs.get(a)
            if s is not None and a in space.names:
                self.steps[a] = max(0.12, 1.0 / max(1, s.max - s.min)) if s.type == "int" else 0.12
        self.pool: List[Move] = []
        self.tried: Set[str] = {centre.key}
        self.inflight = 0
        self.phase = "climb"                             # climb | final
        self.mode = "coarse"                             # coarse | fine
        self.fine_scale = 1.0
        self.flat_axes: Set[str] = set()                 # fine mode: axes whose step changed nothing at this centre
        self.fresh_centre = True                         # the sweep of this centre is not sent yet
        self.sweep_sent = self.sweep_back = 0            # sweep heads of the current centre
        self.generation = 0
        self.equals: List[Config] = []                   # equal-score neighbours of the current centre
        self.sweep: List[Tuple[Move, Config, float]] = []   # results beating the centre, waiting for the sweep
        self.merges: List[Move] = []                     # the last sweep's other winners / inherited transplants
        self.sideways = 0
        self.recent: List[str] = []                      # axes of accepted moves, most recent last
        self.axis_gain: Dict[str, float] = {}
        self.shake: Optional[Shake] = None
        self.shakes = 0
        self.completed = 0                               # heads back since the batch started
        self.batch = 0
        self.batch_start_best = score
        self.gains: List[float] = []                     # gain per closed batch
        self.bracket: Optional[Tuple[int, int]] = None   # (sent, back) of the running bracket
        self.bracket_generation = -1                     # the early bracket ran for this centre
        self.final_scale = 1.0
        self.history: List[Tuple[float, float]] = [(now, score)]
        self.retired = False
        self.frozen = False
        self.real_gain_at = now                          # the last time the gains since then reached RESOLUTION
        self.small_gains = 0.0
        self.exhausted_at_floor = False

    # ---- state ----------------------------------------------------------------------------------------------
    def far(self, other: Config, names: Sequence[str]) -> int:
        return far_axes(self.space, names, self.centre, other)

    def sweep_back_enough(self) -> bool:
        return self.sweep_back >= BACK_FRAC * max(1, self.sweep_sent)

    def gain_recent(self) -> float:
        """Gain over the last two batches (or since birth while younger than two batches): the racing criterion."""
        if len(self.gains) >= 2:
            return self.gains[-1] + self.gains[-2]
        return self.best - self.history[0][1]

    def strongest_axes(self, live: Sequence[str], n: int) -> List[str]:
        ranked = sorted((a for a in live if a in self.axis_gain), key=lambda a: -self.axis_gain[a])
        out = ranked[:n]
        for a in reversed(self.recent):
            if len(out) >= n:
                break
            if a in live and a not in out:
                out.append(a)
        for a in live:
            if len(out) >= n:
                break
            if a not in out:
                out.append(a)
        return out

    # ---- values ---------------------------------------------------------------------------------------------
    def _step_floor(self, axis: str) -> float:
        s = self.space.specs[axis]
        return max(0.02, 1.0 / max(1, (s.max - s.min))) if s.type == "int" else 0.02

    def _allowed(self, axis: str, value) -> bool:
        """Past a learned cliff on a structural axis: no. The genotyping axes are not in the book: always allowed."""
        return axis not in self.book.axes or self.book.allowed(axis, value)

    def step_values(self, axis: str, base: Config, scale: float = 1.0) -> List[object]:
        spec = self.space.specs[axis]
        cur = base.params.get(axis, spec.default)
        if spec.allowed or spec.type == "bool":
            return [v for v in (spec.allowed or (True, False)) if v != cur]
        u = self.space._to_unit(spec, cur)
        h = max(self._step_floor(axis), self.steps[axis] * scale)
        out = []
        for sgn in (1, -1):
            v = self.space._from_unit(spec, min(1.0, max(0.0, u + sgn * h)))
            if v != cur and v not in out and self._allowed(axis, v):
                out.append(v)
        return out

    def bracket_values(self, axis: str, scale: float = 1.0) -> List[object]:
        spec = self.space.specs[axis]
        if spec.allowed or spec.type == "bool":
            return []
        cur = self.centre.params.get(axis, spec.default)
        u = self.space._to_unit(spec, cur)
        out = []
        for d in BRACKET_OFFSETS:
            v = self.space._from_unit(spec, min(1.0, max(0.0, u + d * scale)))
            if v != cur and v not in out and self._allowed(axis, v):
                out.append(v)
        return out

    def probe_values(self, axis: str) -> List[object]:
        spec = self.space.specs[axis]
        if spec.allowed or spec.type == "bool":
            return []
        cur = self.centre.params.get(axis, spec.default)
        out = []
        for i in range(PROBE_POINTS):
            v = self.space._from_unit(spec, i / (PROBE_POINTS - 1))
            if v != cur and v not in out and self._allowed(axis, v):
                out.append(v)
        return out

    def _variant(self, base: Config, changes: Dict[str, object]) -> Optional[Config]:
        params = dict(base.params)
        params.update(changes)
        c = self.space._assemble(params)
        if c is None or c.key in self.tried or c.key == self.centre.key:
            return None
        if self.owner is not None and self.owner(c) not in (self.id, -1):
            return None                                   # another domain's territory
        return c

    def _add(self, kind: str, base: Config, base_score: float, changes: Dict[str, object], priority: int,
             out: List[Move], shake_id: int = -1) -> bool:
        c = self._variant(base, changes)
        if c is None:
            return False
        self.tried.add(c.key)
        out.append(Move(kind, c, tuple(sorted(changes)), base.key, base_score, priority, shake_id))
        return True

    # ---- neighbourhoods ----------------------------------------------------------------------------------------
    def neighbourhood1(self, base: Config, base_score: float, live: Sequence[str], kind: str = "n1",
                       scale: float = 1.0, shake_id: int = -1, skip: Sequence[str] = ()) -> List[Move]:
        out: List[Move] = []
        for axis in live:
            if axis in skip:
                continue
            for v in self.step_values(axis, base, scale):
                self._add(kind, base, base_score, {axis: v}, 0, out, shake_id)
        return out

    def neighbourhood2(self, live: Sequence[str], rng: random.Random, n_axes: int = PAIR_AXES,
                       scale: float = 1.0, kind: str = "n2") -> List[Move]:
        axes = [a for a in self.strongest_axes(live, n_axes) if self.step_values(a, self.centre, scale)]
        out: List[Move] = []
        for a, b in itertools.combinations(axes, 2):
            for x in self.step_values(a, self.centre, scale):
                for y in self.step_values(b, self.centre, scale):
                    self._add(kind, self.centre, self.best, {a: x, b: y}, 1, out)
        rng.shuffle(out)
        return out

    def neighbourhood3(self, live: Sequence[str], rng: random.Random, scale: float = 1.0) -> List[Move]:
        axes = [a for a in self.strongest_axes(live, TRIPLE_AXES) if self.step_values(a, self.centre, scale)]
        out: List[Move] = []
        combos = list(itertools.combinations(axes, 3))
        rng.shuffle(combos)
        for a, b, c in combos:
            for _ in range(4):
                ch = {x: rng.choice(self.step_values(x, self.centre, scale)) for x in (a, b, c)}
                self._add("n3", self.centre, self.best, ch, 2, out)
            if len(out) >= TRIPLES_MAX:
                break
        return out

    def shake_points(self, live: Sequence[str], rng: random.Random, k: int, shake_id: int, n: int,
                     scale: float = SHAKE_SCALE) -> List[Move]:
        out: List[Move] = []
        axes = list(live)
        tries = 0
        while len(out) < n and tries < 10 * n and axes:
            tries += 1
            changes: Dict[str, object] = {}
            for a in rng.sample(axes, min(k, len(axes))):
                s = self.space.specs[a]
                cur = self.centre.params.get(a, s.default)
                if s.allowed or s.type == "bool":
                    others = [v for v in (s.allowed or (True, False)) if v != cur]
                    if others:
                        changes[a] = rng.choice(others)
                else:
                    u = self.space._to_unit(s, cur) + rng.gauss(0, scale)
                    v = self.space._from_unit(s, min(1.0, max(0.0, u)))
                    if v != cur and self._allowed(a, v):
                        changes[a] = v
            if len(changes) >= max(1, k // 2):
                self._add("shake", self.centre, self.best, changes, 2, out, shake_id)
        return out

    def bracket_moves(self, live: Sequence[str]) -> List[Move]:
        scale = 1.0 if self.mode == "coarse" else self.fine_scale * 0.5
        out: List[Move] = []
        for axis in live:
            for v in self.bracket_values(axis, scale):
                self._add("bracket", self.centre, self.best, {axis: v}, 1, out)
        strong = [a for a in self.strongest_axes(live, BRACKET_PAIR_AXES) if self.step_values(a, self.centre, scale)]
        for a, b in itertools.combinations(strong, 2):
            for x in self.step_values(a, self.centre, scale):
                for y in self.step_values(b, self.centre, scale):
                    self._add("bracket", self.centre, self.best, {a: x, b: y}, 1, out)
        return out

    def siblings(self, move: Move, rng: random.Random, n: int = SIBLINGS) -> List[Config]:
        """`n` cheap variations of a structural move: one genotyping step and one PairHMM step of the moved config.
        None for a move that is itself on the cheap tiers only."""
        if n <= 0 or all(a in CHEAP for a in move.axes):
            return []
        out: List[Config] = []
        tiers = ([a for a in GENOTYPING_PARAMS if a in self.steps], [a for a in PAIRHMM_PARAMS if a in self.space.names])
        for i in range(n):
            pool = tiers[i % 2]
            if not pool:
                continue
            ranked = sorted(pool, key=lambda a: -self.axis_gain.get(a, 0.0))
            for axis in ranked[:2] + [rng.choice(pool)]:
                vals = self.step_values(axis, move.config) if axis in self.steps else \
                    [v for v in (self.space.specs[axis].allowed or ()) if v != move.config.params.get(axis)]
                if not vals:
                    continue
                params = dict(move.config.params)
                params[axis] = rng.choice(vals)
                c = self.space._assemble(params)
                if c is not None and c.key not in self.tried and c.key != move.config.key:
                    self.tried.add(c.key)
                    out.append(c)
                    break
        return out

    # ---- the state machine --------------------------------------------------------------------------------------
    def _sweep(self, live: Sequence[str], rng: random.Random) -> List[Move]:
        out: List[Move] = list(self.merges)
        self.merges = []
        if self.mode == "coarse":
            for sc in COARSE_SCALES:
                out += self.neighbourhood1(self.centre, self.best, live, scale=sc)
            out += self.neighbourhood2(live, rng)
        else:
            sc = self.fine_scale
            out += self.neighbourhood1(self.centre, self.best, live, scale=sc, skip=self.flat_axes)
            out += self.neighbourhood1(self.centre, self.best, live, scale=sc * 0.5, skip=self.flat_axes)
            out += self.neighbourhood2(live, rng, scale=sc)
            out += self.neighbourhood3(live, rng, scale=sc)
        return out

    def refill(self, live: Sequence[str], rng: random.Random, target: int = 0) -> List[Move]:
        """More moves for this domain: the sweep of a new centre; once it is back without a gain a sideways move,
        else (coarse) a shake or (fine) the sweep again at half the scale; an early bracket on spare workers."""
        if self.phase == "final":
            return self._refill_final(live, rng, target)
        out: List[Move] = []
        if self.fresh_centre:
            self.fresh_centre = False
            out = self._sweep(live, rng)
            self.sweep_sent, self.sweep_back = len(out), 0
            self.pool.extend(out)
            return out
        if not self.sweep_back_enough():
            if self.bracket is None and self.bracket_generation != self.generation and self.inflight < target:
                self.bracket_generation = self.generation       # spare workers: this centre's bracket, early
                return self.start_bracket_moves(live)
            return out
        if self.equals and not self.sweep and self.sideways < SIDEWAYS_MAX:
            self.sideways += 1
            self._recentre(self.equals[-1], self.best, self.history[-1][0], (), sideways=True)
            return self.refill(live, rng, target)
        if self.mode == "fine":                                   # exhausted: the same sweep at half the scale
            if any(self.fine_scale * st > self._step_floor(a) for a, st in self.steps.items()):
                self.fine_scale *= 0.5
                self.fresh_centre = True
                return self.refill(live, rng, target)
            self.exhausted_at_floor = True                        # nothing left at any scale: frozen
            self.frozen = True
            return out
        # coarse and exhausted: shake
        room = max(SHAKE_POINTS, target - self.inflight)
        sh = self.shake
        if sh is not None and not sh.done:
            if sh.local_sent < 0:
                if sh.back >= BACK_FRAC * max(1, sh.sent) and sh.best is not None:
                    cfg, sc = sh.best
                    local = self.neighbourhood1(cfg, sc, live, kind="local", shake_id=sh.id)
                    sh.local_sent, sh.local_back = len(local), 0
                    out += local
                    if not local:
                        sh.done = True
                elif sh.sent < max(SHAKE_POINTS, target):
                    pts = self.shake_points(live, rng, sh.k, sh.id, min(room, max(SHAKE_POINTS, target) - sh.sent))
                    sh.sent += len(pts)
                    out += pts
                    if not pts and sh.back >= sh.sent:
                        sh.done = True
            elif sh.local_back >= BACK_FRAC * max(1, sh.local_sent):
                sh.done = True
        if self.shake is None or self.shake.done:
            k = SHAKE_K[self.shakes % len(SHAKE_K)]
            self.shakes += 1
            self.shake = Shake(self.shakes, k)
            pts = self.shake_points(live, rng, k, self.shake.id, room)
            self.shake.sent = len(pts)
            out += pts
        self.pool.extend(out)
        return out

    def _refill_final(self, live: Sequence[str], rng: random.Random, target: int) -> List[Move]:
        out: List[Move] = []
        if self.fresh_centre or self.sweep_back_enough():
            self.fresh_centre = False
            self.final_scale = self.final_scale * 0.5 if self.sweep_sent else 0.5
            out += self.neighbourhood1(self.centre, self.best, live, kind="final", scale=self.final_scale)
            out += self.neighbourhood2(live, rng, scale=self.final_scale, kind="final")
            if out:
                self.sweep_sent, self.sweep_back = len(out), 0
        self.pool.extend(out)                                     # no filler: spare workers belong to the scout
        return out

    def start_bracket_moves(self, live: Sequence[str]) -> List[Move]:
        moves = self.bracket_moves(live)
        if moves:
            self.pool[:0] = moves
            self.bracket = (len(moves), 0)
        return moves

    def start_bracket(self, live: Sequence[str]) -> int:
        moves = self.start_bracket_moves(live)
        if not moves:
            self.close_batch()
        return len(moves)

    def close_batch(self) -> float:
        gain = self.best - self.batch_start_best
        self.gains.append(gain)
        self.batch += 1
        self.batch_start_best = self.best
        self.completed = 0
        self.bracket = None
        if self.mode == "coarse" and gain < FINE_GAIN:            # the slope has flattened: finer moves from here
            self.mode = "fine"
            self.shake = None
            for a in self.steps:
                self.steps[a] = max(self._step_floor(a), self.steps[a] * 0.5)
            self.fine_scale = 1.0
            self.fresh_centre = True
            self.pool = [m for m in self.pool if m.kind not in ("shake", "local")]
        return gain

    def stalled(self, now: float) -> bool:
        """No gain the truth can see for STALL_S: the climb is over on this hill."""
        return now - self.real_gain_at >= STALL_S

    def start_final(self, live: Sequence[str], rng: random.Random, target: int = 0) -> int:
        """The last minutes on this one point: neighbourhood 1 at half steps, all pairs of the top axes, small shakes."""
        self.phase = "final"
        self.frozen = False
        self.pool.clear()
        self.shake = None
        self.fresh_centre = True
        self.final_scale = 1.0
        self.sweep_sent = 0
        return len(self._refill_final(live, rng, target))

    def inherit(self, other: "Domain", live: Sequence[str]) -> int:
        """At a merge: the loser's centre values on the axes where the centres differ become transplant moves on this
        centre (each alone and pairs of the first four); its tested set and its finer steps carry over."""
        self.tried |= other.tried
        for a, st in other.steps.items():
            if a in self.steps:
                self.steps[a] = min(self.steps[a], st)
        diff = [(a, other.centre.params.get(a)) for a in live
                if not self.space._same(self.space.specs[a], self.centre.params.get(a), other.centre.params.get(a))
                and self._allowed(a, other.centre.params.get(a))]
        out: List[Move] = []
        for a, v in diff:
            self._add("transplant", self.centre, self.best, {a: v}, 0, out)
        for (a, va), (b, vb) in itertools.combinations(diff[:4], 2):
            self._add("transplant", self.centre, self.best, {a: va, b: vb}, 0, out)
        self.pool[:0] = out
        self.sweep_sent += len(out)
        return len(out)

    def take(self) -> Optional[Move]:
        if not self.pool:
            return None
        m = self.pool.pop(0)
        self.inflight += 1
        return m

    # ---- results ---------------------------------------------------------------------------------------------
    def consider(self, move: Move, is_head: bool, config: Config, score: Optional[float], failed: bool,
                 now: float) -> bool:
        """A result for this domain (a move's head, or one of its siblings). Returns True when the centre moved."""
        bad = failed or score is None or score <= 0
        current = move.centre_key == self.centre.key
        if is_head:
            self.inflight = max(0, self.inflight - 1)
            self.completed += 1
            if move.kind in ("n1", "n2", "n3", "probe", "final", "merge", "transplant") and current:
                self.sweep_back += 1
            if move.kind == "bracket" and self.bracket is not None:
                self.bracket = (self.bracket[0], self.bracket[1] + 1)
            sh = self.shake
            if sh is not None and move.shake_id == sh.id:
                if move.kind == "shake":
                    sh.back += 1
                    if not bad and (sh.best is None or score > sh.best[1]):
                        sh.best = (config, score)
                elif move.kind == "local":
                    sh.local_back += 1
            if move.kind in ("n1", "local", "final", "bracket", "transplant", "probe") and len(move.axes) == 1:
                self._learn(move.axes[0], config, move, score, bad, now)
        if bad:
            return self._maybe_recentre(now)
        if score > self.best + ACCEPT_MIN:
            for a in move.axes:
                if a in self.steps:
                    self.axis_gain[a] = self.axis_gain.get(a, 0.0) + (score - move.centre_score)
            self.sweep.append((move, config, score))
            return self._maybe_recentre(now, decisive=score - self.best >= DECISIVE)
        if abs(score - self.best) <= EQUAL and is_head and move.kind in ("n1", "n2") and current:
            self.equals.append(config)
        return self._maybe_recentre(now)

    def _maybe_recentre(self, now: float, decisive: bool = False) -> bool:
        """The centre moves to the best of the sweep once most of the sweep is back, or at once for a decisive gain."""
        if not self.sweep:
            return False
        if not (decisive or self.sweep_back_enough()):
            return False
        move, config, score = max(self.sweep, key=lambda t: t[2])
        others = [(m, c, s) for m, c, s in self.sweep if c.key != config.key and m.centre_key == move.centre_key]
        self._recentre(config, score, now, move.axes if move.centre_key == self.centre.key else ())
        self._queue_merges(others)
        return True

    def _queue_merges(self, others: List[Tuple[Move, Config, float]]) -> None:
        """The sweep's other winners re-applied on the new centre: each alone, acc2..accN by gain, pairs of the top 4."""
        ranked = sorted(others, key=lambda t: -t[2])
        changes: List[Dict[str, object]] = []
        for m, c, s in ranked:
            ch = {a: c.params.get(a) for a in m.axes
                  if not self.space._same(self.space.specs[a], self.centre.params.get(a), c.params.get(a))}
            if ch:
                changes.append(ch)
        out: List[Move] = []
        acc: Dict[str, object] = {}
        for i, ch in enumerate(changes[:8]):
            self._add("merge", self.centre, self.best, dict(ch), 0, out)
            acc.update(ch)
            if i >= 1:
                self._add("merge", self.centre, self.best, dict(acc), 0, out)
        for a, b in itertools.combinations(changes[:4], 2):
            self._add("merge", self.centre, self.best, {**a, **b}, 0, out)
        self.merges = out

    def _learn(self, axis: str, config: Config, move: Move, score: Optional[float], bad: bool, now: float) -> None:
        spec = self.space.specs[axis]
        base_value = self.centre.params.get(axis, spec.default) if move.centre_key == self.centre.key else None
        delta = None if bad else score - move.centre_score
        if base_value is not None and axis in self.book.axes:
            self.book.note_move(axis, config.params.get(axis), base_value, delta, bad, now)
        if bad or spec.allowed or spec.type == "bool":
            return
        if delta is not None and abs(delta) <= EQUAL:        # the step was below the axis's resolution
            if self.mode == "coarse":
                self.steps[axis] = min(STEP_MAX, self.steps[axis] * 2)
            elif base_value is not None:
                self.flat_axes.add(axis)                      # fine: leave the axis out until the centre moves
        elif delta is not None and delta < 0:
            self.steps[axis] = max(self._step_floor(axis), self.steps[axis] * 0.5)

    def _recentre(self, config: Config, score: float, now: float, axes: Sequence[str], sideways: bool = False) -> None:
        if not sideways:
            self.small_gains += score - self.best
            if self.small_gains >= RESOLUTION - 1e-9:        # gains below the resolution only count once they add up
                self.real_gain_at, self.small_gains = now, 0.0
        self.centre, self.best = config, score
        if not sideways:
            self.improved_at = now
            self.sideways = 0
        self.generation += 1
        if self.bracket is not None:                     # the bracket's unsent moves die with the old centre: its
            dropped = sum(1 for m in self.pool if m.kind == "bracket")      # count shrinks so the batch still closes
            self.bracket = (max(self.bracket[1], self.bracket[0] - dropped), self.bracket[1])
        self.pool.clear()
        self.equals.clear()
        self.sweep.clear()
        self.flat_axes.clear()
        self.shake = None
        self.fresh_centre = True
        self.fine_scale = 1.0
        self.sweep_sent = self.sweep_back = 0
        self.final_scale = 1.0
        for a in axes:
            if a in self.recent:
                self.recent.remove(a)
            self.recent.append(a)
        self.recent = self.recent[-PAIR_AXES:]
        self.history.append((now, score))
