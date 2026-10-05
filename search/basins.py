"""A basin: one region of the config space climbed asynchronously. Owner: search.

The basin keeps a centre (its best config so far) and a pool of moves around it - a sweep. Moves are dispatched
whenever the engine has a free worker for this basin. The centre moves to the BEST of the sweep once most of it
is back (SWEEP_FRAC), never to the first result that happens to beat it, and only for a gain of ACCEPT_MIN or more
- the answer set's measured resolution (its ranking is exact for gaps >= 0.005; the local score itself is
deterministic). The sweep's other winners are merged onto the new centre first (the legacy project's accN).
Results of moves made from an older centre are not wasted either: a winner among them is a merge candidate too,
and one that beats the current centre outright is taken at once.

Move kinds, in the order they are tried:
  step        one structural axis, +-1x, +-0.5x, +-2x the axis's current step (numeric) or another value (enum/bool)
  merge       the centre plus the last winning single-axis moves combined (prefixes by gain, and pairs of the top four)
  transplant  one or two structural axes copied from another basin's best config
  probe       a full-range grid of one axis (at birth, and again when the basin stalls): cliffs out of step reach
  pair        two single steps at once (model-ranked when a ranker is available)
"""
from __future__ import annotations

import itertools
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from common.types import Config
from search.axes import AxisBook, far_axes, rms_distance
from search.space import Space

PROBE_POINTS = 7
ACCEPT_MIN = 0.005                        # a gain below the answer set's resolution is not an improvement
SWEEP_FRAC = 0.8                          # the centre moves once this share of the sweep is back (or a decisive gain)
DECISIVE = 0.05                           # a gain this large does not wait for the rest of the sweep


@dataclass
class Move:
    kind: str
    config: Config
    axes: Tuple[str, ...]                    # structural axes changed from the centre it was made from
    centre_key: str
    centre_score: float
    priority: int                            # lower runs first


@dataclass
class Basin:
    id: int
    space: Space
    book: AxisBook
    centre: Config
    best: float
    born: float
    improved_at: float
    generation: int = 0
    pool: List[Move] = field(default_factory=list)
    inflight: int = 0                                     # heads in flight for this basin
    winners: List[Tuple[str, object, float, int]] = field(default_factory=list)   # (axis, value, gain, generation)
    tried: set = field(default_factory=set)
    probed_at_generation: int = -1
    refilled_generation: int = -1                         # refill() ran for this centre
    stalls: int = 0                                       # refills that produced nothing new
    retired: bool = False
    history: List[Tuple[float, float]] = field(default_factory=list)              # (t, best)
    sources: List[Config] = field(default_factory=list)   # configs on this hill (valley-tested): transplant donors
    gen_sent: int = 0                                     # moves dispatched from the current centre
    gen_back: int = 0                                     # ... whose results are back
    sweep: List[Tuple[Optional[Move], Config, float]] = field(default_factory=list)   # results beating the centre

    # ---- state ------------------------------------------------------------------------------------------------
    def distance(self, other: Config) -> float:
        return rms_distance(self.space, self.book.names, self.centre, other)

    def far(self, other: Config, names: Sequence[str]) -> int:
        """Axes among `names` (the live structural axes) on which `other` is far from this centre."""
        return far_axes(self.space, names, self.centre, other)

    def exhausted(self) -> bool:
        return not self.pool and self.inflight == 0 and self.stalls >= 2

    def wants_refill(self) -> bool:
        """The pool is empty and either the centre is new or nothing is in flight that could still move it."""
        return not self.pool and (self.refilled_generation != self.generation or self.inflight == 0)

    # ---- results --------------------------------------------------------------------------------------------
    def consider(self, move: Optional[Move], config: Config, score: Optional[float], failed: bool, now: float) -> bool:
        """A result for this basin (a move's head, or one of its riders). Returns True when the centre moved."""
        bad = failed or score is None or score <= 0
        current = move is None or move.centre_key == self.centre.key      # part of the sweep around this centre
        if move is not None and len(move.axes) == 1 and move.kind in ("step", "probe", "transplant"):
            axis = move.axes[0]
            spec = self.space.specs[axis]
            centre_value = self.centre.params.get(axis, spec.default)
            delta = None if bad else score - move.centre_score
            self.book.note_move(axis, config.params.get(axis), centre_value, delta, bad, now)
            if delta is not None and delta > ACCEPT_MIN:
                self.winners.append((axis, config.params.get(axis), delta, self.generation))
                self.winners = self.winners[-24:]
        if current and move is not None:
            self.gen_back += 1
        if bad or score <= self.best + ACCEPT_MIN:
            return self._maybe_recentre(now)
        if not current:                                      # an older centre's move, better than where we are now
            return self._recentre(config, score, now)
        self.sweep.append((move, config, score))
        return self._maybe_recentre(now, decisive=score - self.best >= DECISIVE)

    def adopt(self, config: Config, score: float, now: float) -> bool:
        """A scout on this hill that is higher than the centre: the basin moves there (its sweep in flight becomes
        merge material, as after any re-centre)."""
        if config.key == self.centre.key or score <= self.best + ACCEPT_MIN:
            return False
        return self._recentre(config, score, now)

    def _maybe_recentre(self, now: float, decisive: bool = False) -> bool:
        if not self.sweep:
            return False
        if not (decisive or self.gen_back >= SWEEP_FRAC * max(1, self.gen_sent) or self.gen_back >= self.gen_sent):
            return False
        move, config, score = max(self.sweep, key=lambda t: t[2])
        return self._recentre(config, score, now)

    def _recentre(self, config: Config, score: float, now: float) -> bool:
        self.centre, self.best, self.improved_at = config, score, now
        self.generation += 1
        self.pool.clear()                                    # moves of the old centre are stale
        self.sweep.clear()
        self.gen_sent = self.gen_back = 0
        self.stalls = 0
        self.history.append((now, score))
        return True

    # ---- moves ------------------------------------------------------------------------------------------------
    def _fresh(self, c: Optional[Config]) -> bool:
        return c is not None and c.key not in self.tried and c.key != self.centre.key

    def _variant(self, changes: Dict[str, object]) -> Optional[Config]:
        params = dict(self.centre.params)
        params.update(changes)
        c = self.space._assemble(params)
        return c if self._fresh(c) else None

    def _add(self, kind: str, changes: Dict[str, object], priority: int, out: List[Move]) -> None:
        c = self._variant(changes)
        if c is None:
            return
        self.tried.add(c.key)
        out.append(Move(kind, c, tuple(sorted(changes)), self.centre.key, self.best, priority))

    def step_values(self, axis: str, scales: Sequence[float] = (1.0, 0.5, 2.0)) -> List[object]:
        """Values one step away from the centre on `axis`, nearest scale first."""
        spec = self.space.specs[axis]
        cur = self.centre.params.get(axis, spec.default)
        if spec.allowed or spec.type == "bool":
            return [v for v in (spec.allowed or (True, False)) if v != cur]
        u = self.space._to_unit(spec, cur)
        h = self.book.axes[axis].step
        out, seen = [], set()
        for sc in scales:
            for sgn in (1, -1):
                v = self.space._from_unit(spec, u + sgn * sc * h)
                if v != cur and v not in seen and self.book.allowed(axis, v):
                    seen.add(v)
                    out.append(v)
        return out

    def refill(self, rng: random.Random, others: Sequence[Config], now: float, ranker=None,
               want_pairs: int = 40) -> int:
        """Rebuild the move pool around the current centre. Returns how many moves were added."""
        out: List[Move] = []
        live = self.book.live(now)
        for axis in live:                                   # steps: nearest scale first
            for i, v in enumerate(self.step_values(axis)):
                self._add("step", {axis: v}, 0 if i < 2 else 1, out)
        for changes in self._merges():
            self._add("merge", changes, 0, out)
        for changes in self._transplants(list(others) + self.sources[-6:], live):
            self._add("transplant", changes, 1, out)
        if self.probed_at_generation != self.generation and (self.generation == 0 or self.stalls > 0):
            self.probed_at_generation = self.generation
            for axis in live:
                for v in self._probe_values(axis):
                    self._add("probe", {axis: v}, 2, out)
        pairs: List[Dict[str, object]] = []
        axes = [a for a in live if self.step_values(a)]
        for _ in range(want_pairs * 3):
            if len(axes) < 2:
                break
            a, b = rng.sample(axes, 2)
            va, vb = self.step_values(a)[:2], self.step_values(b)[:2]
            pairs.append({a: rng.choice(va), b: rng.choice(vb)})
        if pairs:
            cands = [(ch, self._variant(ch)) for ch in pairs]
            cands = [(ch, c) for ch, c in cands if c is not None]
            if ranker is not None and ranker.ready and cands:
                order = ranker.rank([c for _, c in cands], best=self.best)
                cands = [cands[i] for i in order]
            for ch, c in cands[:want_pairs]:
                if self._fresh(c):
                    self.tried.add(c.key)
                    out.append(Move("pair", c, tuple(sorted(ch)), self.centre.key, self.best, 3))
        out.sort(key=lambda m: m.priority)
        self.pool.extend(out)
        self.stalls = 0 if out else self.stalls + 1
        self.refilled_generation = self.generation
        return len(out)

    def _merges(self) -> List[Dict[str, object]]:
        wins = [(a, v, g) for a, v, g, gen in self.winners if gen >= self.generation - 1]
        best_per_axis: Dict[str, Tuple[object, float]] = {}
        for a, v, g in wins:
            if a not in best_per_axis or g > best_per_axis[a][1]:
                best_per_axis[a] = (v, g)
        ranked = sorted(best_per_axis.items(), key=lambda kv: -kv[1][1])
        ranked = [(a, v) for a, (v, _) in ranked
                  if not self.space._same(self.space.specs[a], self.centre.params.get(a), v)]
        out: List[Dict[str, object]] = []
        for n in range(1, len(ranked) + 1):                # a winner re-applied on the new centre, then acc2 .. accN
            out.append(dict(ranked[:n]))
        for (a, va), (b, vb) in itertools.combinations(ranked[:4], 2):
            out.append({a: va, b: vb})
        return out

    def _transplants(self, others: Sequence[Config], live: Sequence[str]) -> List[Dict[str, object]]:
        out: List[Dict[str, object]] = []
        for other in others:
            diff = [(a, other.params.get(a)) for a in live
                    if not self.space._same(self.space.specs[a], self.centre.params.get(a), other.params.get(a))
                    and self.book.allowed(a, other.params.get(a))]
            for a, v in diff:
                out.append({a: v})
            for (a, va), (b, vb) in itertools.combinations(diff[:4], 2):
                out.append({a: va, b: vb})
        return out

    def _probe_values(self, axis: str) -> List[object]:
        spec = self.space.specs[axis]
        if spec.allowed or spec.type == "bool":
            return []
        cur = self.centre.params.get(axis, spec.default)
        out, seen = [], set()
        for i in range(PROBE_POINTS):
            v = self.space._from_unit(spec, i / (PROBE_POINTS - 1))
            if v != cur and v not in seen and self.book.allowed(axis, v):
                seen.add(v)
                out.append(v)
        return out

    def take(self) -> Optional[Move]:
        if not self.pool:
            return None
        m = self.pool.pop(0)
        self.inflight += 1
        self.gen_sent += 1
        return m
