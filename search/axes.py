"""What the engine knows about the axes of the config space: which ones move the score, their step sizes, where
they crash, and the designed opening wave over the dominant ones. Owner: search.

Structural axes (everything but the genotyping tier) are what a basin climbs; genotyping axes are swept as riders
on every traversal (families.py). Evidence is per round: an axis whose moves keep changing nothing is parked, a
value that keeps failing is a cliff nothing is proposed past.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from common.types import Config
from search.space import GENOTYPING_PARAMS, Space

# The axes that moved the score most in measured rounds (single-axis effects on true scores, 2026-10-03), with the
# levels the opening wave tries. Values outside a round's published range are dropped. Everything else starts at
# the GATK default; what is dead this round is learned, not assumed.
DOMINANT_LEVELS: Dict[str, Sequence] = {
    "min_pruning": (2, 3),
    "pcr_indel_model": ("AGGRESSIVE", "CONSERVATIVE", "NONE"),
    "min_base_quality_score": (15, 25, 35),
    "assembly_region_padding": (13, 50, 100, 200),
    "min_dangling_branch_length": (4, 12, 20),
    "phred_scaled_global_read_mismapping_rate": (30, 45, 60),
    "base_quality_score_threshold": (18, 25, 37),
    "min_mapping_quality_score": (10, 30, 55),
    "active_probability_threshold": (0.001, 0.005, 0.02),
    "max_assembly_region_size": (100, 300),
}
DEAD_AFTER = 8          # single-axis moves in a row that changed nothing: the axis is parked
PARK_S = 600.0          # for this long, then one retest move
STEP0 = 0.12            # first unit step of a numeric axis (12% of its range)
STEP_MIN_FRAC = 0.02    # never below this (ints also never below one unit)
STEP_MAX = 0.4


def structural_names(space: Space) -> List[str]:
    return [n for n in space.names if n not in GENOTYPING_PARAMS]


def far_axes(space: Space, names: Sequence[str], a: Config, b: Config, frac: float = 0.10) -> int:
    """How many of `names` differ by more than `frac` of the axis range (enums and bools: any difference) - the
    legacy project's basin distance: how many knobs separate two centres, which a single RMS number blurs."""
    n = 0
    for name in names:
        s = space.specs[name]
        x, y = a.params.get(name, s.default), b.params.get(name, s.default)
        if s.allowed or s.type == "bool":
            n += x != y
        else:
            n += abs(space._to_unit(s, x) - space._to_unit(s, y)) > frac
    return n


def rms_distance(space: Space, names: Sequence[str], a: Config, b: Config) -> float:
    """Unit-space RMS distance over `names` (0 = same structure)."""
    if not names:
        return 0.0
    tot = 0.0
    for n in names:
        s = space.specs[n]
        tot += (space._to_unit(s, a.params.get(n, s.default)) - space._to_unit(s, b.params.get(n, s.default))) ** 2
    return math.sqrt(tot / len(names))


@dataclass
class Axis:
    name: str
    step: float = STEP0                     # unit step for numeric axes
    fail_streak: int = 0                    # single-axis moves since the last win on this axis
    zero_streak: int = 0                    # moves in a row whose score was identical to the centre's
    parked_until: float = 0.0
    bad_above: Optional[float] = None       # unit value from which moves failed (crash / no calls)
    bad_below: Optional[float] = None


class AxisBook:
    """Per-axis evidence shared by every basin of a round."""

    def __init__(self, space: Space):
        self.space = space
        self.names = structural_names(space)
        self.axes: Dict[str, Axis] = {n: Axis(n) for n in self.names}
        for n in self.names:
            s = space.specs[n]
            if s.type == "int" and s.max is not None and s.min is not None:
                self.axes[n].step = max(STEP0, 1.0 / max(1, s.max - s.min))

    def live(self, now: float) -> List[str]:
        return [n for n in self.names if self.axes[n].parked_until <= now]

    def unit(self, name: str, value) -> float:
        return self.space._to_unit(self.space.specs[name], value)

    def allowed(self, name: str, value) -> bool:
        """False when `value` lies past a learned cliff on this axis."""
        a = self.axes[name]
        s = self.space.specs[name]
        if s.allowed or s.type == "bool":
            return True
        u = self.unit(name, value)
        if a.bad_above is not None and u >= a.bad_above:
            return False
        if a.bad_below is not None and u <= a.bad_below:
            return False
        return True

    # ---- evidence ------------------------------------------------------------------------------------------
    def note_move(self, name: str, value, centre_value, delta: Optional[float], failed: bool, now: float) -> None:
        """One single-axis move's outcome: delta = score - centre score (None when it failed / scored 0)."""
        a = self.axes[name]
        s = self.space.specs[name]
        numeric = not (s.allowed or s.type == "bool")
        if failed:
            if numeric:
                u, uc = self.unit(name, value), self.unit(name, centre_value)
                if u > uc:
                    a.bad_above = u if a.bad_above is None else min(a.bad_above, u)
                elif u < uc:
                    a.bad_below = u if a.bad_below is None else max(a.bad_below, u)
            a.fail_streak += 1
            return
        if delta is not None and delta > 1e-9:
            a.fail_streak, a.zero_streak = 0, 0
            if numeric:
                a.step = min(STEP_MAX, a.step * 1.5)
            return
        a.fail_streak += 1
        if delta is not None and abs(delta) <= 1e-9:
            a.zero_streak += 1
            if a.zero_streak >= DEAD_AFTER:
                a.parked_until = now + PARK_S
                a.zero_streak = 0
        else:
            a.zero_streak = 0
        if numeric and a.fail_streak and a.fail_streak % 6 == 0:
            floor = max(STEP_MIN_FRAC, 1.0 / max(1, (s.max - s.min)) if s.type == "int" else STEP_MIN_FRAC)
            a.step = max(floor, a.step * 0.5)

    def unpark_all(self) -> None:
        for a in self.axes.values():
            a.parked_until = 0.0


# ---- the designed opening wave ------------------------------------------------------------------------------
def design_wave(space: Space, rng: random.Random, n: int, backgrounds: Sequence[Config]) -> List[Config]:
    """`n` configs that spread the dominant axes' levels evenly, in combination, over the round's space: a balanced
    random design (every level of every axis appears about equally often, pairs mixed) rather than a full
    factorial. Other axes come from one of `backgrounds` (the GATK default, a seed), alternating."""
    axes: List[Tuple[str, List]] = []
    for name, levels in DOMINANT_LEVELS.items():
        if name not in space.names:
            continue
        spec = space.specs[name]
        ok = [v for v in levels if (v in spec.allowed) if spec.allowed] if spec.allowed else \
             [v for v in levels if spec.min is None or (spec.min <= v <= spec.max)]
        if len(ok) >= 2:
            axes.append((name, list(ok)))
    if n <= 0 or not axes:
        return []
    out: List[Config] = []
    seen = set()
    # balanced columns: each axis gets its levels repeated and shuffled, so counts per level differ by at most one
    columns = {}
    for name, levels in axes:
        col = [levels[i % len(levels)] for i in range(n)]
        rng.shuffle(col)
        columns[name] = col
    bgs = list(backgrounds) or [space._assemble({})]
    for i in range(n):
        bg = bgs[i % len(bgs)]
        params = dict(bg.params) if bg is not None else {}
        for name, _ in axes:
            params[name] = columns[name][i]
        c = space._assemble(params)
        if c is not None and c.key not in seen:
            seen.add(c.key)
            out.append(c)
    return out


def latin_hypercube(space: Space, rng: random.Random, n: int, names: Optional[Sequence[str]] = None) -> List[Config]:
    """`n` configs that stratify EVERY axis across its whole range (each axis's range cut into n strata, one point per
    stratum, strata paired at random between axes), so every axis is seen at n distinct levels and every pair of axes
    in all its low/mid/high combinations about n/9 times. Categorical axes are stratified over their values. The
    global opening wave: what the dominant-axes design (design_wave) cannot see, a region that needs the other axes
    away from a seed's values, this sees."""
    names = list(names) if names is not None else list(space.names)
    if n <= 0 or not names:
        return []
    columns: Dict[str, List[float]] = {}
    for name in names:
        strata = list(range(n))
        rng.shuffle(strata)
        columns[name] = [(k + rng.random()) / n for k in strata]
    out: List[Config] = []
    seen = set()
    for i in range(n):
        params = {name: space._from_unit(space.specs[name], columns[name][i]) for name in names}
        c = space._assemble(params)
        if c is not None and c.key not in seen:
            seen.add(c.key)
            out.append(c)
    return out
