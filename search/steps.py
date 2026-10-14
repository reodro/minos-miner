"""Each setting's ladder: its values in order at the agreed resolution (search/DESIGN.md, B1). A move of k δ is k
rungs up or down the ladder.

    integers                       every value of the round's range (floors and caps of the Space applied)
    min_base_quality_score,        the round's BAM groups (search/bamgroups.py): every rung is a different GATK run
    base_quality_score_threshold
    enums                          the allowed values in their published order; bools False, True
    active_probability_threshold   0.001 .. 0.05, step 0.001
    contamination_fraction_...     0 .. 0.05, step 0.0025; rungs that remove the same number of reads are one rung
    standard_min_confidence_...    30 .. 100, step 1
    heterozygosity                 1e-4 .. 1e-2, 10 per decade (log)
    indel_heterozygosity           1e-5 .. 1e-3, 10 per decade (log)

A value off the ladder (a float from the portfolio, say) sits between two rungs: moving from it starts at the nearest
rung, so a move of +1 δ never stays in place.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

from common.types import Config

FLOAT_STEP = {"active_probability_threshold": 0.001, "contamination_fraction_to_filter": 0.0025,
              "standard_min_confidence_threshold_for_calling": 1.0}
FLOAT_LOG_PER_DECADE = {"heterozygosity": 10, "indel_heterozygosity": 10}


def ladder(space, name: str) -> list:
    s = space.specs[name]
    if s.allowed:
        return list(s.allowed)
    if s.type == "bool":
        return [False, True]
    lo, hi = s.min, s.max
    if s.type == "int":
        return list(range(int(lo), int(hi) + 1))
    if name in FLOAT_STEP:
        step = FLOAT_STEP[name]
        rungs = [round(lo + i * step, 6) for i in range(int(round((hi - lo) / step)) + 1)]
        bounds = getattr(space, "bounds", None)
        facts = getattr(bounds, "facts", None)
        if name == "contamination_fraction_to_filter" and facts is not None and facts.depths:
            from search.bamgroups import contamination_groups
            reps = contamination_groups(facts.depths, float(lo), float(hi), step)
            if reps:
                return list(reps)
        return rungs
    if name in FLOAT_LOG_PER_DECADE:
        k = FLOAT_LOG_PER_DECADE[name]
        return [float(f"{lo * 10 ** (i / k):.4g}") for i in range(int(round(math.log10(hi / lo) * k)) + 1)]
    n = 50                                          # a float with no agreed resolution: 50 rungs
    return [round(lo + (hi - lo) * i / (n - 1), 6) for i in range(n)]


class Ladders:
    def __init__(self, space):
        self.space = space
        self.rungs: Dict[str, list] = {n: ladder(space, n) for n in space.names}

    def movable(self, name: str) -> bool:
        return len(self.rungs.get(name, [])) > 1

    def index(self, name: str, value) -> float:
        """The value's position: exact for a rung; between two rungs, fractional."""
        r = self.rungs[name]
        if value in r:
            return float(r.index(value))
        s = self.space.specs[name]
        if s.allowed or s.type == "bool":
            return 0.0
        u = self.space._to_unit(s, value)
        units = [self.space._to_unit(s, x) for x in r]
        for i in range(len(r) - 1):
            if units[i] <= u <= units[i + 1]:
                span = units[i + 1] - units[i]
                return i + ((u - units[i]) / span if span else 0.0)
        return 0.0 if u < units[0] else float(len(r) - 1)

    def shift(self, name: str, value, k: int):
        """The value k rungs away (k > 0 up, k < 0 down), or None past the end of the ladder."""
        r = self.rungs[name]
        pos = self.index(name, value)
        if pos != int(pos):                         # between rungs: +1 is the rung above, -1 the rung below
            target = (math.ceil(pos) + k - 1) if k > 0 else (math.floor(pos) + k + 1)
        else:
            target = int(pos) + k
        if 0 <= target < len(r):
            return r[target]
        return None

    def rungs_moved(self, name: str, a, b) -> float:
        """Signed rungs from value a to value b."""
        return self.index(name, b) - self.index(name, a)


def moved(space, centre: Config, other: Config) -> List[str]:
    return [n for n in space.names if other.params.get(n) != centre.params.get(n)]
