"""The config space as numbers: features for the surrogate model, random sampling, and local moves.

Every tunable parameter maps to one feature in [0, 1]: numeric ranges linearly, or on a log scale when they span
more than two orders of magnitude (rates, heterozygosities); bools 0/1; enums by position. Parameters with a
single allowed value are fixed and not features, and so are the values in PINNED: settings that can never score.
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional

from common.params import ParamRanges, ParamSpec
from common.types import Config


# GVCF / BP_RESOLUTION make GATK write reference blocks (one record per position); hap.py rejects that output on the
# validator too, so such a config scores nothing (checked on 20260928T092900Z with the validator-exact path).
PINNED = {"gatk": {"emit_ref_confidence": "NONE"}}


class Space:
    def __init__(self, ranges: ParamRanges, tool: str = "gatk", pinned: Optional[Dict[str, object]] = None):
        self.ranges, self.tool = ranges, tool
        self.specs: Dict[str, ParamSpec] = ranges.tools[tool]
        self.pinned = dict(PINNED.get(tool, {}) if pinned is None else pinned)
        self.names: List[str] = [n for n, s in sorted(self.specs.items()) if not self._fixed(s) and n not in self.pinned]

    @staticmethod
    def _fixed(s: ParamSpec) -> bool:
        return (s.allowed and len(s.allowed) == 1) or (s.type in ("int", "float") and s.min == s.max)

    def _log(self, s: ParamSpec) -> bool:
        return s.type == "float" and s.min is not None and s.min > 0 and s.max / s.min > 100

    # ---- config <-> features -------------------------------------------------------------------------------
    def encode(self, config: Config) -> List[float]:
        return [self._to_unit(self.specs[n], config.params.get(n, self.specs[n].default)) for n in self.names]

    def _to_unit(self, s: ParamSpec, v) -> float:
        if s.allowed:
            return s.allowed.index(v) / max(len(s.allowed) - 1, 1) if v in s.allowed else 0.0
        if s.type == "bool":
            return 1.0 if v else 0.0
        lo, hi, x = float(s.min), float(s.max), float(v)
        if self._log(s):
            return (math.log(x) - math.log(lo)) / (math.log(hi) - math.log(lo))
        return (x - lo) / (hi - lo)

    def _from_unit(self, s: ParamSpec, u: float):
        u = min(1.0, max(0.0, u))
        if s.allowed:
            return s.allowed[min(int(u * len(s.allowed)), len(s.allowed) - 1)]
        if s.type == "bool":
            return u >= 0.5
        lo, hi = float(s.min), float(s.max)
        x = math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo))) if self._log(s) else lo + u * (hi - lo)
        return int(round(x)) if s.type == "int" else round(x, 6)

    # ---- sampling ----------------------------------------------------------------------------------------------
    def random(self, rng: random.Random) -> Config:
        """A uniform random config that passes the platform's ranges and cross-parameter constraints."""
        for _ in range(1000):
            params = {n: self._from_unit(s, rng.random()) if n in self.names
                      else self.pinned.get(n, s.allowed[0] if s.allowed else s.default) for n, s in self.specs.items()}
            c = Config(self.tool, params)
            if not self.ranges.errors(c):
                return c
        raise RuntimeError("could not sample a valid config")

    def neighbour(self, config: Config, rng: random.Random, n_params: int = 2, scale: float = 0.15,
                  only: Optional[List[str]] = None) -> Optional[Config]:
        """Move 1..n_params parameters (from `only` if given): numeric by a Gaussian step of `scale` in unit space,
        bool/enum to another value. None if no valid neighbour was found."""
        pool = [n for n in (only or self.names) if n in self.names]
        if not pool:
            return None
        for _ in range(50):
            params = dict(config.params)
            for n in rng.sample(pool, min(len(pool), rng.randint(1, n_params))):
                s = self.specs[n]
                if s.allowed or s.type == "bool":
                    params[n] = self._from_unit(s, rng.random())
                else:
                    params[n] = self._from_unit(s, self._to_unit(s, params.get(n, s.default)) + rng.gauss(0, scale))
            c = Config(self.tool, params)
            if c.key != config.key and not self.ranges.errors(c):
                return c
        return None
