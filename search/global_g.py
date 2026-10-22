"""Generalized Step A (offline test G, 2026-10-22): one wave that covers all searched settings, not portfolio lines.

The domain is cut into sub-domains along the settings that need their own GATK run and matter most - min_pruning
(4 bins), assembly_region_padding (3) and active_probability_threshold (3): 36 sub-domains. Each gets the same number
of runs. Inside a sub-domain the other run-level settings (activity, assembly and full-run settings) are spread with a
Latin hypercube, so every part of each setting's range gets runs and their combinations are mixed. Each run carries
PER_RUN riders that vary only the cheap settings (genotyping and PairHMM: they ride inside one GATK run, ~0.8-13 s
each), also spread with a Latin hypercube. A few portfolio anchors keep the proven regions in the scan.

Today's Step A (search/global_a.py) sweeps 8 settings one at a time around 24 portfolio starts; in round 04:58 it saw
1.2% of assembly_region_padding's values and 2.5-9% of the region-size and haplotype settings.
"""
from __future__ import annotations

import asyncio
import itertools
import logging
import math
import random
from typing import Dict, List, Optional, Sequence, Tuple

from common.tuning import tune
from common.types import Config
from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS
from search.steps import Ladders

log = logging.getLogger("search.global_a")

PARTS: Dict[str, int] = {"min_pruning": 4, "assembly_region_padding": 3, "active_probability_threshold": 3}
# Fleet-size scaling: the finest split that still gives every sub-domain MIN_PER_CELL runs (one wave of runs).
SPLITS = ({"min_pruning": 4, "assembly_region_padding": 3, "active_probability_threshold": 3},    # 36 sub-domains
          {"min_pruning": 2, "assembly_region_padding": 3, "active_probability_threshold": 2},    # 12
          {"min_pruning": 2, "assembly_region_padding": 2, "active_probability_threshold": 1})    # 4
MIN_PER_CELL = 8
MIN_RUNS = 40            # fewer runs than this: the generalized scan is too thin (the caller uses today's Step A)
RIDERS = GENOTYPING_PARAMS + PAIRHMM_PARAMS
PER_RUN = 7
TAG = "G"


def lhs(n: int, k: int, rng: random.Random) -> List[Tuple[float, ...]]:
    """n points in [0, 1)^k, one per stratum of every axis (Latin hypercube)."""
    cols = []
    for _ in range(k):
        perm = list(range(n))
        rng.shuffle(perm)
        cols.append([(p + rng.random()) / n for p in perm])
    return list(zip(*cols)) if k else [()] * n


class GeneralizedA:
    def __init__(self, client, space, length_s: float, runs: int, seed: int = 0, anchors: Sequence[Config] = ()):
        self.client, self.space = client, space
        self.rng = random.Random(seed)
        self.t0 = client.clock()
        self.stop_at = self.t0 + length_s
        self.runs = max(1, int(runs))
        self.anchors = list(anchors)
        lad = Ladders(space)
        self.values: Dict[str, list] = {}
        for n in space.names:
            s = space.specs[n]
            vals = list(s.allowed) if s.allowed else list(lad.rungs.get(n) or [])
            if not vals and s.type == "bool":
                vals = [False, True]
            if len(vals) > 1:
                self.values[n] = vals
        split = next((sp for sp in SPLITS if self.runs >= MIN_PER_CELL * math.prod(sp.values())), SPLITS[-1])
        self.split = {n: b for n, b in split.items() if n in self.values}
        self.parts = list(self.split)
        self.riders = [n for n in RIDERS if n in self.values]
        self.base = [n for n in self.values if n not in self.parts and n not in self.riders]
        self.scores: Dict[str, Tuple[Config, float]] = {}
        self.cell_of: Dict[str, tuple] = {}
        self.report: Dict[str, object] = {}

    def assemble(self, params) -> Optional[Config]:
        c = self.space._assemble(params)
        return None if c is None or self.client.ranges.errors(c) else c

    @staticmethod
    def pick(vals: list, u: float):
        return vals[min(len(vals) - 1, int(u * len(vals)))]

    def bin_values(self, n: str, b: int) -> list:
        vals, nb = self.values[n], self.split[n]
        lo = b * len(vals) // nb
        return vals[lo:max(lo + 1, (b + 1) * len(vals) // nb)]

    def riders_for(self, base: Dict[str, object], points: Sequence[Tuple[float, ...]]) -> List[Config]:
        out = []
        for rr in points:
            c = self.assemble(dict(base, **{n: self.pick(self.values[n], u) for n, u in zip(self.riders, rr)}))
            if c is not None:
                out.append(c)
        return out

    def plan(self) -> List[Tuple[tuple, Config]]:
        cells = list(itertools.product(*[range(self.split[n]) for n in self.parts]))
        per_cell = max(1, (self.runs - len(self.anchors)) // max(1, len(cells)))
        seed_params = dict(self.anchors[0].params) if self.anchors else {}
        out: List[Tuple[tuple, Config]] = []
        for cell in cells:
            U = lhs(per_cell, len(self.parts) + len(self.base), self.rng)
            R = lhs(per_cell * PER_RUN, len(self.riders), self.rng)
            for i, u in enumerate(U):
                p = dict(seed_params)
                for j, n in enumerate(self.parts):
                    p[n] = self.pick(self.bin_values(n, cell[j]), u[j])
                for j, n in enumerate(self.base):
                    p[n] = self.pick(self.values[n], u[len(self.parts) + j])
                for c in self.riders_for(p, R[i * PER_RUN:(i + 1) * PER_RUN]):
                    out.append((cell, c))
        for a in self.anchors:                              # portfolio anchors, each with its own riders
            for c in [a] + self.riders_for(a.params, lhs(PER_RUN - 1, len(self.riders), self.rng)):
                out.append((("anchor",), c))
        return out

    def cell(self, params: Dict[str, object]) -> tuple:
        """The sub-domain a config falls in (for any config, also one from another Step A)."""
        key = []
        for n in self.parts:
            vals, v = self.values[n], params.get(n)
            i = vals.index(v) if v in vals else min(range(len(vals)), key=lambda k: abs(float(vals[k]) - float(v)))
            key.append(min(self.split[n] - 1, i * self.split[n] // len(vals)))
        return tuple(key)

    def emptying(self, back: int, sent: int) -> bool:
        """islands.global_a_tail (the idle fix of today's Step A, search/global_a.py GlobalA.emptying): the scan ends
        once 60% of it is back and that share of the workers is idle; its last runs still reach the store."""
        from search.global_a import GlobalA
        return GlobalA.emptying(self, back, sent)

    async def run(self) -> Dict[str, float]:
        planned = self.plan()
        uniq: Dict[str, Tuple[tuple, Config]] = {}
        for cell, c in planned:
            uniq.setdefault(c.key, (cell, c))
        configs = [c for _, c in uniq.values()]
        self.cell_of = {k: cell for k, (cell, _) in uniq.items()}
        futs = self.client.propose_many(configs, tag=TAG)
        pending = set(futs)
        quorum = float(tune("islands.global_a_quorum"))
        need = len(futs) - math.ceil(quorum * len(futs))
        left = self.stop_at - self.client.clock()
        while len(pending) > need and left > 0 and not self.emptying(len(futs) - len(pending), len(futs)):
            _, pending = await asyncio.wait(pending, timeout=min(left, 5.0), return_when=asyncio.FIRST_COMPLETED)
            left = self.stop_at - self.client.clock()
        out: Dict[str, float] = {}
        for c, f in zip(configs, futs):
            if f.done() and not f.cancelled() and f.exception() is None:
                o = f.result()
                if o is not None and o.ok:
                    out[c.key] = o.score
                    self.scores[c.key] = (c, o.score)
        self.report = {"split": dict(self.split), "planned": len(planned), "sent": len(configs), "scored": len(out),
                       "cells": len({cell for cell, _ in uniq.values()}), "cut": bool(pending),
                       "minute": round((self.client.clock() - self.t0) / 60, 1),
                       "best": round(max(out.values()), 4) if out else None}
        log.info("generalized A: %d planned, %d sent, %d scored over %d sub-domains%s; best %s", len(planned),
                 len(configs), len(out), self.report["cells"], ", cut at the stop time" if pending else "",
                 self.report["best"])
        return out
