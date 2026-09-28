"""Deterministic fake evaluator for developing and testing main without GATK.

Scores come from a smooth toy landscape over normalised numeric parameters with a
hidden optimum, turned into plausible per-class counts and passed through the real
formula. It says nothing about real configs - it only exercises the plumbing.
"""
from __future__ import annotations

import hashlib
import math
import random

from common.params import ParamRanges
from common.scoring import finalize
from common.types import CLASSES, AnswerSet, CallSet, ClassCounts, Config, Round, ScoreResult

_TRUTH_PER_CLASS = {"snp_hom": 5000, "snp_het": 11000, "indel_1bp": 1600,
                    "indel_2_3bp": 800, "indel_4_7bp": 600, "indel_8bp": 250}


class FakeEvaluator:
    def __init__(self, ranges: ParamRanges, seed: int = 0, run_seconds: float = 0.0):
        self.ranges = ranges
        self.seed = seed
        self.run_seconds = run_seconds
        self._configs = {}

    def _quality(self, config: Config) -> float:
        specs = self.ranges.tools[config.tool]
        rng = random.Random(f"{self.seed}:{config.tool}")
        d2 = 0.0
        for name in sorted(specs):
            s = specs[name]
            if s.type not in ("int", "float") or s.min is None or s.max is None or s.max == s.min:
                continue
            opt = rng.random()
            v = config.params.get(name, s.default)
            d2 += ((float(v) - s.min) / (s.max - s.min) - opt) ** 2
        return math.exp(-d2)

    def run(self, round: Round, config: Config) -> CallSet:
        if self.run_seconds:
            import time
            time.sleep(self.run_seconds)
        self._configs[config.key] = config
        digest = hashlib.sha256(f"{round.round_id}:{config.key}".encode()).hexdigest()
        return CallSet(round.round_id, config.key, vcf=f"fake://{config.key}", sha256=digest, runtime_s=self.run_seconds)

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        q = self._quality(self._configs[callset.config_key])
        counts = {}
        for k in CLASSES:
            t = _TRUTH_PER_CLASS[k]
            tp = int(t * (0.85 + 0.14 * q))
            counts[k] = ClassCounts(tp=tp, fn=t - tp, fp=int(40 * (1 - q)))
        region_fp = 300 + sum(c.fp for c in counts.values())
        return finalize(ScoreResult(callset.round_id, callset.config_key, answer_set.version, counts,
                                    region_fp=region_fp, n_targets=answer_set.n_targets or 600))
