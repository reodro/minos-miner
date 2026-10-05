"""Deterministic fake runner + scorer for developing and testing main without GATK.

Scores come from a smooth toy landscape over normalised numeric parameters with a hidden optimum, turned into
plausible per-class counts and passed through the real formula. It says nothing about real configs.
"""
from __future__ import annotations

import hashlib
import math
import random

from common.params import ParamRanges
from common.scoring import finalize
from common.types import CLASSES, FULL, AnswerSet, CallSet, ClassCounts, Config, Round, ScoreResult

_TRUTH_PER_CLASS = {"snp_hom": 5000, "snp_het": 11000, "indel_1bp": 1600,
                    "indel_2_3bp": 800, "indel_4_7bp": 600, "indel_8bp": 250}


class FakeRunner:
    def __init__(self, ranges: ParamRanges, seed: int = 0, batching: bool = False, hills: int = 1):
        """batching: expose GATK's branch key and variant costs, so an executor batches riders as the fleet does
        (the engine's tests); off, every config is its own run (the older strategies' tests). hills: the landscape
        has this many peaks (the second and later ones a little lower), so births, valleys and merges are tested."""
        self.ranges, self.seed, self.batching, self.hills = ranges, seed, batching, max(1, hills)
        self.configs = {}

    def quality(self, config: Config) -> float:
        specs = self.ranges.tools[config.tool]
        best = 0.0
        for h in range(self.hills):
            rng = random.Random(f"{self.seed}:{config.tool}" + (f":hill{h}" if h else ""))
            d2 = 0.0
            for name in sorted(specs):
                s = specs[name]
                if s.type not in ("int", "float") or s.min is None or s.max is None or s.max == s.min:
                    continue
                opt = rng.random()
                d2 += ((float(config.params.get(name, s.default)) - s.min) / (s.max - s.min) - opt) ** 2
            best = max(best, (1.0 - 0.03 * h) * math.exp(-d2))
        return best

    def branch_key(self, config: Config, fidelity: str = FULL):
        """Like GatkRunner: configs differing only in branchable (genotyping/PairHMM/assembly/activity) params
        share a traversal in the simulator. None (nothing shares) unless `batching`."""
        if not self.batching:
            return None
        from evaluator.gatk import BRANCH_PARAMS
        return Config(config.tool, {k: v for k, v in config.params.items() if k not in BRANCH_PARAMS}).key + fidelity

    def variant_cost(self, main: Config, variant: Config, fidelity: str = FULL):
        from evaluator.gatk import BRANCH_PARAMS, VARIANT_COST, _TIER_OF
        if not self.batching or self.branch_key(main, fidelity) != self.branch_key(variant, fidelity):
            return None
        changed = [p for p in BRANCH_PARAMS if main.params.get(p) != variant.params.get(p)]
        return max((VARIANT_COST[_TIER_OF[p]] for p in changed), default=0.0)

    def run(self, round: Round, config: Config, fidelity: str = FULL) -> CallSet:
        self.configs[config.key] = config
        digest = hashlib.sha256(f"{round.round_id}:{config.key}:{fidelity}".encode()).hexdigest()
        return CallSet(round.round_id, config.key, f"fake://{config.key}", digest, 0.0, fidelity=fidelity)


class FakeScorer:
    def __init__(self, runner: FakeRunner):
        self.runner = runner

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        q = self.runner.quality(self.runner.configs[callset.config_key])
        # a later answer-set version shifts every score slightly, so rescoring is observable in tests
        q = q * (1 - 0.001 * (answer_set.version - 1))
        counts = {}
        for k in CLASSES:
            t = _TRUTH_PER_CLASS[k]
            tp = int(t * (0.85 + 0.14 * q))
            counts[k] = ClassCounts(tp=tp, fn=t - tp, fp=int(40 * (1 - q)))
        return finalize(ScoreResult(callset.round_id, callset.config_key, answer_set.version, counts,
                                    region_fp=300 + sum(c.fp for c in counts.values()),
                                    n_targets=answer_set.n_targets or 600, fidelity=callset.fidelity,
                                    job_id=callset.job_id))
