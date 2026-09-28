"""The validator's v2 formula, imported from the upstream repo - never reimplemented.

Legacy lesson: a reimplementation of the scoring formula carried 5 errors and a
-0.025 bias. Everything that turns counts into a score goes through here.
Matching calls against an answer set (producing the counts) is evaluator's job.
"""
from __future__ import annotations

import sys
from typing import Dict, Optional, Tuple

from common.settings import SETTINGS
from common.types import CLASSES, ClassCounts, ScoreResult

if SETTINGS.minos_subnet_dir not in sys.path:
    sys.path.insert(0, SETTINGS.minos_subnet_dir)

from utils.scoring import (  # noqa: E402  (upstream minos_subnet)
    DIFFICULTY_WEIGHTS,
    AdvancedScorer,
    difficulty_weighted_f1,
    plausibility_gate,
)

assert tuple(DIFFICULTY_WEIGHTS) == CLASSES, "upstream difficulty classes changed - update common.types.CLASSES"


def v2_score(counts: Dict[str, ClassCounts], region_fp: int, n_targets: int,
             ratios: Optional[Dict[str, float]] = None) -> Tuple[Optional[float], Optional[float], bool]:
    """Return (core, score_0_100, gate_passed) exactly as a validator computes them.

    `ratios` carries titv_truth_snp / titv_query_snp / hethom_truth_snp / hethom_query_snp.
    Without them upstream skips the gate (truth ratio unmeasurable) - fine for ranking
    during search, never for claims about the platform score.
    """
    cd = {k: {"tp": c.tp, "fp": c.fp, "fn": c.fn} for k, c in counts.items()}
    metrics = {"fp_per_target": region_fp / max(float(n_targets), 1.0), **(ratios or {})}
    return difficulty_weighted_f1(cd), AdvancedScorer.compute_score_v2(metrics, cd), plausibility_gate(metrics)


def finalize(result: ScoreResult) -> ScoreResult:
    """Fill core/score/gate_passed on a ScoreResult from its counts."""
    result.core, result.score, result.gate_passed = v2_score(
        result.counts, result.region_fp, result.n_targets, result.ratios)
    return result
