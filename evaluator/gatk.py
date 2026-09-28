"""The real evaluator. Owner: evaluator. Not implemented yet.

Plan:
  run():   GATK 4.5.0.0 HaplotypeCaller with exactly the validator's command line
           (minos_subnet/templates/gatk.py). Downstream-only axes
           (stand_call_conf, heterozygosity, indel_heterozygosity,
           max_alternate_alleles) via -ERC GVCF once + GenotypeGVCFs per config
           (legacy: exact, 7.7x cheaper per variant) - must be re-proven here.
  score(): fast in-memory matcher producing hap.py-equivalent TP/FP/FN per class,
           validated against real hap.py (vcfeval engine) on released rounds;
           counts go through common.scoring.finalize, never a local formula.
  KPIs:    seconds/eval, score exactness vs hap.py path (<1e-9), evals/s/core.
"""
from __future__ import annotations

from common.types import AnswerSet, CallSet, Config, Round, ScoreResult


class GatkEvaluator:
    def run(self, round: Round, config: Config) -> CallSet:
        raise NotImplementedError("evaluator.GatkEvaluator.run")

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        raise NotImplementedError("evaluator.GatkEvaluator.score")
