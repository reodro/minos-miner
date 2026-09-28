"""The real runner and scorer. Owner: evaluator. Not implemented yet.

Verified so far (test/, 2026-09-28; 3 rounds x 3 configs, validator-identical baseline):
  - skip annotations, SAFE set only (SAFE_ANNOTATION_EXCLUSIONS): bit-identical calls, x1.135.
    Excluding QualByDepth or DepthPerAlleleBySample changes QUAL; RMSMappingQuality is required.
  - GVCF re-genotyping: NOT exact (heterozygosity changes calls; QUAL shifts flip stand_call_conf
    thresholds in ~13% of cases). Usable only as an approximate screen.
  - QUAL projection: NOT exact (1e-4..1e-3 off).
"""
from __future__ import annotations

from common.types import AnswerSet, CallSet, Config, Round, ScoreResult

SAFE_ANNOTATION_EXCLUSIONS = (
    "Coverage", "DepthPerSampleHC", "ChromosomeCounts", "BaseQualityRankSumTest", "ExcessHet", "FisherStrand",
    "InbreedingCoeff", "MappingQualityRankSumTest", "ReadPosRankSumTest", "StrandOddsRatio", "StrandBiasBySample",
)


class GatkRunner:
    """Validator command (minos_subnet/templates/gatk.py) + `-AX` for each SAFE_ANNOTATION_EXCLUSIONS entry."""

    def run(self, round: Round, config: Config, fidelity: str) -> CallSet:
        raise NotImplementedError("evaluator.GatkRunner.run")


class HappyScorer:
    """Exact: upstream HappyScorer.score_vcf + parse_happy_vcf + compute_score_v2 (to be replaced by the
    per-cluster decision cache once it is proven identical)."""

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        raise NotImplementedError("evaluator.HappyScorer.score")
