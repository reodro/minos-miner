"""The real runner and scorer. Owner: evaluator.

Verified in test/ (2026-09-28; 3 rounds x 3 configs against the validator-identical baseline):
  used     - SAFE_ANNOTATION_EXCLUSIONS (-AX): bit-identical calls, x1.135.
             Excluding QualByDepth or DepthPerAlleleBySample changes QUAL; RMSMappingQuality is required.
           - 1 thread per job: x2.06 throughput vs 4 threads per job on a saturated box; calls identical.
  rejected - GVCF re-genotyping (heterozygosity changes calls; QUAL shifts flip stand_call_conf thresholds in ~13%
             of cases), QUAL projection (1e-4..1e-3 off), native java outside the image (different calls and slower:
             GKL PairHMM missing), interval restriction (~2%), scatter-gather (x0.49 throughput).
Workers must run the validator's docker image; a worker that differs silently changes calls.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Optional

from common.scoring import finalize
from common.settings import SETTINGS
from common.types import CLASSES, FULL, AnswerSet, CallSet, ClassCounts, Config, Round, ScoreResult

if SETTINGS.minos_subnet_dir not in sys.path:
    sys.path.insert(0, SETTINGS.minos_subnet_dir)

GATK_IMAGE = "broadinstitute/gatk:4.5.0.0"
SAFE_ANNOTATION_EXCLUSIONS = (
    "Coverage", "DepthPerSampleHC", "ChromosomeCounts", "BaseQualityRankSumTest", "ExcessHet", "FisherStrand",
    "InbreedingCoeff", "MappingQualityRankSumTest", "ReadPosRankSumTest", "StrandOddsRatio", "StrandBiasBySample",
)


def records_digest(vcf: Path) -> str:
    """Hash of the VCF body only: identical calls give identical digests whatever the header says."""
    h = hashlib.sha256()
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                h.update(line.encode())
    return h.hexdigest()


class GatkRunner:
    """The validator's HaplotypeCaller command (minos_subnet/templates/gatk.py) plus the verified speed-ups."""

    def __init__(self, work_dir: Optional[str] = None, threads: int = 1, memory_gb: int = 4, timeout_s: int = 1800):
        self.work_dir = Path(work_dir or os.path.join(SETTINGS.cache_dir, "callsets"))
        self.threads, self.memory_gb, self.timeout_s = threads, memory_gb, timeout_s

    def command(self, round: Round, config: Config, fidelity: str, out_vcf: Path) -> list:
        from templates.tool_params import validate_and_build_flags   # upstream
        flags = validate_and_build_flags(config.tool, config.params)
        if not flags["valid"]:
            raise ValueError("; ".join(flags["errors"]))
        bam, ref = Path(round.bam).resolve(), Path(round.reference).resolve()
        cmd = ["docker", "run", "--rm", "--name", f"minos-run-{uuid.uuid4().hex[:10]}",
               f"--cpus={self.threads}", f"--memory={self.memory_gb}g",
               "-v", f"{bam.parent}:/data/bams:ro", "-v", f"{ref.parent}:/data/reference:ro",
               "-v", f"{out_vcf.parent}:/data/output", GATK_IMAGE,
               "gatk", "--java-options", f"-Xmx{max(512, int(self.memory_gb * 1024 * 0.8))}m", "HaplotypeCaller",
               "-R", f"/data/reference/{ref.name}", "-I", f"/data/bams/{bam.name}",
               "-O", f"/data/output/{out_vcf.name}", "--native-pair-hmm-threads", str(self.threads),
               "-L", round.region if fidelity == FULL else fidelity]
        for f in flags["flags"]:
            cmd.extend(str(f).split())
        for a in SAFE_ANNOTATION_EXCLUSIONS:
            cmd += ["-AX", a]
        return cmd

    def run(self, round: Round, config: Config, fidelity: str = FULL) -> CallSet:
        if config.tool != "gatk":
            raise ValueError(f"GatkRunner cannot run {config.tool!r}")
        fid = "full" if fidelity == FULL else hashlib.sha1(fidelity.encode()).hexdigest()[:8]
        out = self.work_dir / round.tag / f"{config.key}_{fid}.vcf.gz"
        out.parent.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        p = subprocess.run(self.command(round, config, fidelity, out), capture_output=True, text=True,
                           timeout=self.timeout_s)
        dt = time.time() - t0
        if p.returncode != 0 or not out.exists():
            err = (p.stderr or "").strip()
            return CallSet(round.round_id, config.key, "", "", dt, ok=False, error=f"rc={p.returncode} {err[-600:]}",
                           fidelity=fidelity)
        return CallSet(round.round_id, config.key, str(out), records_digest(out), dt, fidelity=fidelity)


class HappyScorer:
    """Exact: upstream HappyScorer.score_vcf + parse_happy_vcf, counts through common.scoring.finalize.

    Truth comes from the answer set: `notes['truth_vcf']` / `notes['mutations_vcf']` when the provider has files
    (the oracle), otherwise a VCF written from its records. Slow (~8 s); to be replaced by the per-cluster
    decision cache once that is proven identical to this.
    """

    def __init__(self, work_dir: Optional[str] = None):
        self.work_dir = Path(work_dir or os.path.join(SETTINGS.cache_dir, "happy"))

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        from utils.scoring import HappyScorer as UpstreamHappy, difficulty_class_counts, parse_happy_vcf
        truth = answer_set.notes.get("truth_vcf")
        mutations = answer_set.notes.get("mutations_vcf")
        region = answer_set.notes.get("region")
        reference = answer_set.notes.get("reference")
        if not (truth and mutations and region and reference):
            raise NotImplementedError("HappyScorer needs truth/mutations VCFs, region and reference in answer_set.notes "
                                      "(writing them from reconstructed records is not implemented yet)")
        chrom = region.split(":")[0]
        sdf = str(Path(reference).parent / f"{chrom}.sdf")
        out = self.work_dir / uuid.uuid4().hex[:12]
        out.mkdir(parents=True, exist_ok=True)
        m = UpstreamHappy().score_vcf(truth_vcf=truth, query_vcf=callset.vcf, reference_fasta=reference,
                                      confident_bed=None, region=region, reference_sdf=sdf, mutations_vcf=mutations)
        if not m or "happy_vcf_path" not in m:
            raise RuntimeError("hap.py produced no result")
        recs = parse_happy_vcf(m["happy_vcf_path"], truth_vcf_path=truth)
        cc = difficulty_class_counts(recs)
        ratios = {k: m[k] for k in ("titv_truth_snp", "titv_query_snp", "hethom_truth_snp", "hethom_query_snp") if k in m}
        n_targets = answer_set.n_targets or 1
        result = ScoreResult(callset.round_id, callset.config_key, answer_set.version,
                             {k: ClassCounts(**cc[k]) for k in CLASSES}, int(m.get("region_fp_total", 0)), n_targets,
                             ratios, fidelity=callset.fidelity, job_id=callset.job_id)
        return finalize(result)
