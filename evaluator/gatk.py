"""The real runner and scorer. Owner: evaluator.

Verified in test/ (2026-09-28; 3 rounds x 3 configs against the validator-identical baseline):
  used     - SAFE_ANNOTATION_EXCLUSIONS (-AX): bit-identical calls, x1.135.
           - batching (gatk_patch): variants differing only in genotyping params cost ~0.8 s extra each; in PairHMM
             params ~13 s per distinct PairHMM set (vs ~73 s full). 54/54 + 62/62 bit-identical to stock runs;
             variants whose own run would draw random numbers differently come back DIVERGED and are re-run.
           - assembly params: ~41 s extra per distinct assembly set (only traversal is shared), 86/86 bit-identical.
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
from typing import List, Optional

from common.scoring import finalize
from common.settings import SETTINGS
from common.types import CLASSES, FULL, AnswerSet, CallSet, ClassCounts, Config, Round, ScoreResult

if SETTINGS.minos_subnet_dir not in sys.path:
    sys.path.insert(0, SETTINGS.minos_subnet_dir)

GATK_IMAGE = "broadinstitute/gatk:4.5.0.0"
GATK_JAR_IN_IMAGE = "/gatk/gatk-package-4.5.0.0-local.jar"
# JVM flags the image's `gatk` wrapper adds (`gatk --dry-run`); used when we call java directly.
WRAPPER_JVM = ["-Dsamjdk.use_async_io_read_samtools=false", "-Dsamjdk.use_async_io_write_samtools=true",
               "-Dsamjdk.use_async_io_write_tribble=false", "-Dsamjdk.compression_level=2"]
# Read only by HaplotypeCaller's genotyping engine (verified in the GATK 4.5.0.0 source and by 54/54 bit-identical
# batch outputs): configs differing only in these can share one patched run.
GENOTYPING_BRANCH_PARAMS = ("standard_min_confidence_threshold_for_calling", "indel_heterozygosity",
                            "max_alternate_alleles", "contamination_fraction_to_filter")
GENOTYPING_DEFAULTS = {"standard_min_confidence_threshold_for_calling": 30.0, "indel_heterozygosity": 0.000125,
                       "max_alternate_alleles": 6, "contamination_fraction_to_filter": 0.0}
# Read only when HaplotypeCaller builds its PairHMM likelihood engine (verified in the source): variants differing
# in these get one extra likelihood computation per distinct set, still sharing traversal and assembly.
PAIRHMM_BRANCH_PARAMS = ("base_quality_score_threshold", "pcr_indel_model", "pair_hmm_gap_continuation_penalty",
                         "phred_scaled_global_read_mismapping_rate")
PAIRHMM_DEFAULTS = {"base_quality_score_threshold": 18, "pcr_indel_model": "CONSERVATIVE",
                    "pair_hmm_gap_continuation_penalty": 10, "phred_scaled_global_read_mismapping_rate": 45}
# Read only when HaplotypeCaller builds its assembler: variants differing in these get their own assembly (on the
# same finalized region) and everything after it, still sharing traversal and region finding.
ASSEMBLY_BRANCH_PARAMS = ("min_pruning", "min_dangling_branch_length", "recover_all_dangling_branches",
                          "max_num_haplotypes_in_population", "adaptive_pruning_initial_error_rate",
                          "pruning_lod_threshold")
ASSEMBLY_DEFAULTS = {"min_pruning": 2, "min_dangling_branch_length": 4, "recover_all_dangling_branches": False,
                     "max_num_haplotypes_in_population": 128, "adaptive_pruning_initial_error_rate": 0.001,
                     "pruning_lod_threshold": 2.302585}
# Read only when the traversal cuts the activity profile into regions (assembly_region_padding is not: it also sets
# how far each read shard reaches): variants differing in these get their own region stream from the shared per-locus
# activity states and reads, then assembly and everything after it.
ACTIVITY_BRANCH_PARAMS = ("active_probability_threshold", "min_assembly_region_size", "max_assembly_region_size")
ACTIVITY_DEFAULTS = {"active_probability_threshold": 0.002, "min_assembly_region_size": 50, "max_assembly_region_size": 300}
BRANCH_PARAMS = GENOTYPING_BRANCH_PARAMS + PAIRHMM_BRANCH_PARAMS + ASSEMBLY_BRANCH_PARAMS + ACTIVITY_BRANCH_PARAMS
# Extra cost of one batch variant, as a fraction of a solo run, by the most expensive stage it changes
# (measured at production load by test/timing_tiers.py). A variant that differs in several tiers pays the highest.
VARIANT_COST = {"genotyping": 0.011, "pairhmm": 0.12, "assembly": 0.55, "activity": 0.59}   # 1.4 / 16.3 / 72 / 78 s of 132 s (14 workers, 16 cores)
_TIER_OF = {**{p: "genotyping" for p in GENOTYPING_BRANCH_PARAMS}, **{p: "pairhmm" for p in PAIRHMM_BRANCH_PARAMS},
            **{p: "assembly" for p in ASSEMBLY_BRANCH_PARAMS}, **{p: "activity" for p in ACTIVITY_BRANCH_PARAMS}}
def variant_line(config: Config, output: str) -> str:
    """One line of the patch's variants file: output path, then every branch parameter (config value or default)."""
    def fmt(v):
        return str(v).lower() if isinstance(v, bool) else str(v)
    defaults = {**GENOTYPING_DEFAULTS, **PAIRHMM_DEFAULTS, **ASSEMBLY_DEFAULTS, **ACTIVITY_DEFAULTS}
    return "\t".join([output] + [fmt(config.params.get(p, defaults[p])) for p in BRANCH_PARAMS]) + "\n"


DIVERGED = "minos:diverged"      # batch variant whose own run would differ; the executor re-runs it alone
DEFAULT_PATCH_JAR = os.path.join(os.path.dirname(__file__), "gatk_patch", "build", "minos-patch.jar")

SAFE_ANNOTATION_EXCLUSIONS = (
    "Coverage", "DepthPerSampleHC", "ChromosomeCounts", "BaseQualityRankSumTest", "ExcessHet", "FisherStrand",
    "InbreedingCoeff", "MappingQualityRankSumTest", "ReadPosRankSumTest", "StrandOddsRatio", "StrandBiasBySample",
)


def gatk_error(stderr: str) -> str:
    """GATK's actual error, not the command line it echoes: the USER ERROR block or the exception line."""
    lines = (stderr or "").splitlines()
    for i, line in enumerate(lines):
        if "A USER ERROR has occurred" in line:
            return " ".join(l.strip() for l in lines[i:i + 3])
    exc = [l.strip() for l in lines if "Exception" in l or l.startswith("ERROR")]
    return exc[0] if exc else " ".join(lines[-3:])[-600:]


def records_digest(vcf: Path) -> str:
    """Hash of the VCF body only: identical calls give identical digests whatever the header says."""
    h = hashlib.sha256()
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                h.update(line.encode())
    return h.hexdigest()


def splice_repair(batch_vcf: Path, regions: Path, repair_vcf: Path, repair_log: Path, out: Path) -> None:
    """Final VCF of a repaired variant: the repair run's header and records (its computed regions) plus the batch
    records of every region the repair run skipped. Regions are GATK assembly regions; each record belongs to the
    region whose span contains its POS."""
    import bisect
    import pysam
    spans = []                                        # (contig, start, end, key) of the batch's regions
    for line in regions.read_text().splitlines():
        key = line.split("\t", 1)[0]
        contig, rng = key.split("|", 1)[0].rsplit(":", 1)
        a, b = rng.split("-")
        spans.append((contig, int(a), int(b), key))
    spans.sort()
    starts = [(c, a) for c, a, _, _ in spans]
    skipped = {l.split("\t", 1)[1] for l in repair_log.read_text().splitlines() if l.startswith("S\t")}

    def region_of(contig, pos):
        i = bisect.bisect_right(starts, (contig, pos)) - 1
        if i < 0 or spans[i][0] != contig or not spans[i][1] <= pos <= spans[i][2]:
            raise ValueError(f"batch record {contig}:{pos} lies in no logged region")
        return spans[i][3]

    header, records = [], []
    with gzip.open(repair_vcf, "rt") as f:
        for line in f:
            (header if line.startswith("#") else records).append(line)
    contig_order = {}
    for line in header:
        if line.startswith("##contig=<ID="):
            contig_order[line[13:].split(",", 1)[0].rstrip(">\n")] = len(contig_order)
    with gzip.open(batch_vcf, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                c, pos = line.split("\t", 2)[:2]
                if region_of(c, int(pos)) in skipped:
                    records.append(line)
    records.sort(key=lambda l: (contig_order.get(l.split("\t", 1)[0], 1 << 30), int(l.split("\t", 2)[1])))
    plain = out.with_name(out.name[:-3])              # .vcf
    plain.write_text("".join(header + records))
    out.unlink(missing_ok=True)
    pysam.tabix_compress(str(plain), str(out), force=True)
    plain.unlink()
    pysam.tabix_index(str(out), preset="vcf", force=True)


class GatkRunner:
    """The validator's HaplotypeCaller command (minos_subnet/templates/gatk.py) plus the verified speed-ups."""

    def __init__(self, work_dir: Optional[str] = None, threads: int = 1, memory_gb: int = 4, timeout_s: int = 1800,
                 patch_jar: Optional[str] = None, container_prefix: str = "minos-run"):
        """patch_jar: our patched GATK classes (evaluator/gatk_patch), put before the GATK jar on the classpath,
        still inside the validator's image. None = the stock `gatk` wrapper."""
        self.work_dir = Path(work_dir or os.path.join(SETTINGS.cache_dir, "callsets"))
        self.threads, self.memory_gb, self.timeout_s = threads, memory_gb, timeout_s
        self.patch_jar = Path(patch_jar).resolve() if patch_jar else None
        self.container_prefix = container_prefix     # names this runner's containers, so its owner can stop them

    def command(self, round: Round, config: Config, fidelity: str, out_vcf: Path, jvm_props: Optional[dict] = None) -> list:
        from templates.tool_params import validate_and_build_flags   # upstream
        flags = validate_and_build_flags(config.tool, config.params)
        if not flags["valid"]:
            raise ValueError("; ".join(flags["errors"]))
        bam, ref = Path(round.bam).resolve(), Path(round.reference).resolve()
        heap = f"-Xmx{max(512, int(self.memory_gb * 1024 * 0.8))}m"
        cmd = ["docker", "run", "--rm", "--name", f"{self.container_prefix}-{uuid.uuid4().hex[:10]}",
               f"--cpus={self.threads}", f"--memory={self.memory_gb}g",
               "-v", f"{bam.parent}:/data/bams:ro", "-v", f"{ref.parent}:/data/reference:ro",
               "-v", f"{out_vcf.parent}:/data/output"]
        if self.patch_jar is None:
            cmd += [GATK_IMAGE, "gatk", "--java-options", heap, "HaplotypeCaller"]
        else:
            props = [f"-D{k}={v}" for k, v in (jvm_props or {}).items()]
            cmd += ["-v", f"{self.patch_jar.parent}:/patch:ro", GATK_IMAGE, "java", *WRAPPER_JVM, heap, *props,
                    "-cp", f"/patch/{self.patch_jar.name}:{GATK_JAR_IN_IMAGE}", "org.broadinstitute.hellbender.Main",
                    "HaplotypeCaller"]
        cmd += ["-R", f"/data/reference/{ref.name}", "-I", f"/data/bams/{bam.name}",
               "-O", f"/data/output/{out_vcf.name}", "--native-pair-hmm-threads", str(self.threads),
               "-L", round.region if fidelity == FULL else fidelity]
        for f in flags["flags"]:
            cmd.extend(str(f).split())
        for a in SAFE_ANNOTATION_EXCLUSIONS:
            cmd += ["-AX", a]
        return cmd

    def _out(self, round: Round, config: Config, fidelity: str) -> Path:
        fid = "full" if fidelity == FULL else hashlib.sha1(fidelity.encode()).hexdigest()[:8]
        out = self.work_dir / round.tag / f"{config.key}_{fid}.vcf.gz"
        out.parent.mkdir(parents=True, exist_ok=True)
        return out

    def run(self, round: Round, config: Config, fidelity: str = FULL) -> CallSet:
        if config.tool != "gatk":
            raise ValueError(f"GatkRunner cannot run {config.tool!r}")
        out = self._out(round, config, fidelity)
        t0 = time.time()
        p = subprocess.run(self.command(round, config, fidelity, out), capture_output=True, text=True,
                           timeout=self.timeout_s)
        dt = time.time() - t0
        if p.returncode != 0 or not out.exists():
            return CallSet(round.round_id, config.key, "", "", dt, ok=False,
                           error=f"rc={p.returncode} {gatk_error(p.stderr)}", fidelity=fidelity)
        return CallSet(round.round_id, config.key, str(out), records_digest(out), dt, fidelity=fidelity)

    # ---- batching: configs differing only in genotyping/PairHMM parameters share one run (needs patch_jar) ----
    def branch_key(self, config: Config, fidelity: str = FULL) -> Optional[str]:
        """Configs with the same key can be evaluated in one run_batch. None = cannot be batched."""
        if self.patch_jar is None or config.tool != "gatk":
            return None
        if config.params.get("emit_ref_confidence", "NONE") != "NONE":
            return None                   # the patch refuses GVCF modes: such a config always runs alone
        shared = {k: v for k, v in config.params.items() if k not in BRANCH_PARAMS}
        return Config("gatk", shared).key + ":" + fidelity

    def variant_cost(self, main: Config, variant: Config, fidelity: str = FULL) -> Optional[float]:
        """Extra cost of evaluating `variant` inside `main`'s batch, as a fraction of a solo run. None = not batchable."""
        key = self.branch_key(main, fidelity)
        if key is None or key != self.branch_key(variant, fidelity):
            return None
        defaults = {**GENOTYPING_DEFAULTS, **PAIRHMM_DEFAULTS, **ASSEMBLY_DEFAULTS, **ACTIVITY_DEFAULTS}
        changed = [p for p in BRANCH_PARAMS if main.params.get(p, defaults[p]) != variant.params.get(p, defaults[p])]
        return max((VARIANT_COST[_TIER_OF[p]] for p in changed), default=0.0)

    def run_batch(self, round: Round, configs: List[Config], fidelity: str = FULL) -> List[CallSet]:
        """One patched run: configs[0] is the main output, the rest are genotyping variants.

        A variant whose own run would have downsampled reads differently (random-stream divergence, see
        gatk_patch) comes back as CallSet(ok=False, error=DIVERGED): its output here is not exact and it must be
        run on its own.
        """
        if len(configs) == 1 or self.patch_jar is None:
            return [self.run(round, c, fidelity) for c in configs]
        keys = {self.branch_key(c, fidelity) for c in configs}
        if len(keys) != 1 or None in keys:
            raise ValueError("run_batch needs configs that differ only in branchable parameters")
        main = self._out(round, configs[0], fidelity)
        outs = [self._out(round, c, fidelity) for c in configs[1:]]
        spec = main.parent / f"variants_{uuid.uuid4().hex[:10]}.tsv"
        spec.write_text("".join(variant_line(c, f"/data/output/{o.name}") for c, o in zip(configs[1:], outs)))
        t0 = time.time()
        try:
            p = subprocess.run(self.command(round, configs[0], fidelity, main,
                                            jvm_props={"minos.gt.variants": f"/data/output/{spec.name}"}),
                               capture_output=True, text=True, timeout=self.timeout_s)
        finally:
            spec.unlink(missing_ok=True)
        dt = time.time() - t0
        if p.returncode != 0 or not main.exists():
            err = f"rc={p.returncode} {gatk_error(p.stderr)}"
            return [CallSet(round.round_id, c.key, "", "", dt, ok=False, error=err, fidelity=fidelity) for c in configs]
        out = [CallSet(round.round_id, configs[0].key, str(main), records_digest(main), dt, fidelity=fidelity)]
        for c, o in zip(configs[1:], outs):
            status = Path(str(o) + ".minos")
            state = status.read_text().strip() if status.exists() else "missing"
            status.unlink(missing_ok=True)
            if state == "exact" and o.exists():
                out.append(CallSet(round.round_id, c.key, str(o), records_digest(o), dt, fidelity=fidelity))
            elif state == "diverged" and o.exists() and Path(str(o) + ".minos.regions").exists():
                # kept for repair(): the batch output is right everywhere except near deep spots
                kept = o.with_name(o.name.replace(".vcf.gz", ".batch.vcf.gz"))
                o.rename(kept)
                Path(str(o) + ".minos.regions").rename(str(kept) + ".minos.regions")
                Path(str(o) + ".tbi").unlink(missing_ok=True)
                out.append(CallSet(round.round_id, c.key, str(kept), "", dt, ok=False, error=DIVERGED, fidelity=fidelity))
            else:
                o.unlink(missing_ok=True)
                out.append(CallSet(round.round_id, c.key, "", "", dt, ok=False,
                                   error=DIVERGED if state == "diverged" else f"variant status {state}",
                                   fidelity=fidelity))
        return out

    def repair(self, round: Round, config: Config, diverged: CallSet, fidelity: str = FULL) -> CallSet:
        """Exact callset for a diverged batch variant, from its batch output plus a repair run.

        The repair run is the variant's own run, so its read pass keeps exactly the reads its real run keeps; every
        region whose reads equal the batch's is skipped (no assembly/likelihoods/genotyping) and taken from the batch
        output. Only regions touched by deep-spot read selection are computed. Falls back to a plain run if the
        batch files are missing.
        """
        batch = Path(diverged.vcf) if diverged.vcf else None
        regions = Path(str(batch) + ".minos.regions") if batch else None
        if batch is None or not batch.exists() or not regions.exists():
            return self.run(round, config, fidelity)
        out = self._out(round, config, fidelity)
        rep_vcf = out.with_name(out.name.replace(".vcf.gz", ".repair.vcf.gz"))
        rep_log = Path(str(rep_vcf) + ".minos.repair")
        t0 = time.time()
        p = subprocess.run(self.command(round, config, fidelity, rep_vcf,
                                        jvm_props={"minos.repair.log": f"/data/output/{regions.name}",
                                                   "minos.repair.out": f"/data/output/{rep_log.name}"}),
                           capture_output=True, text=True, timeout=self.timeout_s)
        dt = time.time() - t0
        if p.returncode != 0 or not rep_vcf.exists() or not rep_log.exists():
            return CallSet(round.round_id, config.key, "", "", dt, ok=False, fidelity=fidelity,
                           error=f"repair rc={p.returncode} {gatk_error(p.stderr)}")
        splice_repair(batch, regions, rep_vcf, rep_log, out)
        for f in (batch, regions, rep_vcf, Path(str(rep_vcf) + ".tbi"), rep_log):
            f.unlink(missing_ok=True)
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
