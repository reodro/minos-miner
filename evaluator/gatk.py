"""The real runner and scorer. Owner: evaluator.

Verified (2026-09-28; 3 rounds x 3 configs against the validator-identical baseline):
  used     - SAFE_ANNOTATION_EXCLUSIONS (-AX): bit-identical calls, x1.135.
           - batching (gatk_patch): at full load a variant costs the most expensive stage it changes
             (VARIANT_COST: genotyping 1.4 s, PairHMM 16 s, assembly 72 s, activity 78 s, of a 132 s solo).
             54/54 + 62/62 bit-identical to stock runs; variants whose own run would draw random numbers
             differently come back DIVERGED and are repaired.
           - assembly params share traversal and region finding only, 86/86 bit-identical.
             Excluding QualByDepth or DepthPerAlleleBySample changes QUAL; RMSMappingQuality is required.
           - 1 thread per job: x2.06 throughput vs 4 threads per job on a saturated box; calls identical.
           - warm JVM: the process stays up after a traversal. A later genotyping-only variant of that
             family is genotyped from the stored likelihoods. Anything else, or a JVM that ran out of
             memory, is a cold run. Bit-identical to cold runs on chr18:57771001-57871001 (229 calls).
  rejected - GVCF re-genotyping (heterozygosity changes calls; QUAL shifts flip stand_call_conf thresholds in ~13%
             of cases), QUAL projection (1e-4..1e-3 off), native java outside the image (different calls and slower:
             GKL PairHMM missing), interval restriction (~2%), scatter-gather (x0.49 throughput).
Workers must run the validator's image files. MINOS_RUNTIME=sandbox (the default) runs them
without the Docker daemon; MINOS_RUNTIME=docker is the container path.
"""
from __future__ import annotations

import gzip
import hashlib
import logging
import os
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from common.scoring import finalize
from common.settings import SETTINGS
from common.types import CLASSES, FULL, AnswerSet, CallSet, ClassCounts, Config, Round, ScoreResult

log = logging.getLogger(__name__)

if SETTINGS.minos_subnet_dir not in sys.path:
    sys.path.append(SETTINGS.minos_subnet_dir)

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
                     # GATK's unset value is MathUtils.log10ToLog(1), not the 6-digit truncation. A shorter
                     # value makes the warm JVM treat a genotyping-only request as an assembly change and refuse it.
                     "pruning_lod_threshold": 2.302585092994046}
# Read only when the traversal cuts the activity profile into regions (assembly_region_padding is not: it also sets
# how far each read shard reaches): variants differing in these get their own region stream from the shared per-locus
# activity states and reads, then assembly and everything after it.
ACTIVITY_BRANCH_PARAMS = ("active_probability_threshold", "min_assembly_region_size", "max_assembly_region_size")
ACTIVITY_DEFAULTS = {"active_probability_threshold": 0.002, "min_assembly_region_size": 50, "max_assembly_region_size": 300}
BRANCH_PARAMS = GENOTYPING_BRANCH_PARAMS + PAIRHMM_BRANCH_PARAMS + ASSEMBLY_BRANCH_PARAMS + ACTIVITY_BRANCH_PARAMS
# Extra cost of one batch variant, as a fraction of a solo run, by the most expensive stage it changes
# (measured at production load). A variant that differs in several tiers pays the highest.
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
MAX_WARM = 2                     # JVMs kept after a run so a later genotyping variant can reuse it

SAFE_ANNOTATION_EXCLUSIONS = (
    "Coverage", "DepthPerSampleHC", "ChromosomeCounts", "BaseQualityRankSumTest", "ExcessHet", "FisherStrand",
    "InbreedingCoeff", "MappingQualityRankSumTest", "ReadPosRankSumTest", "StrandOddsRatio", "StrandBiasBySample",
)


LOG4J_NOISE = ("StatusConsoleListener", "StatusLogger")   # log4j's own complaints (no hostname in the sandbox ...)


def _without_log4j_noise(stderr: str) -> List[str]:
    """stderr lines minus log4j's own status reports and their stack traces (no hostname in the sandbox, ...)."""
    out, skipping = [], False
    for line in (stderr or "").splitlines():
        if any(n in line for n in LOG4J_NOISE):
            skipping = True
            continue
        if skipping and (line.startswith(("\t", " ", "Caused by:", "...")) or not line.strip()):
            continue
        skipping = False
        out.append(line)
    return out


def gatk_error(stderr: str) -> str:
    """GATK's actual error, not the command line it echoes: the USER ERROR block, else the exception line (the
    first "Caused by" is more telling than the wrapper), never log4j's own status noise."""
    lines = _without_log4j_noise(stderr)
    for i, line in enumerate(lines):
        if "A USER ERROR has occurred" in line:
            return " ".join(l.strip() for l in lines[i:i + 3])
    caused = [l.strip() for l in lines if l.strip().startswith("Caused by:")]
    if caused:
        return caused[0]
    exc = [l.strip() for l in lines if "Exception" in l or "Error" in l.split(":")[0] or l.startswith("ERROR")]
    return exc[0] if exc else " ".join(lines[-3:])[-600:]


def keep_stderr(out: Path, p: subprocess.CompletedProcess) -> None:
    """A failed run's stderr (its tail), next to where its VCF would have been: the reason stays diagnosable after
    the round, the result line carries only the extracted error."""
    try:
        Path(str(out) + ".stderr").write_text(f"rc={p.returncode}\n" + (p.stderr or "")[-20000:])
    except OSError:
        pass


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
    for earlier, later in zip(spans, spans[1:]):
        # Closed intervals. A shared base belongs to two regions, and the bisect below would keep the later one.
        if earlier[0] == later[0] and later[1] <= earlier[2]:
            raise ValueError(f"overlapping assembly regions {earlier[3]} and {later[3]}")
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


class _WarmSession:
    """One GATK process that has finished its traversal and is waiting for genotyping variants."""

    def __init__(self, proc, log_file, directory: Path, main: Config, timeout_s: int):
        self.proc, self.log_file, self.directory, self.main = proc, log_file, directory, main
        self.timeout_s = timeout_s
        self.n = 0
        self.lock = threading.Lock()

    def request(self, line: str) -> str:
        """Send one variant line. Returns exact, diverged, refused, or an error token."""
        with self.lock:
            if self.proc.poll() is not None:
                return "dead"
            self.n += 1
            name = f"req-{self.n:05d}"
            tmp = self.directory / (name + ".tmp")
            tmp.write_text(line if line.endswith("\n") else line + "\n")
            os.replace(tmp, self.directory / (name + ".tsv"))
            done = self.directory / (name + ".done")
            deadline = time.time() + self.timeout_s
            while time.time() < deadline:
                if done.exists():
                    return done.read_text().strip() or "error"
                if self.proc.poll() is not None:
                    return "dead"
                time.sleep(0.02)
            return "timeout"

    def stop(self) -> None:
        (self.directory / "stop").write_text("1\n")
        try:
            self.proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            killer = getattr(self.proc, "minos_kill", None)
            if killer:
                killer()
            else:
                args = self.proc.args if isinstance(self.proc.args, (list, tuple)) else []
                if "--name" in args:
                    subprocess.run(["docker", "kill", args[args.index("--name") + 1]], capture_output=True)
                self.proc.kill()
        self.log_file.close()


class GatkRunner:
    """The validator's HaplotypeCaller command (minos_subnet/templates/gatk.py) plus the verified speed-ups."""

    def __init__(self, work_dir: Optional[str] = None, threads: int = 1, memory_gb: int = 4, timeout_s: int = 1800,
                 patch_jar: Optional[str] = None, container_prefix: str = "minos-run", runtime: Optional[str] = None):
        """patch_jar: our patched GATK classes (evaluator/gatk_patch), put before the GATK jar on the classpath,
        still inside the validator's image. None = the stock `gatk` wrapper."""
        self.work_dir = Path(work_dir or os.path.join(SETTINGS.cache_dir, "callsets"))
        self.threads, self.memory_gb, self.timeout_s = threads, memory_gb, timeout_s
        self.patch_jar = Path(patch_jar).resolve() if patch_jar else None
        self.container_prefix = container_prefix     # containers are named <prefix>-<round tag>-<id>, so the owner can
                                                     # stop all of them, or exactly one round's (kill_round)
        self.runtime = runtime or os.environ.get("MINOS_RUNTIME", "sandbox")
        self._warm = {}                              # (round_id, branch_key) -> _WarmSession
        self._warm_lock = threading.Lock()
        self._warm_reserved = 0
        self._live: Dict[str, tuple] = {}            # container name -> (round_id, config keys it computes)
        self._live_lock = threading.Lock()

    def command(self, round: Round, config: Config, fidelity: str, out_vcf: Path, jvm_props: Optional[dict] = None) -> list:
        from templates.tool_params import validate_and_build_flags   # upstream
        flags = validate_and_build_flags(config.tool, config.params)
        if not flags["valid"]:
            raise ValueError("; ".join(flags["errors"]))
        bam, ref = Path(round.bam).resolve(), Path(round.reference).resolve()
        heap = f"-Xmx{max(512, int(self.memory_gb * 1024 * 0.8))}m"
        cmd = ["docker", "run", "--rm", "--name", f"{self.container_prefix}-{round.tag}-{uuid.uuid4().hex[:10]}",
               f"--cpus={self.threads}", f"--memory={self.memory_gb}g",
               "--network", "none",     # GATK needs no network; without one docker creates no interface per container
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

    def kill_round(self, round_id: str) -> int:
        """Kill this runner's running containers of one round (the round is over: their results are useless)."""
        with self._warm_lock:
            for key in [k for k in self._warm if k[0] == round_id]:
                self._warm.pop(key).stop()
        tag = Round.tag_of(round_id)
        if self.runtime != "docker":
            from evaluator.sandbox import kill_prefix
            return kill_prefix(f"{self.container_prefix}-{tag}-")
        ids = subprocess.run(["docker", "ps", "-q", "--filter", f"name={self.container_prefix}-{tag}-"],
                             capture_output=True, text=True).stdout.split()
        if ids:
            subprocess.run(["docker", "kill", *ids], capture_output=True)
        return len(ids)

    def kill_configs(self, round_id: str, keys: Iterable[str]) -> int:
        """Kill the running containers computing any of these configs for the round (their results are not wanted:
        a duplicate run elsewhere delivered first). Returns how many were killed."""
        want = set(keys)
        with self._live_lock:
            names = [n for n, (rid, ks) in self._live.items() if rid == round_id and ks & want]
        for name in names:
            self._kill(name)
        return len(names)

    def _docker(self, cmd: list, warm_dir: Optional[Path] = None, warm_main: Optional[tuple] = None,
                owner: Optional[tuple] = None, reserved: bool = True) -> subprocess.CompletedProcess:
        """Run a GATK container; on timeout also kill the container (killing the docker client leaves it running).
        `owner` = (round_id, config keys): while it runs, kill_configs can find it.

        warm_dir set: the process stays up after the traversal when it writes ready=ok, so a later genotyping
        variant can be submitted to it (warm_main set: it is kept as that family's session; None: it is told to
        stop at once - a warm-style run). A ready of "no", or a process that exits, is an ordinary run.
        `reserved`: warm_dir holds one of the MAX_WARM slots, given back here.
        """
        name = cmd[cmd.index("--name") + 1]
        if owner is not None:
            with self._live_lock:
                self._live[name] = owner
        try:
            return self._docker_run(cmd, name, warm_dir, warm_main, reserved)
        finally:
            with self._live_lock:
                self._live.pop(name, None)

    def _docker_run(self, cmd: list, name: str, warm_dir: Optional[Path], warm_main: Optional[tuple],
                    reserved: bool = True) -> subprocess.CompletedProcess:
        if warm_dir is None:
            try:
                if self.runtime == "docker":
                    return subprocess.run(cmd, capture_output=True, text=True, timeout=self.timeout_s)
                from evaluator.sandbox import run as sandbox_run
                return sandbox_run(cmd, timeout=self.timeout_s)
            except subprocess.TimeoutExpired:
                self._kill(name)
                return subprocess.CompletedProcess(cmd, 124, "", f"timed out after {self.timeout_s}s")
        log_file = open(warm_dir / "jvm.log", "w")
        if self.runtime == "docker":
            proc = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        else:
            from evaluator.sandbox import start as sandbox_start
            proc = sandbox_start(cmd, stdout=log_file)
        state = self._wait_ready(proc, warm_dir)
        if state == "ok" and warm_main is not None:
            round, config, fidelity = warm_main
            session = _WarmSession(proc, log_file, warm_dir, config, self.timeout_s)
            with self._warm_lock:
                self._warm_reserved = max(0, self._warm_reserved - 1)
                old = self._warm.get((round.round_id, self.branch_key(config, fidelity)))
                self._warm[(round.round_id, self.branch_key(config, fidelity))] = session
            if old is not None:
                old.stop()
            return subprocess.CompletedProcess(cmd, 0, "", "")
        if reserved:
            with self._warm_lock:
                self._warm_reserved = max(0, self._warm_reserved - 1)
        if proc.poll() is None:
            (warm_dir / "stop").write_text("1\n")
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self._kill(name, proc)
        log_file.close()
        rc = proc.returncode if proc.returncode is not None else 1
        err = ""
        try:
            err = (warm_dir / "jvm.log").read_text()[-20000:]   # enough to hold the exception line, not just its stack
        except OSError:
            pass
        if self.runtime != "docker":
            from evaluator.sandbox import release
            release(proc)
        return subprocess.CompletedProcess(cmd, rc, "", err)

    def _wait_ready(self, proc: subprocess.Popen, warm_dir: Path) -> str:
        ready = warm_dir / "ready"
        deadline = time.time() + self.timeout_s
        while time.time() < deadline:
            if ready.exists():
                return ready.read_text().strip().split()[0] if ready.read_text().strip() else "no"
            if proc.poll() is not None:
                return "exited"
            time.sleep(0.05)
        self._kill(None, proc)
        return "timeout"

    def _kill(self, name: Optional[str], proc=None) -> None:
        killer = getattr(proc, "minos_kill", None)
        if killer:
            killer()
            return
        if name:
            if self.runtime == "docker":
                subprocess.run(["docker", "kill", name], capture_output=True)
            else:
                from evaluator.sandbox import kill_name
                kill_name(name)

    def _release_warm_reservation(self, warm_dir: Optional[Path]) -> None:
        """command() failed before the container started: give the reserved slot back."""
        if warm_dir is None:
            return
        with self._warm_lock:
            self._warm_reserved = max(0, self._warm_reserved - 1)

    def _reserve_warm(self, parent: Path, config: Config, fidelity: str) -> Optional[Path]:
        """A directory for a warm session, or None when this run should exit as usual."""
        if self.patch_jar is None or self.branch_key(config, fidelity) is None:
            return None
        with self._warm_lock:
            alive = sum(1 for s in self._warm.values() if s.proc.poll() is None)
            if alive + self._warm_reserved >= MAX_WARM:
                return None
            self._warm_reserved += 1
        directory = parent / f"warm_{uuid.uuid4().hex[:10]}"
        directory.mkdir(parents=True)
        return directory

    def _out(self, round: Round, config: Config, fidelity: str) -> Path:
        fid = "full" if fidelity == FULL else hashlib.sha1(fidelity.encode()).hexdigest()[:8]
        out = self.work_dir / round.tag / f"{config.key}_{fid}.vcf.gz"
        out.parent.mkdir(parents=True, exist_ok=True)
        return out

    def run(self, round: Round, config: Config, fidelity: str = FULL) -> CallSet:
        if config.tool != "gatk":
            raise ValueError(f"GatkRunner cannot run {config.tool!r}")
        out = self._out(round, config, fidelity)
        warm_dir = self._reserve_warm(out.parent, config, fidelity)
        props = {"minos.warm.dir": f"/data/output/{warm_dir.name}"} if warm_dir else None
        t0 = time.time()
        try:
            p = self._docker(self.command(round, config, fidelity, out, jvm_props=props),
                             warm_dir=warm_dir, warm_main=(round, config, fidelity) if warm_dir else None,
                             owner=(round.round_id, frozenset([config.key])))
        except Exception:
            self._release_warm_reservation(warm_dir)
            raise
        dt = time.time() - t0
        if p.returncode != 0 or not out.exists():
            keep_stderr(out, p)
            return CallSet(round.round_id, config.key, "", "", dt, ok=False,
                           error=f"rc={p.returncode} {gatk_error(p.stderr)}", fidelity=fidelity)
        return CallSet(round.round_id, config.key, str(out), records_digest(out), dt, fidelity=fidelity)

    # ---- batching: configs differing only in branchable parameters (BRANCH_PARAMS) share one run (needs patch_jar) ----
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

    def run_batch(self, round: Round, configs: List[Config], fidelity: str = FULL, warm: bool = False) -> List[CallSet]:
        """One patched run: configs[0] is the main output, the rest are variants of it (BRANCH_PARAMS only).

        A variant whose own run would have downsampled reads differently (random-stream divergence, see
        gatk_patch) comes back as CallSet(ok=False, error=DIVERGED, vcf=<kept batch output>): its output is not
        exact everywhere; repair() fixes it.

        `warm`: run warm-style even when no warm slot is free (the JVM is told to stop right after the traversal).
        The plain batch path dies in the main's assembly ("Graph must have ref source and sink vertices") for some
        heads - 6 of 6 runs with assembly_region_padding 0-1 - and the same batch passes 6 of 6 times warm-style;
        the executor retries a crashed batch this way.
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
        warm_dir = self._reserve_warm(main.parent, configs[0], fidelity)
        kept = warm_dir is not None
        if warm_dir is None and warm:
            warm_dir = main.parent / f"warm_{uuid.uuid4().hex[:10]}"
            warm_dir.mkdir(parents=True)
        props = {"minos.gt.variants": f"/data/output/{spec.name}"}
        if warm_dir is not None:
            props["minos.warm.dir"] = f"/data/output/{warm_dir.name}"
        t0 = time.time()
        try:
            try:
                p = self._docker(self.command(round, configs[0], fidelity, main, jvm_props=props),
                                 warm_dir=warm_dir, warm_main=(round, configs[0], fidelity) if kept else None,
                                 owner=(round.round_id, frozenset(c.key for c in configs)), reserved=kept)
            except Exception:
                self._release_warm_reservation(warm_dir if kept else None)
                raise
        finally:
            spec.unlink(missing_ok=True)
        dt = time.time() - t0
        if p.returncode != 0 or not main.exists():
            keep_stderr(main, p)
            err = f"rc={p.returncode} {gatk_error(p.stderr)}"
            return [CallSet(round.round_id, c.key, "", "", dt, ok=False, error=err, fidelity=fidelity) for c in configs]
        out = [CallSet(round.round_id, configs[0].key, str(main), records_digest(main), dt, fidelity=fidelity)]
        for c, o in zip(configs[1:], outs):
            out.append(self._adopt_variant(round.round_id, c, o, dt, fidelity))
        return out

    def _genotyping_only(self, main: Config, variant: Config, fidelity: str) -> bool:
        """True when variant can be genotyped from main's stored likelihoods."""
        cost = self.variant_cost(main, variant, fidelity)
        return cost is not None and cost <= VARIANT_COST["genotyping"]

    def _adopt_variant(self, round_id: str, config: Config, out: Path, dt: float, fidelity: str) -> CallSet:
        """Turn a variant's sidecar status into a CallSet. Diverged output is kept for repair()."""
        status = Path(str(out) + ".minos")
        state = status.read_text().strip() if status.exists() else "missing"
        status.unlink(missing_ok=True)
        if state == "exact" and out.exists():
            return CallSet(round_id, config.key, str(out), records_digest(out), dt, fidelity=fidelity)
        if state == "diverged" and out.exists() and Path(str(out) + ".minos.regions").exists():
            kept = out.with_name(out.name.replace(".vcf.gz", ".batch.vcf.gz"))
            out.rename(kept)
            Path(str(out) + ".minos.regions").rename(str(kept) + ".minos.regions")
            Path(str(out) + ".tbi").unlink(missing_ok=True)
            return CallSet(round_id, config.key, str(kept), "", dt, ok=False, error=DIVERGED, fidelity=fidelity)
        out.unlink(missing_ok=True)
        return CallSet(round_id, config.key, "", "", dt, ok=False,
                       error=DIVERGED if state == "diverged" else f"variant status {state}", fidelity=fidelity)

    def warm_ready(self, round: Round, config: Config, fidelity: str = FULL) -> bool:
        """True when a live JVM has this family's likelihoods and config only changes genotyping."""
        key = self.branch_key(config, fidelity)
        if key is None:
            return False
        with self._warm_lock:
            session = self._warm.get((round.round_id, key))
        if session is None or session.proc.poll() is not None:
            return False
        return self._genotyping_only(session.main, config, fidelity)

    def warm_submit(self, round: Round, config: Config, fidelity: str = FULL) -> Optional[CallSet]:
        """Genotype config inside the JVM that already walked this family. None = run it from scratch."""
        key = self.branch_key(config, fidelity)
        if key is None:
            return None
        with self._warm_lock:
            session = self._warm.get((round.round_id, key))
        if session is None or session.proc.poll() is not None:
            with self._warm_lock:
                self._warm.pop((round.round_id, key), None)
            return None
        if not self._genotyping_only(session.main, config, fidelity):
            return None
        out = self._out(round, config, fidelity)
        t0 = time.time()
        state = session.request(variant_line(config, f"/data/output/{out.name}"))
        dt = time.time() - t0
        if state not in ("exact", "diverged"):
            log.info("warm JVM refused %s (%s); running it on its own", config.key, state)
            return None
        cs = self._adopt_variant(round.round_id, config, out, dt, fidelity)
        if cs.error and cs.error != DIVERGED:
            return None
        return cs

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
        p = self._docker(self.command(round, config, fidelity, rep_vcf,
                                      jvm_props={"minos.repair.log": f"/data/output/{regions.name}",
                                                 "minos.repair.out": f"/data/output/{rep_log.name}"}),
                         owner=(round.round_id, frozenset([config.key])))
        dt = time.time() - t0
        spliced = False
        if p.returncode == 0 and rep_vcf.exists() and rep_log.exists():
            try:
                splice_repair(batch, regions, rep_vcf, rep_log, out)
                spliced = True
            except Exception as e:  # noqa: BLE001 - e.g. a record outside every logged region: run it plainly
                log.warning("splicing the repair of %s failed (%s); running it on its own", config.key, e)
        for f in (batch, regions, rep_vcf, Path(str(rep_vcf) + ".tbi"), rep_log):
            f.unlink(missing_ok=True)
        if not spliced:                           # still exact: a plain solo run
            return self.run(round, config, fidelity)
        return CallSet(round.round_id, config.key, str(out), records_digest(out), dt, fidelity=fidelity)


class HappyScorer:
    """Exact: upstream HappyScorer.score_vcf + parse_happy_vcf, counts through common.scoring.finalize.

    One upstream hap.py run per callset, exactly what a validator runs. Truth comes from the answer set:
    `notes['truth_vcf']` / `notes['mutations_vcf']` (writing them from reconstructed records is not implemented).
    hap.py writes its outputs next to the callset's VCF.
    """

    def __init__(self, work_dir: Optional[str] = None):
        # Accepted so a caller can pass a score directory. hap.py writes beside the query VCF
        # (upstream score_vcf); this path is not an output directory and is not read.
        self.work_dir = Path(work_dir or os.path.join(SETTINGS.cache_dir, "happy"))

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        from evaluator import truth_slices
        from evaluator.sandbox import install_scoring_shim
        from utils.scoring import HappyScorer as UpstreamHappy, difficulty_class_counts, parse_happy_vcf
        install_scoring_shim()                  # hap.py and bcftools: the image files, no container, when sandbox is on
        truth_slices.install()                 # the truth is sliced once per round, not on every hap.py run
        truth = answer_set.notes.get("truth_vcf")
        mutations = answer_set.notes.get("mutations_vcf")
        region = answer_set.notes.get("region")
        reference = answer_set.notes.get("reference")
        if not (truth and mutations and region and reference):
            raise NotImplementedError("HappyScorer needs truth/mutations VCFs, region and reference in answer_set.notes "
                                      "(writing them from reconstructed records is not implemented yet)")
        n_targets = answer_set.n_targets
        if n_targets is None or n_targets < 1:
            raise ValueError("HappyScorer needs a positive answer_set.n_targets "
                             "(the fp_per_target denominator); refusing to invent one")
        chrom = region.split(":")[0]
        sdf = str(Path(reference).parent / f"{chrom}.sdf")
        m = UpstreamHappy().score_vcf(truth_vcf=truth, query_vcf=callset.vcf, reference_fasta=reference,
                                      confident_bed=None, region=region, reference_sdf=sdf, mutations_vcf=mutations)
        if not m or "happy_vcf_path" not in m:
            raise RuntimeError("hap.py produced no result")
        recs = parse_happy_vcf(m["happy_vcf_path"], truth_vcf_path=truth)
        cc = difficulty_class_counts(recs)
        ratios = {k: m[k] for k in ("titv_truth_snp", "titv_query_snp", "hethom_truth_snp", "hethom_query_snp") if k in m}
        result = ScoreResult(callset.round_id, callset.config_key, answer_set.version,
                             {k: ClassCounts(**cc[k]) for k in CLASSES}, int(m.get("region_fp_total", 0)), n_targets,
                             ratios, fidelity=callset.fidelity, job_id=callset.job_id)
        return finalize(result)
