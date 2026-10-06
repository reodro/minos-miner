import shutil
import subprocess
import tempfile
import threading
import time
import unittest
from pathlib import Path

import evaluator.gatk as G
from common.types import Config, Round


class _FakeProc:
    def __init__(self):
        self.alive = True
        self.args = []

    def poll(self):
        return None if self.alive else 0

    def wait(self, timeout=None):
        self.alive = False
        return 0

    def kill(self):
        self.alive = False


def _answer(directory: Path, answers: list) -> threading.Thread:
    """A stand-in MinosServer: answers job-00001, job-00002, ... with `answers`, recording each job's lines."""
    def serve():
        for n, a in enumerate(answers, 1):
            job = directory / f"job-{n:05d}.tsv"
            while not job.exists():
                time.sleep(0.01)
            (directory / f"job-{n:05d}.seen").write_text(job.read_text())
            (directory / f"job-{n:05d}.done").write_text(a + "\n")
    t = threading.Thread(target=serve, daemon=True)
    t.start()
    return t


class GatkServerTest(unittest.TestCase):
    """The runner's side of the persistent GATK JVMs, without GATK: paths, protocol, fallbacks."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        (self.tmp / "patch").mkdir()
        (self.tmp / "patch" / "minos-patch.jar").write_text("")
        self.work = self.tmp / "work"
        (self.work / "R").mkdir(parents=True)
        bam, ref = self.work / "in" / "input.bam", self.work / "ref" / "chr1.fa"
        for p in (bam, ref):
            p.parent.mkdir(parents=True)
            p.write_text("")
        self.round = Round("2026-01-01T00:00:00+00:00", "chr1", 1, 1000, str(bam), str(ref))
        self.runner = G.GatkRunner(work_dir=str(self.work), patch_jar=str(self.tmp / "patch" / "minos-patch.jar"),
                                   runtime="docker")
        self.addCleanup(self.runner.stop_servers)

    def cmd(self, props=None):
        out = self.work / "R" / "abc_full.vcf.gz"
        return G.GatkRunner.command(self.runner, self.round, Config("gatk", {}), "full", out, jvm_props=props)

    def test_paths_move_into_the_server_mounts(self):
        spec = self.work / "R" / "variants_x.tsv"
        spec.write_text("/data/output/v1_full.vcf.gz\t30.0\n")
        try:
            cmd = self.cmd({"minos.gt.variants": "/data/output/variants_x.tsv", "minos.kbest.out": "/data/output/k"})
        except Exception as e:  # noqa: BLE001 - upstream templates missing on this machine
            self.skipTest(f"upstream command builder unavailable: {e}")
        props, args, out = self.runner._server_job(cmd)
        self.assertEqual(out, self.work / "R" / "abc_full.vcf.gz")
        self.assertIn("/srv/work/in/input.bam", args)
        self.assertIn("/srv/work/ref/chr1.fa", args)
        self.assertIn("/srv/work/R/abc_full.vcf.gz", args)
        self.assertEqual(props["minos.kbest.out"], "/srv/work/R/k")
        self.assertEqual(props["minos.gt.variants"], "/srv/work/R/variants_x.srv.tsv")
        self.assertEqual((self.work / "R" / "variants_x.srv.tsv").read_text(), "/srv/work/R/v1_full.vcf.gz\t30.0\n")
        self.assertFalse(any(a.startswith("/data/") for a in args))

    def test_a_path_outside_the_mounts_stays_cold(self):
        cmd = ["docker", "run", "--name", "minos-run-T-1", "-v", "/nowhere/x:/data/bams:ro", "img", "java",
               "-cp", "x", "org.broadinstitute.hellbender.Main", "HaplotypeCaller", "-I", "/data/bams/in.bam",
               "-O", "/data/bams/o.vcf.gz"]
        self.assertIsNone(self.runner._server_job(cmd))

    def test_jobs_go_to_one_server_and_a_failing_server_is_replaced_and_the_job_runs_cold(self):
        d = self.work / "servers" / "s1"
        d.mkdir(parents=True)
        srv = G._GatkServer(_FakeProc(), open(d / "jvm.log", "w"), d, "minos-run-srv-s1", 30)
        started = []

        def start():
            started.append(1)
            return srv if len(started) == 1 else None
        self.runner._start_server = start
        _answer(d, ["ok 1.0", "error boom"])
        cmd = ["docker", "run", "--name", "minos-run-T-1", "-v", f"{self.work}:/data/output", "img", "java",
               "-Dminos.kbest.out=/data/output/k", "-cp", "x", "org.broadinstitute.hellbender.Main", "HaplotypeCaller",
               "-O", "/data/output/R/o.vcf.gz"]
        first = self.runner._run_on_server(cmd)
        self.assertIsInstance(first, subprocess.CompletedProcess)
        self.assertEqual(first.returncode, 0)
        seen = (d / "job-00001.seen").read_text().splitlines()
        self.assertEqual(seen[0], "minos.kbest.out=/srv/work/k")
        self.assertEqual(seen[1:], ["HaplotypeCaller", "-O", "/srv/work/R/o.vcf.gz"])
        self.assertIsNone(self.runner._run_on_server(cmd))      # the server failed: run it cold
        self.assertEqual(len(started), 1)                       # the same server took the second job
        self.assertEqual(self.runner.server_stats["fallbacks"], 1)
        self.assertEqual(self.runner._servers_idle, [])         # and it was dropped
        self.assertTrue((d / "stop").exists())


class ActivityCacheTest(unittest.TestCase):
    """GatkRunner names one activity-cache directory per skeleton (what a batch shares) for every patched run."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        (self.tmp / "minos-patch.jar").write_text("")
        self.runner = G.GatkRunner(work_dir=str(self.tmp), patch_jar=str(self.tmp / "minos-patch.jar"), runtime="docker")
        self.addCleanup(self.runner.stop_servers)

    def test_one_directory_per_skeleton_and_never_for_repairs_or_unbatchable_configs(self):
        from tests.helpers import ranges
        d = ranges().defaults("gatk").params
        dir_of = lambda c, props=None: self.runner._activity_cache(c, "full", props).get("minos.activity.dir")  # noqa: E731
        a = Config("gatk", {**d, "min_pruning": 2})
        late = Config("gatk", {**d, "min_pruning": 3, "standard_min_confidence_threshold_for_calling": 50.0,
                               "pcr_indel_model": "NONE", "active_probability_threshold": 0.005})
        early = Config("gatk", {**d, "min_pruning": 2, "min_mapping_quality_score": d["min_mapping_quality_score"] + 5})
        self.assertTrue(dir_of(a).startswith("/data/output/activity/"))
        self.assertEqual(dir_of(a), dir_of(late))                  # late stages only: the same skeleton
        self.assertNotEqual(dir_of(a), dir_of(early))              # a read setting: another skeleton
        self.assertIsNone(dir_of(a, {"minos.repair.log": "/data/output/r"}))
        self.assertIsNone(dir_of(Config("gatk", {**d, "emit_ref_confidence": "GVCF"})))
        self.addCleanup(setattr, G, "ACTIVITY_CACHE", G.ACTIVITY_CACHE)
        G.ACTIVITY_CACHE = False
        self.assertIsNone(dir_of(a))


if __name__ == "__main__":
    unittest.main()
