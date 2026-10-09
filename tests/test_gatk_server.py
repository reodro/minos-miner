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

    def _crashing(self, name):
        d = self.work / "servers" / name
        d.mkdir(parents=True)
        srv = G._GatkServer(_FakeProc(), open(d / "jvm.log", "w"), d, f"minos-run-srv-{name}", 30)
        self.runner._start_server = lambda: srv
        _answer(d, ["error java.lang.IllegalStateException: Graph must have ref source and sink vertices"])
        return d

    def test_a_batch_that_crashes_in_gatk_is_not_rerun_cold(self):
        self._crashing("s3")
        spec = self.work / "R" / "variants_y.tsv"
        spec.parent.mkdir(parents=True, exist_ok=True)
        spec.write_text("/data/output/R/v1_full.vcf.gz\t30.0\n")
        cmd = ["docker", "run", "--name", "minos-run-T-1", "-v", f"{self.work}:/data/output", "img", "java",
               "-Dminos.gt.variants=/data/output/R/variants_y.tsv", "-cp", "x", "org.broadinstitute.hellbender.Main",
               "HaplotypeCaller", "-O", "/data/output/R/o.vcf.gz"]
        p = self.runner._run_on_server(cmd)
        self.assertIsNotNone(p)                                 # not None: no cold re-run of the batch
        self.assertEqual(p.returncode, 3)                       # rc=3: the executor re-runs the head alone
        self.assertIn("ref source and sink", G.gatk_error(p.stderr))
        self.assertEqual(self.runner.server_stats["batch_crashes"], 1)
        self.assertEqual(self.runner.server_stats["fallbacks"], 0)

    def test_a_solo_job_that_crashes_still_runs_cold(self):
        self._crashing("s4")
        cmd = ["docker", "run", "--name", "minos-run-T-1", "-v", f"{self.work}:/data/output", "img", "java",
               "-cp", "x", "org.broadinstitute.hellbender.Main", "HaplotypeCaller", "-O", "/data/output/R/o.vcf.gz"]
        self.assertIsNone(self.runner._run_on_server(cmd))      # the solo run is the stock-equal check: cold re-run
        self.assertEqual(self.runner.server_stats["fallbacks"], 1)

    def test_an_ended_round_runs_nothing_more_not_even_cold(self):
        d = self.work / "servers" / "s2"
        d.mkdir(parents=True)
        srv = G._GatkServer(_FakeProc(), open(d / "jvm.log", "w"), d, "minos-run-srv-s2", 30)
        self.runner._start_server = lambda: srv
        _answer(d, ["dead server exited"])                      # kill_round stopped the server under the job
        cmd = ["docker", "run", "--name", "minos-run-T-1", "-v", f"{self.work}:/data/output", "img", "java",
               "-cp", "x", "org.broadinstitute.hellbender.Main", "HaplotypeCaller", "-O", "/data/output/R/o.vcf.gz"]
        self.runner._ended_tags.add("T")
        p = self.runner._run_on_server(cmd)
        self.assertEqual((p.returncode, p.stderr), (137, "round ended"))   # no cold rerun of the dead job
        self.assertEqual(self.runner.server_stats["fallbacks"], 0)
        started = []
        self.runner._docker_run = lambda *a, **k: started.append(1)
        self.assertEqual(self.runner._docker(cmd).returncode, 137)        # and no new job of that round starts
        self.assertEqual(started, [])
        other = [a.replace("minos-run-T-1", "minos-run-T2-1") for a in cmd]
        self.runner.server = False
        self.runner._docker(other)                                         # another round still runs
        self.assertEqual(started, [1])


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
        self.addCleanup(setattr, G, "ACTIVITY_CACHE", G.ACTIVITY_CACHE)
        G.ACTIVITY_CACHE = True                                                   # off by default (MINOS_ACTIVITY_CACHE=1)
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
        G.ACTIVITY_CACHE = False
        self.assertIsNone(dir_of(a))

    def test_the_assembly_cache_is_the_rounds_and_carries_the_skeleton(self):
        from tests.helpers import ranges
        self.addCleanup(setattr, G, "ASSEMBLY_CACHE", G.ASSEMBLY_CACHE)
        G.ASSEMBLY_CACHE = True                                                   # off by default (MINOS_ASSEMBLY_CACHE=1)
        d = ranges().defaults("gatk").params
        props = lambda c, p=None: self.runner._assembly_cache(c, "full", p)                                        # noqa: E731
        a = Config("gatk", {**d, "min_pruning": 2})
        late = Config("gatk", {**d, "min_pruning": 3, "pcr_indel_model": "NONE", "active_probability_threshold": 0.005})
        early = Config("gatk", {**d, "min_pruning": 2, "min_base_quality_score": d["min_base_quality_score"] + 3})
        self.assertEqual(props(a)["minos.assembly.dir"], "/data/output/assembly")                # one directory per round
        self.assertEqual(props(a)["minos.assembly.skeleton"], props(late)["minos.assembly.skeleton"])
        self.assertNotEqual(props(a)["minos.assembly.skeleton"], props(early)["minos.assembly.skeleton"])
        self.assertEqual(props(a, {"minos.repair.log": "/data/output/r"}), {})
        self.assertEqual(props(Config("gatk", {**d, "emit_ref_confidence": "GVCF"})), {})
        G.ASSEMBLY_CACHE = False
        self.assertEqual(props(a), {})


class RoundCleanupTest(unittest.TestCase):
    """finish_round: a round's caches go when it ends, its VCFs once KEEP_ROUNDS newer rounds ended on this machine."""

    def setUp(self):
        from types import SimpleNamespace
        from unittest import mock
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        for p in (mock.patch.object(G, "SETTINGS", SimpleNamespace(cache_dir=str(self.tmp))),
                  mock.patch.object(G, "KEEP_ROUNDS", 3)):
            p.start()
            self.addCleanup(p.stop)
        self.work = self.tmp / "callsets"
        self.ids = [f"2026-01-0{d}T00:00:00+00:00" for d in range(1, 7)]
        for rid in self.ids + ["2025-12-31T00:00:00+00:00"]:      # the last one: a stored round never mined here
            t = Round.tag_of(rid)
            for d in (self.work / t / "activity" / "sk", self.work / t / "assembly" / "ab", self.tmp / "dcache" / t):
                d.mkdir(parents=True)
                (d / "f").write_text("x")
            (self.work / t / "k_full.vcf.gz").write_text("x")

    def runner(self):
        return G.GatkRunner(work_dir=str(self.work), patch_jar=None, runtime="docker")

    def kept(self):
        return sorted(t for t in [Round.tag_of(r) for r in self.ids] if (self.work / t / "k_full.vcf.gz").exists())

    def test_caches_go_at_once_and_vcfs_after_three_newer_rounds_also_across_restarts(self):
        tags = [Round.tag_of(r) for r in self.ids]
        r = self.runner()
        r.finish_round(self.ids[0])
        self.assertFalse((self.work / tags[0] / "activity").exists())
        self.assertFalse((self.work / tags[0] / "assembly").exists())
        self.assertTrue((self.work / tags[0] / "k_full.vcf.gz").exists())      # kept for the review
        for rid in self.ids[1:3]:
            r.finish_round(rid)
        self.assertEqual(self.kept(), tags)                                     # 3 ended: nothing older to delete
        self.runner().finish_round(self.ids[3])                                 # a restarted miner
        self.assertFalse((self.work / tags[0]).exists())
        self.assertFalse((self.tmp / "dcache" / tags[0]).exists())
        self.assertEqual(self.kept(), tags[1:])
        r.finish_round(self.ids[1])                                             # ended again: the newest now
        r.finish_round(self.ids[4])
        self.assertEqual(self.kept(), [tags[1], tags[3], tags[4], tags[5]])     # tags[2] was the oldest
        self.assertTrue((self.work / "20251231T000000Z" / "k_full.vcf.gz").exists())   # never ended here: untouched
        self.assertTrue((self.work / "20251231T000000Z" / "activity").exists())

    def test_keep_rounds_0_deletes_only_the_caches(self):
        G.KEEP_ROUNDS = 0
        r = self.runner()
        for rid in self.ids:
            r.finish_round(rid)
        self.assertEqual(len(self.kept()), len(self.ids))
        self.assertFalse(any((self.work / Round.tag_of(rid) / "activity").exists() for rid in self.ids))


if __name__ == "__main__":
    unittest.main()
