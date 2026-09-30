"""FleetExecutor over the real node agent (subprocesses in --fake mode, no ssh, no GATK): results come back scored,
a family is one group, a node that dies mid-work has its jobs re-queued, bad nodes are rejected."""
import json
import os
import shutil
import tempfile
import time
import unittest
from pathlib import Path

from common.types import AnswerSet, Config, EvalJob, Regime, Round
from common.settings import SETTINGS
from evaluator.fleet import selftest
from evaluator.fleet.coordinator import FleetExecutor
from evaluator.fleet.inventory import REPO
from tests.helpers import ranges

GT = "standard_min_confidence_threshold_for_calling"


class Planner:
    """Families = same config except the genotyping threshold (like GatkRunner.branch_key for that tier)."""

    def branch_key(self, config, fidelity):
        return Config("gatk", {k: v for k, v in config.params.items() if k != GT}).key + fidelity

    def variant_cost(self, main, variant, fidelity):
        return 0.01 if self.branch_key(main, fidelity) == self.branch_key(variant, fidelity) else None


SELFTEST_SOURCE = Path(SETTINGS.rounds_dir) / selftest.spec()["tag"] / "input.bam"


@unittest.skipUnless(SELFTEST_SOURCE.exists(), "needs main's stored round for the self-test bundle")
class FleetTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        main = self.tmp / "main"
        (main / "rounds" / "r").mkdir(parents=True)
        for f in ("input.bam", "input.bam.bai"):
            (main / "rounds" / "r" / f).write_text("x")
        for f in ("truth.vcf.gz", "mutations.vcf.gz"):
            (main / f).write_text("x")
        self.round = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, str(main / "rounds" / "r" / "input.bam"), "")
        self.aset = AnswerSet(self.round.round_id, 1, Regime.KG, "S", [], [], 600,
                              notes={"truth_vcf": str(main / "truth.vcf.gz"), "mutations_vcf": str(main / "mutations.vcf.gz")})
        d = ranges().defaults("gatk").params
        self.family = [Config("gatk", {**d, GT: c}) for c in (30.0, 40.0, 50.0)]
        self.others = [Config("gatk", {**d, "min_pruning": p}) for p in (3, 4, 5)]
        self.ex = None

    def tearDown(self):
        if self.ex:
            self.ex.close()

    def fleet(self, *nodes, reference=True):
        """nodes: (name, workers, env) -> local hosts running the real agent in --fake mode."""
        inv = []
        for name, workers, env in nodes:
            data = self.tmp / name
            (data / "reference" / "chr20").mkdir(parents=True, exist_ok=True)
            if reference:
                (data / "reference" / "chr20" / "chr20.fa").write_text("x")
            inv.append({"host": "local", "port": 2200 + len(inv), "repo": str(REPO), "data_dir": str(data),
                        "workers": workers, "node_args": "--fake", "env": env})
        path = self.tmp / "fleet.json"
        path.write_text(json.dumps(inv))
        self.ex = FleetExecutor(Planner(), inventory=path, reconnect_s=3600).start()
        self.ex.publish_round(self.round, self.aset)
        return self.ex

    def wait(self, n, timeout=60):
        got, t0 = {}, time.time()
        while len(got) < n and time.time() - t0 < timeout:
            for cs in self.ex.poll(0.5):
                got[cs.job_id] = cs
        return got

    def wait_ready(self, n, timeout=30):
        t0 = time.time()
        while time.time() - t0 < timeout:
            nodes = self.ex.nodes()
            if sum(1 for x in nodes if x["ready"] and x["rounds"]) >= n or any(x["rejected"] for x in nodes):
                return nodes
            time.sleep(0.2)
        return self.ex.nodes()

    def test_results_come_back_scored_and_a_family_is_one_group(self):
        ex = self.fleet(("a", 2, {}), ("b", 2, {}))
        self.wait_ready(2)
        jobs = [EvalJob(f"j{i}", self.round, c) for i, c in enumerate(self.family + self.others)]
        ex.submit(jobs)
        got = self.wait(len(jobs))
        self.assertEqual(sorted(got), sorted(j.job_id for j in jobs))
        self.assertTrue(all(cs.ok and cs.score is not None for cs in got.values()))
        with ex._lock:                                           # grouping is the scheduler's decision
            ex._queue.clear()
        ex.submit([EvalJob(f"f{i}", self.round, Config("gatk", {**c.params, "min_mapping_quality_score": 30}))
                   for i, c in enumerate(self.family)])
        time.sleep(0.5)
        self.assertEqual(len(self.wait(3)), 3)

    def test_dead_node_jobs_are_requeued_to_the_other(self):
        ex = self.fleet(("slow", 3, {"MINOS_FAKE_DELAY": "30"}), ("fast", 1, {"MINOS_FAKE_DELAY": "0.2"}))
        self.wait_ready(2)
        jobs = [EvalJob(f"j{i}", self.round, c) for i, c in enumerate(self.others)]
        ex.submit(jobs)
        time.sleep(1.0)
        slow = next(s for s in ex._sessions.values() if s.host.data_dir.endswith("slow"))
        self.assertGreater(len(slow.jobs), 0)
        slow.proc.kill()                                         # the ssh session dies
        got = self.wait(len(jobs), timeout=20)
        self.assertEqual(len(got), len(jobs))
        self.assertTrue(all(cs.ok for cs in got.values()))

    def test_end_round_frees_every_node_at_once(self):
        ex = self.fleet(("a", 2, {"MINOS_FAKE_DELAY": "30"}), ("b", 2, {"MINOS_FAKE_DELAY": "30"}))
        self.wait_ready(2)
        ex.submit([EvalJob(f"j{i}", self.round, c) for i, c in enumerate(self.others + self.family)])
        time.sleep(1.0)
        self.assertGreater(sum(n["open_groups"] for n in ex.nodes()), 0)
        ex.end_round(self.round.round_id)
        self.assertEqual(sum(n["open_groups"] for n in ex.nodes()), 0)     # workers free for the next round
        self.assertEqual(ex.idle_slots(), 4)
        self.assertEqual(self.wait(1, timeout=3), {})                      # nothing of the ended round arrives

    def test_a_round_no_node_can_prepare_fails_its_jobs(self):
        ex = self.fleet(("noref", 2, {}), reference=False)          # the node lacks the round's reference genome
        t0 = time.time()
        while time.time() - t0 < 20 and not any(x["ready"] for x in ex.nodes()):
            time.sleep(0.2)
        ex.submit([EvalJob(f"j{i}", self.round, c) for i, c in enumerate(self.others)])
        got = self.wait(len(self.others), timeout=20)
        self.assertEqual(len(got), len(self.others))                 # they come back instead of waiting forever
        self.assertTrue(all(not cs.ok for cs in got.values()))

    def test_bad_nodes_are_rejected(self):
        ex = self.fleet(("badtest", 1, {"MINOS_FAKE_SELFTEST": "fail"}), ("oldcode", 1, {"MINOS_FAKE_CODE": "0" * 16}))
        time.sleep(3)
        rejected = {n["host"]: n["rejected"] for n in ex.nodes()}
        self.assertTrue(all(rejected.values()), rejected)
        self.assertEqual(ex.idle_slots(), 0)


class GroupingTest(unittest.TestCase):
    def test_family_is_taken_as_one_group(self):
        ex = FleetExecutor(Planner(), inventory=Path("/nonexistent"))
        r = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")
        d = ranges().defaults("gatk").params
        fam = [Config("gatk", {**d, GT: c}) for c in (30.0, 40.0)] + [Config("gatk", {**d, "min_pruning": 4})]
        ex.submit([EvalJob(f"j{i}", r, c) for i, c in enumerate(fam)])
        with ex._lock:
            g1 = ex._take_groups({r.round_id})
            g2 = ex._take_groups({r.round_id})
            none = ex._take_groups({"other-round"})
        self.assertEqual(sorted(len(g[0]) for g in (g1, g2)), [1, 2])
        self.assertEqual(none, [])


class CodecTest(unittest.TestCase):
    def test_round_trip(self):
        from dataclasses import replace
        from common.types import CallSet, ClassCounts, ScoreResult
        from evaluator.fleet import codec
        r = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "/x/input.bam", "/x/chr20.fa")
        job = EvalJob("j", r, Config("gatk", {"x": 1}), "full", 3, "t")
        self.assertEqual(codec.job_from(codec.job_to(job)), job)
        s = ScoreResult(r.round_id, "k", 1, {"snp_het": ClassCounts(1, 0, 0)}, 0, 600, score=90.0, job_id="j")
        cs = CallSet(r.round_id, "k", "/n/x.vcf.gz", "d", 1.0, job_id="j", score=s)
        self.assertEqual(codec.callset_from(codec.callset_to(cs)), cs)
        self.assertEqual(codec.callset_from(codec.callset_to(replace(cs, score=None))), replace(cs, score=None))


class InventoryTest(unittest.TestCase):
    def test_new_hosts_inherit_user_and_key(self):
        from evaluator.fleet.inventory import load
        p = Path(tempfile.mkdtemp()) / "fleet.json"
        p.write_text(json.dumps([{"host": "1.1.1.1", "user": "ubuntu", "key_filename": "~/.ssh/k"},
                                 {"host": "2.2.2.2", "port": 2222},
                                 {"host": "3.3.3.3", "enabled": False}]))
        hosts = load(p)
        self.assertEqual([(h.name, h.user, h.key_filename) for h in hosts],
                         [("1.1.1.1", "ubuntu", "~/.ssh/k"), ("2.2.2.2:2222", "ubuntu", "~/.ssh/k")])

    def test_a_host_with_its_own_key_does_not_get_another_keys_passphrase(self):
        from evaluator.fleet.inventory import load
        p = Path(tempfile.mkdtemp()) / "fleet.json"
        p.write_text(json.dumps([{"host": "1.1.1.1", "key_filename": "~/.ssh/a", "key_passphrase": "pa"},
                                 {"host": "2.2.2.2"},
                                 {"host": "3.3.3.3", "key_filename": "~/.ssh/b"}]))
        h = {x.host: x for x in load(p)}
        self.assertEqual(h["2.2.2.2"].key_passphrase, "pa")
        self.assertIsNone(h["3.3.3.3"].key_passphrase)
        self.assertEqual(len(load(p, enabled_only=False)), 3)


class PassphraseTest(unittest.TestCase):
    def test_agent_unlocks_an_encrypted_key(self):
        import subprocess
        from evaluator.fleet.inventory import ensure_agent, load
        d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d, True)
        key = d / "k"
        subprocess.run(["ssh-keygen", "-q", "-t", "ed25519", "-N", "pw-123", "-f", str(key)], check=True)
        inv = d / "fleet.json"
        inv.write_text(json.dumps([{"host": "1.1.1.1", "key_filename": str(key), "key_passphrase": "pw-123"}]))
        hosts = load(inv)
        self.assertNotIn("pw-123", repr(hosts))
        ensure_agent(hosts)
        listed = subprocess.run(["ssh-add", "-l"], capture_output=True, text=True).stdout
        fingerprint = subprocess.run(["ssh-keygen", "-lf", str(key) + ".pub"], capture_output=True, text=True).stdout.split()[1]
        self.assertIn(fingerprint, listed)
