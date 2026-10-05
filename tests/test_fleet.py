"""FleetExecutor over the real node agent (subprocesses in --fake mode, no ssh, no GATK): results come back scored,
a family is one group, GATK runs before the answer set exists, a late latency run is duplicated on a faster node,
a node that dies mid-work has its jobs re-queued, bad nodes are rejected."""
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


PCR = "pcr_indel_model"


class Planner:
    """Families = same config except the genotyping threshold (cheap rider, 0.01 of a run) or the PCR model (a
    PairHMM rider, 0.12) - like GatkRunner.branch_key/variant_cost for those tiers."""

    def branch_key(self, config, fidelity):
        return Config("gatk", {k: v for k, v in config.params.items() if k not in (GT, PCR)}).key + fidelity

    def variant_cost(self, main, variant, fidelity):
        if self.branch_key(main, fidelity) != self.branch_key(variant, fidelity):
            return None
        return 0.12 if main.params.get(PCR) != variant.params.get(PCR) else 0.01


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

    def fleet(self, *nodes, reference=True, answer_set=True):
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
        self.ex.publish_round(self.round)
        if answer_set:
            self.ex.publish_answer_set(self.round.round_id, self.aset)
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

    def test_latency_lane_runs_on_the_fastest_node_with_cheap_riders_only(self):
        ex = self.fleet(("slow", 2, {"MINOS_FAKE_DELAY": "1.2"}), ("fast", 2, {"MINOS_FAKE_DELAY": "0.15"}))
        self.wait_ready(2)
        d = ranges().defaults("gatk").params
        warm = [Config("gatk", {**d, "min_pruning": p, GT: 30.0 + i}) for i, p in enumerate((3, 4, 5, 6, 7, 8))]
        ex.submit([EvalJob(f"w{i}", self.round, c) for i, c in enumerate(warm)])
        self.wait(len(warm), timeout=30)                        # both nodes ran solo jobs: their speeds are measured
        speeds = {n["host"].split("/")[-1] if "/" in n["host"] else n["host"]: n["solo_s"] for n in ex.nodes()}
        fast = min(speeds, key=speeds.get)
        self.assertLess(min(speeds.values()) * 2, max(speeds.values()))
        # a latency head with a genotyping family + a PCR variant: the family rides, the PCR variant does not; the
        # whole group lands on the fast node although the slow one is free too
        head = Config("gatk", {**d, "min_pruning": 9})
        kids = [Config("gatk", {**d, "min_pruning": 9, GT: v}) for v in (40.0, 50.0)]
        pcr = Config("gatk", {**d, "min_pruning": 9, PCR: "NONE"})
        ex.submit([EvalJob("L", self.round, head, lane="latency")]
                  + [EvalJob(f"k{i}", self.round, c) for i, c in enumerate(kids)] + [EvalJob("P", self.round, pcr)])
        with ex._lock:                                          # submit dispatched synchronously: inspect at once
            fast_s = next(s for n, s in ex._sessions.items() if n == fast)
            gid = fast_s.job_group.get("L")
            self.assertIsNotNone(gid, "the latency job was not sent to the fastest node")
            with_head = {j for j, g in fast_s.job_group.items() if g == gid}
        self.assertEqual(with_head, {"L", "k0", "k1"})         # cheap riders packed, the PairHMM variant left out
        self.assertEqual(len(self.wait(4, timeout=30)), 4)

    def test_gatk_runs_before_the_answer_set_and_results_wait_for_it(self):
        ex = self.fleet(("a", 2, {}), answer_set=False)
        self.wait_ready(1)                                        # the round (BAM, reference) is enough to run GATK
        jobs = [EvalJob(f"j{i}", self.round, c) for i, c in enumerate(self.others)]
        ex.submit(jobs)
        self.assertEqual(self.wait(1, timeout=3), {})             # computed, but nothing can be scored yet
        node = next(iter(ex._sessions.values()))
        self.assertEqual((len(node.jobs), ex.idle_slots()), (2, 0))   # both slots busy until their results come back
        ex.publish_answer_set(self.round.round_id, self.aset)
        got = self.wait(len(jobs), timeout=20)                    # the held results are scored the moment it lands
        self.assertEqual(sorted(got), sorted(j.job_id for j in jobs))
        self.assertTrue(all(cs.ok and cs.score is not None for cs in got.values()))

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
        from dataclasses import replace
        self.round = replace(self.round, chrom="chr0")              # no such reference anywhere: a node lacking it
        ex = self.fleet(("noref", 2, {}), reference=False)          # cannot fetch it either
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
        from evaluator.fleet.coordinator import _Session
        node = _Session(host=None, workers=4, rounds_ready={r.round_id})
        other = _Session(host=None, workers=4, rounds_ready={"other-round"})
        with ex._lock:
            g1 = ex._take_groups(node)
            g2 = ex._take_groups(node)
            none = ex._take_groups(other)
        self.assertEqual(sorted(len(g[0][0]) for g in (g1, g2)), [1, 2])
        self.assertAlmostEqual(max(g[0][1] for g in (g1, g2)), 1.01)   # the family: one run + one cheap rider
        self.assertEqual(none, [])


class SpeculationTest(unittest.TestCase):
    def test_an_overdue_latency_run_is_duplicated_on_a_free_fast_node_and_the_first_result_wins(self):
        from dataclasses import replace
        from common.types import CallSet
        from evaluator.fleet.coordinator import _Session
        ex = FleetExecutor(Planner(), inventory=Path("/nonexistent"))
        ex._now = lambda: 1000.0
        r = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")
        d = ranges().defaults("gatk").params
        a = _Session(host=None, workers=2, ready=True, rounds_ready={r.round_id}, solo_s=150.0)
        b = _Session(host=None, workers=2, ready=True, rounds_ready={r.round_id}, solo_s=140.0)
        ex._sessions = {"a": a, "b": b}
        ex.submit([EvalJob("L", r, Config("gatk", {**d, "min_pruning": 3}), lane="latency")])
        self.assertIn("L", b.jobs)                                 # the fastest free node
        ex._speculate()
        self.assertEqual(ex.speculated, 0)                         # fresh start: 150 s there is no worthwhile gain
        ex._now = lambda: 1000.0 + 1.6 * 140.0                    # overdue (> 1.5x its expected time)
        ex._speculate()
        self.assertEqual(ex.speculated, 1)
        self.assertIn("L", a.jobs)                                 # the duplicate runs on the other node
        g_b, g_a = next(iter(b.groups.values())), next(iter(a.groups.values()))
        self.assertTrue(g_a.duplicate and g_b.twin is not None)
        cs = CallSet(r.round_id, Config("gatk", {**d, "min_pruning": 3}).key, "/x.vcf.gz", "d", 1.0, job_id="L")
        ex._result(a, cs)                                          # the duplicate delivers first
        self.assertEqual(ex.speculation_won, 1)
        self.assertEqual((a.open_groups(), b.open_groups(), list(b.jobs)), (0, 0, []))   # the loser was told to abort
        self.assertEqual([c.job_id for c in ex.poll(0)], ["L"])


class LatencyLaneTest(unittest.TestCase):
    def test_a_latency_job_waits_for_a_fast_node_while_slow_nodes_take_width(self):
        from evaluator.fleet.coordinator import _Session
        ex = FleetExecutor(Planner(), inventory=Path("/nonexistent"))
        r = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")
        d = ranges().defaults("gatk").params
        fast = _Session(host=None, workers=1, ready=True, rounds_ready={r.round_id}, solo_s=100.0)
        slow = _Session(host=None, workers=4, ready=True, rounds_ready={r.round_id}, solo_s=300.0)
        ex._sessions = {"fast": fast, "slow": slow}
        fast.groups[0] = None                                     # the fast node is busy
        fast.groups.__setitem__(0, type("G", (), {"jobs": {"x": None}})())
        ex.submit([EvalJob("L", r, Config("gatk", {**d, "min_pruning": 3}), lane="latency", priority=3),
                   EvalJob("T", r, Config("gatk", {**d, "min_pruning": 4}))])
        self.assertIn("T", slow.jobs)                              # the throughput job ran on the slow node
        self.assertNotIn("L", slow.jobs)                           # the move did not: it waits for the fast node
        self.assertEqual([it[2].job_id for it in ex._queue], ["L"])
        del fast.groups[0]
        ex._pump()
        self.assertIn("L", fast.jobs)


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


class FleetScheduleTest(unittest.TestCase):
    """Grouping rules, without SSH: a long family is split, and a job that cannot finish stays queued."""

    def test_splits_a_long_family_and_holds_a_job_that_cannot_finish(self):
        rnd = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")
        d = ranges().defaults("gatk").params
        conf = [Config("gatk", {**d, GT: c}) for c in (30.0, 40.0, 50.0, 60.0)]
        ex = FleetExecutor(Planner(), inventory=Path("/tmp/does-not-exist-fleet.json"))
        ex.runner.variant_cost = lambda main, variant, fidelity: 0.6
        ex.submit([EvalJob(f"c{i}", rnd, c) for i, c in enumerate(conf)])
        from evaluator.fleet.coordinator import _Session
        node = _Session(host=None, workers=4, rounds_ready={rnd.round_id})
        first = ex._take_groups(node)
        second = ex._take_groups(node)
        self.assertEqual([len(first[0][0]), len(second[0][0])], [2, 2])
        slow = _Session(host=None, workers=4, rounds_ready={rnd.round_id}, solo_s=200.0)   # a slower node packs the
        more = [Config("gatk", {**d, GT: c}) for c in (70.0, 80.0, 90.0, 95.0)]              # same shapes: its cap and
        ex.submit([EvalJob(f"d{i}", rnd, c) for i, c in enumerate(more)])                 # rider costs both scale
        self.assertEqual([len(g[0][0]) for g in (ex._take_groups(slow), ex._take_groups(slow))], [2, 2])

        held = FleetExecutor(Planner(), inventory=Path("/tmp/does-not-exist-fleet.json"))
        held._now = lambda: 0.0
        held.set_finish_by(10.0)                              # a solo run is 132 s
        held.submit([EvalJob("a", rnd, conf[0])])
        self.assertEqual(held._take_groups(node), [])
        self.assertEqual(len(held._queue), 1)


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
