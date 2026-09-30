"""FleetExecutor scheduling over its HTTP API (in-process, no GATK): auth, families as one group, node scores
reaching main, lease expiry re-queueing, late results."""
import time
import unittest
from dataclasses import replace

from common.types import FULL, CallSet, ClassCounts, Config, EvalJob, Round, ScoreResult
from evaluator.fleet import codec
from evaluator.fleet.coordinator import FleetExecutor
from tests.helpers import ranges

ROUND = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "/x/input.bam", "/x/chr20.fa")
GT = "standard_min_confidence_threshold_for_calling"


class Planner:
    """Families = same config except the genotyping threshold (like GatkRunner.branch_key for that tier)."""

    def branch_key(self, config, fidelity):
        return Config("gatk", {k: v for k, v in config.params.items() if k != GT}).key + fidelity

    def variant_cost(self, main, variant, fidelity):
        return 0.01 if self.branch_key(main, fidelity) == self.branch_key(variant, fidelity) else None


def scored(job, score=90.0):
    r = ScoreResult(job.round.round_id, job.config.key, 1, {"snp_het": ClassCounts(1, 0, 0)}, 0, 600, score=score,
                    job_id=job.job_id)
    return CallSet(job.round.round_id, job.config.key, "/node/x.vcf.gz", "d", 1.0, job_id=job.job_id, score=r)


class FleetTest(unittest.TestCase):
    def setUp(self):
        from fastapi.testclient import TestClient
        self.ex = FleetExecutor(Planner(), token="t0k", lease_s=0.5)
        self.api = TestClient(self.ex.app(), headers={"Authorization": "Bearer t0k"})
        d = ranges().defaults("gatk").params
        self.fam = [Config("gatk", {**d, GT: c}) for c in (30.0, 40.0, 50.0)]
        self.other = Config("gatk", {**d, "min_pruning": 4})
        self.node = self.api.post("/v1/register", json={"name": "n1", "workers": 2}).json()["node_id"]

    def jobs(self):
        return [EvalJob(f"j{i}", ROUND, c) for i, c in enumerate(self.fam + [self.other])]

    def test_token_required(self):
        from fastapi.testclient import TestClient
        bad = TestClient(self.ex.app(), headers={"Authorization": "Bearer nope"})
        self.assertEqual(bad.post("/v1/register", json={"name": "x", "workers": 1}).status_code, 401)

    def test_family_is_one_group_and_scores_come_back(self):
        self.ex.submit(self.jobs())
        self.assertEqual(self.ex.cost_hints([EvalJob("k", ROUND, self.fam[0])]), [0.0])     # already in flight
        w = self.api.post("/v1/work", json={"node_id": self.node, "slots": 2}).json()
        groups = [[codec.job_from(d) for d in g] for g in w["groups"]]
        self.assertEqual(sorted(len(g) for g in groups), [1, 3])                              # family together
        css = [codec.callset_to(scored(j)) for g in groups for j in g]
        self.api.post("/v1/results", json={"node_id": self.node, "lease_id": w["lease_id"], "callsets": css})
        got = self.ex.poll(1.0)
        self.assertEqual(sorted(c.job_id for c in got), ["j0", "j1", "j2", "j3"])
        self.assertTrue(all(c.score is not None and c.score.score == 90.0 for c in got))
        self.assertEqual(self.ex.idle_slots(), 2)                                              # workers freed

    def test_expired_lease_is_requeued_and_late_results_still_count(self):
        self.ex.submit(self.jobs()[:1])
        w1 = self.api.post("/v1/work", json={"node_id": self.node, "slots": 1}).json()
        job = codec.job_from(w1["groups"][0][0])
        time.sleep(0.7)                                                                         # node went silent
        self.assertEqual(self.ex.poll(0.1), [])                                                 # expiry happens here
        w2 = self.api.post("/v1/work", json={"node_id": self.node, "slots": 1}).json()
        self.assertEqual(codec.job_from(w2["groups"][0][0]).job_id, job.job_id)                # handed out again
        self.api.post("/v1/results", json={"node_id": self.node, "lease_id": w1["lease_id"],     # the old lease
                                           "callsets": [codec.callset_to(scored(job))]})
        self.assertEqual([c.job_id for c in self.ex.poll(1.0)], [job.job_id])
        self.api.post("/v1/results", json={"node_id": self.node, "lease_id": w2["lease_id"],     # duplicate: ignored
                                           "callsets": [codec.callset_to(scored(job, 1.0))]})
        self.assertEqual(self.ex.poll(0.2), [])

    def test_duplicate_submissions_run_once(self):
        a, b = EvalJob("a", ROUND, self.other), EvalJob("b", ROUND, self.other)
        self.ex.submit([a, b])
        w = self.api.post("/v1/work", json={"node_id": self.node, "slots": 2}).json()
        self.assertEqual(sum(len(g) for g in w["groups"]), 1)
        self.api.post("/v1/results", json={"node_id": self.node, "lease_id": w["lease_id"],
                                           "callsets": [codec.callset_to(scored(a))]})
        self.assertEqual(sorted(c.job_id for c in self.ex.poll(1.0)), ["a", "b"])


class CodecTest(unittest.TestCase):
    def test_round_trip(self):
        job = EvalJob("j", ROUND, Config("gatk", {"x": 1}), FULL, 3, "t")
        self.assertEqual(codec.job_from(codec.job_to(job)), job)
        cs = scored(job)
        self.assertEqual(codec.callset_from(codec.callset_to(cs)), cs)
        self.assertEqual(codec.callset_from(codec.callset_to(replace(cs, score=None))), replace(cs, score=None))
