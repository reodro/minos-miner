import asyncio
import os
import time
import unittest
from dataclasses import replace

from answer_set import OracleProvider
from backtest.kpi import oracle_gap, rank_in_field, would_win
from common.settings import SETTINGS
from common.store import ResultStore
from common.types import FULL, AnswerSet, CallSet, ClassCounts, Config, EvalJob, Regime, Round, ScoreResult
from evaluator import FakeRunner, FakeScorer, LocalExecutor, SimExecutor
from main import DryRunSubmitter, EvalClient, LoopSettings, RandomSampler, random_stream, run_round
from tests.helpers import ranges

ROUND_0817 = os.path.join(SETTINGS.rounds_dir, "20260928T081700Z")
ROUND = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")


class StubProvider:
    """Answer set that upgrades itself on every update() call."""

    def __init__(self):
        self.updates = 0

    def build(self, round):
        return AnswerSet(round.round_id, 1, Regime.UNKNOWN, None, [], [], 600)

    def update(self, answer_set, new_callsets):
        self.updates += 1
        return replace(answer_set, version=answer_set.version + 1)


def barrier_strategy(r, pop=16, quorum=1.0, generations=6, seed=1):
    """GA-shaped: waits for each generation before proposing the next."""
    sampler = RandomSampler(r, seed=seed)

    async def run(client):
        for _ in range(generations):
            if client.exhausted:
                return
            await client.map(sampler.sample(pop), quorum=quorum, tag="ga")
    return run


def sequential_strategy(r, seed=2):
    """Coordinate-descent-shaped: one config at a time."""
    sampler = RandomSampler(r, seed=seed)

    async def run(client):
        while not client.exhausted:
            await client.submit(sampler.sample(1)[0], tag="cd")
    return run


class LoopTest(unittest.TestCase):
    def setUp(self):
        self.r = ranges()
        self.runner = FakeRunner(self.r)

    def run_loop(self, strategies, executor=None, provider=None, max_evals=60, baselines=None, **kw):
        sub = DryRunSubmitter()
        ex = executor or SimExecutor(self.runner, workers=8)
        out = run_round(ROUND, provider or StubProvider(), ex, FakeScorer(self.runner), strategies, sub, self.r,
                        LoopSettings(update_every=10, max_evals=max_evals, **kw),
                        baselines=baselines if baselines is not None else [self.r.defaults("gatk")])
        return out, sub, ex

    def test_baseline_first_scores_only_improve(self):
        out, sub, _ = self.run_loop([random_stream(self.r, seed=1)])
        self.assertEqual(sub.sent[0][1], self.r.defaults("gatk").key)
        self.assertEqual(out.evaluated + out.failed, 60)
        scores = [s for _, _, s in out.submissions]
        self.assertEqual(scores, sorted(scores))
        self.assertEqual(out.best.config_key, sub.sent[-1][1])

    def test_answer_set_upgrade_rescores_everything(self):
        provider = StubProvider()
        out, _, _ = self.run_loop([random_stream(self.r, seed=1)], provider=provider)
        self.assertGreaterEqual(provider.updates, 3)
        self.assertEqual(out.best.answer_set_version, out.answer_set_version)

    def test_crashes_do_not_end_the_round(self):
        ex = SimExecutor(self.runner, workers=8, crashes=lambda job: int(job.job_id[1:]) % 3 == 2)
        out, sub, _ = self.run_loop([random_stream(self.r, seed=1)], executor=ex)
        self.assertGreater(out.failed, 0)
        self.assertTrue(sub.sent)

    def test_invalid_and_crashing_baselines_fall_through(self):
        ex = SimExecutor(self.runner, workers=8, crashes=lambda job: job.config.params.get("min_pruning") == 3)
        baselines = [Config("gatk", {"min_pruning": 99}), Config("gatk", {"min_pruning": 3}), self.r.defaults("gatk")]
        out, sub, _ = self.run_loop([], executor=ex, baselines=baselines)
        self.assertEqual(sub.sent[0][1], self.r.defaults("gatk").key)

    def test_mixed_strategies_keep_workers_busy(self):
        # stragglers: 1 in 8 jobs takes 5x longer. A lone barrier strategy idles at every generation;
        # adding a streaming strategy fills the gaps.
        dur = lambda job: 300.0 if int(job.job_id[1:]) % 8 == 0 else 60.0  # noqa: E731
        lone, _, _ = self.run_loop([barrier_strategy(self.r)], executor=SimExecutor(self.runner, 16, dur),
                                   max_evals=96, baselines=[])
        mixed, _, _ = self.run_loop([barrier_strategy(self.r), random_stream(self.r, in_flight=16, seed=3),
                                     sequential_strategy(self.r)],
                                    executor=SimExecutor(self.runner, 16, dur), max_evals=300, baselines=[])
        self.assertLess(lone.utilization, 0.7)
        self.assertGreater(mixed.utilization, 0.85)

    def test_quorum_returns_early_and_stragglers_still_arrive(self):
        dur = lambda job: 500.0 if job.job_id == "j0" else 60.0  # noqa: E731
        seen = []

        async def strategy(client):
            client.on_result.append(lambda o: seen.append(o.config.key))
            res = await client.map(RandomSampler(self.r, seed=5).sample(8), quorum=0.75)
            self.assertIsNone(res[0])
            self.assertEqual(sum(o is not None for o in res), 7)

        self.run_loop([strategy], executor=SimExecutor(self.runner, 8, dur), baselines=[])
        self.assertEqual(len(seen), 8)       # the straggler's result was still delivered after the barrier

    def test_duplicate_configs_run_once(self):
        cfg = RandomSampler(self.r, seed=9).sample(1)[0]

        async def twice(client):
            await asyncio.gather(client.map([cfg]), client.map([cfg]))

        _, _, ex = self.run_loop([twice], baselines=[])
        self.assertAlmostEqual(ex.busy, 60.0)


class DeadlineTest(unittest.TestCase):
    def test_virtual_deadline(self):
        r = ranges(); runner = FakeRunner(r); ex = SimExecutor(runner, 8)
        out = run_round(ROUND, StubProvider(), ex, FakeScorer(runner), [random_stream(r, seed=1)], DryRunSubmitter(),
                        r, LoopSettings(stop_before_deadline_s=120), baselines=[], deadline_ts=72 * 60)
        self.assertLessEqual(ex.now, 72 * 60 - 120 + 60)
        self.assertGreater(out.evaluated, 8 * 50)


class LocalExecutorTest(unittest.TestCase):
    def test_runs_dedups_and_survives_crashes(self):
        r = ranges(); runner = FakeRunner(r)

        class Flaky(FakeRunner):
            def run(self, round, config, fidelity):
                if config.params.get("min_pruning") == 3:
                    raise RuntimeError("Graph must have ref source and sink vertices")
                time.sleep(0.01)
                return super().run(round, config, fidelity)

        ex = LocalExecutor(Flaky(r), workers=4)
        good, bad = r.defaults("gatk"), Config("gatk", {"min_pruning": 3})
        ex.submit([EvalJob("a", ROUND, good), EvalJob("b", ROUND, good), EvalJob("c", ROUND, bad)])
        got = {}
        deadline = time.time() + 5
        while len(got) < 3 and time.time() < deadline:
            for cs in ex.poll(0.5):
                got[cs.job_id] = cs
        self.assertTrue(got["a"].ok and got["b"].ok)
        self.assertEqual(got["a"].sha256, got["b"].sha256)
        self.assertFalse(got["c"].ok)
        self.assertEqual(ex.cost_hint(EvalJob("d", ROUND, good)), 0.0)


class StoreTest(unittest.TestCase):
    def test_roundtrip(self):
        s = ResultStore()
        res = ScoreResult("r1", "k", 2, {"snp_het": ClassCounts(1, 2, 3)}, 5, 600, score=80.0, gate_passed=True)
        s.add(res, Config("gatk", {"min_pruning": 3}), 12.0, "ga")
        self.assertEqual(s.results("r1"), [res])
        self.assertEqual(s.results("r1", answer_set_version=1), [])
        self.assertEqual(s.configs("r1"), [Config("gatk", {"min_pruning": 3})])


class KpiTest(unittest.TestCase):
    def test_oracle_gap(self):
        self.assertEqual(oracle_gap({"a": 90, "b": 80}, {"a": 88, "b": 89}), 1)
        self.assertEqual(oracle_gap({"a": 90, "b": 80}, {"a": 89, "b": 88}), 0)

    def test_rank(self):
        self.assertEqual(rank_in_field(89.0, [89.3, 89.1, 88.0]), 3)
        self.assertTrue(would_win(90.0, [89.3, 89.1]))


@unittest.skipUnless(os.path.exists(os.path.join(ROUND_0817, "truth.vcf.gz")), "round 08:17 not downloaded")
class OracleTest(unittest.TestCase):
    def test_counts_match_released_files(self):
        r = Round("2026-09-28T08:17:00+00:00", "chr14", 45912001, 53254000, "", "", released_dir=ROUND_0817)
        a = OracleProvider().build(r)
        self.assertEqual((len(a.germline), len(a.synthetic), a.n_targets), (14744, 5056, 606))
        self.assertIs(OracleProvider().update(a, [CallSet(r.round_id, "k", "", "", 0)]), a)


if __name__ == "__main__":
    unittest.main()
