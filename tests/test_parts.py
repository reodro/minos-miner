import asyncio
import os
import time
import unittest
from dataclasses import replace

from answer_set import OracleProvider
from backtest.kpi import oracle_gap, rank_in_field, would_win
from common.settings import SETTINGS
from common.store import ResultStore
from common.types import FULL, AnswerSet, CallSet, ClassCounts, Config, EvalJob, Regime, Round, ScoreResult, Variant
from evaluator import FakeRunner, FakeScorer, LocalExecutor, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search import EvalClient, RandomSampler, random_stream, rf_ucb
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


class BatchingExecutorTest(unittest.TestCase):
    """LocalExecutor batching with a fake runner: grouping, diverged re-run, cost hints."""

    class BatchRunner(FakeRunner):
        GT = "standard_min_confidence_threshold_for_calling"

        def __init__(self, r, diverge=()):
            super().__init__(r)
            self.batches, self.solo, self.diverge = [], [], set(diverge)

        def branch_key(self, config, fidelity):
            return Config("gatk", {k: v for k, v in config.params.items() if k != self.GT}).key + fidelity

        def variant_cost(self, main, variant, fidelity):
            return 0.02 if self.branch_key(main, fidelity) == self.branch_key(variant, fidelity) else None

        def run(self, round, config, fidelity):
            time.sleep(0.05)
            self.solo.append(config.key)
            return super().run(round, config, fidelity)

        def run_batch(self, round, configs, fidelity):
            time.sleep(0.05)
            self.batches.append([c.key for c in configs])
            out = []
            for i, c in enumerate(configs):
                cs = super().run(round, c, fidelity)
                out.append(replace(cs, ok=False, error="minos:diverged") if i and c.key in self.diverge else cs)
            return out

    def drain(self, ex, n):
        got, deadline = {}, time.time() + 10
        while len(got) < n and time.time() < deadline:
            for cs in ex.poll(0.2):
                got[cs.job_id] = cs
        return got

    def test_groups_queued_genotyping_variants_and_reruns_diverged(self):
        r = ranges(); d = r.defaults("gatk").params
        conf = [Config("gatk", {**d, self.BatchRunner.GT: c}) for c in (30.0, 40.0, 50.0, 60.0)]
        other = Config("gatk", {**d, "min_pruning": 4})
        runner = self.BatchRunner(r, diverge={conf[2].key})
        ex = LocalExecutor(runner, workers=1)
        blocker = Config("gatk", {**d, "min_pruning": 7})
        ex.submit([EvalJob("block", ROUND, blocker)])              # occupies the only worker
        jobs = [EvalJob(f"c{i}", ROUND, c) for i, c in enumerate(conf)] + [EvalJob("o", ROUND, other)]
        self.assertEqual(ex.cost_hint(jobs[0]), 1.0)
        ex.submit(jobs[:1])
        self.assertEqual(ex.cost_hint(jobs[1]), 0.02)              # would ride along with c0: the runner's variant cost
        ex.submit(jobs[1:])
        got = self.drain(ex, 6)
        self.assertEqual(len(got), 6)
        self.assertTrue(all(cs.ok for cs in got.values()))
        self.assertIn([c.key for c in conf], runner.batches)       # one batch for the four variants
        self.assertIn(conf[2].key, runner.solo)                    # the diverged one re-ran alone
        self.assertIn(other.key, runner.solo)                      # different branch key: not batched

    def test_no_batching_without_run_batch(self):
        r = ranges(); runner = FakeRunner(r)
        ex = LocalExecutor(runner, workers=2)
        cfgs = RandomSampler(r, seed=4).sample(5)
        ex.submit([EvalJob(f"x{i}", ROUND, c) for i, c in enumerate(cfgs)])
        self.assertEqual(len(self.drain(ex, 5)), 5)
        self.assertEqual(ex.batches_run, 0)


class FamilySubmitTest(unittest.TestCase):
    """submit_many: a family reaches the executor in one call, runs as one batch, and is charged as a family."""

    def test_family_is_one_batch_and_cheap(self):
        r = ranges(); d = r.defaults("gatk").params
        runner = BatchingExecutorTest.BatchRunner(r)
        runner.variant_cost = lambda main, v, fid: (0.02 if runner.branch_key(main, fid) == runner.branch_key(v, fid)
                                                    else None)
        ex = LocalExecutor(runner, workers=1)
        gt = BatchingExecutorTest.BatchRunner.GT
        fam = [Config("gatk", {**d, gt: c}) for c in (30.0, 40.0, 50.0)]

        async def go():
            client = EvalClient(ROUND, ex, FakeScorer(runner), r, ResultStore(), StubProvider().build(ROUND), budget=10)
            futs = client.submit_many(fam)
            while client.in_flight:
                for cs in await asyncio.to_thread(ex.poll, 0.2):
                    client.deliver(cs)
            return client.spent, [f.result().ok for f in futs]

        spent, oks = asyncio.run(go())
        self.assertEqual(runner.batches, [[c.key for c in fam]])
        self.assertTrue(all(oks))
        self.assertAlmostEqual(spent, 1.0 + 2 * 0.02)


class RfUcbTest(unittest.TestCase):
    def test_beats_random_on_the_fake_landscape_and_respects_budget(self):
        r = ranges()
        best = {"random": [], "rf": []}
        for seed in range(3):
            for name, strat in (("random", random_stream(r, in_flight=8, seed=seed)), ("rf", rf_ucb(r, seed=seed))):
                runner = FakeRunner(r, seed=seed)
                out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner), [strat],
                                DryRunSubmitter(), r, LoopSettings(max_evals=120), baselines=[])
                self.assertLessEqual(out.evaluated + out.failed, 120 + 8)     # in-flight ones may finish after
                best[name].append(out.best.score)
        self.assertGreater(sum(best["rf"]), sum(best["random"]))


class ResubmitIntervalTest(unittest.TestCase):
    def test_interval_limits_submissions_and_the_best_still_goes_out_last(self):
        r = ranges(); runner = FakeRunner(r)
        sub = DryRunSubmitter()
        out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner), [random_stream(r, seed=1)],
                        sub, r, LoopSettings(max_evals=120, resubmit_interval_s=600.0), baselines=[r.defaults("gatk")])
        times = [t for t, _, _ in out.submissions]
        self.assertTrue(all(b - a >= 600.0 for a, b in zip(times[:-2], times[1:-1])))    # the final one may be forced
        self.assertEqual(sub.sent[-1][1], out.best.config_key)


class FallbackChainTest(unittest.TestCase):
    """mine_round never ends a round without a submission, and never searches without a truth."""

    def parts(self, provider):
        from main import Parts
        r = ranges(); runner = FakeRunner(r)
        return Parts(provider, SimExecutor(runner, 8), FakeScorer(runner), [random_stream(r, seed=1)],
                     DryRunSubmitter(), r, [Config("gatk", {"min_pruning": 99}), r.defaults("gatk")],
                     LoopSettings(max_evals=40))

    def test_answer_set_failure_submits_the_first_valid_baseline(self):
        from main import mine_round

        class Broken:
            def build(self, round):
                raise RuntimeError("parent fetch down")

        p = self.parts(Broken())
        m = mine_round(ROUND, p, deadline_ts=None)
        self.assertEqual((m.mode, p.submitter.sent), ("static", [(ROUND.round_id, ranges().defaults("gatk").key)]))

    def test_no_donor_means_no_search(self):
        from main import mine_round
        p = self.parts(StubProvider())                  # Regime.UNKNOWN, no germline
        self.assertEqual(mine_round(ROUND, p, deadline_ts=None).mode, "static")

    def test_known_donor_searches(self):
        from main import mine_round

        class Known(StubProvider):
            def build(self, round):
                return AnswerSet(round.round_id, 1, Regime.KG, "NA00001",
                                 [Variant("chr20", 5, "A", "G", "0/1")], [], 600)

        p = self.parts(Known())
        m = mine_round(ROUND, p, deadline_ts=None)
        self.assertEqual(m.mode, "searched")
        self.assertGreater(m.outcome.evaluated, 10)
        self.assertEqual(p.submitter.sent[-1][1], m.outcome.best.config_key)


class SpaceTest(unittest.TestCase):
    def test_unscorable_settings_are_never_proposed(self):
        import random as _r
        from search.space import Space
        sp, rng = Space(ranges()), _r.Random(0)
        cs = [sp.random(rng) for _ in range(200)] + RandomSampler(ranges(), seed=0).sample(200)
        cs += [sp.neighbour(c, rng, n_params=4) for c in cs[:100]]
        self.assertEqual({c.params["emit_ref_confidence"] for c in cs if c}, {"NONE"})
        self.assertNotIn("emit_ref_confidence", sp.names)
