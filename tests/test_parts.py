import asyncio
import os
import time
import unittest
from dataclasses import replace

from answer_set import OracleProvider
from backtest.kpi import oracle_gap, rank_in_field, would_win
from common.settings import SETTINGS
from common.store import ResultStore
from common.types import FULL, CLASSES, AnswerSet, CallSet, ClassCounts, Config, EvalJob, Regime, Round, ScoreResult, Variant
from evaluator import FakeRunner, FakeScorer, LocalExecutor, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search import EvalClient, RandomSampler, random_stream, rf_ucb, tpe
from tests.helpers import ranges, together

ROUND_0817 = os.path.join(SETTINGS.rounds_dir, "20260928T081700Z")
def jobno(job):
    """The running number of an EvalClient job id ("<round tag>-j<n>")."""
    return int(job.job_id.rsplit("-j", 1)[1])


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
            await client.propose_all(sampler.sample(pop), quorum=quorum, tag="ga")
    return run


def sequential_strategy(r, seed=2):
    """Coordinate-descent-shaped: one config at a time."""
    sampler = RandomSampler(r, seed=seed)

    async def run(client):
        while not client.exhausted:
            await client.propose(sampler.sample(1)[0], tag="cd")
    return run


class LoopTest(unittest.TestCase):
    def setUp(self):
        self.r = ranges()
        self.runner = FakeRunner(self.r)

    def run_loop(self, strategies, executor=None, provider=None, max_evals=60, baselines=None, accept=None,
                 scoring=None, **kw):
        sub = DryRunSubmitter()
        ex = executor or SimExecutor(self.runner, workers=8)
        out = run_round(ROUND, provider or StubProvider(), ex, FakeScorer(self.runner), together(*strategies), sub, self.r,
                        LoopSettings(update_every=10, max_evals=max_evals, **kw),
                        baselines=baselines if baselines is not None else [self.r.defaults("gatk")], accept=accept,
                        scoring=scoring)
        return out, sub, ex

    def test_a_scoring_pool_scores_units_in_parallel_and_reports_what_it_cannot_score(self):
        from evaluator import ScoringPool
        ex = LocalExecutor(self.runner, workers=4)
        out, sub, _ = self.run_loop([random_stream(self.r, seed=1)], executor=ex, max_evals=30,
                                    scoring=ScoringPool(3, scorer=FakeScorer(self.runner)))
        self.assertEqual(out.evaluated + out.failed, 30)
        self.assertEqual(out.failed, 0)
        self.assertEqual(sub.sent[-1][1], out.best.config_key)

        class Down:                                          # hap.py failing: every callset is a failure, never a 0.0
            def score(self, cs, aset):
                raise RuntimeError("hap.py down")

        out, sub, _ = self.run_loop([random_stream(self.r, seed=2)], executor=LocalExecutor(self.runner, workers=2),
                                    max_evals=12, scoring=ScoringPool(2, scorer=Down()))
        self.assertEqual((out.evaluated, out.failed, sub.sent), (0, 12, []))
        self.assertEqual(list(out.failures), ["scoring failed: hap.py down"])

    def test_gatk_runs_while_the_answer_set_is_built_and_its_results_wait(self):
        ex = LocalExecutor(self.runner, workers=4)          # real threads: the executor runs while the truth is built
        seen = []

        class Late(StubProvider):
            def build(self, round):
                time.sleep(0.6)
                seen.append(ex.done)                         # callsets the executor had delivered by then
                return AnswerSet(round.round_id, 1, Regime.KG, "S", [], [], 600)

        out, sub, _ = self.run_loop([random_stream(self.r, seed=1)], executor=ex, provider=Late(), max_evals=20)
        self.assertGreater(seen[0], 0)                       # GATK results existed before the answer set did
        self.assertEqual(out.evaluated + out.failed, 20)     # ... and every one of them was scored once it arrived
        self.assertIsNotNone(out.best)
        self.assertEqual(sub.sent[0][1], self.r.defaults("gatk").key)

    def test_a_build_that_never_finishes_ends_the_round_on_time(self):
        from main.loop import AnswerSetUnavailable

        class Hung(StubProvider):
            def build(self, round):
                time.sleep(4.0)
                return AnswerSet(round.round_id, 1, Regime.KG, "S", [], [], 600)

        t0 = time.time()
        with self.assertRaises(AnswerSetUnavailable):
            self.run_loop([random_stream(self.r, seed=1)], executor=LocalExecutor(self.runner, workers=2),
                          provider=Hung(), max_evals=20, build_timeout_s=0.5)
        self.assertLess(time.time() - t0, 3.0)                 # the round ended at the guard, not when the build did

    def test_a_refused_answer_set_ends_the_round(self):
        from main.loop import AnswerSetUnavailable
        with self.assertRaises(AnswerSetUnavailable):
            self.run_loop([random_stream(self.r, seed=1)], accept=lambda a: "no donor identified")

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
        ex = SimExecutor(self.runner, workers=8, crashes=lambda job: jobno(job) % 3 == 2)
        out, sub, _ = self.run_loop([random_stream(self.r, seed=1)], executor=ex)
        self.assertGreater(out.failed, 0)
        self.assertEqual(sum(out.failures.values()), out.failed)        # every failure has a counted reason
        self.assertTrue(sub.sent)

    def test_invalid_and_crashing_baselines_fall_through(self):
        ex = SimExecutor(self.runner, workers=8, crashes=lambda job: job.config.params.get("min_pruning") == 3)
        baselines = [Config("gatk", {"min_pruning": 99}), Config("gatk", {"min_pruning": 3}), self.r.defaults("gatk")]
        out, sub, _ = self.run_loop([], executor=ex, baselines=baselines)
        self.assertEqual(sub.sent[0][1], self.r.defaults("gatk").key)

    def test_mixed_strategies_keep_workers_busy(self):
        # stragglers: 1 in 8 jobs takes 5x longer. A lone barrier strategy idles at every generation;
        # adding a streaming strategy fills the gaps.
        dur = lambda job: 300.0 if jobno(job) % 8 == 0 else 60.0  # noqa: E731
        lone, _, _ = self.run_loop([barrier_strategy(self.r)], executor=SimExecutor(self.runner, 16, dur),
                                   max_evals=96, baselines=[])
        mixed, _, _ = self.run_loop([barrier_strategy(self.r), random_stream(self.r, in_flight=16, seed=3),
                                     sequential_strategy(self.r)],
                                    executor=SimExecutor(self.runner, 16, dur), max_evals=300, baselines=[])
        self.assertLess(lone.utilization, 0.7)
        self.assertGreater(mixed.utilization, 0.85)

    def test_quorum_returns_early_and_stragglers_still_arrive(self):
        dur = lambda job: 500.0 if jobno(job) == 0 else 60.0  # noqa: E731
        seen = []

        got = []

        async def strategy(client):
            client.on_result.append(lambda o: seen.append(o.config.key))
            got.append(await client.propose_all(RandomSampler(self.r, seed=5).sample(8), quorum=0.75))

        self.run_loop([strategy], executor=SimExecutor(self.runner, 8, dur), baselines=[])
        # asserted out here: the round loop logs (does not raise) exceptions from a strategy
        self.assertIsNone(got[0][0])                              # the straggler was not waited for
        self.assertEqual(sum(o is not None for o in got[0]), 7)
        self.assertEqual(len(seen), 8)       # the straggler's result was still delivered after the barrier

    def test_duplicate_configs_run_once(self):
        cfg = RandomSampler(self.r, seed=9).sample(1)[0]

        async def twice(client):
            await asyncio.gather(client.propose_all([cfg]), client.propose_all([cfg]))

        _, _, ex = self.run_loop([twice], baselines=[])
        self.assertAlmostEqual(ex.busy, 60.0)


class DeadlineTest(unittest.TestCase):
    def test_virtual_deadline(self):
        r = ranges(); runner = FakeRunner(r); ex = SimExecutor(runner, 8)
        out = run_round(ROUND, StubProvider(), ex, FakeScorer(runner), random_stream(r, seed=1), DryRunSubmitter(),
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

    def test_splits_a_family_that_would_run_for_minutes(self):
        r = ranges(); d = r.defaults("gatk").params
        conf = [Config("gatk", {**d, self.BatchRunner.GT: c}) for c in (30.0, 40.0, 50.0, 60.0)]
        runner = self.BatchRunner(r)
        runner.variant_cost = lambda main, variant, fidelity: 0.6   # 132 + 79 s; a third variant passes 4.5 min
        ex = LocalExecutor(runner, workers=2)
        ex.submit([EvalJob(f"c{i}", ROUND, c) for i, c in enumerate(conf)])
        got = self.drain(ex, 4)
        self.assertEqual(len(got), 4)
        self.assertEqual(sorted(len(b) for b in runner.batches), [2, 2])

    def test_a_later_genotyping_variant_uses_the_warm_jvm(self):
        r = ranges(); d = r.defaults("gatk").params
        runner = self.BatchRunner(r)
        runner.runs = []

        def run(round, config, fidelity):
            runner.runs.append(config.key)
            return BatchingExecutorTest.BatchRunner.run(runner, round, config, fidelity)

        def warm_submit(round, config, fidelity):
            if runner.runs and config.params.get(self.BatchRunner.GT) == 40.0:
                runner.warmed = config.key
                return BatchingExecutorTest.BatchRunner.run(runner, round, config, fidelity)
            return None

        runner.run = run
        runner.warm_submit = warm_submit
        ex = LocalExecutor(runner, workers=1)
        first = Config("gatk", {**d, self.BatchRunner.GT: 30.0})
        second = Config("gatk", {**d, self.BatchRunner.GT: 40.0})
        ex.submit([EvalJob("a", ROUND, first)])
        got = self.drain(ex, 1)
        ex.submit([EvalJob("b", ROUND, second)])
        got.update(self.drain(ex, 1))
        self.assertEqual(len(got), 2)
        self.assertEqual(runner.warmed, second.key)
        self.assertEqual(runner.runs, [first.key])          # the second config did not start its own GATK run

    def test_queued_genotyping_variants_reuse_a_warm_jvm(self):
        r = ranges(); d = r.defaults("gatk").params
        runner = self.BatchRunner(r)
        runner.warmed = []
        runner.warm_ready = lambda round, config, fidelity: True

        def warm_submit(round, config, fidelity):
            runner.warmed.append(config.key)
            return FakeRunner.run(runner, round, config, fidelity)

        runner.warm_submit = warm_submit
        ex = LocalExecutor(runner, workers=1, now=lambda: 0.0)
        ex.set_finish_by(10.0)                             # a solo run is 132 s; a warm genotype is about 1.5 s
        conf = [Config("gatk", {**d, self.BatchRunner.GT: c}) for c in (40.0, 50.0)]
        ex.submit([EvalJob("a", ROUND, conf[0]), EvalJob("b", ROUND, conf[1])])
        got = self.drain(ex, 2)
        self.assertEqual(sorted(got), ["a", "b"])
        self.assertEqual(runner.warmed, [c.key for c in conf])
        self.assertEqual(runner.batches, [])
        self.assertEqual(runner.solo, [])

    def test_does_not_start_a_job_that_cannot_finish_and_still_starts_a_repair(self):
        r = ranges(); d = r.defaults("gatk").params
        runner = self.BatchRunner(r)
        runner.repairs = []
        runner.repair = lambda round, config, diverged, fidelity: runner.repairs.append(config.key) or FakeRunner.run(
            runner, round, config, fidelity)
        now = [0.0]
        ex = LocalExecutor(runner, workers=1, now=lambda: now[0])
        ex.set_finish_by(40.0)                                 # a solo run is 132 s; a repair is 34 s
        normal = EvalJob("n", ROUND, Config("gatk", {**d, "min_pruning": 4}))
        repair_job = EvalJob("r", ROUND, Config("gatk", {**d, "min_pruning": 5}))
        ex._solo.add(repair_job.job_id)
        ex._repair[repair_job.job_id] = CallSet(ROUND.round_id, repair_job.config.key, "batch.vcf.gz", "", 1.0,
                                                ok=False, error="minos:diverged", job_id=repair_job.job_id)
        ex.submit([normal, repair_job])
        got = self.drain(ex, 1)
        self.assertEqual(list(got), ["r"])
        self.assertEqual(runner.repairs, [repair_job.config.key])
        self.assertEqual(runner.solo, [])
        self.assertEqual(ex.stats()["queued"], 1)

    def test_measured_speed_holds_a_run_the_table_says_fits_and_the_timeline_says_why(self):
        r = ranges(); d = r.defaults("gatk").params
        now = [0.0]
        runner = self.BatchRunner(r)
        plain = runner.run

        def slow(round, config, fidelity):
            now[0] += 264.0                                    # this machine: twice the 132 s table
            return plain(round, config, fidelity)

        runner.run = slow
        ex = LocalExecutor(runner, workers=1, now=lambda: now[0])
        ex.submit([EvalJob("a", ROUND, Config("gatk", {**d, "min_pruning": 4}))])
        self.assertEqual(list(self.drain(ex, 1)), ["a"])
        self.assertEqual((ex.stats()["solo_s"], ex.stats()["scale"]), (264, 2.0))
        ex.set_finish_by(now[0] + 200.0)                       # the table (132 s) fits, this machine (264 s) does not
        ex.submit([EvalJob("b", ROUND, Config("gatk", {**d, "min_pruning": 5}))])
        self.assertEqual((ex.stats()["running"], ex.stats()["queued"]), (0, 1))
        now[0] += 60.0
        ex.poll(0)
        last = ex.stats()["timeline"][-1]
        self.assertEqual((last["busy"], last["idle"], last["queued"], last["held"]), (0, 1, 1, 1))
        self.assertEqual(last["scale"], 2.0)

    def test_speed_ignores_warm_and_instant_runs(self):
        r = ranges(); d = r.defaults("gatk").params
        ex = LocalExecutor(self.BatchRunner(r), workers=1)    # fake runs take 0.05 s: not this machine's speed
        ex.submit([EvalJob("a", ROUND, Config("gatk", {**d, "min_pruning": 4}))])
        self.drain(ex, 1)
        self.assertEqual((ex.stats()["solo_s"], ex.stats()["scale"]), (132, 1.0))

    def test_a_haplotype_cap_the_run_never_reached_reuses_that_run(self):
        r = ranges(); d = r.defaults("gatk").params
        MNH = "max_num_haplotypes_in_population"
        runner = self.BatchRunner(r)
        runner.kbest_needed = lambda round_id, key, fidelity: 40   # every run here needed a cap of 40
        ex = LocalExecutor(runner, workers=1)
        blocker = Config("gatk", {**d, "min_pruning": 7})
        first, higher, lower = (Config("gatk", {**d, MNH: v}) for v in (128, 64, 30))
        ex.submit([EvalJob("block", ROUND, blocker), EvalJob("first", ROUND, first), EvalJob("higher", ROUND, higher)])
        got = self.drain(ex, 3)                                 # "higher" was queued when "first" taught the class
        ex.submit([EvalJob("again", ROUND, Config("gatk", {**d, MNH: 50})), EvalJob("lower", ROUND, lower)])
        got.update(self.drain(ex, 2))
        self.assertEqual(sorted(got), ["again", "block", "first", "higher", "lower"])
        self.assertEqual(runner.solo, [blocker.key, first.key, lower.key])   # 64 and 50 never ran: 30 < 40 did
        self.assertEqual(got["higher"].config_key, higher.key)               # each under the config it asked for
        self.assertEqual(got["higher"].sha256, got["first"].sha256)
        self.assertEqual(got["again"].config_key, Config("gatk", {**d, MNH: 50}).key)
        self.assertNotEqual(got["lower"].sha256, got["first"].sha256)

    def test_configs_the_bam_cannot_tell_apart_run_once(self):
        from evaluator.equivalence import RoundFacts
        r = ranges(); d = r.defaults("gatk").params
        runner = self.BatchRunner(r)
        ex = LocalExecutor(runner, workers=1)
        ex._facts[ROUND.round_id] = RoundFacts(frozenset({2, 14, 30, 40}), frozenset({60}), 12)
        a, b, c = (Config("gatk", {**d, "min_base_quality_score": v}) for v in (10, 13, 14))   # 3..13 | 14
        ex.submit([EvalJob("a", ROUND, a), EvalJob("b", ROUND, b), EvalJob("c", ROUND, c)])
        got = self.drain(ex, 3)
        self.assertEqual(len(runner.solo), 2)                  # 10 and 13 are one run, 14 another
        self.assertEqual([got[j].config_key for j in "abc"], [a.key, b.key, c.key])
        self.assertEqual(got["a"].sha256, got["b"].sha256)
        self.assertEqual(ex.stats()["canonicalized"], 3)       # each ran as its class's representative

    def test_no_batching_without_run_batch(self):
        r = ranges(); runner = FakeRunner(r)
        ex = LocalExecutor(runner, workers=2)
        cfgs = RandomSampler(r, seed=4).sample(5)
        ex.submit([EvalJob(f"x{i}", ROUND, c) for i, c in enumerate(cfgs)])
        self.assertEqual(len(self.drain(ex, 5)), 5)
        self.assertEqual(ex.batches_run, 0)


class FamilySubmitTest(unittest.TestCase):
    """propose_many: a family reaches the executor in one call, runs as one batch, and is charged as a family."""

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
            futs = client.propose_many(fam)
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
                out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner), strat,
                                DryRunSubmitter(), r, LoopSettings(max_evals=200), baselines=[])
                self.assertLessEqual(out.evaluated + out.failed, 200 + 8)     # in-flight ones may finish after
                best[name].append(out.best.score)
        self.assertGreater(sum(best["rf"]), sum(best["random"]))


class ResubmitIntervalTest(unittest.TestCase):
    def test_interval_limits_submissions_and_the_best_still_goes_out_last(self):
        r = ranges(); runner = FakeRunner(r)
        sub = DryRunSubmitter()
        out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner), random_stream(r, seed=1),
                        sub, r, LoopSettings(max_evals=120, resubmit_interval_s=600.0), baselines=[r.defaults("gatk")])
        times = [t for t, _, _ in out.submissions]
        self.assertTrue(all(b - a >= 600.0 for a, b in zip(times[:-2], times[1:-1])))    # the final one may be forced
        self.assertEqual(sub.sent[-1][1], out.best.config_key)


class FallbackChainTest(unittest.TestCase):
    """mine_round never ends a round without a submission, and never searches without a truth."""

    def parts(self, provider):
        from main import Parts
        r = ranges(); runner = FakeRunner(r)
        return Parts(provider, SimExecutor(runner, 8), FakeScorer(runner), random_stream(r, seed=1),
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
        own = [sp.random(rng) for _ in range(200)]
        own += [sp.neighbour(c, rng, n_params=4) for c in own[:100]] + sp.initial_design(rng, 32)
        cs = own + RandomSampler(ranges(), seed=0).sample(200)
        self.assertEqual({c.params["emit_ref_confidence"] for c in cs if c}, {"NONE"})
        self.assertNotIn("emit_ref_confidence", sp.names)
        self.assertEqual({c.params["adaptive_pruning_initial_error_rate"] for c in cs if c}, {0.001})
        self.assertEqual({c.params["pruning_lod_threshold"] for c in cs if c}, {2.302585})
        self.assertNotIn("adaptive_pruning_initial_error_rate", sp.names)
        self.assertNotIn("pruning_lod_threshold", sp.names)
        self.assertGreaterEqual(min(c.params["assembly_region_padding"] for c in own if c), 5)   # the space's floor
        self.assertEqual(sp._assemble({"assembly_region_padding": 0}).params["assembly_region_padding"], 5)
        capped = Space(ranges(), caps={"min_base_quality_score": 30})          # the round's BAM tops out at 30
        rng = _r.Random(1)
        got = [capped.random(rng) for _ in range(100)] + capped.initial_design(rng, 16)
        self.assertLessEqual(max(c.params["min_base_quality_score"] for c in got if c), 30)
        self.assertEqual(capped._assemble({"min_base_quality_score": 45}).params["min_base_quality_score"], 30)
        from search.space import bam_caps
        self.assertEqual(bam_caps("/nonexistent.bam", "chr20", 1, 10), {})

    def test_tiers_partition_the_free_parameters(self):
        from search.space import Space
        sp = Space(ranges())
        self.assertEqual(set(sp.genotyping) | set(sp.pairhmm) | set(sp.expensive), set(sp.names))
        self.assertFalse(set(sp.genotyping) & set(sp.pairhmm))
        self.assertFalse(set(sp.expensive) & (set(sp.genotyping) | set(sp.pairhmm)))
        self.assertIn("heterozygosity", sp.expensive)
        self.assertIn("min_pruning", sp.expensive)
        self.assertIn("pcr_indel_model", sp.pairhmm)
        self.assertIn("standard_min_confidence_threshold_for_calling", sp.genotyping)

    def test_sobol_init_keeps_cheap_parameters_at_defaults(self):
        import random as _r
        from search.space import Space
        sp, rng = Space(ranges()), _r.Random(1)
        design = sp.initial_design(rng, 16)
        self.assertEqual(len(design), 16)
        self.assertEqual(len({c.key for c in design}), 16)
        self.assertEqual(design[0].key, ranges().defaults("gatk").key)
        changed = [n for n in sp.names if design[1].params[n] != design[0].params[n]]
        self.assertEqual(changed, ["min_pruning"])
        self.assertEqual(design[1].params["min_pruning"], 3)
        for c in design:
            self.assertEqual(ranges().errors(c), [])
            for name in sp.genotyping + sp.pairhmm:
                self.assertEqual(c.params[name], sp.specs[name].default, name)

    def test_designed_variants_touch_one_cheap_axis(self):
        from search.space import Space
        sp = Space(ranges())
        base = ranges().defaults("gatk")
        for c in sp.genotyping_variants(base):
            changed = [n for n in sp.names if c.params[n] != base.params[n]]
            self.assertEqual(len(changed), 1)
            self.assertIn(changed[0], sp.genotyping)
        models = {c.params["pcr_indel_model"] for c in sp.pairhmm_variants(base)}
        self.assertEqual(models, set(sp.specs["pcr_indel_model"].allowed) - {base.params["pcr_indel_model"]})


class LineupTest(unittest.TestCase):
    def test_search_hands_main_one_searcher(self):
        from search import searcher
        self.assertTrue(callable(searcher(ranges(), 8)))       # the basin race; main chooses nothing


class GradientRegionTest(unittest.TestCase):
    def test_beats_random_on_the_fake_landscape_and_keeps_regions(self):
        from search import gradient_regions
        r = ranges()
        best = {"regions": [], "random": []}
        born = []
        alive = []
        for seed in range(3):
            for name, strat in (("random", random_stream(r, in_flight=8, seed=seed)),
                                ("regions", gradient_regions(r, in_flight=8, seed=seed))):
                runner = FakeRunner(r, seed=seed)
                out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner), strat,
                                DryRunSubmitter(), r, LoopSettings(max_evals=800), baselines=[])
                self.assertLessEqual(out.evaluated + out.failed, 800 + 8)
                best[name].append(out.best.score)
                if name == "regions":
                    born.append(strat.state["born"])
                    alive.append(strat.state["alive"])
        self.assertGreater(sum(best["regions"]), sum(best["random"]), msg=str(best))
        self.assertTrue(all(n >= 2 for n in born), msg=str(born))
        self.assertTrue(all(n >= 1 for n in alive), msg=str(list(zip(born, alive))))

    def test_restarts_leave_a_local_basin(self):
        import math
        from search import gradient_regions
        from search.space import Space
        r = ranges()

        class TwoBasin(FakeRunner):
            def quality(self, config):
                specs = self.ranges.tools[config.tool]

                def unit(name):
                    spec = specs[name]
                    v = float(config.params.get(name, spec.default))
                    return (v - spec.min) / (spec.max - spec.min)

                up, uh = unit("min_pruning"), unit("heterozygosity")
                local = math.exp(-((up - 0.125) ** 2 + (uh - 0.0909) ** 2) / 0.08)
                glob = math.exp(-((up - 0.875) ** 2 + (uh - 0.8485) ** 2) / 0.08)
                return 0.50 * local + 0.99 * glob

        params = dict(r.defaults("gatk").params)
        params["min_pruning"] = 3
        params["heterozygosity"] = 0.001
        local = Space(r)._assemble(params)
        self.assertIsNotNone(local)

        def once(strat, evals):
            runner = TwoBasin(r)
            out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner), strat,
                            DryRunSubmitter(), r, LoopSettings(max_evals=evals), baselines=[])
            return out.best.score

        trapped = once(gradient_regions(r, n_regions=1, restart_frac=0.0, pair_frac=0.0, start=[local], seed=0), 80)
        escaped = once(gradient_regions(r, in_flight=8, seed=1), 240)
        self.assertGreater(escaped, trapped)


class LateWindowTest(unittest.TestCase):
    def test_last_minutes_do_not_start_full_runs_the_executor_cannot_batch(self):
        r = ranges()
        runner = FakeRunner(r)
        out = run_round(ROUND, StubProvider(), SimExecutor(runner, 4), FakeScorer(runner),
                        rf_ucb(r, late_s=1e9, seed=0), DryRunSubmitter(), r,
                        LoopSettings(max_evals=40), baselines=[], deadline_ts=10_000)
        self.assertEqual(out.evaluated, 0)


class TpeTest(unittest.TestCase):
    def test_runs_inside_the_budget(self):
        r = ranges()
        runner = FakeRunner(r)
        out = run_round(ROUND, StubProvider(), SimExecutor(runner, 4), FakeScorer(runner),
                        tpe(r, seed=0, in_flight=4), DryRunSubmitter(), r,
                        LoopSettings(max_evals=40), baselines=[])
        self.assertGreater(out.evaluated, 10)
        self.assertLessEqual(out.evaluated + out.failed, 40 + 8)


class PortfolioTest(unittest.TestCase):
    def test_only_earlier_rounds_configs_are_reused(self):
        import os
        import tempfile
        from common.scoring import finalize
        d = tempfile.mkdtemp()
        past = ResultStore(os.path.join(d, "20261001T000000Z.sqlite"))
        current = ResultStore(os.path.join(d, "20261002T000000Z.sqlite"))
        later = ResultStore(os.path.join(d, "20261003T000000Z.sqlite"))
        cfg = ranges().defaults("gatk")
        future = Config("gatk", {**cfg.params, "min_pruning": 4})
        res = lambda c: finalize(ScoreResult("r", c.key, 1, {k: ClassCounts(1, 0, 0) for k in CLASSES}, 0, 600))
        past.add(res(cfg), cfg)
        later.add(res(future), future)
        self.assertEqual([c.key for c in current.portfolio()], [cfg.key])   # the later round never leaks back
        self.assertEqual(past.portfolio(), [])
        self.assertEqual(ResultStore(os.path.join(d, "scratch.sqlite")).portfolio(), [])


class RoundSafetyTest(unittest.TestCase):
    def test_a_late_loop_failure_does_not_replace_a_submitted_config_with_the_baseline(self):
        from main import Parts, mine_round
        r = ranges(); runner = FakeRunner(r)

        class Known(StubProvider):
            def build(self, round):
                return AnswerSet(round.round_id, 1, Regime.KG, "S", [Variant("chr20", 5, "A", "G", "0/1")], [], 600)

        async def explode_later(client):
            await client.propose_all(RandomSampler(r, seed=1).sample(8))
            await asyncio.sleep(0)
            raise RuntimeError("strategy bug")

        sub = DryRunSubmitter()
        base = r.defaults("gatk")
        parts = Parts(Known(), SimExecutor(runner, 8), FakeScorer(runner), explode_later, sub, r, [base],
                      LoopSettings(max_evals=40))
        m = mine_round(ROUND, parts, deadline_ts=None)
        self.assertEqual(m.mode, "searched")
        self.assertNotEqual(sub.sent[-1][1], base.key)          # the best config found, not the baseline again

    def test_results_of_another_round_are_ignored(self):
        r = ranges(); runner = FakeRunner(r)

        async def go():
            client = EvalClient(ROUND, SimExecutor(runner, 2), FakeScorer(runner), r, ResultStore(),
                                StubProvider().build(ROUND))
            fut = client.propose(r.defaults("gatk"))
            job_id = next(iter(client._pending))
            other = CallSet("another-round", r.defaults("gatk").key, "fake://x", "d", 0.0, job_id=job_id)
            client.deliver(other)                                # same job id, other round: must not resolve it
            return fut.done(), client.in_flight
        self.assertEqual(asyncio.run(go()), (False, 1))


if __name__ == "__main__":
    unittest.main()


class CrashRetryTest(unittest.TestCase):
    def test_a_batch_that_died_in_gatk_reruns_its_head_alone(self):
        r = ranges(); d = r.defaults("gatk").params
        calls = []

        class BatchCrashes(FakeRunner):                    # the batch path dies; a config on its own passes
            def run_batch(self, round, configs, fidelity=FULL):
                calls.append(("batch", len(configs)))
                return [CallSet(round.round_id, c.key, "", "", 0.1, ok=False, fidelity=fidelity,
                                error="rc=3 java.lang.IllegalStateException: Graph must have ref source and sink vertices")
                        for c in configs]

            def run(self, round, config, fidelity=FULL):
                calls.append(("solo", 1))
                return super().run(round, config, fidelity)

        ex = LocalExecutor(BatchCrashes(r, batching=True), workers=1)
        conf = "standard_min_confidence_threshold_for_calling"
        jobs = [EvalJob(f"j{i}", ROUND, Config("gatk", {**d, conf: c})) for i, c in enumerate((30.0, 40.0, 50.0))]
        ex.submit(jobs)                                    # one family: one batch
        got = {}
        for _ in range(30):
            for cs in ex.poll(0.2):
                got[cs.job_id] = cs
            if len(got) == 3:
                break
        self.assertEqual(calls, [("batch", 3), ("solo", 1)])
        self.assertTrue(got["j0"].ok)                                        # the head, rerun alone
        self.assertEqual([got[j].ok for j in ("j1", "j2")], [False, False])  # the riders fail with the batch
        self.assertIn("batch crashed; its head runs alone", got["j1"].error)
        self.assertEqual(ex.stats()["retries"], 1)

        class Always(FakeRunner):                          # a config that crashes on its own is a failure
            def run(self, round, config, fidelity=FULL):
                return CallSet(round.round_id, config.key, "", "", 0.1, ok=False, error="rc=3 boom", fidelity=fidelity)

        ex = LocalExecutor(Always(r), workers=1)
        ex.submit([EvalJob("k", ROUND, Config("gatk", {**d, "min_pruning": 5}))])
        got = []
        for _ in range(30):
            got += ex.poll(0.2)
            if got:
                break
        self.assertEqual([(cs.ok, cs.error) for cs in got], [(False, "rc=3 boom")])


class EndRoundTest(unittest.TestCase):
    def test_local_executor_drops_and_kills_the_rounds_work(self):
        r = ranges(); d = r.defaults("gatk").params
        killed = []

        class Slow(FakeRunner):
            def run(self, round, config, fidelity=FULL):
                time.sleep(0.5)
                return super().run(round, config, fidelity)

            def kill_round(self, round_id):
                killed.append(round_id)

        ex = LocalExecutor(Slow(r), workers=1)
        ex.submit([EvalJob(f"j{i}", ROUND, Config("gatk", {**d, "min_pruning": p})) for i, p in enumerate((3, 4, 5))])
        time.sleep(0.1)                                   # j0 is running, j1/j2 are queued
        ex.end_round(ROUND.round_id)
        self.assertEqual(killed, [ROUND.round_id])
        got = []
        for _ in range(4):
            got += ex.poll(0.4)
        self.assertEqual(got, [])                          # neither the running job nor the queued ones report
        self.assertEqual(ex.idle_slots(), 1)

    def test_cancel_running_kills_the_batch_and_delivers_nothing_for_it(self):
        r = ranges(); d = r.defaults("gatk").params
        killed = []

        class Slow(FakeRunner):
            def run(self, round, config, fidelity=FULL):
                time.sleep(0.5)
                return super().run(round, config, fidelity)

            def kill_configs(self, round_id, keys):
                killed.append((round_id, sorted(keys)))

        ex = LocalExecutor(Slow(r), workers=1)
        jobs = [EvalJob(f"j{i}", ROUND, Config("gatk", {**d, "min_pruning": p})) for i, p in enumerate((3, 4))]
        ex.submit(jobs)
        time.sleep(0.1)                                   # j0 is running, j1 is queued
        ex.cancel_running(["j0"])                         # its result came from elsewhere first
        self.assertEqual(killed, [(ROUND.round_id, [jobs[0].config.key])])
        got = []
        for _ in range(4):
            got += ex.poll(0.4)
        self.assertEqual([cs.job_id for cs in got], ["j1"])   # j0 never reports; j1 runs as usual
