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
from search import EvalClient, RandomSampler, random_stream
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
        self.assertEqual((out.evaluated, out.failed), (0, 12))
        self.assertEqual(sub.sent, [(ROUND.round_id, ranges().defaults("gatk").key)])   # nothing scored: the baseline
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
        runner.variant_cost = lambda main, variant, fidelity: 0.7   # 132 + 92 s; a third passes 2.3 runs (304 s)
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
        from search.axes import latin_hypercube
        own += latin_hypercube(sp, rng, 64)
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
        got = [capped.random(rng) for _ in range(100)] + latin_hypercube(capped, rng, 16)
        self.assertLessEqual(max(c.params["min_base_quality_score"] for c in got if c), 30)
        self.assertEqual(capped._assemble({"min_base_quality_score": 45}).params["min_base_quality_score"], 30)
        from search.space import bam_caps
        self.assertEqual(bam_caps("/nonexistent.bam", "chr20", 1, 10), {})

    @unittest.skip("built on validator-edge values outside the published ranges, removed 2026-10-18 (search/space.py VALIDATOR_EDGES)")

    def test_validator_edges_assemble_and_stay_off_the_dense_ladder(self):
        from search.space import Space
        from search.steps import Ladders
        from search.candidates import in_indel_family
        sp = Space(ranges())
        edge = sp._assemble({"min_pruning": 1, "active_probability_threshold": 0.0001,
                             "max_num_haplotypes_in_population": 256, "dont_use_soft_clipped_bases": True})
        self.assertIsNotNone(edge)
        self.assertEqual(edge.params["min_pruning"], 1)
        self.assertEqual(edge.params["max_num_haplotypes_in_population"], 256)
        self.assertTrue(in_indel_family(edge))
        ladders = Ladders(sp)
        self.assertIn(1, ladders.rungs["min_pruning"])
        self.assertIn(256, ladders.rungs["max_num_haplotypes_in_population"])
        self.assertNotIn(129, ladders.rungs["max_num_haplotypes_in_population"])
        self.assertIn(0.0001, ladders.rungs["active_probability_threshold"])
        self.assertTrue(edge.params["dont_use_soft_clipped_bases"])
        self.assertNotIn("sample_ploidy", sp.names)
        self.assertNotIn("dont_use_soft_clipped_bases", sp.names)

    def test_a_round_bam_replaces_three_ranges_with_the_groups_it_can_tell_apart(self):
        import random as _r
        from search.bamgroups import BQST, MBQ, MRPAS, BamFacts, bam_bounds, bounds_from
        from search.space import Space
        # 20261006T010500Z: eight base qualities, 74 reads at the busiest start.
        facts = BamFacts(frozenset({2, 3, 4, 5, 6, 10, 20, 30}), frozenset(range(61)), 74)
        bounds = bounds_from(facts, ranges())
        self.assertEqual(bounds.groups[MBQ], (10, 11, 20, 21, 30))
        self.assertEqual(bounds.count[BQST], 26)          # finest list, at mapping quality 0
        self.assertEqual(bounds.count[MRPAS], 50)         # 25..74; everything above 74 is 74
        self.assertIn("1/78", bounds.report)
        sp = Space(ranges(), bounds=bounds)
        rng = _r.Random(0)
        got = [c for c in (sp.random(rng) for _ in range(80)) if c]
        self.assertTrue(got)
        self.assertLessEqual(set(c.params[MBQ] for c in got), set(bounds.groups[MBQ]))
        self.assertLessEqual(set(c.params[BQST] for c in got), set(bounds.groups[BQST]))
        self.assertLessEqual(max(c.params[MRPAS] for c in got), 74)
        snapped = sp._assemble({MBQ: 16, BQST: 13, "min_mapping_quality_score": 59, MRPAS: 200})
        self.assertEqual(snapped.params[MBQ], 11)
        self.assertEqual(snapped.params[BQST], 11)        # at mapping quality 59 the 26 collapse toward 5
        self.assertEqual(snapped.params[MRPAS], 74)
        self.assertEqual(sp._assemble({MBQ: 45}).params[MBQ], 30)
        nothing = bounds_from(BamFacts(frozenset({2, 10, 20, 30}), frozenset(range(61)), 13), ranges())
        fixed = Space(ranges(), bounds=nothing)
        self.assertNotIn(MRPAS, fixed.names)              # no start has 25 reads, so the axis is one value
        self.assertEqual(fixed._assemble({MRPAS: 116}).params[MRPAS], 25)
        self.assertIsNone(bam_bounds("/nonexistent.bam", ranges()))

    @unittest.skip("built on validator-edge values outside the published ranges, removed 2026-10-18 (search/space.py VALIDATOR_EDGES)")

    def test_contamination_rungs_that_remove_the_same_number_of_reads_are_one_rung(self):
        from search.bamgroups import CONT, BamFacts, bounds_from, contamination_groups, snap
        from search.space import Space
        from search.steps import Ladders
        depths = frozenset(i for i in range(1, 72) if i not in (64, 65))
        facts = BamFacts(frozenset({2, 10, 20, 30}), frozenset(range(61)), 13, depths)
        reps = contamination_groups(depths, 0.0, 0.05)
        self.assertEqual(reps[0], 0.0)
        self.assertNotIn(0.0025, reps)
        self.assertNotIn(0.0125, reps)
        self.assertIn(0.015, reps)
        self.assertIn(0.045, reps)
        self.assertIn(0.0475, reps)
        self.assertLess(len(reps), 21)
        bounds = bounds_from(facts, ranges())
        self.assertNotIn(CONT, bounds.groups)                 # allowed would encode an off-ladder fraction as rung 0
        self.assertIn(f"{CONT} 21 -> {len(reps)}", bounds.report)
        from search.space import EDGE_RUNGS
        sp = Space(ranges(), bounds=bounds)
        dense = [r for r in Ladders(sp).rungs[CONT] if r <= 0.05 + 1e-12]
        self.assertEqual(dense, list(reps))
        self.assertTrue(set(EDGE_RUNGS[CONT]) <= set(Ladders(sp).rungs[CONT]))
        self.assertEqual(sp._assemble({CONT: 0.01}).params[CONT], 0.0)
        self.assertEqual(sp._assemble({CONT: 0.01875}).params[CONT], 0.01875)
        self.assertEqual(snap(facts, CONT, 0.005, mmq=0, lo=0.0, hi=0.05), 0.0)
        shallow = BamFacts(frozenset({2, 10, 20, 30}), frozenset(range(61)), 13, frozenset(range(1, 11)))
        one = bounds_from(shallow, ranges())
        fixed = Space(ranges(), bounds=one)
        self.assertEqual([r for r in Ladders(fixed).rungs[CONT] if r <= 0.05 + 1e-12], [0.0])
        self.assertEqual(fixed._assemble({CONT: 0.04}).params[CONT], 0.0)

    def test_snapping_keeps_every_class_the_bam_can_tell_apart(self):
        """A representative has to be the same GATK run as the value it replaces, and every coarser class has to
        stay reachable from the finest list. Checked against the executor's own canonical() on four BAM shapes."""
        from evaluator.equivalence import RoundFacts, canonical
        from search.bamgroups import BQST, MBQ, MMQ, MRPAS, BamFacts, bounds_from, bqst_groups, bqst_values, snap
        shapes = (
            BamFacts(frozenset({2, 3, 4, 5, 6, 10, 20, 30}), frozenset(range(61)), 74),
            BamFacts(frozenset(range(1, 41)), frozenset({0, 10, 20, 60}), 13),
            BamFacts(frozenset({4, 8, 15}), frozenset(range(61)), 33),
            BamFacts(frozenset(range(0, 50)), frozenset(range(61)), 400),
        )
        specs = ranges().tools["gatk"]
        mbq_lo, mbq_hi = int(specs[MBQ].min), int(specs[MBQ].max)
        bq_lo, bq_hi = int(specs[BQST].min), int(specs[BQST].max)
        mr_lo, mr_hi = int(specs[MRPAS].min), int(specs[MRPAS].max)
        mmq_lo, mmq_hi = int(specs[MMQ].min), int(specs[MMQ].max)

        def same_closed(a, b, values):
            if a == b:
                return True
            lo, hi = (a, b) if a < b else (b, a)
            return not any(lo <= q <= hi for q in values)

        def same_half(a, b, values):
            if a == b:
                return True
            lo, hi = (a, b) if a < b else (b, a)
            return not any(lo <= q < hi for q in values)

        for facts in shapes:
            bounds = bounds_from(facts, ranges())
            cap = bounds.caps.get(MBQ, mbq_hi)
            self.assertFalse(any(q > cap for q in facts.quals))
            compared = facts.compared_quals
            for v in range(mbq_lo, cap + 1):
                got = snap(facts, MBQ, v, mmq=0, lo=mbq_lo, hi=cap)
                self.assertTrue(same_closed(got, v, compared), (v, got))
                self.assertIn(got, bounds.groups.get(MBQ, (got,)))
            fine = bounds.groups.get(BQST) or bqst_groups(facts, mmq_lo, bq_lo, bq_hi)
            for mmq in range(mmq_lo, mmq_hi + 1):
                capped = bqst_values(facts, mmq)
                coarse = bqst_groups(facts, mmq, bq_lo, bq_hi)
                reached = set()
                for v in range(6, bq_hi + 1):
                    got = snap(facts, BQST, v, mmq=mmq, lo=bq_lo, hi=bq_hi)
                    self.assertTrue(same_half(got, v, capped), (mmq, v, got))
                    reached.add(got)
                for rep in coarse:
                    self.assertTrue(any(same_half(rep, g, capped) for g in reached), (mmq, rep))
                for rep in fine:
                    got = snap(facts, BQST, rep, mmq=mmq, lo=bq_lo, hi=bq_hi)
                    self.assertIn(got, fine)
                    ask = {MBQ: bounds.groups.get(MBQ, (mbq_lo,))[0], BQST: rep, MMQ: mmq, MRPAS: mr_lo}
                    out = canonical(Config("gatk", ask), RoundFacts(facts.quals, facts.mapqs, facts.max_reads_per_start))
                    self.assertTrue(same_half(out.params[BQST], got, capped), (mmq, rep, got, out.params[BQST]))
                    self.assertTrue(same_closed(out.params[MBQ], ask[MBQ], compared))
            start_cap = facts.max_reads_per_start
            for v in range(mr_lo, mr_hi + 1):
                got = snap(facts, MRPAS, v, mmq=0, lo=mr_lo, hi=mr_hi)
                same = (got >= start_cap and v >= start_cap) or got == min(max(v, mr_lo), mr_hi)
                self.assertTrue(same, (v, got, start_cap))


class LineupTest(unittest.TestCase):
    def test_search_hands_main_one_searcher(self):
        from search import searcher
        self.assertTrue(callable(searcher(ranges(), 8)))       # the engine; main chooses nothing


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

    def test_clearing_one_round_leaves_the_others(self):
        import os
        import tempfile
        from common.store import clear_round
        d = tempfile.mkdtemp()
        tag = "20261011T141800Z"
        kept = "20261010T141800Z"
        os.close(os.open(os.path.join(d, tag + ".sqlite"), os.O_CREAT | os.O_WRONLY, 0o644))
        os.close(os.open(os.path.join(d, tag + ".sqlite-wal"), os.O_CREAT | os.O_WRONLY, 0o644))
        os.close(os.open(os.path.join(d, kept + ".sqlite"), os.O_CREAT | os.O_WRONLY, 0o644))
        self.assertEqual(clear_round(d, tag), 2)
        self.assertFalse(os.path.exists(os.path.join(d, tag + ".sqlite")))
        self.assertFalse(os.path.exists(os.path.join(d, tag + ".sqlite-wal")))
        self.assertTrue(os.path.exists(os.path.join(d, kept + ".sqlite")))


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
        # the head reruns alone; the riders are batched again without it (that batch dies too here), then its new
        # head runs alone and the rider lost twice runs alone: every config gets a result of its own
        self.assertEqual(calls, [("batch", 3), ("solo", 1), ("batch", 2), ("solo", 1), ("solo", 1)])
        self.assertTrue(all(got[j].ok for j in ("j0", "j1", "j2")))
        self.assertEqual(ex.stats()["retries"], 5)

        class OneBad(FakeRunner):                          # one rider crashes GATK wherever it runs
            def run_batch(self, round, configs, fidelity=FULL):
                calls2.append(len(configs))
                if any(c.params[conf] == 50.0 for c in configs):
                    return [CallSet(round.round_id, c.key, "", "", 0.1, ok=False, fidelity=fidelity,
                                    error="rc=3 java.lang.IllegalStateException: boom") for c in configs]
                return super().run_batch(round, configs, fidelity)

            def run(self, round, config, fidelity=FULL):
                calls2.append(1)
                if config.params[conf] == 50.0:
                    return CallSet(round.round_id, config.key, "", "", 0.1, ok=False, error="rc=3 boom",
                                   fidelity=fidelity)
                return super().run(round, config, fidelity)

        calls2 = []
        ex = LocalExecutor(OneBad(r, batching=True), workers=1)
        jobs = [EvalJob(f"b{i}", ROUND, Config("gatk", {**d, conf: c})) for i, c in enumerate((30.0, 40.0, 50.0, 60.0))]
        ex.submit(jobs)
        got = {}
        for _ in range(60):
            for cs in ex.poll(0.2):
                got[cs.job_id] = cs
            if len(got) == 4:
                break
        self.assertEqual(sorted(got), ["b0", "b1", "b2", "b3"])                # every job answered, no endless retry
        self.assertEqual([got[j].ok for j in ("b0", "b1", "b2", "b3")], [True, True, False, True])
        self.assertLessEqual(len(calls2), 8)

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
