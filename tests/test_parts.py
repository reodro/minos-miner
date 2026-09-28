import os
import unittest
from dataclasses import replace

from answer_set import OracleProvider
from backtest.kpi import oracle_gap, rank_in_field, would_win
from common.settings import SETTINGS
from common.types import AnswerSet, CallSet, Config, Regime, Round
from evaluator import FakeEvaluator
from main import DryRunSubmitter, LoopSettings, RandomSearcher, run_round
from tests.helpers import ranges

ROUND_0817 = os.path.join(SETTINGS.rounds_dir, "20260928T081700Z")


class StubProvider:
    """Answer set that upgrades itself on every update() call."""

    def __init__(self):
        self.updates = 0

    def build(self, round):
        return AnswerSet(round.round_id, 1, Regime.UNKNOWN, None, [], [], 600)

    def update(self, answer_set, new_callsets):
        self.updates += 1
        return replace(answer_set, version=answer_set.version + 1)


class LoopTest(unittest.TestCase):
    def setUp(self):
        self.r = ranges()
        self.round = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")

    def run_loop(self, evaluator, provider=None, max_evals=40, baselines=None):
        sub = DryRunSubmitter()
        out = run_round(self.round, provider or StubProvider(), evaluator, RandomSearcher(self.r, seed=1), sub,
                        self.r, LoopSettings(workers=4, update_every=10, max_evals=max_evals),
                        baselines=baselines if baselines is not None else [self.r.defaults("gatk")])
        return out, sub

    def test_baseline_submitted_first_and_scores_only_improve(self):
        out, sub = self.run_loop(FakeEvaluator(self.r))
        self.assertEqual(sub.sent[0][1], self.r.defaults("gatk").key)
        self.assertEqual(out.evaluated, 40)
        scores = [s for _, _, s in out.submissions]
        self.assertEqual(scores, sorted(scores))
        self.assertEqual(out.best.config_key, sub.sent[-1][1])

    def test_answer_set_upgrades_trigger_rescoring(self):
        provider = StubProvider()
        out, _ = self.run_loop(FakeEvaluator(self.r), provider)
        self.assertGreaterEqual(provider.updates, 3)
        self.assertEqual(out.best.answer_set_version, out.answer_set_version)

    def test_crashing_runs_do_not_end_the_round(self):
        class Flaky(FakeEvaluator):
            n = 0

            def run(self, round, config):
                Flaky.n += 1
                if Flaky.n % 3 == 0:
                    raise RuntimeError("Graph must have ref source and sink vertices")
                return super().run(round, config)

        out, sub = self.run_loop(Flaky(self.r))
        self.assertGreater(out.failed, 0)
        self.assertEqual(out.evaluated + out.failed, 40)
        self.assertTrue(sub.sent)

    def test_crashing_baseline_falls_through_to_next(self):
        class CrashFirst(FakeEvaluator):
            def run(self, round, config):
                if config.params.get("min_pruning") == 3:
                    raise RuntimeError("boom")
                return super().run(round, config)

        first = Config("gatk", {"min_pruning": 3})
        out, sub = self.run_loop(CrashFirst(self.r), baselines=[first, self.r.defaults("gatk")])
        self.assertEqual(sub.sent[0][1], self.r.defaults("gatk").key)

    def test_invalid_baseline_is_skipped(self):
        bad = Config("gatk", {"min_pruning": 99})
        out, sub = self.run_loop(FakeEvaluator(self.r), baselines=[bad, self.r.defaults("gatk")])
        self.assertEqual(sub.sent[0][1], self.r.defaults("gatk").key)


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
