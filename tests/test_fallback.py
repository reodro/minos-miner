"""The unscored baseline fallback (main.loop): a round in which nothing is scored still gets a submission."""
import unittest

from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search import engine
from tests.helpers import ranges
from tests.test_engine import ROUND, Stub


class FallbackTest(unittest.TestCase):
    def round(self, ex, runner, r, **settings):
        sub = DryRunSubmitter()
        with self.assertLogs("main.loop", level="INFO") as logs:
            run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), sub, r,
                      LoopSettings(stop_before_deadline_s=60, **settings), baselines=[r.defaults("gatk")],
                      deadline_ts=72 * 60)
        return sub, "\n".join(logs.output)

    def test_nothing_scored_submits_the_baseline(self):
        r = ranges()
        runner = FakeRunner(r, seed=0)
        sub, text = self.round(SimExecutor(runner, 16, crashes=lambda job: True), runner, r)   # every run fails
        self.assertEqual(sub.sent, [(ROUND.round_id, r.defaults("gatk").key)])
        self.assertIn("submitted unscored", text)

    def test_a_scored_best_replaces_the_fallback(self):
        r = ranges()
        runner = FakeRunner(r, seed=0)
        base = r.defaults("gatk").key
        # the baseline run fails at once (fallback due: nothing submitted), the other configs score 10 minutes later
        ex = SimExecutor(runner, 16, duration=lambda job: 1.0 if job.config.key == base else 600.0,
                         crashes=lambda job: job.config.key == base, batching=False)
        sub, text = self.round(ex, runner, r, fallback_before_deadline_s=72 * 60)
        self.assertIn("submitted unscored", text)
        self.assertEqual(sub.sent[0][1], base)
        self.assertGreaterEqual(len(sub.sent), 2)
        self.assertNotEqual(sub.sent[-1][1], base)

    def test_a_normal_round_never_falls_back(self):
        r = ranges()
        runner = FakeRunner(r, seed=0)
        sub, text = self.round(SimExecutor(runner, 16), runner, r)
        self.assertNotIn("submitted unscored", text)
        self.assertTrue(sub.sent)


class BaselinesTest(unittest.TestCase):
    def test_generalist_first_then_defaults(self):
        import unittest.mock
        from common.tuning import TUNING
        from main.__main__ import baselines
        r = ranges()
        with unittest.mock.patch.dict(TUNING.values, {"main.generalist_first": 1}):
            b = baselines(r)
        self.assertEqual(b[-1].key, r.defaults("gatk").key)
        if len(b) == 2:                                    # g0 present and in range on this machine
            self.assertNotEqual(b[0].key, b[1].key)
            self.assertEqual(r.errors(b[0]), [])
        with unittest.mock.patch.dict(TUNING.values, {"main.generalist_first": 0}):
            self.assertEqual([c.key for c in baselines(r)], [r.defaults("gatk").key])


if __name__ == "__main__":
    unittest.main()
