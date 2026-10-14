"""search/excess.py: the core-above-curve measure (band-edge rule, lane results excluded) and an excess lane (switched
on here; off by default) running inside a whole virtual round."""
import unittest
from types import SimpleNamespace

from common.types import FULL
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search import engine, excess
from tests.helpers import ranges
from tests.test_engine import ROUND, Stub


def _r(score, core, fp):
    return SimpleNamespace(score=score, core=core, region_fp=fp, gate_passed=True)


class _Client:
    def __init__(self, rows):
        self.latest = {(k, FULL): r for k, r in rows.items()}

    def best(self, fid=FULL):
        return max(self.latest.values(), key=lambda r: r.score)


class _Ctx:
    def __init__(self, rows):
        self.client = _Client(rows)
        self.said = []

    def say(self, msg, *args):
        self.said.append(msg % args)


class CurveTest(unittest.TestCase):
    def setUp(self):
        # the best score sits at FP 1000; the curve rises with FP
        self.ctx = _Ctx({"best": _r(82.0, 0.860, 1000), "lo": _r(80.0, 0.850, 900),
                         "hi": _r(81.0, 0.866, 1150), "far": _r(70.0, 0.600, 600)})
        self.cv = excess.Curve(self.ctx)
        self.cv.refresh(force=True)

    def test_inside_the_band_it_is_core_above_the_curve_at_the_same_fp(self):
        self.assertAlmostEqual(self.cv.j(0.862, 1000), 0.002, places=6)
        self.assertAlmostEqual(self.cv.j(0.860, 1000), 0.0, places=6)

    def test_below_the_band_it_is_compared_at_the_band_edge_not_the_thin_curve(self):
        # at FP 700 the curve holds only "far" (0.600): a naive comparison would give +0.25
        self.assertAlmostEqual(self.cv.j(0.850, 700), 0.850 - 0.850, places=6)     # compared at 0.95 x 1000 = 950
        self.assertLess(self.cv.j(0.700, 700), -0.1)                              # trading core for FP is not rewarded

    def test_above_the_band_it_pays_the_scores_own_price_of_fp(self):
        # band_high 1.05 x 1000 = 1050: compared with the curve there (0.860), and each % of FP above it costs what
        # the score charges for it in core (N 548 when the answer set gives none): 20 % here
        import math
        per_fp = 30.0 / (8 * 548) * math.exp(-1050 / (8 * 548))
        price = per_fp * 1000 / 100 / 70
        self.assertAlmostEqual(self.cv.over_cost(1050), price, places=9)
        self.assertAlmostEqual(self.cv.j(0.870, 1250), 0.870 - 0.860 - 20 * price, places=6)
        # the price is the score's exchange rate: same score change from 1 % more FP or `price` less core
        s = lambda core, fp: 100 * (0.70 * core + 0.30 * math.exp(-fp / (8 * 548)))
        self.assertAlmostEqual(s(0.86, 1050) - s(0.86, 1060), s(0.86, 1050) - s(0.86 - price, 1050), places=3)

    def test_inside_the_band_fp_is_free_but_only_within_five_percent(self):
        self.assertAlmostEqual(self.cv.j(0.862, 1040), 0.862 - 0.860, places=6)
        self.assertLess(self.cv.j(0.866, 1150), self.cv.j(0.866, 1050))         # the old band's edge now pays

    def test_lane_results_are_not_part_of_the_curve(self):
        self.ctx.client.latest[("lane", FULL)] = _r(81.5, 0.880, 1000)
        self.cv.lane_keys.add("lane")
        self.cv.refresh(force=True)
        self.assertAlmostEqual(self.cv.of(SimpleNamespace(key="lane")), 0.880 - 0.860, places=6)


class LaneOffByDefaultTest(unittest.TestCase):
    def test_the_lane_is_off_by_default(self):
        from common.tuning import tune
        self.assertEqual(int(tune("excess.lanes")), 0)


class LaneInRoundTest(unittest.TestCase):
    def test_an_excess_lane_runs_inside_a_whole_round(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        searcher = engine.searcher(r, 300)
        from unittest import mock
        from common.tuning import tune as real_tune
        on = lambda name: 1 if name == "excess.lanes" else real_tune(name)          # the lane is off by default
        with mock.patch("search.excess.tune", on), self.assertLogs("search.engine", level="INFO") as logs:
            out = run_round(ROUND, Stub(), ex, FakeScorer(runner), searcher, DryRunSubmitter(), r,
                            LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=72 * 60)
        lane = [l for l in logs.output if "excess lane" in l]
        self.assertTrue(lane, "no excess lane message")
        self.assertTrue(any("starts at" in l or "not started" in l for l in lane))
        self.assertIsNotNone(out.best)


if __name__ == "__main__":
    unittest.main()
