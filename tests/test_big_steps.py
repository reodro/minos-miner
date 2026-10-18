"""climb.big_steps (search/climb.py step_plan, Climber.step_size / adapt_reach) and search.leaps.summary (2026-10-18)."""
import random
import unittest
import unittest.mock

from common.tuning import TUNING
from search import engine, leaps
from search.climb import Climber, step_plan
from search.space import Space
from tests.helpers import ranges
from tests.test_islands import FakeBase, JumpCtx

MMQ = "min_mapping_quality_score"


class StepPlanTest(unittest.TestCase):
    def test_three_plans_by_ladder_length(self):
        self.assertEqual(step_plan(30), ((1, 2, 3, 5), (1,)))
        self.assertEqual(step_plan(61), ((3, 5, 7, 9), (3, 1)))
        self.assertEqual(step_plan(601), ((3, 9, 18, 27), (9, 3, 1)))
        with unittest.mock.patch.dict(TUNING.values, {"climb.big_step_large_max": 40}):
            self.assertEqual(step_plan(601)[0][-1], 40)


class ReachTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.ctx = JumpCtx(self.space, lambda c: 80.0)
        self.centre = self.space._assemble({MMQ: 30})
        self.cl = Climber(self.ctx, 0, self.centre, 80.0)
        lad = self.ctx.ladders.rungs
        self.big = max(self.space.names, key=lambda n: len(lad.get(n) or ()))           # > 150 rungs
        self.mid = MMQ                                                                    # 61 rungs

    def steps_without_moving(self, name, n):
        out = [self.cl.step_size(name)]
        for _ in range(n):
            self.cl.adapt_reach([name], self.centre, None)
            out.append(self.cl.step_size(name))
        return out

    def test_a_long_ladder_widens_to_its_maximum_then_narrows_to_one(self):
        self.assertEqual(self.steps_without_moving(self.big, 8), [3, 9, 18, 27, 9, 3, 1, 1, 1])

    def test_a_medium_ladder_widens_to_nine_then_narrows(self):
        self.assertEqual(self.steps_without_moving(self.mid, 7), [3, 5, 7, 9, 3, 1, 1, 1])

    def test_offsets_follow_the_step_and_a_moved_setting_keeps_it(self):
        self.assertEqual(self.cl.offsets(self.mid), [0, 3, -3])
        self.cl.adapt_reach([self.mid], self.centre, None)                                # 3 -> 5
        moved = self.space._assemble({MMQ: 35})
        self.cl.adapt_reach([self.mid], self.centre, moved)                               # moved: stays 5
        self.assertEqual(self.cl.step_size(self.mid), 5)

    def test_off_means_one_rung(self):
        with unittest.mock.patch.dict(TUNING.values, {"climb.big_steps": 0}):
            self.assertEqual(self.cl.offsets(self.mid), [0, 1, -1])

    def test_a_new_centre_starts_the_steps_again(self):
        self.cl.adapt_reach([self.mid], self.centre, None)
        self.cl.reach.clear()
        self.assertEqual(self.cl.step_size(self.mid), 3)


class LeapSummaryTest(unittest.TestCase):
    def test_counts_tries_wins_and_two_setting_jumps(self):
        ctx = engine.Context(FakeBase(), Space(ranges()), random.Random(0), None)
        arms = leaps.arms_of(ctx)
        arms.learn((MMQ, 1, -1), 0.02, 0.005)
        arms.learn((MMQ, 1, -1), 0.0, 0.005)
        ctx._leap2 = [3, 1, 0.01]
        s = leaps.summary(ctx)
        self.assertEqual((s["tries"], s["wins"]), (2, 1))
        self.assertEqual(s["by_setting"][MMQ]["tries"], 2)
        self.assertEqual(s["two_setting"], {"tries": 3, "wins": 1, "gain": 0.01})


if __name__ == "__main__":
    unittest.main()
