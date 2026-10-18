"""climb.after_peak (2026-10-18): a peaked climber keeps sending centre jumps, and kicks after three dry jump waves."""
import unittest
import unittest.mock

from common.tuning import TUNING
from search import climb
from search.climb import Climber
from search.space import Space
from tests.helpers import ranges
from tests.test_islands import JumpCtx

MMQ = "min_mapping_quality_score"


class AfterPeakTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.ctx = JumpCtx(self.space, lambda c: 80.0)
        self.cl = Climber(self.ctx, 0, self.space._assemble({MMQ: 30}), 82.0)
        self.until = self.ctx.clock() + 3600.0
        self.jump = self.space._assemble({MMQ: 45})

    def begin(self, leaps):
        """begin_jumps with no untried neighbour: the certificate (n None) is empty, a numbered pick returns `leaps`."""
        def plan(centre, n, taken=()):
            return [] if n is None else list(leaps)
        with unittest.mock.patch.object(climb, "one_rung_neighbors", return_value=[]), \
                unittest.mock.patch.object(self.cl, "plan_leaps", side_effect=plan):
            return self.cl.begin_jumps(self.until)

    def test_off_names_the_peak_and_stops(self):
        with unittest.mock.patch.dict(TUNING.values, {"climb.after_peak": 0}):
            self.assertFalse(self.begin([self.jump]))
        self.assertEqual(self.cl.state, "peak")

    def test_on_keeps_jumping_after_the_peak(self):
        self.assertTrue(self.begin([self.jump]))
        self.assertEqual(self.cl.state, "climbing")
        self.assertTrue(self.cl.after_peak_wave)
        self.assertFalse(self.cl.kicking)
        self.assertEqual(self.cl.kind, "jumps")

    def test_three_dry_waves_kick(self):
        self.cl.dry_jumps = climb.KICK_AFTER
        self.assertTrue(self.begin([self.jump]))
        self.assertTrue(self.cl.kicking)

    def test_no_kick_late_in_the_climb(self):
        self.cl.dry_jumps = climb.KICK_AFTER
        self.until = self.ctx.clock() + 60.0 * climb.KICK_MIN_LEFT_MIN - 1
        self.assertTrue(self.begin([self.jump]))
        self.assertFalse(self.cl.kicking)
        self.assertTrue(self.cl.after_peak_wave)

    def test_kick_configs_move_four_to_six_settings_inside_the_ranges(self):
        out = self.cl.kick_configs(8)
        self.assertGreaterEqual(len(out), 4)
        for c in out:
            moved = [n for n in self.space.names if c.params.get(n) != self.cl.centre.params.get(n)]
            self.assertTrue(1 <= len(moved) <= 6, moved)            # repairs can merge a move, never add one
            self.assertIsNotNone(self.space._assemble(dict(c.params)))


if __name__ == "__main__":
    unittest.main()
