"""The trade-off sweep (search/engine.py Context.send_sweep): cheap riders over the genotyping and PairHMM settings
around the round's newest best, off unless climb.tradeoff_sweep is 1."""
import asyncio
import random
import unittest
import unittest.mock

from common.tuning import TUNING
from search import engine
from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS, Space
from tests.helpers import ranges
from tests.test_islands import FakeBase

SWEPT = set(GENOTYPING_PARAMS) | set(PAIRHMM_PARAMS)


class SweepTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        self.knobs = unittest.mock.patch.dict(TUNING.values, {"climb.tradeoff_sweep": 1, "dups.enabled": 0})
        self.knobs.start()
        self.base = FakeBase()
        self.ctx = engine.Context(self.base, self.space, random.Random(0), None)
        self.ctx.phase = "B: climbing"
        self.best = self.space._assemble({"assembly_region_padding": 50})
        self.ctx.propose([self.best], "G", ["big"])
        self.base.finish(self.best.key, 85.0)                    # the round's best so far

    def tearDown(self):
        self.knobs.stop()
        self.loop.close()

    def sent(self):
        return [self.base.configs[k] for k in self.base.pending if k != self.best.key]

    def test_riders_change_one_genotyping_or_pairhmm_setting_of_the_best(self):
        self.ctx.send_sweep()
        sent = self.sent()
        self.assertTrue(0 < len(sent) <= engine.SWEEP_MAX)
        for c in sent:
            diff = [n for n in self.space.names if c.params.get(n) != self.best.params.get(n)]
            self.assertEqual(len(diff), 1, diff)
            self.assertIn(diff[0], SWEPT)
            self.assertIn(c.key, self.ctx.ga_keys)              # a higher one can move a climber
        self.assertIn("standard_min_confidence_threshold_for_calling",
                      {n for c in sent for n in self.space.names if c.params.get(n) != self.best.params.get(n)})

    def test_off_by_default_and_outside_steps_b_and_c(self):
        with unittest.mock.patch.dict(TUNING.values, {"climb.tradeoff_sweep": 0}):
            self.ctx.send_sweep()
        self.assertEqual(self.sent(), [])
        self.ctx.phase = "A: wave 0"
        self.ctx.send_sweep()
        self.assertEqual(self.sent(), [])

    def test_one_sweep_per_best_and_at_most_one_a_minute(self):
        self.ctx.send_sweep()
        self.ctx.send_sweep()                                    # same best: nothing new
        self.assertEqual(self.ctx._swept, {self.best.key})
        higher = self.space._assemble({"assembly_region_padding": 75})
        self.ctx.propose([higher], "G", ["big"])
        self.base.finish(higher.key, 86.0)
        self.ctx.send_sweep()                                    # a new best, but within the minute
        self.assertEqual(self.ctx._swept, {self.best.key})
        self.base.now += engine.SWEEP_EVERY_S
        self.ctx.send_sweep()
        self.assertEqual(self.ctx._swept, {self.best.key, higher.key})

    def test_nothing_already_scored_or_running_is_sent_again(self):
        first = self.ctx.sweep_riders(self.best)
        self.ctx.propose(first[:5], "G", ["big"] * 5)
        again = {c.key for c in self.ctx.sweep_riders(self.best)}
        self.assertFalse(again & {c.key for c in first[:5]})


if __name__ == "__main__":
    unittest.main()
