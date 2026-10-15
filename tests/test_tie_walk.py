"""search/climb.py tie walk (v9.18.3): a grid that gains nothing walks each line that tied the centre exactly, every
untried rung out to climb.tie_walk rungs, in one more wave; a rung that scored differently ends its line."""
import asyncio
import unittest
from unittest import mock

from search.climb import Climber
from search.space import Space
from tests.helpers import ranges
from tests.test_engine import FakeCtx

MP = "min_pruning"                        # ladder 2..10: index = value - 2


def landscape(c):
    """A flat from min_pruning 5 to 7 (95.0), the exit at 8 (96.0), lower beyond and below; nothing else matters."""
    v = c.params[MP]
    if 5 <= v <= 7:
        return 95.0
    if v == 8:
        return 96.0
    return 94.0


class TieWalkTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.ctx = FakeCtx(self.space, landscape)
        self.ctx.cap = 1                  # no line points beside the compass: only the compass and the tie walk
        self.centre = self.space._assemble({MP: 5})
        self.cl = Climber(self.ctx, 0, self.centre, 95.0)
        self.loop = asyncio.new_event_loop()          # FakeCtx.propose needs a current loop; an earlier asyncio.run
        asyncio.set_event_loop(self.loop)             # in the suite leaves none

    def tearDown(self):
        asyncio.set_event_loop(None)
        self.loop.close()

    def one(self, v):
        return self.space._assemble({**self.centre.params, MP: v})

    def test_a_tied_line_walks_every_untried_rung_out_to_the_reach(self):
        results = [(self.one(6), 95.0), (self.one(4), 94.0)]
        walk = self.cl.tie_walk_heads(self.centre, [MP], results, 95.0)
        self.assertEqual([c.params[MP] for c in walk], [7, 8, 9])        # 6 tied already; 4 rungs out = 9

    def test_a_losing_side_is_not_walked(self):
        results = [(self.one(6), 95.0), (self.one(4), 94.0)]
        walk = self.cl.tie_walk_heads(self.centre, [MP], results, 95.0)
        self.assertFalse([c for c in walk if c.params[MP] < 5])

    def test_a_rung_scored_differently_ends_the_line(self):
        self.ctx.sent = [[self.one(7)]]                                  # 7 already scored (the landscape says 95)
        with mock.patch.object(self.ctx, "f", lambda c: 93.0 if c.params[MP] == 7 else landscape(c)):
            walk = self.cl.tie_walk_heads(self.centre, [MP], [(self.one(6), 95.0)], 95.0)
        self.assertEqual(walk, [])

    def test_a_grid_with_a_gain_does_not_walk(self):
        results = [(self.one(6), 95.0), (self.one(4), 95.5)]
        self.assertEqual(self.cl.tie_walk_heads(self.centre, [MP], results, 95.0), [])

    def test_switched_off_it_does_not_walk(self):
        with mock.patch("search.climb.tune", lambda name: 0 if name == "climb.tie_walk" else __import__(
                "common.tuning", fromlist=["tune"]).tune(name)):
            self.assertEqual(self.cl.tie_walk_heads(self.centre, [MP], [(self.one(6), 95.0)], 95.0), [])

    def test_the_interval_side_goes_first(self):
        # a second tied line on another setting, on the side away from its interval: it goes after MP's line
        other = "min_dangling_branch_length"
        down = self.space._assemble({**self.centre.params, other: self.centre.params[other] - 1})
        results = [(self.one(6), 95.0), (down, 95.0)]
        with mock.patch("search.plateau.intervals", return_value={other: {len(self.ctx.ladders.rungs[other]) - 1}}):
            walk = self.cl.tie_walk_heads(self.centre, [MP, other], results, 95.0)
        names = [MP if c.params[MP] != 5 else other for c in walk]
        self.assertEqual(names[0], MP)
        self.assertEqual(names[-1], other)

    def live_grid(self):
        """One grid the way the live engine runs it (begin_grid, send every head, close_wave until the wave ends)."""
        cl = self.cl
        cl.next_group = lambda: ("assembly", [MP])
        self.assertTrue(cl.begin_grid(1e9))
        for _ in range(10):
            while cl.queued():
                cl.submit_one()
            self.assertTrue(cl.gate_done())
            cl.close_wave()
            if not cl.batch_live:
                return
        self.fail("the grid never closed")

    def test_the_live_grid_crosses_the_flat_to_the_exit(self):
        self.live_grid()
        self.assertEqual(self.cl.centre.params[MP], 8)
        self.assertEqual(self.cl.score, 96.0)

    def test_without_the_walk_the_live_grid_stays_on_the_flat(self):
        with mock.patch("search.climb.tune", lambda name: 0 if name == "climb.tie_walk" else __import__(
                "common.tuning", fromlist=["tune"]).tune(name)):
            self.live_grid()
        self.assertEqual(self.cl.score, 95.0)
        self.assertIn(MP, self.cl.flat)


if __name__ == "__main__":
    unittest.main()
