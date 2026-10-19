"""climb.grid_fill (search/climb.py fill_grids) and RF-UCB in Step B (search/rf_ucb.py step_b, climb.submit_rf),
2026-10-19."""
import asyncio
import random
import unittest
import unittest.mock

from common.tuning import TUNING
from search import climb, rf_ucb
from search.climb import fill_grids, submit_rf
from search.space import Space
from tests.helpers import ranges
from tests.test_islands import JumpCtx

MMQ = "min_mapping_quality_score"
BQ = "min_base_quality_score"


class FakeGridClimber:
    def __init__(self, kind="read", pending=()):
        self.state, self.batch_live, self.kind = "climbing", True, kind
        self.pending, self.extra, self.sent, self.asked = list(pending), [], 0, 0
        self.leaps = []

    def extend(self, n):
        self.asked += n
        self.extra += [object()] * n
        return n

    def submit_one(self):
        self.extra.pop()
        self.sent += 1


class FillGridsTest(unittest.TestCase):
    def test_off_sends_nothing(self):
        cl = FakeGridClimber()
        with unittest.mock.patch.dict(TUNING.values, {"climb.grid_fill": 0.0}):
            self.assertEqual(fill_grids([cl], 40), 0)
        self.assertEqual(cl.asked, 0)

    def test_half_of_the_free_workers_go_to_grid_waves_only(self):
        grid, jumps, busy = FakeGridClimber(), FakeGridClimber(kind="jumps"), FakeGridClimber(pending=[1])
        with unittest.mock.patch.dict(TUNING.values, {"climb.grid_fill": 0.5}):
            self.assertEqual(fill_grids([grid, jumps, busy], 40), 20)
        self.assertEqual((grid.sent, jumps.sent, busy.sent), (20, 0, 0))

    def test_queued_centre_jumps_go_first(self):
        grid, jumper = FakeGridClimber(), FakeGridClimber()
        jumper.leaps = [object()]
        with unittest.mock.patch.dict(TUNING.values, {"climb.grid_fill": 1.0}):
            self.assertEqual(fill_grids([grid, jumper], 10), 0)
        self.assertEqual(grid.asked, 0)

    def test_shared_between_grid_waves(self):
        a, b = FakeGridClimber(), FakeGridClimber()
        with unittest.mock.patch.dict(TUNING.values, {"climb.grid_fill": 1.0}):
            self.assertEqual(fill_grids([a, b], 10), 10)
        self.assertEqual((a.sent, b.sent), (5, 5))


class Keyed:
    def __init__(self, key):
        self.key = key


class QuotaClimber:
    """Queues of heads (grid, or jumps when listed in leap_meta), centre jumps and line points; each send is one run."""

    def __init__(self, heads=0, jump_heads=0, leaps=0, lines=0, kind="read"):
        self.state, self.kind, self.batch_live = "climbing", kind, True
        self.pending = [Keyed(f"h{i}") for i in range(heads)] + [Keyed(f"j{i}") for i in range(jump_heads)]
        self.leap_meta = {f"j{i}": 1 for i in range(jump_heads)}
        self.leaps, self.extra, self.later = [object()] * leaps, [object()] * lines, []
        self.evaluated, self.live = 0, {"grid": 0, "jump": 0}
        self.wave_centre = self.centre = None
        self.refill, self.planned = 0, 0

    def plan_leaps(self, centre, n, taken=()):
        k = min(n, self.refill)
        self.refill -= k
        self.planned += k
        return [object()] * k

    def head_part(self):
        return "jump" if self.pending and self.pending[0].key in self.leap_meta else "grid"

    def runs_in_flight(self):
        return dict(self.live)

    def _send(self, queue, part):
        queue.pop(0)
        self.evaluated += 1
        self.live[part] += 1

    def _submit_head(self):
        self._send(self.pending, self.head_part())

    def _submit_leap(self):
        self._send(self.leaps, "jump")

    def _submit_extra(self):
        self._send(self.extra, "grid")

    def extend(self, n):
        return 0


class QuotaCtx:
    def __init__(self, cap, idle):
        self.cap = cap
        self.client = type("C", (), {"idle_slots": lambda _self: idle})()


class QuotaTest(unittest.TestCase):
    def quotas(self, cap=100, g=0.4, j=0.3):
        with unittest.mock.patch.dict(TUNING.values, {"climb.grid_cores": g, "climb.jump_cores": j}):
            return climb.core_quotas(QuotaCtx(cap, cap))

    def test_forty_and_thirty_percent_of_the_island(self):
        self.assertEqual(self.quotas(419), {"grid": 168, "jump": 126})
        self.assertIsNone(self.quotas(419, g=0.0))

    def test_each_part_stops_at_its_quota(self):
        a, b = QuotaClimber(heads=30, leaps=40, lines=10), QuotaClimber(heads=30, leaps=40)
        out = climb.submit_by_quota(QuotaCtx(100, 100), [a, b], {"grid": 40, "jump": 30})
        self.assertEqual(out, {"grid": 40, "jump": 30})
        self.assertEqual(a.live["grid"] + b.live["grid"], 40)
        self.assertEqual(a.live["jump"] + b.live["jump"], 30)

    def test_a_jump_wave_head_counts_as_a_jump(self):
        cl = QuotaClimber(jump_heads=50, kind="jumps")
        out = climb.submit_by_quota(QuotaCtx(100, 100), [cl], {"grid": 40, "jump": 30})
        self.assertEqual(out, {"grid": 0, "jump": 30})

    def test_every_climber_gets_an_equal_share(self):
        big, small = QuotaClimber(heads=60, leaps=60), QuotaClimber(heads=5, leaps=60)
        out = climb.submit_by_quota(QuotaCtx(100, 100), [big, small], {"grid": 40, "jump": 30})
        self.assertEqual((big.live["grid"], small.live["grid"]), (20, 5))          # small's unused grid share stays free
        self.assertEqual((big.live["jump"], small.live["jump"]), (15, 15))
        self.assertEqual(out, {"grid": 25, "jump": 30})

    def test_jumps_are_topped_up_to_the_quota(self):
        cl = QuotaClimber(leaps=10)
        cl.refill = 100
        out = climb.submit_by_quota(QuotaCtx(100, 100), [cl], {"grid": 40, "jump": 30})
        self.assertEqual(out["jump"], 30)
        self.assertEqual(cl.planned, 20)

    def test_no_new_jumps_are_asked_again_in_the_same_wave(self):
        cl = QuotaClimber(leaps=0)
        climb.submit_by_quota(QuotaCtx(100, 100), [cl], {"grid": 40, "jump": 30})
        self.assertEqual(cl._leaps_dry, 0)

    def test_heads_first_go_past_the_grid_quota_but_extras_do_not(self):
        cl = QuotaClimber(heads=60, lines=30)
        with unittest.mock.patch.dict(TUNING.values, {"climb.heads_first": 1}):
            out = climb.submit_by_quota(QuotaCtx(100, 100), [cl], {"grid": 40, "jump": 30})
        self.assertEqual(out["grid"], 60)                       # every head; no line point once past the quota
        self.assertEqual(len(cl.extra), 30)

    def test_a_flat_climber_moves_half_its_grid_share_to_jumps(self):
        flat, live = QuotaClimber(heads=60, leaps=60), QuotaClimber(heads=60, leaps=60)
        flat.last_wave_flat = True
        with unittest.mock.patch.dict(TUNING.values, {"climb.flat_shift": 0.5, "climb.heads_first": 0}):
            climb.submit_by_quota(QuotaCtx(200, 200), [flat, live], {"grid": 80, "jump": 60})
        self.assertEqual((flat.live["grid"], flat.live["jump"]), (20, 50))
        self.assertEqual((live.live["grid"], live.live["jump"]), (40, 30))

    def test_runs_already_out_count_against_the_quota(self):
        cl = QuotaClimber(heads=50, leaps=50)
        cl.live = {"grid": 35, "jump": 30}
        out = climb.submit_by_quota(QuotaCtx(100, 100), [cl], {"grid": 40, "jump": 30})
        self.assertEqual(out, {"grid": 5, "jump": 0})


class RfCtx(JumpCtx):
    """Scored configs score 80 + MMQ / 10: the forest should rank high MMQ first."""

    def __init__(self, space, n=80):
        super().__init__(space, lambda c: 80.0 + float(c.params.get(MMQ, 0)) / 10.0)
        self.rng = random.Random(0)
        self.proposed, self.sent = set(), []
        lad, other = self.ladders.rungs[MMQ], self.ladders.rungs[BQ]
        for i in range(n):
            c = space._assemble({MMQ: lad[self.rng.randrange(len(lad) // 2)], BQ: self.rng.choice(other)})
            if c is not None:
                self.done[c.key] = (c, self.f(c))

    def scored(self):
        return list(self.done.values())

    def fresh(self, c):
        return c is not None and c.key not in self.done and c.key not in self.proposed

    def propose(self, configs, tag, batches, priority=0):
        loop = asyncio.get_event_loop()
        out = []
        for c in configs:
            self.proposed.add(c.key)
            self.sent.append((c, tag))
            out.append(loop.create_future())
        return out


class StepBTest(unittest.TestCase):
    def setUp(self):
        self.ctx = RfCtx(Space(ranges()))

    async def two_calls(self, n):
        async def accept(c):
            return self.ctx.fresh(c)
        first = await rf_ucb.step_b(self.ctx, n, accept)
        await self.ctx._rf["task"]
        return first, await rf_ucb.step_b(self.ctx, n, accept)

    def test_nothing_until_the_first_fit_then_a_ranked_batch(self):
        first, second = asyncio.run(self.two_calls(12))
        self.assertEqual(first, [])
        self.assertTrue(1 <= len(second) <= 12)
        self.assertEqual(self.ctx._rf["fits"], 1)
        self.assertEqual(self.ctx._rf["sent"], len(second))
        self.assertTrue(all(self.ctx.fresh(c) for c in second))

    def test_too_few_scores_propose_nothing(self):
        self.ctx.done = dict(list(self.ctx.done.items())[:rf_ucb.MIN_TRAIN - 1])
        async def accept(c):
            return True
        self.assertEqual(asyncio.run(rf_ucb.step_b(self.ctx, 5, accept)), [])
        self.assertIsNone(self.ctx._rf["task"])

    def test_submit_rf_tags_r_and_counts_as_ga(self):
        async def run():
            await submit_rf(self.ctx, 6)                        # starts the fit
            await self.ctx._rf["task"]
            return await submit_rf(self.ctx, 6)
        futs = asyncio.run(run())
        self.assertEqual(len(futs), len(self.ctx.sent))
        self.assertTrue(futs)
        self.assertTrue(all(tag == "R" for _, tag in self.ctx.sent))
        self.assertTrue({c.key for c, _ in self.ctx.sent} <= self.ctx.ga_keys)

    def test_share_zero_sends_no_rf(self):
        self.assertEqual(asyncio.run(submit_rf(self.ctx, 0)), [])
        self.assertFalse(hasattr(self.ctx, "_rf"))


if __name__ == "__main__":
    unittest.main()


class RebalanceTest(unittest.TestCase):
    class View:
        def __init__(self, climbing, quota):
            self.quota = quota
            self.notes = {"engine": {"basins": [{"status": "climbing"}] * climbing + [{"status": "peak"}]}}

    def test_weighted_quotas_add_up(self):
        from search.islands import weighted_quotas
        self.assertEqual(weighted_quotas(838, [3, 1]), [629, 209])
        self.assertEqual(sum(weighted_quotas(839, [2, 2, 3])), 839)

    def test_quotas_follow_the_climbing_climbers(self):
        from search.islands import rebalance
        views = [self.View(3, 419), self.View(1, 419)]
        with unittest.mock.patch.dict(TUNING.values, {"islands.by_climber": 1}):
            rebalance(views, 838)
        self.assertEqual([v.quota for v in views], [629, 209])

    def test_an_island_without_climbers_keeps_one_share_and_off_changes_nothing(self):
        from search.islands import rebalance
        views = [self.View(2, 419), self.View(0, 419)]
        with unittest.mock.patch.dict(TUNING.values, {"islands.by_climber": 1}):
            rebalance(views, 900)
        self.assertEqual([v.quota for v in views], [600, 300])
        views = [self.View(2, 419), self.View(0, 419)]
        with unittest.mock.patch.dict(TUNING.values, {"islands.by_climber": 0}):
            rebalance(views, 838)
        self.assertEqual([v.quota for v in views], [419, 419])
