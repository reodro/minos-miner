import random
import unittest

from common.types import AnswerSet, Config, Regime, Round, Variant
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search.axes import DOMINANT_LEVELS, AxisBook, design_wave, structural_names
from search.basins import Basin
from search.families import genotyping_riders, pairhmm_riders
from search.ranker import Ranker
from search.space import GENOTYPING_PARAMS, Space
import search.engine as engine
from tests.helpers import ranges

ROUND = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")


class Stub:
    def build(self, rnd):
        return AnswerSet(rnd.round_id, 1, Regime.GIAB, "HG002", [Variant("chr20", 5, "A", "G", "0/1")], [], 600, {})

    def update(self, a, c):
        return a


def virtual_round(cap, minutes=72, seed=0, duration=150.0, crashes=None, hills=1):
    r = ranges()
    runner = FakeRunner(r, seed=seed, batching=True, hills=hills)
    ex = SimExecutor(runner, cap, duration=lambda job: duration, crashes=crashes or (lambda job: False))
    out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, cap), DryRunSubmitter(), r,
                    LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=minutes * 60)
    return out, ex


class AxesTest(unittest.TestCase):
    def test_design_wave_balances_levels_and_stays_valid(self):
        space = Space(ranges())
        cells = design_wave(space, random.Random(0), 60, [space._assemble({})])
        self.assertGreaterEqual(len(cells), 55)
        counts = {}
        for c in cells:
            counts[c.params["min_pruning"]] = counts.get(c.params["min_pruning"], 0) + 1
        self.assertEqual(set(counts), set(DOMINANT_LEVELS["min_pruning"]))
        self.assertLessEqual(max(counts.values()) - min(counts.values()), 5)
        self.assertEqual(len({c.key for c in cells}), len(cells))

    def test_axis_book_learns_cliffs_and_parks_dead_axes(self):
        space = Space(ranges())
        book = AxisBook(space)
        # moves above the base-quality cliff fail: nothing is proposed past the lowest failing value
        book.note_move("min_base_quality_score", 45, 30, None, failed=True, now=0.0)
        book.note_move("min_base_quality_score", 42, 30, None, failed=True, now=0.0)
        self.assertFalse(book.allowed("min_base_quality_score", 43))
        self.assertTrue(book.allowed("min_base_quality_score", 40))
        # an axis whose moves change nothing is parked, then retested later
        for _ in range(8):
            book.note_move("max_reads_per_alignment_start", 60, 50, 0.0, failed=False, now=100.0)
        self.assertNotIn("max_reads_per_alignment_start", book.live(now=100.0))
        self.assertIn("max_reads_per_alignment_start", book.live(now=100.0 + 601))
        # a winning move widens the step, repeated failures narrow it
        step0 = book.axes["assembly_region_padding"].step
        book.note_move("assembly_region_padding", 60, 50, 0.05, failed=False, now=0.0)
        self.assertGreater(book.axes["assembly_region_padding"].step, step0)


class FarAxesTest(unittest.TestCase):
    def test_counts_axes_more_than_a_tenth_of_their_range_apart(self):
        from search.axes import far_axes
        space = Space(ranges())
        a = space._assemble({})
        b = space._assemble({"min_pruning": 6, "assembly_region_padding": 13, "min_mapping_quality_score": 21,
                             "pcr_indel_model": "NONE"})
        names = ["min_pruning", "assembly_region_padding", "min_mapping_quality_score", "pcr_indel_model",
                 "min_dangling_branch_length"]
        self.assertEqual(far_axes(space, names, a, b), 3)       # pruning 2->6 (4 of 8), padding 100->13, the enum; mapq 20->21 is near
        self.assertEqual(far_axes(space, names, a, a), 0)


class WidthTest(unittest.TestCase):
    def test_structural_riders_and_screens_change_the_right_axes(self):
        from search.families import STRUCTURAL, screen_moves, structural_riders
        space = Space(ranges())
        base = space._assemble({})
        rng = random.Random(0)
        kids = structural_riders(space, base, 4, rng)
        self.assertEqual(len(kids), 4)
        for k in kids:
            changed = [n for n in space.names if k.params.get(n) != base.params.get(n)]
            self.assertTrue(1 <= len(changed) <= 2 and all(n in STRUCTURAL for n in changed), changed)
        names = [n for n in space.names if n not in GENOTYPING_PARAMS]
        screens = screen_moves(space, base, rng, 6, names)
        self.assertGreaterEqual(len(screens), 4)
        for c in screens:
            changed = [n for n in names if c.params.get(n) != base.params.get(n)]
            self.assertTrue(3 <= len(changed) <= 8, changed)      # several axes at once


class FamiliesTest(unittest.TestCase):
    def test_riders_differ_from_the_base_in_one_cheap_axis(self):
        space = Space(ranges())
        base = space._assemble({})
        kids = genotyping_riders(space, base, 24)
        self.assertGreaterEqual(len(kids), 12)
        for k in kids:
            changed = [n for n in space.names if k.params.get(n) != base.params.get(n)]
            self.assertEqual(len(changed), 1)
            self.assertIn(changed[0], GENOTYPING_PARAMS)
        self.assertEqual(len({k.key for k in kids}), len(kids))
        pair = pairhmm_riders(space, base, 6)
        self.assertGreaterEqual(len(pair), 4)
        self.assertTrue(all("pcr_indel_model" in c.params for c in pair))


class BasinTest(unittest.TestCase):
    def test_climb_accepts_better_moves_and_merges_winners(self):
        space = Space(ranges())
        book = AxisBook(space)
        centre = space._assemble({})
        b = Basin(0, space, book, centre, 80.0, 0.0, 0.0)
        n = b.refill(random.Random(0), [], now=0.0, want_pairs=0)   # no pairs: a pair tried first would own the merge's config
        self.assertGreater(n, 50)
        kinds = {m.kind for m in b.pool}
        self.assertTrue({"step", "probe"} <= kinds)
        m1 = next(m for m in b.pool if m.kind == "step" and m.axes == ("assembly_region_padding",))
        m2 = next(m for m in b.pool if m.kind == "step" and m.axes == ("min_mapping_quality_score",))
        b.pool.remove(m1); b.pool.remove(m2)
        b.inflight, b.gen_sent = 2, 2                                             # the sweep: these two moves
        self.assertFalse(b.consider(m2, m2.config, 80.03, False, now=1.0))     # a winner, but half the sweep is out
        self.assertTrue(b.consider(m1, m1.config, 80.045, False, now=1.5))     # sweep back: the best of it is the centre
        self.assertEqual(b.centre.key, m1.config.key)
        self.assertEqual(b.pool, [])
        b.refill(random.Random(1), [], now=2.0, want_pairs=0)
        merges = [m for m in b.pool if m.kind == "merge"]                        # the other winner rides on the new centre
        self.assertTrue(merges)
        self.assertIn("min_mapping_quality_score", merges[0].axes)
        self.assertTrue(all(m.priority == 0 for m in merges))
        # a gain below the answer set's resolution is not an improvement; a decisive one does not wait
        m3 = b.take()
        self.assertFalse(b.consider(m3, m3.config, b.best + 0.004, False, now=3.0))
        m4 = b.take()
        self.assertTrue(b.consider(m4, m4.config, b.best + 0.06, False, now=3.5))

    def test_transplants_copy_axes_from_another_basin(self):
        space = Space(ranges())
        book = AxisBook(space)
        a = Basin(0, space, book, space._assemble({}), 80.0, 0.0, 0.0)
        other = space._assemble({"min_mapping_quality_score": 59, "assembly_region_padding": 13})
        a.refill(random.Random(0), [other], now=0.0)
        tr = [m for m in a.pool if m.kind == "transplant"]
        self.assertTrue(any(m.axes == ("min_mapping_quality_score",) for m in tr))
        self.assertTrue(any(set(m.axes) == {"assembly_region_padding", "min_mapping_quality_score"} for m in tr))


class RankerTest(unittest.TestCase):
    def test_ranks_known_good_region_first_once_fitted(self):
        r = ranges()
        space = Space(r)
        runner = FakeRunner(r, batching=True)
        rk = Ranker(space, n_trees=20, min_points=20)
        rng = random.Random(0)
        for _ in range(120):
            c = space.random(rng)
            q = runner.quality(c)
            rk.add(c, 0.5 + 0.4 * q, 5.0 * (1 - q))
        rk.fit()
        self.assertTrue(rk.ready)
        cands = [space.random(rng) for _ in range(40)]
        order = rk.rank(cands, best=0.0)
        top, bottom = [runner.quality(cands[i]) for i in order[:10]], [runner.quality(cands[i]) for i in order[-10:]]
        self.assertGreater(sum(top) / 10, sum(bottom) / 10)


class CapacityTest(unittest.TestCase):
    def test_a_fleets_capacity_is_read_once_its_nodes_have_joined(self):
        import asyncio
        from unittest import mock

        class Fleet:                                  # nodes join over the first polls; the last one is slow
            def __init__(self):
                self.calls = 0

            def nodes(self):                          # a node still joining holds the wait even when the count is still
                joining = self.calls < 9
                return [{"ready": True, "rejected": ""}, {"ready": not joining, "rejected": ""}, {"ready": False, "rejected": "old code"}]

            def capacity(self):
                self.calls += 1
                return [0, 0, 300, 300, 300, 300, 300, 300, 300, 700, 700, 700][min(self.calls, 11)]

        fleet = Fleet()
        client = mock.Mock(executor=fleet, capacity=fleet.capacity)
        with mock.patch.multiple(engine, CAPACITY_WAIT_S=5.0, CAPACITY_SETTLE_S=0.25, CAPACITY_POLL_S=0.05):
            cap = asyncio.run(engine.settled_capacity(client, 14))
        self.assertEqual(cap, 700)
        local = mock.Mock(executor=object(), capacity=lambda: 14)       # no nodes(): a local executor answers at once
        self.assertEqual(asyncio.run(engine.settled_capacity(local, 8)), 14)


class EngineRoundTest(unittest.TestCase):
    def test_small_fleet_climbs_all_round_and_fills_the_workers(self):
        out, ex = virtual_round(16)
        self.assertEqual(out.failed, 0)
        self.assertGreater(ex.utilization(), 0.9)
        opening = max(r.score for (k, fid), r in out.callsets.items() if False) if False else None
        self.assertGreater(out.best.score, 91.0)                # the fake landscape's seeds score ~90
        self.assertGreater(out.evaluated, 1500)                 # riders: many configs per traversal

    def test_big_fleet_races_several_basins(self):
        with self.assertLogs("search.engine", level="INFO") as logs:
            out, ex = virtual_round(950, minutes=30)
        text = "\n".join(logs.output)
        self.assertIn("4 basins max", text)
        self.assertGreaterEqual(text.count(" born at "), 3)
        self.assertGreater(ex.utilization(), 0.85)
        self.assertGreater(out.evaluated, 20000)

    def test_two_hills_give_two_basins_after_a_valley_test(self):
        with self.assertLogs("search.engine", level="INFO") as logs:
            out, ex = virtual_round(950, minutes=40, hills=2)
        text = "\n".join(logs.output)
        self.assertIn("a valley between", text)
        self.assertGreaterEqual(text.count(" born at "), 2)
        self.assertGreater(out.best.score, 91.0)
        last = [l for l in logs.output if " basins " in l][-1]                 # the leader holds the best found
        leader = max(float(x.split(":")[1].split("/")[0]) for x in last.split("basins ")[1].split(" | ")[0].split())
        self.assertGreaterEqual(leader, out.best.score - 0.3)

    def test_crashes_are_counted_and_not_fatal(self):
        bad = lambda job: int(job.config.key[:2], 16) % 8 == 0           # ~12% of configs crash, at random
        out, ex = virtual_round(16, minutes=40, crashes=bad)
        self.assertGreater(out.failed, 0)
        self.assertEqual(sum(out.failures.values()), out.failed)
        self.assertGreater(out.best.score, 90.5)                       # the climb goes on around the crashes


if __name__ == "__main__":
    unittest.main()
