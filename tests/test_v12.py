"""v12: trust-region rules (box weights, steps, clock cap, restarts, duplicates, endgame) and a whole virtual round."""
import random
import unittest

import numpy as np

from common.tuning import DEFAULTS, TUNING
from search.space import Space
from search.v12.engine import (Region, TRLoop, box_weights, clock_cap, coarse_to_fine, endgame_centres, indel_part,
                               loop_spent, region_lines, screen_keep, step_update)
from tests.helpers import ranges
from tests.test_v10 import ROUND, Stub

PAD = "assembly_region_padding"


class RulesTest(unittest.TestCase):
    def test_box_weights_narrow_what_matters(self):
        w = box_weights(np.array([0.5, 0.05, 0.005]), clip=5.0)
        self.assertLess(w[0], w[1])
        self.assertLess(w[1], w[2])
        self.assertAlmostEqual(float(np.exp(np.mean(np.log(np.array([0.2, 2.0, 2.5]))))), 1.0, places=0)
        self.assertTrue(np.all(w <= 5.0) and np.all(w >= 0.2))

    def test_steps_grow_after_successes_and_shrink_after_failures(self):
        c = Space(ranges())._assemble({})
        r = Region(0, c, 90.0, 0.8, 0.0, step_start=90.0)
        for i in range(3):
            r.best += 0.01
            out = step_update(r, gain=0.005, succ_tol=3, fail_tol=2, l_max=1.6)
        self.assertEqual(out, "grow")
        self.assertAlmostEqual(r.length, 1.6)
        self.assertEqual(step_update(r, 0.005, 3, 2, 1.6), "failure")
        self.assertEqual(step_update(r, 0.005, 3, 2, 1.6), "shrink")
        self.assertAlmostEqual(r.length, 0.8)
        r.best += 0.001                                   # below the gain: still a failure
        self.assertEqual(step_update(r, 0.005, 3, 2, 1.6), "failure")

    def test_the_clock_cap_falls_to_the_minimum_at_the_endgame(self):
        self.assertAlmostEqual(clock_cap(0.0, 1.6, 0.0078125), 1.6)
        self.assertAlmostEqual(clock_cap(1.0, 1.6, 0.0078125), 0.0078125)
        self.assertGreater(clock_cap(0.5, 1.6, 0.0078125), clock_cap(0.6, 1.6, 0.0078125))

    def test_endgame_takes_the_centres_near_the_best(self):
        c = Space(ranges())._assemble({})
        rs = [Region(i, c, b, 0.1, 0.0) for i, b in enumerate((82.40, 82.65, 82.36, 82.10))]
        self.assertEqual([r.id for r in endgame_centres(rs, 4, 0.3)], [1, 0, 2])


class LinesTest(unittest.TestCase):
    def test_coarse_to_fine_spans_the_ladder_first_and_covers_it(self):
        order = coarse_to_fine(51, 16)
        self.assertEqual(sorted(order), [i for i in range(51) if i != 16])
        self.assertIn(32, order[:3])                     # the winning +17 move (16 -> 33) is next to an early point
        self.assertEqual({0, 50} & set(order[:3]), {0, 50})

    def test_region_lines_cover_cheap_settings_far_and_early(self):
        from search.v9.steps import Ladders
        space = Space(ranges())

        class Snap:
            pass
        snap = Snap()
        snap.space, snap.ladders = space, Ladders(space)
        snap.scored = lambda: iter([])
        centre = space._assemble({})
        lines = region_lines(snap, centre, coarse=8)
        BQST = "base_quality_score_threshold"
        n = len(snap.ladders.rungs[BQST])
        p = snap.ladders.index(BQST, centre.params[BQST])
        early = [c for c in lines[:80] if c.params[BQST] != centre.params[BQST]]
        self.assertTrue(early)                                            # a cheap (PairHMM) setting has a line
        self.assertTrue(any(abs(snap.ladders.index(BQST, c.params[BQST]) - p) >= n / 4 for c in early))   # far, early
        moved = {next(k for k in c.params if c.params[k] != centre.params[k]) for c in lines[:60]}
        self.assertGreater(len(moved), 10)                                # interleaved across settings
        self.assertTrue(all(sum(c.params[k] != centre.params[k] for k in c.params) == 1 for c in lines))


class SpentTest(unittest.TestCase):
    def test_spent_when_every_region_is_quiet(self):
        c = Space(ranges())._assemble({})
        a, b = Region(0, c, 90, 0.1, 0.0, last_gain=100.0), Region(1, c, 89, 0.1, 0.0, last_gain=500.0)
        self.assertFalse(loop_spent([a, b], now=700.0, quiet_s=360.0))
        self.assertTrue(loop_spent([a, b], now=900.0, quiet_s=360.0))
        self.assertFalse(loop_spent([], now=900.0, quiet_s=360.0))


class StubLoop:
    """A TRLoop with a stub context: settings distance only (no sketches)."""

    def loop(self, now=600.0):
        space = Space(ranges())

        class Dist:
            def between(self, a, b):
                from search.v9.distance import setting_distance
                return setting_distance(space, a, b), "settings"

        class Ctx:
            pass
        ctx = Ctx()
        ctx.space, ctx.dist, ctx.rng = space, Dist(), random.Random(0)
        ctx.elapsed = lambda: now
        ctx.T = lambda s: s
        ctx.fresh = lambda c: True
        loop = TRLoop.__new__(TRLoop)
        loop.ctx, loop.regions, loop.finished, loop.next_id = ctx, {}, [], 0
        loop.key_region, loop.events, loop.loop_start, loop.best = {}, [], 0.0, float("-inf")
        loop.reopened, loop.names = set(), [n for n in space.names if space.specs[n].type != "bool"]
        loop.screening, loop.screen_opened, loop.last_halving = False, False, 0.0
        loop.lines, loop.by_region, loop.order = {}, {}, None
        return space, loop

    @staticmethod
    def far(space, base, k):
        p = dict(base.params)
        names = [n for n in space.names if space.specs[n].type != "bool" and not space.specs[n].allowed][:k]
        for n in names:
            s = space.specs[n]
            p[n] = s.max if p[n] != s.max else s.min
        return space._assemble(p)

class ManageTest(StubLoop, unittest.TestCase):
    def test_regions_start_far_apart_and_a_collapsed_one_restarts_elsewhere(self):
        space, loop = self.loop()
        a = space._assemble({})
        a2 = space._assemble({**a.params, PAD: a.params[PAD] + 1})          # near a: never a second start
        b, c = self.far(space, a, 5), self.far(space, a, 9)
        pts = [(a, 90.0, 0.9, 1), (a2, 89.9, 0.9, 1), (b, 89.5, 0.9, 1), (c, 89.0, 0.9, 1)]
        TUNING.values = {**TUNING.values, "tr.regions": 2, "tr.sh_regions": 1}     # normal mode (no screening)
        try:
            loop.manage(pts)
            self.assertEqual(sorted(r.centre.key for r in loop.regions.values()), sorted([a.key, b.key]))
            ra = next(r for r in loop.regions.values() if r.centre.key == a.key)
            ra.length = 1e-6                                                 # collapsed
            loop.manage(pts)
            centres = sorted(r.centre.key for r in loop.regions.values())
            self.assertEqual(centres, sorted([b.key, c.key]))               # a is tabu now; c is the next far start
            self.assertEqual([r.centre.key for r in loop.finished], [a.key])
        finally:
            TUNING.values = dict(DEFAULTS)

    def test_two_regions_on_one_hill_keep_the_better(self):
        space, loop = self.loop()
        a = space._assemble({})
        a2 = space._assemble({**a.params, PAD: a.params[PAD] + 1})
        loop.open_region((a, 90.0, 0.9, 1), "test")
        loop.open_region((a2, 89.0, 0.9, 1), "test")
        TUNING.values = {**TUNING.values, "tr.regions": 2}
        try:
            loop.manage([(a, 90.0, 0.9, 1), (a2, 89.0, 0.9, 1)])
        finally:
            TUNING.values = dict(DEFAULTS)
        centres = [r.centre.key for r in loop.regions.values()]
        self.assertIn(a.key, centres)
        self.assertNotIn(a2.key, centres)                                   # the weaker one on the same hill closed
        self.assertIn("same hill", loop.finished[0].why)                    # (its slot may take a fresh start)


class FixesTest(StubLoop, unittest.TestCase):
    def test_the_best_in_no_region_gets_a_small_region_once(self):
        space, loop = self.loop()
        a = space._assemble({})
        b = self.far(space, a, 9)
        loop.open_region((b, 89.0, 0.9, 1), "test")
        TUNING.values = {**TUNING.values, "tr.regions": 1}
        try:
            loop.manage([(a, 90.0, 0.9, 1), (b, 89.0, 0.9, 1)])
            got = {r.centre.key: r for r in loop.regions.values()}
            self.assertIn(a.key, got)
            self.assertAlmostEqual(got[a.key].length, DEFAULTS["tr.reopen_length"])
            loop.close_region(got[a.key], "collapsed")
            loop.manage([(a, 90.0, 0.9, 1), (b, 89.0, 0.9, 1)])
            self.assertNotIn(a.key, {r.centre.key for r in loop.regions.values()})   # once per config
        finally:
            TUNING.values = dict(DEFAULTS)

    def test_a_fresh_start_when_nothing_evaluated_is_far_enough(self):
        space, loop = self.loop()
        a = space._assemble({})
        loop.points = lambda: [(a, 90.0, 0.9, 1)]
        TUNING.values = {**TUNING.values, "tr.regions": 2, "tr.sh_regions": 1}
        try:
            loop.manage([(a, 90.0, 0.9, 1)])
        finally:
            TUNING.values = dict(DEFAULTS)
        fresh = [r for r in loop.regions.values() if r.centre.key != a.key]
        self.assertEqual(len(fresh), 1)
        self.assertEqual(fresh[0].best, float("-inf"))
        self.assertEqual(list(fresh[0].seeds), [fresh[0].centre])        # its own centre is evaluated first
        self.assertTrue(any("fresh start" in e for e in loop.events))

    def test_lines_are_picked_and_counted_for_their_region(self):
        space, loop = self.loop()
        a = space._assemble({})
        r = loop.open_region((a, 90.0, 0.9, 1), "test")
        from collections import deque
        l1 = space._assemble({**a.params, PAD: a.params[PAD] + 7})
        l2 = space._assemble({**a.params, PAD: a.params[PAD] + 9})
        loop.lines = {r.id: deque([l1, l2])}
        loop.inject_q = deque()
        loop.next_inject = lambda n: None
        TUNING.values = {**TUNING.values, "tr.share_lines": 1.0}
        try:
            got = loop.take(2)
        finally:
            TUNING.values = dict(DEFAULTS)
        self.assertEqual([(c.key, src) for c, src in got], [(l1.key, "line"), (l2.key, "line")])
        self.assertEqual(loop.key_region[l1.key], r.id)


class ScreenTest(StubLoop, unittest.TestCase):
    def test_keep_ranks_by_best_plus_gain(self):
        c = Space(ranges())._assemble({})
        hi_flat = Region(0, c, 82.40, 0.8, 0.0, start_best=82.39)      # high, barely climbing
        lo_fast = Region(1, c, 82.35, 0.8, 0.0, start_best=82.00)      # lower, climbing fast
        lo_flat = Region(2, c, 82.10, 0.8, 0.0, start_best=82.09)
        fresh = Region(3, c, float("-inf"), 0.8, 0.0)                 # nothing back yet
        self.assertEqual([r.id for r in screen_keep([hi_flat, lo_fast, lo_flat, fresh], 2, 1.0)], [1, 0])
        self.assertEqual([r.id for r in screen_keep([hi_flat, lo_fast, lo_flat, fresh], 2, 0.0)], [0, 1])

    def test_screening_opens_many_then_halves_to_the_target(self):
        space, loop = self.loop()
        base = space._assemble({})
        names = [n for n in space.names if space.specs[n].type != "bool" and not space.specs[n].allowed]
        pts = []
        for i in range(8):                                           # 8 configs far apart in settings
            p = dict(base.params)
            for j, n in enumerate(names):
                s_ = space.specs[n]
                p[n] = s_.max if (i >> (j % 3)) & 1 else s_.min
            p[names[-1 - i % len(names)]] = space.specs[names[-1 - i % len(names)]].max
            c = space._assemble(p)
            if c is not None and all(c.key != q[0].key for q in pts):
                pts.append((c, 90.0 - 0.1 * i, 0.9, 1))
        TUNING.values = {**TUNING.values, "tr.regions": 2, "tr.sh_regions": 8, "tr.sh_step_min": 1.0,
                         "tr.setting_new_far": 0.05}
        try:
            loop.manage(pts)
            n_open = len(loop.regions)
            self.assertGreater(n_open, 2)
            self.assertTrue(loop.screening)
            for r in loop.regions.values():                          # the last-opened climbs fastest
                r.best = r.start_best + 0.05 * r.id
            loop.ctx.elapsed = lambda: 600.0 + 61.0                  # one step later
            loop.manage(pts)
            self.assertEqual(len(loop.regions), max(2, (n_open + 1) // 2))
            self.assertTrue(all("screened out" in r.why for r in loop.finished))
            while loop.screening:
                t = loop.ctx.elapsed() + 61.0
                loop.ctx.elapsed = lambda t=t: t
                loop.manage(pts)
            self.assertEqual(len(loop.regions), 2)
        finally:
            TUNING.values = dict(DEFAULTS)


class IndelFarTest(StubLoop, unittest.TestCase):
    def test_indel_part_is_core_minus_the_snp_terms(self):
        from common.scoring import v2_score
        from common.types import ClassCounts
        counts = {"snp_hom": ClassCounts(90, 2, 10), "snp_het": ClassCounts(200, 5, 20), "indel_1bp": ClassCounts(40, 6, 9),
                  "indel_2_3bp": ClassCounts(12, 1, 5), "indel_4_7bp": ClassCounts(5, 0, 3), "indel_8bp": ClassCounts(0, 0, 0)}
        core = v2_score(counts, 100, 50)[0]
        f1 = lambda c: 2 * c.tp / (2 * c.tp + c.fp + c.fn)
        tw = 0.02 + 0.18 + 0.40 + 0.16 + 0.10                        # indel_8bp has no truth and no calls: dropped
        snp = (0.02 * f1(counts["snp_hom"]) + 0.18 * f1(counts["snp_het"])) / tw
        self.assertAlmostEqual(indel_part(counts), core - snp, places=12)
        counts["indel_8bp"] = ClassCounts(0, 3, 0)                    # calls without truth: kept at F1 0
        self.assertLess(indel_part(counts), core - snp)

    def test_default_screen_does_not_reward_a_low_start_for_climbing_fast(self):
        # fleet 044200Z, 2026-10-11: far starts from ~75 gained +6-7 in 8 minutes and took 3 of the 4 places at gain
        # weight 1; the best normal start (82.03, +0.08) was cut in the first halving
        c = Space(ranges())._assemble({})
        normal = [Region(i, c, b, 0.8, 0.0, start_best=s0) for i, (b, s0) in
                  enumerate([(82.03, 81.94), (82.05, 81.72), (81.92, 81.42), (81.58, 81.19), (80.98, 80.68)])]
        fars = [Region(5 + i, c, 81.7 - 0.3 * i, 0.8, 0.0, start_best=74.7 + 0.4 * i, far=True) for i in range(3)]
        w = float(DEFAULTS["tr.sh_gain_weight"])
        kept = {r.id for r in screen_keep(normal + fars, 4, w, indel_keep=1, far_keep=1)}
        self.assertEqual(kept, {0, 1, 2, 5})                         # the three best normal starts + the best far start

    def test_keep_reserves_the_best_indel_region_and_the_best_far_start(self):
        c = Space(ranges())._assemble({})
        top = [Region(i, c, 82.4 - 0.01 * i, 0.8, 0.0, start_best=82.3, indel=0.660) for i in range(4)]
        indel = Region(4, c, 82.0, 0.8, 0.0, start_best=81.99, indel=0.668)         # lowest score, best indels
        far1 = Region(5, c, 81.0, 0.8, 0.0, start_best=80.0, far=True)
        far2 = Region(6, c, 80.5, 0.8, 0.0, start_best=80.4, far=True)
        regions = top + [indel, far1, far2]
        self.assertEqual([r.id for r in screen_keep(regions, 4, 1.0)], [0, 1, 2, 3])
        kept = [r.id for r in screen_keep(regions, 4, 1.0, indel_keep=1, far_keep=1)]
        self.assertEqual(sorted(kept), [0, 1, 4, 5])                  # far1 climbed faster than far2
        self.assertEqual(sorted(r.id for r in screen_keep(regions, 1, 1.0, indel_keep=1, far_keep=1)), [4])

    def test_screening_opens_far_starts_and_the_halving_keeps_the_indel_region(self):
        from search.v9.distance import setting_distance
        space, loop = self.loop()
        home = space._assemble({})
        names = [n for n in space.names if space.specs[n].type != "bool" and not space.specs[n].allowed]
        near = [space._assemble({**home.params, PAD: home.params[PAD] + k}) for k in range(1, 4)]
        mid = [self.far(space, home, k) for k in (3, 4, 5)]           # other hills, not far from home in settings
        far = self.far(space, home, 12)
        pts = [(home, 90.0, 0.9, 100)] + [(c, 89.9 - 0.01 * i, 0.9, 100) for i, c in enumerate(near)] + \
              [(c, 89.5 - 0.1 * i, 0.9, 100) for i, c in enumerate(mid)] + [(far, 85.0, 0.8, 100)]
        TUNING.values = {**TUNING.values, "tr.regions": 2, "tr.sh_regions": 4, "tr.sh_far_starts": 1, "tr.sh_home": 4,
                         "tr.sh_far_setting": 0.3, "tr.setting_new_far": 0.05, "tr.sh_step_min": 1.0,
                         "tr.sh_indel_keep": 1, "tr.sh_far_keep": 1}
        try:
            loop.manage(pts)
            self.assertEqual(len(loop.regions), 4)
            fars = [r for r in loop.regions.values() if r.far]
            self.assertEqual([r.centre.key for r in fars], [far.key])
            self.assertTrue(all(setting_distance(space, far, h) >= 0.3 for h, *_ in pts[:4]))
            ids = {r.centre.key: r.id for r in loop.regions.values()}
            self.assertEqual(set(ids), {home.key, mid[0].key, mid[1].key, far.key})
            weak = ids[mid[1].key]                                    # the lowest-scoring normal start
            loop.indel_of = lambda c: 0.70 if loop.key_region.get(c.key) == weak else 0.60
            over = space._assemble({**home.params, PAD: home.params[PAD] + 10})
            loop.key_region[over.key] = ids[home.key]                 # best indels, but beyond the reference FP
            loop.indel_of = (lambda f: lambda c: 0.99 if c.key == over.key else f(c))(loop.indel_of)
            pts2 = pts + [(over, 80.0, 0.9, 200)]
            self.assertEqual(loop.region_indel(pts2)[weak], 0.70)
            self.assertEqual(loop.region_indel(pts2)[ids[home.key]], 0.60)
            loop.ctx.elapsed = lambda: 600.0 + 61.0
            loop.manage(pts2)
            self.assertEqual(len(loop.regions), 2)
            self.assertEqual(sorted(loop.regions), sorted([weak, fars[0].id]))
            self.assertFalse(loop.screening)
        finally:
            TUNING.values = dict(DEFAULTS)


class VirtualRoundTest(unittest.TestCase):
    def test_a_whole_round(self):
        from evaluator import FakeRunner, FakeScorer, SimExecutor
        from main import DryRunSubmitter, LoopSettings, run_round
        from search.v12.engine import searcher
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=3)
        ex = SimExecutor(runner, 64, duration=lambda job: 150.0)
        seen = []
        inner = searcher(r, 64, seed=1)

        async def wrapped(client):
            seen.append(client)
            await inner(client)

        TUNING.values = {**TUNING.values, "rf.trees": 40, "tr.candidates": 500}
        try:
            out = run_round(ROUND, Stub(), ex, FakeScorer(runner), wrapped, DryRunSubmitter(), r,
                            LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=72 * 60)
        finally:
            TUNING.values = dict(DEFAULTS)
        eng = seen[0].notes["engine"]
        self.assertEqual(eng["version"], "v12")
        rf = eng["rf"]
        self.assertGreater(rf["refits"], 3)
        self.assertGreater(rf["picks"].get("tr", 0), 100)                  # Thompson-sampled picks
        self.assertTrue(rf["regions"])
        self.assertGreaterEqual(sum(x["picks"] > 0 for x in rf["regions"] + rf["finished"]), 2)  # several regions worked
        self.assertTrue(eng["basins"])                                      # the endgame climbed
        self.assertIsNotNone(out.best)


if __name__ == "__main__":
    unittest.main()
