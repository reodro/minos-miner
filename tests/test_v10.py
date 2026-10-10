"""search engine v10 (search/v10): the forest, UCB, crowding weights, QLock's test, the proposers, the pick quotas, and
a whole virtual round."""
import random
import unittest

import numpy as np

from common.tuning import DEFAULTS, TUNING
from common.types import AnswerSet, Regime, Round, Variant
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search.space import Space
from search.v9.steps import Ladders
from search.v10 import proposers
from search.v10.forest import Forest
from tests.helpers import ranges

ROUND = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")
MP, PAD = "min_pruning", "assembly_region_padding"


class Stub:
    def build(self, rnd):
        return AnswerSet(rnd.round_id, 1, Regime.GIAB, "HG002", [Variant("chr20", 5, "A", "G", "0/1")], [], 600, {})

    def update(self, a, c):
        return a


def landscape(c):
    """A step in min_pruning (jagged, as GATK's thresholds are) plus a slope in padding; nothing else matters."""
    return 90.0 + (1.0 if c.params[MP] >= 5 else 0.0) - abs(c.params[PAD] - 60) / 500.0


class ForestTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.lad = Ladders(self.space)
        rng = random.Random(0)
        self.train = [self.space.random(rng) for _ in range(400)]
        self.forest = Forest(self.space, self.lad, trees=60, min_leaf=2, jobs=1)
        self.forest.fit(self.train, [landscape(c) for c in self.train], radius=0.0)

    def test_the_forest_learns_a_step_and_ucb_ranks_by_it(self):
        lo = self.space._assemble({MP: 2, PAD: 60})
        hi = self.space._assemble({MP: 7, PAD: 60})
        mu, sd = self.forest.mean_sd([lo, hi])
        self.assertGreater(mu[1] - mu[0], 0.7)
        u, _, _, _ = self.forest.ucb([lo, hi], kappa=1.0)
        self.assertGreater(u[1], u[0])

    def test_qlock_finds_a_flat_setting_and_keeps_a_live_one(self):
        centre = self.space._assemble({MP: 7, PAD: 60})
        self.assertFalse(self.forest.flat(centre, MP, 12, 0.003))
        flat = [n for n in self.forest.names if n not in (MP, PAD) and self.forest.flat(centre, n, 12, 0.5)]
        self.assertTrue(flat)                     # some setting the score never reads reads as flat at 0.5

    def test_novelty_is_the_distance_to_the_data(self):
        seen = self.train[0]
        far_ = self.space._assemble({n: self.lad.rungs[n][-1] for n in self.forest.names})
        nov = self.forest.novelty([seen, far_])
        self.assertEqual(nov[0], 0.0)
        self.assertGreater(nov[1], 0.0)
        u0, _, _, _ = self.forest.ucb([seen, far_], kappa=0.0, kappa_dist=0.0)
        u1, _, _, _ = self.forest.ucb([seen, far_], kappa=0.0, kappa_dist=5.0)
        self.assertGreater(u1[1] - u0[1], u1[0] - u0[0])     # the far config gains more from the novelty term

    def test_the_novelty_bonus_is_capped(self):
        from search.v10.forest import NOV_CAP
        far_ = self.space._assemble({n: self.lad.rungs[n][-1] for n in self.forest.names})
        pool = self.train[:50] + [far_]
        u, mu, _, _ = self.forest.ucb(pool, kappa=0.0, kappa_dist=1.0)
        top = np.argsort(-mu)[:max(2, len(mu) // 20)]
        cap = NOV_CAP * float(np.std(mu[top]))
        self.assertLessEqual(float(u[-1] - mu[-1]), cap + 1e-9)

    def test_a_cluster_shares_its_weight(self):
        X = np.array([[0.5, 0.5]] * 9 + [[0.0, 1.0]], dtype=np.float32)
        w = Forest.density_weights(X, radius=0.01)
        self.assertAlmostEqual(w[0], 1 / 3)       # nine identical configs: each 1/sqrt(9)
        self.assertAlmostEqual(w[-1], 1.0)


    def test_a_bad_zone_sinks_without_moving_the_scores(self):
        bad = [c for c in self.train if c.params[PAD] > 200][:60]
        good = [c for c in self.train if c.params[PAD] <= 200]
        f = Forest(self.space, self.lad, trees=60, min_leaf=2, jobs=1)
        f.fit(good, [landscape(c) for c in good], radius=0.0)
        mu0, _ = f.mean_sd(self.train[:100])
        f.fit_bad(good, bad)
        mu1, _ = f.mean_sd(self.train[:100])
        np.testing.assert_allclose(mu0, mu1)       # the score forest is not taught a floor
        inside = self.space._assemble({MP: 7, PAD: 300})
        outside = self.space._assemble({MP: 7, PAD: 60})
        pb = f.p_bad([inside, outside])
        self.assertGreater(pb[0], 0.5)
        self.assertLess(pb[1], 0.5)
        u, mu, _, _ = f.ucb([inside, outside] + self.train[:30], kappa=0.0, kappa_dist=0.0)
        self.assertLess(u[0] - mu[0], u[1] - mu[1])


class ProposerTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.lad = Ladders(self.space)
        self.rng = random.Random(1)
        self.base = self.space._assemble({MP: 7, PAD: 60})

    def test_qga_rotates_toward_the_elites_and_samples_on_the_ladder(self):
        q = proposers.QGA(self.space, self.lad, 32, self.rng)
        for _ in range(30):
            q.rotate([self.base], 0.10, 0.0, 0.01)
        j = q.bin_of(MP, 7)
        self.assertEqual(max(range(len(q.amp[MP])), key=lambda i: q.amp[MP][i]), j)
        for c in q.sample(50, self.base, q.names):
            self.assertIn(c.params[MP], self.lad.rungs[MP])

    def test_frozen_settings_keep_the_base_value(self):
        live = [n for n in self.space.names if self.lad.movable(n) and n != MP]
        q = proposers.QGA(self.space, self.lad, 32, self.rng)
        for c in q.sample(30, self.base, live) + proposers.inject(self.space, self.rng, 30, self.base, live):
            self.assertEqual(c.params[MP], 7)

    def test_ga_children_mix_parents_and_mutate_by_rungs(self):
        a = self.space._assemble({MP: 2, PAD: 20})
        b = self.space._assemble({MP: 9, PAD: 300})
        names = [n for n in self.space.names if self.lad.movable(n)]
        kids = proposers.breed(self.space, self.lad, self.rng, [(a, 1.0), (b, 2.0)], 40, a, names, 0.0, 3)
        self.assertTrue(all(k.params[MP] in (2, 9) and k.params[PAD] in (20, 300) for k in kids))

    def test_the_trade_off_ends(self):
        a, b, c, d = (self.space._assemble({MP: m}) for m in (2, 4, 6, 8))
        pts = [(a, 76.43, 0.864, 1553), (b, 76.28, 0.871, 1653), (c, 76.30, 0.858, 1450), (d, 72.0, 0.95, 900)]
        ends = proposers.trade_off_ends(pts, 0.5)
        self.assertEqual([(x.params[MP], label) for x, label in ends],
                         [(2, "best score"), (4, "best core"), (6, "best FP")])   # d is junk: 4.4 below the best

    def test_pareto_keeps_the_core_and_fp_trade_offs(self):
        c = self.base
        pts = [(c, 90.0, 0.95, 400), (c, 91.0, 0.94, 100), (c, 80.0, 0.90, 300), (c, 85.0, 0.96, 900)]
        self.assertEqual(sorted(s for _, s in proposers.pareto(pts, 20.0, 10)), [85.0, 90.0, 91.0])


class TakeTest(unittest.TestCase):
    """Loop.take on its own: with nothing ranked yet, Inject keeps its share instead of filling every free worker."""

    def test_no_pool_means_inject_keeps_its_share(self):
        import search.v10.engine as eng
        space = Space(ranges())

        class Ctx:
            pass
        ctx = Ctx()
        ctx.space, ctx.ladders, ctx.rng, ctx.cap = space, Ladders(space), random.Random(0), 100
        ctx.fresh = lambda c: True
        ctx.best = lambda: None
        loop = eng.Loop.__new__(eng.Loop)
        loop.ctx, loop.names, loop.frozen = ctx, [n for n in space.names if ctx.ladders.movable(n)], set()
        loop.pool, loop.qga_q, loop.sweep_q, loop.top = eng.deque(), eng.deque(), eng.deque(), []
        loop.inject_q, loop.reinit_left = eng.deque(), 0
        got = loop.take(800)
        self.assertLessEqual(len(got), 80)        # 10% of 800, not all 800 (the first fleet replay's 767)
        self.assertTrue(all(src == "inject" for _, src in got))


class CentresTest(unittest.TestCase):
    def test_the_trade_off_ends_keep_their_slots_however_near(self):
        import search.v10.engine as eng
        space = Space(ranges())
        a, b, c, d = (space._assemble({MP: m}) for m in (2, 4, 6, 8))

        class Dist:
            def far(self, x, y, mc, ms):
                return False                       # every config "near" every other by callset distance

        class Ctx:
            pass
        ctx = Ctx()
        ctx.dist = Dist()
        loop = eng.Loop.__new__(eng.Loop)
        loop.ctx = ctx
        pts = [(a, 76.43, 0.864, 1553), (b, 76.28, 0.8745, 1738), (c, 76.30, 0.858, 1450), (d, 76.1, 0.86, 1600)]
        got = loop.pick_centres(pts, 4)
        self.assertEqual([x.params[MP] for x in got], [2, 4, 6])   # all three ends; the 4th slot needs distance


class CancelTest(unittest.TestCase):
    def test_queued_jobs_are_dropped_and_leave_in_flight(self):
        import asyncio
        from common.types import Config, EvalJob
        r = ranges(); d = r.defaults("gatk").params
        ex = SimExecutor(FakeRunner(r, seed=0), 1, duration=lambda job: 100.0, batching=False)
        jobs = [EvalJob(f"j{i}", ROUND, Config("gatk", {**d, MP: p})) for i, p in enumerate((3, 4, 5))]
        ex.submit(jobs)                            # j0 runs, j1 and j2 wait
        self.assertEqual(sorted(ex.cancel_queued(["j0", "j1", "j2"])), ["j1", "j2"])
        ex.now = 1000.0
        self.assertEqual([cs.job_id for cs in ex.poll(0)], ["j0"])    # the running one still delivers

        from search.client import EvalClient
        loop = asyncio.new_event_loop()
        try:
            client = EvalClient.__new__(EvalClient)
            client.executor = ex
            futs = {t: loop.create_future() for t in ("qga", "climb")}
            client._pending = {"a": (futs["qga"], None, "qga"), "b": (futs["climb"], None, "climb")}
            ex.submit([EvalJob("x", ROUND, Config("gatk", {**d, MP: 6}))])      # takes the one worker
            ex.submit([EvalJob("a", ROUND, jobs[1].config), EvalJob("b", ROUND, jobs[2].config)])
            self.assertEqual(client.drop_queued({"qga"}), 1)
            self.assertEqual(list(client._pending), ["b"])
            self.assertTrue(futs["qga"].cancelled())
        finally:
            loop.close()


class ConvergenceTest(unittest.TestCase):
    def test_no_gain_for_the_window_converges_and_a_restart_starts_a_new_window(self):
        from search.v10.engine import Convergence
        c = Convergence()
        args = dict(after=20 * 60, window=10 * 60, gain=0.02)
        self.assertFalse(c.check(76.10, 0, 5 * 60, **args))
        self.assertFalse(c.check(76.11, 0, 19 * 60, **args))        # +0.01 is no gain, but not before minute 20
        self.assertTrue(c.check(76.11, 0, 20 * 60, **args))         # the window from minute 5 is long over
        c = Convergence()
        c.check(76.10, 0, 30 * 60, **args)
        self.assertFalse(c.check(76.10, 1, 40 * 60, **args))         # the restart at minute 40 gets its own window
        self.assertFalse(c.check(76.10, 1, 49 * 60, **args))
        self.assertTrue(c.check(76.10, 1, 50 * 60, **args))
        c = Convergence()
        c.check(76.10, 0, 30 * 60, **args)
        self.assertFalse(c.check(76.13, 0, 39 * 60, **args))         # a gain resets the window too
        self.assertFalse(c.check(76.13, 0, 48 * 60, **args))
        self.assertTrue(c.check(76.13, 0, 49 * 60, **args))


class SeedTest(unittest.TestCase):
    def test_each_run_has_its_own_seed_unless_one_is_given(self):
        import os
        from unittest import mock
        import search.v10.engine as eng
        seeds, real = [], random.Random
        with mock.patch.object(eng.random, "Random", side_effect=lambda s: seeds.append(s) or real(s)):
            eng.searcher(ranges(), 8)
            eng.searcher(ranges(), 8)
            with mock.patch.dict(os.environ, {"MINOS_SEED": "42"}):
                eng.searcher(ranges(), 8)
        self.assertNotEqual(seeds[0], seeds[1])
        self.assertEqual(seeds[2], 42)


class VirtualRoundTest(unittest.TestCase):
    def test_a_whole_round(self):
        import search
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=1)
        ex = SimExecutor(runner, 64, duration=lambda job: 150.0)
        seen = []
        inner = search.searcher(r, 64)

        async def wrapped(client):
            seen.append(client)
            await inner(client)

        TUNING.values = {**TUNING.values, "rf.trees": 40}
        try:
            out = run_round(ROUND, Stub(), ex, FakeScorer(runner), wrapped, DryRunSubmitter(), r,
                            LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=72 * 60)
        finally:
            TUNING.values = dict(DEFAULTS)
        eng = seen[0].notes["engine"]
        self.assertEqual(eng["version"], "v10")
        self.assertGreater(eng["rf"]["refits"], 3)
        self.assertTrue(set(eng["rf"]["picks"]) >= {"opening", "inject", "mutate", "qga", "sweep"})   # the ranking decides the rest
        if eng["pick"] is not None:                # named from pick time on (a virtual clock can jump past it)
            self.assertGreaterEqual(round(out.best.score, 4), eng["pick"]["score"])
        self.assertTrue(eng["rf"]["centres"])
        self.assertIn("bad", eng["rf"])
        cal = eng["rf"]["calibration"]
        self.assertTrue(cal and all(v["n"] > 0 for v in cal.values()))
        self.assertGreater(ex.utilization(), 0.9)


if __name__ == "__main__":
    unittest.main()
