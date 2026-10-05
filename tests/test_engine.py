import random
import unittest

from common.types import AnswerSet, Regime, Round, Variant
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search.axes import AxisBook, latin_hypercube
from search.space import GENOTYPING_PARAMS, Space
from search import genetic
from common.tuning import tune
from search.climber import Climber
from search.hills import Hills
import search.engine as engine
from tests.helpers import ranges

ACCEPT_MIN, DECISIVE = tune("climb.accept_min"), tune("climb.decisive")

ROUND = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")


class Stub:
    def build(self, rnd):
        return AnswerSet(rnd.round_id, 1, Regime.GIAB, "HG002", [Variant("chr20", 5, "A", "G", "0/1")], [], 600, {})

    def update(self, a, c):
        return a


def virtual_round(cap, minutes=72, hills=1, duration=150.0, crashes=None):
    r = ranges()
    runner = FakeRunner(r, seed=0, batching=True, hills=hills)
    dur = duration if callable(duration) else (lambda job: duration)
    ex = SimExecutor(runner, cap, duration=dur, crashes=crashes or (lambda job: False))
    out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, cap), DryRunSubmitter(), r,
                    LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=minutes * 60)
    return out, ex


def setup():
    space = Space(ranges())
    book = AxisBook(space)
    names = list(space.names)
    return space, book, names


class HillsTest(unittest.TestCase):
    def test_one_rule_for_distance_ownership_and_clusters(self):
        space, book, names = setup()
        hills = Hills(space, lambda: names)
        a = space._assemble({})
        near = space._assemble({"min_pruning": 4, "assembly_region_padding": 20})
        far = space._assemble({"min_pruning": 6, "assembly_region_padding": 300, "min_mapping_quality_score": 50,
                               "pcr_indel_model": "NONE", "max_assembly_region_size": 650})
        self.assertTrue(hills.same(a, near))
        self.assertFalse(hills.same(a, far))
        self.assertEqual(hills.owner(near, {7: a}), 7)
        self.assertIsNone(hills.owner(far, {7: a}))
        picked = hills.cluster([(a, 90.0), (near, 89.9), (far, 89.8)], 6, band=0.5)
        self.assertEqual([c.key for c, _ in picked], [a.key, far.key])        # near is on a's hill


class GeneticTest(unittest.TestCase):
    def test_merge_children_are_every_combination_when_few_settings_differ(self):
        space, book, names = setup()
        a = space._assemble({})
        b = space._assemble({"min_pruning": 4, "assembly_region_padding": 300, "min_mapping_quality_score": 50})
        kids = genetic.merge_children(space, a, b, random.Random(0), names, n_mutated=0)
        self.assertEqual(len(kids), 2 ** 3 - 2)                                 # every proper subset, parents excluded
        for c in kids:
            diff = genetic.differing(space, names, a, c)
            self.assertTrue(set(diff) <= {"min_pruning", "assembly_region_padding", "min_mapping_quality_score"})
        mutated = genetic.merge_children(space, a, b, random.Random(0), names, n_mutated=20)
        self.assertGreater(len(mutated), len(kids))

    def test_merge_children_sample_when_many_settings_differ(self):
        space, book, names = setup()
        rng = random.Random(1)
        cells = latin_hypercube(space, rng, 2)
        kids = genetic.merge_children(space, cells[0], cells[1], rng, names, n_random=40, n_mutated=0)
        self.assertGreater(len(genetic.differing(space, names, cells[0], cells[1])), tune("merges.exhaustive_max"))
        self.assertTrue(10 < len(kids) <= 40)

    def test_a_generation_keeps_several_hills(self):
        space, book, names = setup()
        rng = random.Random(2)
        hills = Hills(space, lambda: names)
        pool = [(c, 90.0 - 0.01 * i) for i, c in enumerate(latin_hypercube(space, rng, 120))]
        kids = genetic.generation(space, pool, rng, 200, hills, names)
        self.assertGreaterEqual(len(kids), 150)
        self.assertEqual(len({c.key for c in kids}), len(kids))


class ClimberTest(unittest.TestCase):
    def make(self):
        space, book, names = setup()
        geno = [a for a in GENOTYPING_PARAMS if a in space.names]
        live = book.live(0.0) + geno
        d = Climber(0, space, book, space._assemble({}), 88.0, 0.0, extra_axes=geno)
        return space, d, live

    def drain(self, d, moves, score):
        for m in moves:
            d.take()
            d.consider(m, True, m.config, score(m), False, 1.0)

    def test_best_of_sweep_with_the_other_winners_reapplied(self):
        space, d, live = self.make()
        rng = random.Random(0)
        moves = d.refill(live, rng, 500)
        self.assertEqual({m.kind for m in moves}, {"n1", "n2"})
        small = next(m for m in moves if m.axes == ("assembly_region_padding",))
        big = next(m for m in moves if m.axes == ("min_mapping_quality_score",))
        for m in moves:
            d.take()
        d.consider(small, True, small.config, 88.0 + 2 * ACCEPT_MIN, False, 1.0)
        d.consider(big, True, big.config, 88.0 + 4 * ACCEPT_MIN, False, 1.0)
        self.assertEqual(d.best, 88.0)                                          # waits for the sweep
        for m in moves:
            if m is not small and m is not big:
                d.consider(m, True, m.config, 87.0, False, 1.0)
        self.assertEqual(d.centre.key, big.config.key)
        self.assertIn(("reapply", {"assembly_region_padding": small.config.params["assembly_region_padding"]}), d.queued)
        self.assertGreater(d.progress(2.0), 0)

    def test_coarse_then_fine_then_perturbation_then_exhausted(self):
        space, d, live = self.make()
        rng = random.Random(0)
        moves = d.refill(live, rng, 500)
        self.drain(d, moves, lambda m: 87.0)                                    # nothing wins anywhere
        moves = d.refill(live, rng, 500)
        self.assertEqual(d.mode, "fine")
        self.assertIn("n3", {m.kind for m in moves})
        seen_perturb = False
        for _ in range(60):
            self.drain(d, moves, lambda m: 87.0)
            moves = d.refill(live, rng, 500)
            seen_perturb = seen_perturb or d.state == "perturb"
            if d.state == "exhausted":
                break
        self.assertTrue(seen_perturb)
        self.assertEqual(d.state, "exhausted")
        self.assertEqual(d.best, 88.0)

    def test_perturbation_accepts_better_or_equal(self):
        space, d, live = self.make()
        rng = random.Random(0)
        d.mode, d.state, d.cycle, d.cycle_sent, d.fresh = "fine", "perturb", 0, 0, False
        pts = d.refill(live, rng, 30)
        self.assertTrue(pts and all(m.kind == "shake" for m in pts))
        self.drain(d, pts, lambda m: 88.0 if m is pts[0] else 87.0)            # one shaken point ties the centre
        local = d.refill(live, rng, 30)
        self.assertTrue(local and all(m.kind == "local" for m in local))
        self.drain(d, local, lambda m: 87.0)
        nxt = d.refill(live, rng, 30)                                          # equal: the centre moves (ParamILS)
        self.assertEqual(d.centre.key, pts[0].config.key)
        self.assertEqual(d.state, "climb")
        self.assertTrue(nxt)

    def test_offer_and_decisive_gains(self):
        space, d, live = self.make()
        rng = random.Random(0)
        other = space._assemble({"min_pruning": 3})
        self.assertFalse(d.offer(other, 88.0 + ACCEPT_MIN / 2, 1.0))
        self.assertTrue(d.offer(other, 88.1, 1.0))
        self.assertEqual(d.progress(2.0), 0.0)                                  # a taken point is not own progress
        d.refill(live, rng, 500)
        m = d.take()
        self.assertTrue(d.consider(m, True, m.config, 88.1 + 2 * DECISIVE, False, 2.0))

    def test_siblings_are_one_genotyping_and_one_pairhmm_step(self):
        space, d, live = self.make()
        moves = d.refill(live, random.Random(0), 500)
        pad = next(m for m in moves if m.axes == ("assembly_region_padding",))
        sibs = d.siblings(pad, random.Random(0))
        self.assertEqual(len(sibs), 2)


class EngineRoundTest(unittest.TestCase):
    def test_big_fleet_two_hills_runs_the_whole_pipeline(self):
        with self.assertLogs("search.engine", level="INFO") as logs:
            out, ex = virtual_round(950, minutes=72, hills=2)
        text = "\n".join(logs.output)
        for i in range(4):
            self.assertIn(f"wave {i} ", text)                                   # wave 0 + three generations
        self.assertIn("founding at", text)
        self.assertIn("merge 0", text)
        self.assertIn("final at", text)
        self.assertIn("plateau-centre pick", text)
        self.assertGreater(ex.utilization(), 0.9)
        self.assertGreater(out.best.score, 94.5)
        pick = [l for l in logs.output if "plateau-centre pick" in l][0].split("pick ")[1].split(" ")[0]
        self.assertEqual(out.submissions[-1][1], pick)                         # main submitted the search's pick
        for l in logs.output:                                                   # at most 3 domains from minute 40
            if " domains " in l and "t+" in l:
                t = float(l.split("t+")[1].split("m")[0])
                if 41 <= t < 61:
                    doms = l.split("domains ")[1].split(" | ")[0].split()
                    self.assertLessEqual(len([x for x in doms if not x.endswith("/m")]), 3, l)

    def test_small_fleet(self):
        with self.assertLogs("search.engine", level="INFO") as logs:
            out, ex = virtual_round(16, minutes=72)
        self.assertGreater(ex.utilization(), 0.9)
        self.assertGreater(out.best.score, 92.0)
        self.assertIn("final at", "\n".join(logs.output))

    def test_a_slow_tail_does_not_hold_the_opening(self):
        heavy = lambda job: 540.0 if int(job.config.key[:2], 16) % 10 < 3 else 150.0
        with self.assertLogs("search.engine", level="INFO") as logs:
            out, ex = virtual_round(950, minutes=30, duration=heavy)
        founding = [l for l in logs.output if "founding at" in l][0]
        self.assertLessEqual(float(founding.split("t+")[1].split("m")[0]), 8.4)  # 14 min scaled to 30 + one run (virtual time)
        self.assertGreater(ex.utilization(), 0.85)

    def test_crashes_are_counted_and_not_fatal(self):
        bad = lambda job: int(job.config.key[:2], 16) % 8 == 0
        out, ex = virtual_round(16, minutes=40, crashes=bad)
        self.assertGreater(out.failed, 0)
        self.assertEqual(sum(out.failures.values()), out.failed)
        self.assertGreater(out.best.score, 90.5)


if __name__ == "__main__":
    unittest.main()
