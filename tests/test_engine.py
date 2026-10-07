import collections
import random
import unittest

from common.types import AnswerSet, Config, Regime, Round, Variant
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search.axes import AxisBook, latin_hypercube
from search.space import GENOTYPING_PARAMS, Space
from search import genetic
from common.tuning import tune
from search.climber import Climber
from search.reductions import MNH, Reductions
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

    def test_coarse_then_wide_then_fine_then_perturbation_that_never_gives_up(self):
        space, d, live = self.make()
        rng = random.Random(0)
        moves = d.refill(live, rng, 500)
        self.drain(d, moves, lambda m: 87.0)                                    # nothing wins anywhere
        moves = d.refill(live, rng, 500)
        self.assertEqual({m.kind for m in moves}, {"wide"})                    # the all-pairs round first
        self.assertTrue(all(len(m.axes) == 2 for m in moves))
        self.drain(d, moves, lambda m: 87.0)
        moves = d.refill(live, rng, 500)
        self.assertEqual(d.mode, "fine")
        self.assertIn("n3", {m.kind for m in moves})
        seen = set()
        for _ in range(80):
            self.drain(d, moves, lambda m: 87.0)
            moves = d.refill(live, rng, 500)
            seen.add(d.state)
        self.assertIn("perturb", seen)
        self.assertNotIn("exhausted", seen)                                     # ParamILS never gives up
        self.assertGreater(d.rounds, 0)
        self.assertTrue(moves)
        self.assertEqual(d.best, 88.0)

    def test_an_exact_tie_is_not_a_move(self):
        space, d, live = self.make()
        rng = random.Random(0)
        moves = d.refill(live, rng, 500)
        centre = d.centre.key
        self.drain(d, moves, lambda m: 88.0)                                    # identical calls everywhere
        self.assertEqual(d.centre.key, centre)
        self.assertEqual({m.kind for m in d.refill(live, rng, 500)}, {"wide"})

    def test_a_new_domain_scans_every_setting_over_its_range(self):
        space, book, names = setup()
        geno = [a for a in GENOTYPING_PARAMS if a in space.names]
        live = book.live(0.0) + geno
        d = Climber(0, space, book, space._assemble({}), 88.0, 0.0, extra_axes=geno, scan=True)
        moves = d.refill(live, random.Random(0), 2000)
        self.assertEqual({m.kind for m in moves}, {"scan"})
        per_axis = {}
        for m in moves:
            per_axis[m.axes[0]] = per_axis.get(m.axes[0], 0) + 1
        self.assertGreater(len(per_axis), 10)
        self.assertLessEqual(max(per_axis.values()), 9)
        self.assertTrue(any("full-range scan" in e for e in d.events))

    def test_chains_fill_the_room_and_a_better_chain_result_moves_the_centre(self):
        space, d, live = self.make()
        rng = random.Random(0)
        sweep = d.refill(live, rng, 50)
        for m in sweep:
            d.take()                                                            # the sweep is out, nothing back
        chains = d.refill(live, rng, 200)                                       # the main line waits: chains start
        self.assertTrue(chains and all(m.kind == "chain" for m in chains))
        self.assertEqual(len(d.chains), len(chains))
        m = d.take()
        d.consider(m, True, m.config, 87.9, False, 1.0)                         # the shaken point is back
        local = [x for x in d.pool if x.kind == "chainl" and x.chain == m.chain]
        self.assertTrue(local)
        x = local[0]
        d.pool.remove(x); d.inflight += 1
        self.assertTrue(d.consider(x, True, x.config, 88.0 + 2 * ACCEPT_MIN, False, 2.0))   # beats the domain
        self.assertEqual(d.centre.key, x.config.key)

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
        self.assertTrue(d.scan_next)                                             # a gain of 0.05 arms the scan again
        self.assertEqual(d.progress(2.0), 0.0)                                  # a taken point is not own progress
        d.refill(live, rng, 500)
        m = d.take()
        self.assertEqual(m.kind, "scan")
        jumped = 88.1 + 2 * DECISIVE
        self.assertTrue(d.consider(m, True, m.config, jumped, False, 2.0))
        self.assertEqual(d.first_scan_score, jumped)                             # the first scan moved the centre
        self.assertFalse(d.first_scan_held)

    def test_a_resolved_scan_gain_moves_the_centre_without_a_rescan(self):
        space, book, names = setup()
        geno = [a for a in GENOTYPING_PARAMS if a in space.names]
        live = book.live(0.0) + geno
        d = Climber(0, space, book, space._assemble({}), 88.0, 0.0, extra_axes=geno, scan=True)
        centre = d.centre.key
        moves = d.refill(live, random.Random(0), 500)
        self.assertTrue(moves and all(m.kind == "scan" for m in moves))
        head = d.take()
        self.assertFalse(d.consider(head, True, head.config, 88.0 + 2 * ACCEPT_MIN, False, 1.0))
        self.assertEqual(d.centre.key, centre)                                  # under 0.05 does not abort the scan
        self.assertFalse(d.offer(space._assemble({"min_pruning": 3}), 88.02, 1.0))
        for m in moves:
            if m is head:
                continue
            d.take()
            d.consider(m, True, m.config, 88.0, False, 1.0)
        self.assertFalse(d.scan_flat)
        self.assertFalse(d.scan_next)                                            # 0.05 is what restarts the full scan
        self.assertEqual(d.centre.key, head.config.key)
        self.assertAlmostEqual(d.first_scan_score, 88.0 + 2 * ACCEPT_MIN)
        self.assertFalse(d.first_scan_held)
        nxt = d.refill(live, random.Random(1), 50)
        self.assertTrue(nxt)
        self.assertNotIn("scan", {m.kind for m in nxt})

    def test_a_scan_under_the_answer_set_resolution_holds_the_peak(self):
        space, book, names = setup()
        geno = [a for a in GENOTYPING_PARAMS if a in space.names]
        live = book.live(0.0) + geno
        d = Climber(0, space, book, space._assemble({}), 88.0, 0.0, extra_axes=geno, scan=True)
        centre = d.centre.key
        moves = d.refill(live, random.Random(0), 500)
        for m in moves:
            d.take()
            d.consider(m, True, m.config, 88.0 + ACCEPT_MIN / 2, False, 1.0)
        self.assertTrue(d.scan_flat)
        self.assertEqual(d.centre.key, centre)
        self.assertTrue(d.first_scan_held)
        self.assertEqual(d.refill(live, random.Random(1), 50), [])

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
        self.assertIn("wave 0 ", text)
        self.assertIn("bred into idle workers", text)                         # opening.steady: no generations
        self.assertIn("founding at", text)
        self.assertIn("6 domains, 6 held", text)                              # six climb; six candidates wait
        self.assertIn("merge (final)", text)                                  # one merge, at the final window
        self.assertNotIn("scheduled 3 to 1", text)
        self.assertNotIn("one hill", text)                                     # merges.on_evidence 0
        self.assertIn("final at", text)
        self.assertIn("pick ", text)
        self.assertNotIn("plateau-centre", text)
        self.assertGreater(ex.utilization(), 0.85)                    # children of the final merge go out before the next sweep
        self.assertGreater(out.best.score, 94.5)
        pick = [l for l in logs.output if "engine: pick " in l][0].split("pick ")[1].split(" ")[0]
        self.assertEqual(out.submissions[-1][1], pick)                         # main submitted the search's pick
        still_separate = False
        for l in logs.output:                                                   # the basins are still separate at minute 51-61
            if " domains " in l and "t+" in l:
                t = float(l.split("t+")[1].split("m")[0])
                if 51 <= t < 61:
                    doms = l.split("domains ")[1].split(" | ")[0].split()
                    if len([x for x in doms if not x.endswith("/m")]) > 1:
                        still_separate = True
        self.assertTrue(still_separate)

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
        self.assertLessEqual(float(founding.split("t+")[1].split("m")[0]), 8.8)  # 15 min scaled to 30 + one run (virtual time)
        self.assertGreater(ex.utilization(), 0.85)

    def test_crashes_are_counted_and_not_fatal(self):
        bad = lambda job: int(job.config.key[:2], 16) % 8 == 0
        out, ex = virtual_round(16, minutes=40, crashes=bad)
        self.assertGreater(out.failed, 0)
        self.assertEqual(sum(out.failures.values()), out.failed)
        self.assertGreater(out.best.score, 90.5)


class knobs:
    """Settings in force inside a `with` block (tests otherwise run on defaults)."""
    def __init__(self, **kv):
        self.kv = {k.replace("__", "."): v for k, v in kv.items()}

    def __enter__(self):
        from common.tuning import TUNING
        self.saved = TUNING.values
        TUNING.values = {**TUNING.values, **self.kv}

    def __exit__(self, *exc):
        from common.tuning import TUNING
        TUNING.values = self.saved


class V41PathTest(unittest.TestCase):
    def test_generations_and_recombination_merges_still_run_with_the_switches_off(self):
        with knobs(opening__steady=0, opening__islands=0, merges__retire=0, merges__on_evidence=1):
            with self.assertLogs("search.engine", level="INFO") as logs:
                out, ex = virtual_round(950, minutes=72, hills=2)
        text = "\n".join(logs.output)
        for i in range(4):
            self.assertIn(f"wave {i} ", text)                                   # wave 0 + three generations
        self.assertIn("merge 0", text)
        self.assertGreater(out.best.score, 94.5)


class V42Test(unittest.TestCase):
    def test_islands_never_mate_across_and_keep_off_each_others_hills(self):
        space, book, names = setup()
        rng = random.Random(3)
        hills = Hills(space, lambda: names)
        pool = [(c, 90.0 - 0.01 * i) for i, c in enumerate(latin_hypercube(space, rng, 200))]
        island_of = {}
        kids = genetic.island_generation(space, pool, island_of, 6, rng, 300, hills, names, min_far=6)
        heads = {i: next(c for c, _ in pool if island_of.get(c.key) == i) for i in set(island_of.values())}
        self.assertEqual(len(heads), 6)
        per = collections.Counter(i for _, i in kids)
        self.assertLessEqual(max(per.values()) - min(per.values()), 5)                 # equal shares
        for c, i in kids:
            self.assertFalse(any(j != i and hills.same(c, h) for j, h in heads.items()))

    def test_spaced_founding_keeps_hills_apart(self):
        space, book, names = setup()
        hills = Hills(space, lambda: names)
        a = space._assemble({})
        five = space._assemble({"min_pruning": 6, "assembly_region_padding": 300, "min_mapping_quality_score": 50,
                                "pcr_indel_model": "NONE", "max_assembly_region_size": 650})
        picked = hills.cluster([(a, 90.0), (five, 89.9)], 6, band=0.5, min_far=hills.distance(a, five) + 1)
        self.assertEqual([c.key for c, _ in picked], [a.key])
        more = hills.take([(five, 89.0), (a, 90.0)], 1, min_far=4, floor=80.0, avoid=[a])
        self.assertEqual([c.key for c, _ in more], [five.key])

    def test_full_reapply_adds_running_sums_and_pairs(self):
        with knobs(climb__reapply_full=1):
            space, d, live = ClimberTest().make()
            moves = d.refill(live, random.Random(0), 500)
            for m in moves:
                d.take()
            singles = [m for m in moves if m.kind == "n1" and len(m.axes) == 1]
            axes, wins = set(), []
            for m in singles:
                if m.axes[0] not in axes and len(wins) < 4:
                    axes.add(m.axes[0]); wins.append(m)
            for k, m in enumerate(wins):
                d.consider(m, True, m.config, 88.0 + (k + 2) * ACCEPT_MIN, False, 1.0)
            for m in moves:
                if m not in wins:
                    d.consider(m, True, m.config, 87.0, False, 1.0)
            kinds = [ch for kind, ch in d.queued if kind == "reapply"]
            self.assertEqual(len(kinds), 3 + 2 + 3)                    # 3 alone, 2 running sums, 3 pairs

    def test_few_domains_send_the_wide_round_with_the_first_coarse_sweep(self):
        space, d, live = ClimberTest().make()
        d.eager_wide = True
        moves = d.refill(live, random.Random(0), 2000)
        kinds = collections.Counter(m.kind for m in moves)
        self.assertGreater(kinds["wide"], 0)
        self.assertGreater(kinds["n1"], 0)

    def test_retire_merge_hands_over_settings_as_transplants(self):
        space, d, live = ClimberTest().make()
        other = space._assemble({"min_pruning": 6, "assembly_region_padding": 300, "min_mapping_quality_score": 50})
        n = d.transplant(other, live)
        self.assertEqual(n, 3 + 3)
        self.assertEqual(sum(1 for kind, _ in d.queued if kind == "transplant"), 6)

    def test_every_variant_finishes_a_virtual_round(self):
        for kv in ({"opening.islands": 1, "domains.found_min_far": 6}, {"opening.steady": 1},
                   {"merges.retire": 1, "climb.reapply_full": 1}):
            with knobs(**{k.replace(".", "__"): v for k, v in kv.items()}):
                out, ex = virtual_round(200, minutes=40, hills=2)
                self.assertGreater(out.best.score, 90.0, kv)


class ReductionsTest(unittest.TestCase):
    def test_a_finished_run_closes_higher_caps_and_a_flat_calling_threshold(self):
        space = Space(ranges())
        red = Reductions()
        space.reductions = red
        base = space._assemble({}).params
        ran = space._assemble({**base, MNH: 128})
        red.note(ran, 80.0, haplotype_cap=40)
        self.assertEqual(space._assemble({**base, MNH: 100}).params[MNH], 128)
        self.assertEqual(space._assemble({**base, MNH: 40}).params[MNH], 40)
        self.assertEqual(space._assemble({**base, MNH: 16}).params[MNH], 16)
        red.note(space._assemble({**base, MNH: 30}), 80.0, haplotype_cap=10)
        self.assertEqual(space._assemble({**base, MNH: 35}).params[MNH], 35)
        # A needed-0 report names a small cap. An in-flight higher cap can still report a real floor afterwards;
        # that report is the config that was submitted, and it must not keep rewriting down to the small cap.
        red2 = Reductions()
        space.reductions = red2
        low = space._assemble({**base, MNH: 16})
        red2.note(low, 0.0, haplotype_cap=0)
        self.assertEqual(space._assemble({**base, MNH: 100}).params[MNH], 16)
        high = Config(low.tool, {**low.params, MNH: 128})
        red2.note(high, 80.0, haplotype_cap=40)
        self.assertEqual(space._assemble({**base, MNH: 100}).params[MNH], 128)
        self.assertEqual(space._assemble({**base, MNH: 40}).params[MNH], 40)
        self.assertEqual(space._assemble({**base, MNH: 16}).params[MNH], 16)
        alt = "max_alternate_alleles"
        self.assertEqual(space._assemble({**base, alt: 2}).params[alt], 2)
        self.assertEqual(space._assemble({**base, alt: 12}).params[alt], 12)
        # Equal scores at 70 and 90 do not close 80: the maximum can sit between two equal ends.
        calling = "standard_min_confidence_threshold_for_calling"
        red.note(space._assemble({**base, calling: 70.0}), 83.0)
        red.note(space._assemble({**base, calling: 90.0}), 83.0)
        self.assertEqual(space._assemble({**base, calling: 80.0}).params[calling], 80.0)

if __name__ == "__main__":
    unittest.main()
