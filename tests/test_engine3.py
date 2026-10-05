import random
import re
import unittest

from common.types import AnswerSet, Config, Regime, Round, Variant
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search.axes import AxisBook
from search.domain import ACCEPT_MIN, DECISIVE, SHAKE_K, SHAKE_POINTS, SIDEWAYS_MAX, Domain, Move
from search.space import GENOTYPING_PARAMS, Space
import search.engine3 as engine3
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
    dur = duration if callable(duration) else (lambda job: duration)
    ex = SimExecutor(runner, cap, duration=dur, crashes=crashes or (lambda job: False))
    out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine3.searcher(r, cap), DryRunSubmitter(), r,
                    LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=minutes * 60)
    return out, ex


def domain(score=88.0):
    space = Space(ranges())
    book = AxisBook(space)
    centre = space._assemble({})
    return space, book, Domain(0, space, book, centre, score, 0.0)


def live_axes(book, now=0.0):
    return book.live(now) + [a for a in GENOTYPING_PARAMS if a in book.space.names]


class DomainTest(unittest.TestCase):
    def test_a_new_centre_sends_its_coarse_sweep_and_no_shake(self):
        space, book, d = domain()
        live = live_axes(book)
        moves = d.refill(live, random.Random(0), target=120)
        self.assertEqual({m.kind for m in moves}, {"n1", "n2"})                 # no injection: the hypercube covers the ranges
        n1 = [m for m in moves if m.kind == "n1"]
        self.assertTrue(all(len(m.axes) == 1 for m in n1))
        pcr = [m for m in n1 if m.axes == ("pcr_indel_model",)]
        self.assertEqual({m.config.params["pcr_indel_model"] for m in pcr}, {"NONE", "HOSTILE", "AGGRESSIVE"})
        self.assertTrue(all(len(m.axes) == 2 for m in moves if m.kind == "n2"))
        self.assertIn(("standard_min_confidence_threshold_for_calling",), [m.axes for m in n1])   # genotyping climbs too
        self.assertEqual(len({m.config.key for m in moves}), len(moves))
        self.assertGreater(len(moves), 100)                        # three scales of single steps plus the pairs

    def test_small_gains_wait_for_the_sweep_and_the_best_wins_with_the_others_merged(self):
        space, book, d = domain()
        rng = random.Random(0)
        live = live_axes(book)
        moves = d.refill(live, rng)
        small = next(m for m in moves if m.axes == ("assembly_region_padding",))
        big = next(m for m in moves if m.axes == ("min_mapping_quality_score",))
        for m in moves:
            d.take()
        self.assertFalse(d.consider(small, True, small.config, 88.0 + 2 * ACCEPT_MIN, False, 1.0))   # waits
        self.assertFalse(d.consider(big, True, big.config, 88.0 + 4 * ACCEPT_MIN, False, 1.0))
        self.assertEqual(d.centre.key, space._assemble({}).key)
        moved = False
        for m in moves:
            if m.config.key in (small.config.key, big.config.key):
                continue
            moved = d.consider(m, True, m.config, 87.0, False, 2.0) or moved
        self.assertTrue(moved)                                                        # 80% back: the best of the sweep
        self.assertEqual(d.centre.key, big.config.key)
        self.assertEqual(d.recent, ["min_mapping_quality_score"])
        self.assertTrue(d.merges and all(m.kind == "merge" for m in d.merges))        # the other winner, re-applied
        self.assertIn(("assembly_region_padding",), [m.axes for m in d.merges])
        nxt = d.refill(live, rng)
        self.assertEqual(nxt[0].kind, "merge")
        d.take()
        self.assertTrue(d.consider(nxt[0], True, nxt[0].config, 88.0 + 2 * DECISIVE, False, 3.0))   # decisive: at once
        self.assertEqual(d.centre.key, nxt[0].config.key)

    def test_steps_double_on_equal_in_coarse_and_halve_on_worse(self):
        space, book, d = domain()
        live = live_axes(book)
        moves = d.refill(live, random.Random(0))
        pad = next(m for m in moves if m.axes == ("assembly_region_padding",))
        step0 = d.steps["assembly_region_padding"]
        d.take()
        d.consider(pad, True, pad.config, 88.0, False, 1.0)
        self.assertAlmostEqual(d.steps["assembly_region_padding"], step0 * 2)
        self.assertEqual([c.key for c in d.equals], [pad.config.key])
        mq = next(m for m in moves if m.axes == ("min_mapping_quality_score",))
        s0 = d.steps["min_mapping_quality_score"]
        d.take()
        d.consider(mq, True, mq.config, 87.0, False, 2.0)
        self.assertLess(d.steps["min_mapping_quality_score"], s0)

    def test_exhausted_sweep_walks_sideways_then_shakes_with_growing_k(self):
        space, book, d = domain()
        rng = random.Random(0)
        live = live_axes(book)

        def drain(moves, score):
            for m in moves:
                d.take()
                d.consider(m, True, m.config, score(m), False, 1.0)

        moves = d.refill(live, rng, target=0)
        drain(moves, lambda m: 88.0 if m.kind in ("n1", "n2") else 87.0)      # a plateau (probes lose)
        moves = d.refill(live, rng)                                            # -> sideways onto an equal neighbour
        self.assertEqual(d.sideways, 1)
        self.assertEqual({m.kind for m in moves}, {"n1", "n2"})
        for _ in range(SIDEWAYS_MAX):
            drain(moves, lambda m: 88.0 if m.kind in ("n1", "n2") else 87.0)
            moves = d.refill(live, rng)
        self.assertEqual(d.sideways, SIDEWAYS_MAX)                             # the plateau budget is spent: shake
        self.assertTrue(moves and all(m.kind == "shake" for m in moves))
        k = d.shake.k
        self.assertEqual(k, SHAKE_K[0])
        drain(moves, lambda m: 87.0 + 0.01 * len(m.axes))                      # the points come back, one is best
        local = d.refill(live, rng)                                            # neighbourhood 1 around that best
        self.assertTrue(local and all(m.kind == "local" for m in local))
        self.assertEqual(local[0].centre_key, d.shake.best[0].key)
        drain(local, lambda m: 87.0)
        nxt = d.refill(live, rng)                                              # nothing beat the best: next k
        self.assertEqual(d.shake.k, SHAKE_K[1])
        self.assertTrue(nxt and all(m.kind == "shake" for m in nxt))
        self.assertEqual(d.best, 88.0)
        win = nxt[0]
        d.take()
        self.assertTrue(d.consider(win, True, win.config, 88.2, False, 5.0))    # a shake point above the best: move there
        self.assertEqual(d.centre.key, win.config.key)

    def test_a_flat_batch_switches_to_fine_mode_and_exhaustion_at_the_floor_freezes(self):
        space, book, d = domain()
        rng = random.Random(0)
        live = live_axes(book)
        steps0 = dict(d.steps)
        d.completed = 300
        d.close_batch()                                                        # gain 0 < FINE_GAIN
        self.assertEqual(d.mode, "fine")
        self.assertFalse(d.frozen)
        self.assertTrue(all(d.steps[a] <= steps0[a] for a in steps0))
        moves = d.refill(live, rng)
        self.assertEqual({m.kind for m in moves}, {"n1", "n2", "n3"})          # fine: singles, pairs and triples
        pad = next(m for m in moves if m.axes == ("assembly_region_padding",))
        st = d.steps["assembly_region_padding"]
        d.take()
        d.consider(pad, True, pad.config, 88.0, False, 1.0)                     # equal in fine mode: no doubling
        self.assertEqual(d.steps["assembly_region_padding"], st)
        self.assertIn("assembly_region_padding", d.flat_axes)
        for m in moves:
            if m is not pad:
                d.take(); d.consider(m, True, m.config, 87.0, False, 2.0)
        d.sideways = SIDEWAYS_MAX                                              # no plateau budget left
        again = d.refill(live, rng)                                            # exhausted: the sweep at half the scale
        self.assertTrue(again)
        self.assertEqual(d.fine_scale, 0.5)
        self.assertNotIn(("assembly_region_padding",), [m.axes for m in again if m.kind == "n1"])
        d.close_batch()
        self.assertFalse(d.frozen)                                             # flat batches alone do not freeze
        for _ in range(12):                                                    # every finer scale, all losing
            for m in again:
                d.take(); d.consider(m, True, m.config, 87.0, False, 3.0)
            d.sideways = SIDEWAYS_MAX
            again = d.refill(live, rng)
            if d.frozen:
                break
        self.assertTrue(d.frozen)                                              # exhausted down to every axis's floor
        self.assertTrue(d.exhausted_at_floor)

    def test_a_domain_without_a_visible_gain_for_12_minutes_is_stalled(self):
        from search.domain import RESOLUTION, STALL_S
        space, book, d = domain()
        self.assertFalse(d.stalled(STALL_S - 1))
        self.assertTrue(d.stalled(STALL_S + 1))
        d._recentre(space._assemble({"min_pruning": 3}), 88.0 + RESOLUTION / 2, 100.0, ("min_pruning",))   # below resolution
        self.assertTrue(d.stalled(100.0 + STALL_S + 1))
        d._recentre(space._assemble({"min_pruning": 4}), 88.0 + RESOLUTION, 200.0, ("min_pruning",))       # adds up: real
        self.assertFalse(d.stalled(200.0 + STALL_S - 1))

    def test_siblings_are_cheap_variations_of_a_structural_move_only(self):
        space, book, d = domain()
        rng = random.Random(0)
        live = live_axes(book)
        moves = d.refill(live, rng)
        pad = next(m for m in moves if m.axes == ("assembly_region_padding",))
        sibs = d.siblings(pad, rng)
        self.assertEqual(len(sibs), 2)
        for c in sibs:
            changed = [a for a in space.names if c.params.get(a) != pad.config.params.get(a)]
            self.assertEqual(len(changed), 1)
            self.assertIn(changed[0], GENOTYPING_PARAMS + ("base_quality_score_threshold", "pcr_indel_model",
                                                            "pair_hmm_gap_continuation_penalty", "phred_scaled_global_read_mismapping_rate"))
        conf = next(m for m in moves if m.axes == ("standard_min_confidence_threshold_for_calling",))
        self.assertEqual(d.siblings(conf, rng), [])

    def test_inherit_queues_the_losers_differing_axes_as_transplants(self):
        space, book, a = domain()
        b = Domain(1, space, book, space._assemble({"assembly_region_padding": 30, "min_mapping_quality_score": 40,
                                                    "min_pruning": 3}), 87.9, 0.0)
        b.tried.add("somekey"); b.steps["assembly_region_padding"] = 0.02
        live = live_axes(book)
        n = a.inherit(b, live)
        self.assertGreaterEqual(n, 5)                                          # three singles and their pairs
        self.assertNotIn(b.centre.key, [m.config.key for m in a.pool])        # the loser's centre itself was tested
        self.assertTrue(all(m.kind == "transplant" for m in a.pool[:n]))
        self.assertIn(("assembly_region_padding",), [m.axes for m in a.pool])
        self.assertIn(("assembly_region_padding", "min_mapping_quality_score"), [m.axes for m in a.pool])
        self.assertIn("somekey", a.tried)
        self.assertEqual(a.steps["assembly_region_padding"], 0.02)

    def test_a_recentre_during_a_bracket_still_closes_the_batch(self):
        space, book, d = domain()
        sent = d.start_bracket(live_axes(book))
        first = d.take()
        d.consider(first, True, first.config, 88.5, False, 1.0)           # a bracket point wins: the pool is cleared
        self.assertEqual(d.centre.key, first.config.key)
        self.assertEqual(d.pool, [])
        self.assertIsNotNone(d.bracket)
        self.assertGreaterEqual(d.bracket[1], 0.8 * d.bracket[0])

    def test_final_phase_keeps_tightening_around_one_point(self):
        space, book, d = domain()
        rng = random.Random(0)
        live = live_axes(book)
        n = d.start_final(live, rng, target=60)
        self.assertGreater(n, 0)
        self.assertEqual({m.kind for m in d.pool}, {"final"})                  # no filler: the scout takes spare workers
        scale0 = d.final_scale
        for m in list(d.pool):
            d.take(); d.consider(m, True, m.config, 87.0, False, 1.0)
        more = d.refill(live, rng, 60)
        self.assertTrue(more)
        self.assertLess(d.final_scale, scale0)


class HypercubeTest(unittest.TestCase):
    def test_every_axis_is_stratified_across_its_range(self):
        from search.axes import latin_hypercube
        space = Space(ranges())
        cells = latin_hypercube(space, random.Random(0), 200)
        self.assertGreaterEqual(len(cells), 190)
        for name in ("assembly_region_padding", "min_mapping_quality_score", "standard_min_confidence_threshold_for_calling"):
            s = space.specs[name]
            us = sorted(space._to_unit(s, c.params[name]) for c in cells)
            self.assertLess(us[0], 0.05); self.assertGreater(us[-1], 0.95)
            tenths = {min(9, int(u * 10)) for u in us}
            self.assertEqual(len(tenths), 10)                                    # every tenth of the range is visited
        self.assertEqual({c.params["pcr_indel_model"] for c in cells}, {"NONE", "HOSTILE", "AGGRESSIVE", "CONSERVATIVE"})
        self.assertEqual(len({c.key for c in cells}), len(cells))


class EngineRoundTest(unittest.TestCase):
    def test_small_fleet_climbs_all_round_and_fills_the_workers(self):
        with self.assertLogs("search.engine3", level="INFO") as logs:
            out, ex = virtual_round(16)
        self.assertEqual(out.failed, 0)
        self.assertGreater(ex.utilization(), 0.9)
        self.assertGreater(out.best.score, 91.0)
        self.assertGreater(out.evaluated, 400)                     # no riders after the opening: every config is a move
        self.assertIn("final at", "\n".join(logs.output))

    def test_big_fleet_opens_six_domains_and_funnels_them(self):
        with self.assertLogs("search.engine3", level="INFO") as logs:
            out, ex = virtual_round(950, minutes=40)
        text = "\n".join(logs.output)
        self.assertIn("up to 6 domains", text)
        self.assertIn("hypercube cells", text)
        self.assertIn("wave 2 at", text)
        self.assertGreaterEqual(text.count(" born at "), 3)
        self.assertIn("(a scout find)", text)                                  # rolling domains from the scout
        self.assertIn("funnel at", text)
        funnel = [l for l in logs.output if "funnel at" in l and " keeps " in l][0]
        kept = [float(x.split(":")[1].split("/")[0]) for x in funnel.split("keeps ")[1].split()]
        best_before = max(float(m) for l in logs.output[:logs.output.index(funnel)]
                          for m in re.findall(r"-> ([\d.]+) via", l))
        self.assertGreaterEqual(max(kept), best_before - 1e-6)                # the holder of the best is kept
        self.assertGreater(ex.utilization(), 0.85)
        self.assertGreater(out.evaluated, 8000)
        self.assertGreater(out.best.score, 91.0)

    def test_two_hills_keep_the_leader_on_the_best(self):
        with self.assertLogs("search.engine3", level="INFO") as logs:
            out, ex = virtual_round(950, minutes=40, hills=2)
        self.assertGreater(out.best.score, 91.0)
        last = [l for l in logs.output if " domains " in l][-1]
        leader = max(float(x.split(":")[1].split("/")[0]) for x in last.split("domains ")[1].split(" | ")[0].split())
        self.assertGreaterEqual(leader, out.best.score - 0.3)

    def test_a_slow_tail_does_not_hold_the_timed_steps(self):
        """30% of the configs run 9 minutes (the hypercube's heavy cells): wave 2 still goes out on time, the
        idle workers are used, the opening comes when due."""
        heavy = lambda job: 540.0 if int(job.config.key[:2], 16) % 10 < 3 else 150.0
        with self.assertLogs("search.engine3", level="INFO") as logs:
            out, ex = virtual_round(950, minutes=20, duration=heavy)
        wave2 = [l for l in logs.output if "wave 2 at" in l][0]
        self.assertLessEqual(float(wave2.split("t+")[1].split("m")[0]), 5.1)
        opening = [l for l in logs.output if "opening done" in l][0]
        self.assertLessEqual(float(opening.split("t+")[1].split("m")[0]), 11.0)   # virtual time moves by result arrivals
        self.assertGreater(ex.utilization(), 0.8)

    def test_crashes_are_counted_and_not_fatal(self):
        bad = lambda job: int(job.config.key[:2], 16) % 8 == 0
        out, ex = virtual_round(16, minutes=40, crashes=bad)
        self.assertGreater(out.failed, 0)
        self.assertEqual(sum(out.failures.values()), out.failed)
        self.assertGreater(out.best.score, 90.5)


if __name__ == "__main__":
    unittest.main()
