"""search engine v9: ladders, callset distance, Step A's groups, Step B's grid and growth, Step C's children, and a
whole virtual round."""
import asyncio
import itertools
import gzip
import os
import random
import tempfile
import unittest
import unittest.mock

from common.types import AnswerSet, Config, Regime, Round, Variant
from evaluator import FakeRunner, FakeScorer, SimExecutor
from evaluator.sketch import callset_sketch
from main import DryRunSubmitter, LoopSettings, run_round
from search import engine, far, lines
from search.candidates import Point, StepA, pareto_layers
from search.climb import (Climber, adopt_ga_best, far_replacement, far_scouts, held_by_climber, held_by_ga,
                          next_probes, submit_ga)
from search.repeats import TiedCallsets
from search.lines import single
from search.distance import callset_distance, setting_distance
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges

ROUND = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")
MDBL = "min_dangling_branch_length"
MP = "min_pruning"


class Stub:
    def build(self, rnd):
        return AnswerSet(rnd.round_id, 1, Regime.GIAB, "HG002", [Variant("chr20", 5, "A", "G", "0/1")], [], 600, {})

    def update(self, a, c):
        return a


class FakeCtx:
    """What a Climber and Step C need, without a client: scores come from a function of the config."""

    def __init__(self, space, f):
        self.space, self.f = space, f
        self.ladders = Ladders(space)
        self.rng = random.Random(0)
        self.sent = []
        self.climbers = []
        self.now = 0.0
        self.cap = 500
        self.c_log = []

    def over(self):
        return False

    def scored(self):
        seen = set()
        for batch in self.sent:
            for c in batch:
                if c.key not in seen:
                    seen.add(c.key)
                    yield c, self.f(c)

    def clock(self):
        return self.now

    def elapsed(self):
        return self.now

    def say(self, *a):
        pass

    def fresh(self, c):
        return c is not None

    def score(self, c):
        return self.f(c)

    def owner(self, c, me):
        return None

    def propose(self, batch, tag, batches, priority=0):
        self.sent.append(list(batch))
        loop = asyncio.get_event_loop()
        out = []
        for _ in batch:
            fut = loop.create_future()
            fut.set_result(None)
            out.append(fut)
        return out

    async def wait_share(self, futs, share, until):
        self.now += 60.0


def vcf(path, records):
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for pos, ref, alt, gt in records:
            f.write(f"chr20\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n")


class LadderTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.lad = Ladders(self.space)

    def test_integers_move_by_one_and_stop_at_the_ends(self):
        n = "min_pruning"
        lo, hi = self.lad.rungs[n][0], self.lad.rungs[n][-1]
        self.assertEqual(self.lad.shift(n, lo + 2, 1), lo + 3)
        self.assertEqual(self.lad.shift(n, lo + 2, -2), lo)
        self.assertIsNone(self.lad.shift(n, hi, 1))

    def test_an_off_ladder_float_moves_to_the_next_rung(self):
        n = "active_probability_threshold"
        self.assertAlmostEqual(self.lad.shift(n, 0.0025, 1), 0.003)
        self.assertAlmostEqual(self.lad.shift(n, 0.0025, -1), 0.002)

    def test_log_settings_have_ten_rungs_a_decade(self):
        r = self.lad.rungs["heterozygosity"]
        self.assertAlmostEqual(r[10] / r[0], 10.0, places=2)


class DistanceTest(unittest.TestCase):
    def test_identical_callsets_are_zero_and_disjoint_ones_one(self):
        with tempfile.TemporaryDirectory() as d:
            a, b, c = (os.path.join(d, f"{x}.vcf.gz") for x in "abc")
            recs = [(100 + i, "A", "G", "0/1") for i in range(50)] + [(500, "AT", "A", "1/1")]
            vcf(a, recs)
            vcf(b, recs)
            vcf(c, [(9000 + i, "C", "T", "0/1") for i in range(50)] + [(9500, "GT", "G", "1/1")])
            sa, sb, sc = callset_sketch(a), callset_sketch(b), callset_sketch(c)
        self.assertEqual(sa["n"]["snp_het"], 50)
        self.assertEqual(sa["n"]["indel_1bp"], 1)
        self.assertEqual(callset_distance(unpacked(sa), unpacked(sb)), 0.0)
        self.assertEqual(callset_distance(unpacked(sa), unpacked(sc)), 1.0)

    def test_setting_distance(self):
        space = Space(ranges())
        d = space._assemble({})
        self.assertEqual(setting_distance(space, d, d), 0.0)
        far_ = space._assemble({n: space.specs[n].max for n in ("min_pruning", "max_reads_per_alignment_start")})
        self.assertGreater(setting_distance(space, d, far_), 0.0)


def unpacked(sk):
    from search.distance import _unpack
    return {"k": sk["k"], "h": {c: _unpack(v) for c, v in sk["h"].items()}}


class CandidatesTest(unittest.TestCase):
    def test_pareto_layers_take_the_front_first(self):
        c = Config("gatk", {})
        pts = [Point(c, 90, 0.95, 400), Point(c, 91, 0.94, 100), Point(c, 80, 0.90, 300), Point(c, 85, 0.96, 900)]
        front = pareto_layers(pts, 3)
        self.assertEqual(sorted(p.score for p in front), [85, 90, 91])      # 80 is dominated by 91


class ClimbTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.target = {"min_pruning": 7, "assembly_region_padding": 60}

    def peak(self, c):
        return 100.0 - sum(abs(c.params[n] - v) for n, v in self.target.items())

    def climber(self):
        ctx = FakeCtx(self.space, self.peak)
        start = self.space._assemble({"min_pruning": 2, "assembly_region_padding": 50})
        return ctx, Climber(ctx, 0, start, self.peak(start))

    def test_a_hot_ga_sends_no_climber_grid(self):
        from search.climb import LINE_COARSE, Climber, ga_hot, leftover_breeds, line_coarse
        ctx, cl = self.climber()
        ctx.T = lambda seconds: seconds
        ctx.ga_moved = 0.0
        ctx.now = 60.0
        self.assertTrue(ga_hot(ctx))
        if leftover_breeds():
            self.assertEqual(line_coarse(ctx), LINE_COARSE)
            self.assertTrue(cl.begin_grid(1e9))
            self.assertEqual(cl.state, "climbing")
            self.assertTrue(cl.pending)
        else:
            self.assertEqual(line_coarse(ctx), 0)
            self.assertFalse(cl.begin_grid(1e9))
            self.assertEqual(cl.state, "climbing")
            self.assertEqual(cl.pending, [])
            self.assertEqual(cl.iteration, 0)
        ctx.now = 4 * 60
        self.assertFalse(ga_hot(ctx))
        self.assertEqual(line_coarse(ctx), LINE_COARSE)
        settled = Climber(ctx, 1, cl.centre, cl.score)
        self.assertTrue(settled.begin_grid(1e9))
        self.assertTrue(settled.later)
        fine = Climber(ctx, 2, cl.centre, cl.score, origin="step C")
        ctx.now = 60.0
        self.assertTrue(ga_hot(ctx))
        self.assertTrue(fine.begin_grid(1e9))
        self.assertEqual(fine.later, [])
        self.assertTrue(fine.pending)

    def test_a_free_worker_walks_further_up_the_line(self):
        ctx, cl = self.climber()
        ctx.ga_moved = None
        self.assertTrue(cl.begin_grid(1e9))
        cl.pending = []
        added = cl.extend(40)
        self.assertGreater(added, 0)
        cl.submit_one()
        sent = [c for batch in ctx.sent for c in batch]
        self.assertEqual(len(sent), 1)
        self.assertEqual(len(cl.extra), added - 1)

    def test_step_b_runs_until_the_fine_tune(self):
        from search.far import fine_tune_at
        ctx, _cl = self.climber()
        ctx.t0 = 0.0
        ctx.length = 72 * 60
        ctx.scale = 1.0
        ctx.T = lambda seconds: seconds * ctx.scale
        ctx.pick_at = lambda: ctx.t0 + ctx.length - 90.0
        self.assertAlmostEqual(fine_tune_at(ctx), ctx.pick_at() - 8 * 60)

    def test_step_c_fine_tunes_the_best_and_does_not_breed(self):
        ctx, cl = self.climber()
        bred = []

        class Step:
            async def breed(self, *a, **k):
                bred.append("breed")
                return [cl.centre]

            async def wave(self, *a, **k):
                bred.append("wave")

        ctx.step_a = Step()
        ctx.t0 = 0.0
        ctx.length = 72 * 60
        ctx.scale = 1.0
        ctx.T = lambda seconds: seconds * ctx.scale
        ctx.pick_at = lambda: ctx.t0 + ctx.length - 90.0
        ctx.now = 54 * 60
        ctx.best = lambda: (cl.centre, cl.score)
        ctx.phase = ""
        ctx.publish = lambda force=False: None
        ran = []

        async def fake_run(self, until, iterations=None, jumps=True):
            ran.append((self.origin, self.centre.key, until))
            self.state = "peak"

        with unittest.mock.patch.object(Climber, "run", fake_run):
            asyncio.run(far.refine(ctx))
        self.assertEqual(bred, [])
        self.assertEqual(ctx.phase, "C: fine-tune")
        self.assertEqual(ran, [("step C", cl.centre.key, ctx.pick_at())])

    def test_the_fine_tune_restarts_from_a_higher_score(self):
        ctx, cl = self.climber()
        cl.origin = "step C"
        ctx.climbers = [cl]
        start = cl.centre
        higher = self.space._assemble({"min_pruning": 9, "assembly_region_padding": 90})
        offered = self.peak(start) + 5
        ctx.client = type("C", (), {"idle_slots": staticmethod(lambda: 400)})()
        cl.flat.update(["min_pruning"])
        low = self.space._assemble({"min_pruning": 3})
        ctx.best = lambda: (low, cl.score - 1)
        self.assertFalse(cl.follow_best())
        self.assertEqual(cl.centre.key, start.key)
        ctx.best = lambda: (cl.centre, cl.score) if cl.history else (higher, offered)
        asyncio.run(cl.run(1e9))
        self.assertEqual(cl.history[0], (0.0, round(offered, 4), "best"))
        self.assertTrue(ctx.sent)
        moved = [n for n in cl.all_names if ctx.sent[0][0].params.get(n) != higher.params.get(n)]
        self.assertTrue(moved)
        self.assertLessEqual(len(moved), 5)

    def test_a_live_grid_keeps_the_centre_when_a_higher_score_arrives(self):
        ctx, cl = self.climber()
        cl.origin = "step C"
        cl.batch_live = True
        cl.pending = [cl.centre]
        start = cl.centre.key
        higher = self.space._assemble({"min_pruning": 9, "assembly_region_padding": 90})
        ctx.best = lambda: (higher, cl.score + 5)
        ctx.client = type("C", (), {"idle_slots": staticmethod(lambda: 0)})()
        ctx.climbers = [cl]
        called = []
        cl.follow_best = lambda: called.append(1) or False
        ticks = {"n": 0}

        def over():
            ticks["n"] += 1
            return ticks["n"] > 4

        ctx.over = over
        asyncio.run(cl.run_filled(1e9, None, False))
        self.assertEqual(called, [])
        self.assertEqual(cl.centre.key, start)

    def test_the_fine_tune_does_not_breed_on_idle_cores(self):
        ctx, cl = self.climber()
        ctx.client = type("C", (), {"idle_slots": staticmethod(lambda: 400)})()
        ctx.climbers = [cl]
        bred = []

        async def boom(*a, **k):
            bred.append(1)
            return []

        with unittest.mock.patch("search.climb.submit_ga", boom):
            asyncio.run(cl.run(1e9))
        self.assertEqual(bred, [])
        self.assertIn(cl.state, ("peak", "stopped"))

    def test_a_dead_read_grid_is_not_repeated(self):
        ctx, cl = self.climber()
        kind, group = cl.next_group()
        self.assertEqual(kind, "read")
        cl.flat.update(group)
        kind, group = cl.next_group()
        self.assertEqual(kind, "assembly")
        kind, again = cl.next_group()
        self.assertEqual(kind, "assembly")
        self.assertNotEqual(again, group)

    def test_a_read_pass_grid_is_the_compass_plus_lines_not_the_full_cube(self):
        ctx, cl = self.climber()
        kind, group = cl.next_group()
        self.assertEqual((kind, len(group)), ("read", 5))
        centre = cl.centre
        compass = {c.key for c in cl.compass_heads(centre, group)}
        full = {c.key for c in (cl.with_(centre, dict(zip(group, v)))
                                for v in itertools.product(*[cl.values(n, centre.params[n]) for n in group])) if c}
        two_plus = full - compass - {centre.key}
        asyncio.run(cl.grid(kind, group, until=1e9))
        heads = [b[0] for b in ctx.sent]
        sent = {h.key for h in heads}
        self.assertTrue(compass <= sent)                               # each setting at -d and +d
        self.assertFalse(two_plus and two_plus <= sent)                # not the 3^5 cube
        line = []
        for h in heads:
            if h.key in compass or h.key == centre.key:
                continue
            moved = [n for n in cl.all_names if h.params.get(n) != centre.params.get(n)]
            if len(moved) == 1 and moved[0] in group:
                line.append(h)
        self.assertTrue(line)                                          # a line point moves one grid setting only
        het = "heterozygosity"                                         # 21 rungs: the whole line fits
        seen = {h.params[het] for h in heads}
        self.assertEqual(seen, set(ctx.ladders.rungs[het]))
        compass_batch = next(b for b in ctx.sent if b[0].key in compass)
        riders = compass_batch[1:]                                     # a compass head: genotyping only
        from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS
        for r in riders:
            moved = [n for n in cl.all_names if r.params.get(n) != compass_batch[0].params.get(n)]
            self.assertTrue(moved)
            self.assertTrue(set(moved) <= set(GENOTYPING_PARAMS))
            self.assertFalse(set(moved) & set(PAIRHMM_PARAMS))
        later = [c for batch in ctx.sent[1:] for c in batch]
        self.assertEqual({c.params["pcr_indel_model"] for c in later} | {centre.params["pcr_indel_model"]},
                         set(ctx.ladders.rungs["pcr_indel_model"]))

    def test_two_live_axes_keep_the_cube_until_it_scores(self):
        ctx, cl = self.climber()
        cl.wave_centre = cl.centre
        cl.wave_results = []
        cl.kind = "read"
        cl.iteration = 1
        self.assertTrue(cl.begin_confirm(["min_pruning", "assembly_region_padding"]))
        self.assertEqual(cl.phase, "confirm")
        self.assertTrue(cl.batch_live)
        self.assertGreater(len(cl.pending), 0)

    def test_two_live_axes_send_the_joint_cube(self):
        space = self.space
        target = {"min_pruning": 7, "min_dangling_branch_length": 8}

        def peak(c):
            return 100.0 - sum(abs(c.params[n] - v) for n, v in target.items())

        ctx = FakeCtx(space, peak)
        start = space._assemble({"min_pruning": 2, "min_dangling_branch_length": 4})
        cl = Climber(ctx, 0, start, peak(start))
        cl.turn = 1
        kind, group = cl.next_group()
        self.assertEqual(kind, "assembly")
        self.assertTrue({"min_pruning", "min_dangling_branch_length"} <= set(group))
        centre = cl.centre
        asyncio.run(cl.grid(kind, group, until=1e9))
        cube = cl.confirm_heads(centre, ["min_pruning", "min_dangling_branch_length"], {centre.key})
        sent = {h.key for h in (b[0] for b in ctx.sent)}
        self.assertTrue({c.key for c in cube} <= sent)

    def test_the_climb_reaches_the_peak_and_grows_its_step(self):
        ctx, cl = self.climber()
        asyncio.run(cl.run(until=1e9))
        self.assertEqual(cl.score, 100.0)
        self.assertEqual({n: cl.centre.params[n] for n in self.target}, self.target)
        self.assertEqual(cl.state, "peak")

    def test_a_long_ladder_grows_one_two_three_five(self):
        ctx, cl = self.climber()
        wide, short = "assembly_region_padding", "min_pruning"
        self.assertGreaterEqual(len(ctx.ladders.rungs[wide]), 30)
        self.assertLess(len(ctx.ladders.rungs[short]), 30)

        def won(name, value, step):
            cl.growth[name] = (1, step)
            old = ctx.space._assemble({name: value})
            new = ctx.space._assemble({name: value + step})
            cl.grow([name], old, new)
            return cl.growth[name]

        self.assertEqual(won(wide, 100, 1), (1, 2))
        self.assertEqual(won(wide, 100, 2), (1, 3))
        self.assertEqual(won(wide, 100, 3), (1, 5))
        self.assertEqual(won(wide, 100, 5), (1, 5))
        self.assertEqual(won(short, 2, 1), (1, 2))
        self.assertEqual(won(short, 2, 2), (1, 3))
        self.assertEqual(won(short, 2, 3), (1, 4))

    def test_the_lowest_climber_takes_a_ga_best(self):
        ctx, low = self.climber()
        high_start = self.space._assemble({"min_pruning": 6, "assembly_region_padding": 55})
        high = Climber(ctx, 1, high_start, 90.0)
        far = dict(high_start.params)
        for name in ("min_pruning", "assembly_region_padding", "heterozygosity", "min_dangling_branch_length"):
            far[name] = ctx.ladders.shift(name, far[name], 4)
        ga = self.space._assemble(far)
        ctx.climbers = [low, high]
        ctx.ga_keys = {ga.key}
        ctx.best = lambda: (ga, 100.0)
        low.score = low.found = 10.0
        high.score = high.found = 90.0
        low.pending = [low.centre]
        self.assertTrue(adopt_ga_best(ctx, [low, high]))
        self.assertEqual(low.centre.key, ga.key)
        self.assertEqual(low.score, 100.0)
        self.assertEqual(high.centre.key, high_start.key)
        self.assertTrue(low.pending)                                 # the new grid is ready
        self.assertFalse(any(c.key == low.start.key for c in low.pending[:1]))
        self.assertFalse(adopt_ga_best(ctx, [low, high]))            # it already sits on that score
        ctx.ga_keys.clear()
        ctx.best = lambda: (self.space._assemble({"min_pruning": 3}), 200.0)
        self.assertFalse(adopt_ga_best(ctx, [low, high]))            # a grid point is not a GA handoff

    def test_a_ga_best_replaces_the_climber_with_the_lowest_best_found(self):
        ctx, low = self.climber()
        high_start = self.space._assemble({"min_pruning": 6, "assembly_region_padding": 55})
        high = Climber(ctx, 1, high_start, 90.0)
        returned = self.space._assemble({"min_pruning": 3, "assembly_region_padding": 40})
        far = dict(high_start.params)
        for name in ("min_pruning", "assembly_region_padding", "heterozygosity", "min_dangling_branch_length"):
            far[name] = ctx.ladders.shift(name, far[name], 4)
        ga = self.space._assemble(far)
        ctx.climbers = [low, high]
        ctx.ga_keys = {ga.key}
        low.score, low.found = 10.0, 10.0
        high.score = high.found = 90.0
        low.waiting = {returned.key: returned}
        ctx.score = lambda c: 96.0 if c.key == returned.key else None
        ctx.best = lambda: (ga, 95.0)
        self.assertTrue(adopt_ga_best(ctx, [low, high]))             # high's best is 90, under the GA's 95
        self.assertEqual(high.centre.key, ga.key)                   # low already found 96, so low stays
        self.assertEqual(low.centre.key, low.start.key)
        self.assertEqual(high.found, 95.0)
        self.assertFalse(adopt_ga_best(ctx, [low, high]))            # that GA config is already a centre

    def test_a_ga_below_the_round_best_still_replaces_the_lowest_climber(self):
        ctx, low = self.climber()
        high_start = self.space._assemble({"min_pruning": 6, "assembly_region_padding": 55})
        high = Climber(ctx, 1, high_start, 99.0)
        far = dict(high_start.params)
        for name in ("min_pruning", "assembly_region_padding", "heterozygosity", "min_dangling_branch_length"):
            far[name] = ctx.ladders.shift(name, far[name], 4)
        ga = self.space._assemble(far)
        ctx.climbers = [low, high]
        ctx.ga_keys = {ga.key}
        low.score = low.found = 90.0
        high.score = high.found = 99.0
        ctx.best = lambda: (high.centre, 99.0)                      # the round's best is a climber, not the GA
        ctx.scored = lambda: [(ga, 95.0), (high.centre, 99.0)]
        self.assertTrue(adopt_ga_best(ctx, [low, high]))
        self.assertEqual(low.centre.key, ga.key)
        self.assertEqual(high.centre.key, high_start.key)
        self.assertFalse(adopt_ga_best(ctx, [low, high]))            # that config is already this climber's centre

    def test_the_next_free_ga_replaces_the_lowest_when_the_best_is_seated(self):
        ctx, low = self.climber()
        high_start = self.space._assemble({"min_pruning": 6, "assembly_region_padding": 55})
        high = Climber(ctx, 1, high_start, 99.0)
        far = dict(high_start.params)
        farther = dict(high_start.params)
        names = ("min_pruning", "assembly_region_padding", "heterozygosity", "min_dangling_branch_length")
        other = ("min_base_quality_score", "min_mapping_quality_score", "max_num_haplotypes_in_population",
                 "min_assembly_region_size")
        for name in names:
            far[name] = ctx.ladders.shift(name, far[name], 4) or far[name]
        for name in other:
            farther[name] = ctx.ladders.shift(name, farther[name], 4) or farther[name]
        seated_ga = self.space._assemble(far)
        free_ga = self.space._assemble(farther)
        self.assertIsNotNone(seated_ga)
        self.assertIsNotNone(free_ga)
        high.centre = high.start = seated_ga
        high.score = high.found = 99.0
        ctx.climbers = [low, high]
        ctx.ga_keys = {seated_ga.key, free_ga.key}
        low.score = low.found = 90.0
        ctx.best = lambda: (seated_ga, 99.0)
        ctx.scored = lambda: [(seated_ga, 99.0), (free_ga, 95.0)]
        self.assertTrue(adopt_ga_best(ctx, [low, high]))
        self.assertEqual(low.centre.key, free_ga.key)
        self.assertEqual(high.centre.key, seated_ga.key)

    def test_adopt_during_the_ga_opening_does_not_start_a_grid(self):
        ctx, low = self.climber()
        high_start = self.space._assemble({"min_pruning": 6, "assembly_region_padding": 55})
        high = Climber(ctx, 1, high_start, 90.0)
        far = dict(high_start.params)
        for name in ("min_pruning", "assembly_region_padding", "heterozygosity", "min_dangling_branch_length"):
            far[name] = ctx.ladders.shift(name, far[name], 4)
        ga = self.space._assemble(far)
        ctx.climbers = [low, high]
        ctx.ga_keys = {ga.key}
        ctx.best = lambda: (ga, 100.0)
        ctx.ga_opening = True
        low.score = low.found = 10.0
        high.score = high.found = 90.0
        self.assertTrue(adopt_ga_best(ctx, [low, high]))
        self.assertEqual(low.centre.key, ga.key)
        self.assertFalse(low.pending)
        self.assertFalse(low.batch_live)

    def test_a_flat_compass_step_lands_on_the_next_untried_rung(self):
        ctx, cl = self.climber()
        centre = cl.centre
        name = "assembly_region_padding"
        tied_v = ctx.ladders.shift(name, centre.params[name], 1)
        further = ctx.ladders.shift(name, centre.params[name], 2)
        down = ctx.ladders.shift(name, centre.params[name], -1)
        tied_c = cl.with_(centre, {name: tied_v})
        ctx.scored = lambda: [(centre, cl.score), (tied_c, cl.score)]
        vals = {h.params[name] for h in cl.compass_heads(centre, [name])}
        self.assertIn(further, vals)
        self.assertNotIn(tied_v, vals)
        self.assertIn(down, vals)
        # a rung that already scored differently is measured; the compass does not walk past it
        changed = cl.with_(centre, {name: further})
        ctx.scored = lambda: [(centre, cl.score), (tied_c, cl.score), (changed, cl.score + 0.2)]
        vals = {h.params[name] for h in cl.compass_heads(centre, [name])}
        self.assertNotIn(further, vals)
        self.assertNotIn(ctx.ladders.shift(name, further, 1), vals)

    def test_a_compass_step_inside_one_bam_class_lands_on_the_next_class(self):
        from search.bamgroups import BQST, BamFacts, bounds_from
        facts = BamFacts(frozenset({2, 3, 4, 5, 6, 10, 20, 30}), frozenset(range(61)), 74)
        bounds = bounds_from(facts, ranges())
        space = Space(ranges(), bounds=bounds)
        ctx = FakeCtx(space, lambda c: 80.0)
        centre = space._assemble({BQST: 11, "min_mapping_quality_score": 59})
        cl = Climber(ctx, 0, centre, 80.0)
        nxt = ctx.ladders.shift(BQST, centre.params[BQST], 1)
        self.assertEqual(space._assemble({**centre.params, BQST: nxt}).params[BQST], centre.params[BQST])
        landed = cl.land(centre, BQST, 1)
        self.assertIsNotNone(landed)
        self.assertNotEqual(space._assemble({**centre.params, BQST: landed}).params[BQST], centre.params[BQST])

    def test_a_value_that_already_tied_is_not_sent_again(self):
        space = self.space
        tied = TiedCallsets(space.names)
        base = space._assemble({"min_pruning": 4, "assembly_region_padding": 40})
        pruning = space._assemble({"min_pruning": 5, "assembly_region_padding": 40})
        padding = space._assemble({"min_pruning": 4, "assembly_region_padding": 45})
        both = space._assemble({"min_pruning": 5, "assembly_region_padding": 45})
        untried = space._assemble({"min_pruning": 6, "assembly_region_padding": 40})
        other_setting = space._assemble({"min_pruning": 5, "assembly_region_padding": 40, "heterozygosity": 0.01})
        tied.note(base.key, base.params, 50.0)
        tied.note(pruning.key, pruning.params, 50.0)
        self.assertFalse(tied.holds(both.params))                    # padding has not tied yet
        tied.note(padding.key, padding.params, 50.0)
        self.assertTrue(tied.holds(both.params))                     # both changed values have tied
        self.assertFalse(tied.holds(untried.params))                 # pruning 6 has not tied
        self.assertFalse(tied.holds(other_setting.params))           # heterozygosity has not tied
        lost = space._assemble({"min_pruning": 6, "assembly_region_padding": 40})
        tied.note(lost.key, lost.params, 40.0)                      # a lower score does not tie the value at 50
        self.assertFalse(tied.holds(untried.params))

    def test_a_config_within_two_settings_of_a_climber_is_not_climbed_again(self):
        ctx, low = self.climber()
        low.centre = low.start = self.space._assemble({
            "min_pruning": 2, "assembly_region_padding": 50,
            "min_base_quality_score": 30, "min_mapping_quality_score": 40, "heterozygosity": 0.008})
        high_start = self.space._assemble({"min_pruning": 6, "assembly_region_padding": 55, "heterozygosity": 0.001})
        high = Climber(ctx, 1, high_start, 90.0)
        near = self.space._assemble({"min_pruning": 7, "assembly_region_padding": 55, "heterozygosity": 0.001})
        ctx.climbers = [low, high]
        self.assertIs(held_by_climber(ctx, near, me=low), high)
        self.assertIsNone(held_by_climber(ctx, near, me=high))       # the owner still grids its own neighbourhood
        three = dict(high.centre.params)
        four = dict(high.centre.params)
        axes = ("min_pruning", "min_dangling_branch_length", "assembly_region_padding", "heterozygosity")
        for i, name in enumerate(axes):
            step = ctx.ladders.shift(name, high.centre.params[name], 4) or ctx.ladders.shift(name, high.centre.params[name], -4)
            self.assertNotEqual(step, high.centre.params[name])
            if i < 3:
                three[name] = step
            four[name] = step
        three_cfg = self.space._assemble(three)
        four_cfg = self.space._assemble(four)
        self.assertIs(held_by_climber(ctx, three_cfg, me=low), high)   # three settings is still this hill
        self.assertIsNone(held_by_climber(ctx, four_cfg, me=low))      # four settings is a different hill
        ctx.ga_keys = {near.key}
        ctx.best = lambda: (near, 100.0)
        low.score, high.score = 10.0, 90.0
        self.assertFalse(adopt_ga_best(ctx, [low, high]))
        self.assertEqual(low.centre.key, low.start.key)
        low.pending = [near]
        low._submit_head()
        self.assertEqual(ctx.sent, [])
        self.assertEqual(low.pending, [])
        high.pending = [near]
        high._submit_head()
        self.assertEqual(ctx.sent[0][0].key, near.key)

    def test_a_far_seat_does_not_take_a_near_hill_ga(self):
        ctx, high = self.climber()
        high.id = 0
        far_start = dict(high.centre.params)
        names = ("min_pruning", "assembly_region_padding", "heterozygosity", "min_dangling_branch_length")
        for name in names:
            far_start[name] = ctx.ladders.shift(name, far_start[name], 4) or far_start[name]
        far_cfg = self.space._assemble(far_start)
        far = Climber(ctx, 1, far_cfg, 10.0, origin="far")
        near_ga = self.space._assemble({**high.centre.params, "min_pruning": high.centre.params["min_pruning"] + 1})
        ctx.climbers = [high, far]
        ctx.ga_keys = {near_ga.key}
        ctx.best = lambda: (near_ga, 100.0)
        high.score = high.found = 80.0
        far.score = far.found = 10.0
        self.assertTrue(adopt_ga_best(ctx, [high, far]))
        self.assertEqual(high.centre.key, near_ga.key)
        self.assertEqual(far.centre.key, far_cfg.key)

    def test_a_far_seat_takes_a_ga_that_is_far_from_the_near_centres(self):
        ctx, high = self.climber()
        high.id = 0

        def shove(params, names):
            out = dict(params)
            for name in names:
                v = ctx.ladders.shift(name, out[name], 4) or ctx.ladders.shift(name, out[name], -4)
                if v is not None:
                    out[name] = v
            return out

        far_cfg = self.space._assemble(shove(high.centre.params, (
            "min_pruning", "assembly_region_padding", "heterozygosity", "min_dangling_branch_length")))
        ga = self.space._assemble(shove(high.centre.params, (
            "min_base_quality_score", "min_mapping_quality_score", "max_num_haplotypes_in_population",
            "min_assembly_region_size", "active_probability_threshold")))
        self.assertGreaterEqual(sum(ga.params.get(n) != high.centre.params.get(n) for n in ctx.space.names), 4)
        far = Climber(ctx, 1, far_cfg, 10.0, origin="far")
        ctx.climbers = [high, far]
        ctx.ga_keys = {ga.key}
        ctx.best = lambda: (ga, 100.0)
        high.score = high.found = 80.0
        far.score = far.found = 10.0
        self.assertTrue(adopt_ga_best(ctx, [high, far]))
        self.assertEqual(far.centre.key, ga.key)
        self.assertEqual(high.centre.key, high.start.key)

    def test_the_parents_skip_a_higher_score_that_is_close(self):
        ctx, cl = self.climber()
        centre = cl.centre
        near = self.space._assemble({**centre.params, "min_pruning": centre.params["min_pruning"] + 1})
        moved = dict(centre.params)
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            moved[name] = ctx.ladders.shift(name, centre.params[name], 4) or ctx.ladders.shift(name, centre.params[name], -4)
        far_cfg = self.space._assemble(moved)
        step = StepA(ctx)
        step.pool = lambda: [Point(centre, 50.0, 0.9, 10.0), Point(near, 49.9, 0.9, 11.0),
                             Point(far_cfg, 49.0, 0.88, 12.0)]
        space = self.space

        class _Dist:
            def far(self, a, b, min_callset, min_settings):
                return setting_distance(space, a, b) >= min_settings

        ctx.dist = _Dist()
        got = asyncio.run(step.far_parents(2))
        self.assertEqual([p.config.key for p in got], [centre.key, far_cfg.key])

    def test_the_parents_are_reused_until_the_best_score_moves(self):
        ctx, cl = self.climber()
        centre = cl.centre
        moved = dict(centre.params)
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            moved[name] = ctx.ladders.shift(name, centre.params[name], 4) or ctx.ladders.shift(name, centre.params[name], -4)
        far_cfg = self.space._assemble(moved)
        extra = self.space._assemble({**centre.params, "min_pruning": centre.params["min_pruning"] + 2})
        pool = [Point(centre, 50.0, 0.9, 10.0), Point(far_cfg, 49.0, 0.88, 12.0)]
        step = StepA(ctx)
        step.pool = lambda: list(pool)
        space = self.space

        class _Dist:
            def far(self, a, b, min_callset, min_settings):
                return setting_distance(space, a, b) >= min_settings

        ctx.dist = _Dist()
        first = asyncio.run(step.far_parents(2))
        pool.append(Point(extra, 49.5, 0.89, 11.0))
        again = asyncio.run(step.far_parents(2))
        self.assertIs(first, again)
        pool[0] = Point(centre, 50.1, 0.9, 10.0)
        moved_best = asyncio.run(step.far_parents(2))
        self.assertIsNot(moved_best, first)

    def test_the_last_opening_seat_is_the_best_far_config(self):
        ctx, cl = self.climber()
        ctx.cap = 500
        centre = cl.centre
        near = self.space._assemble({**centre.params, "min_pruning": centre.params["min_pruning"] + 1})
        moved = dict(centre.params)
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            moved[name] = ctx.ladders.shift(name, centre.params[name], 4) or ctx.ladders.shift(name, centre.params[name], -4)
        far_cfg = self.space._assemble(moved)
        step = StepA(ctx)
        step.pool = lambda: [Point(centre, 50.0, 0.9, 10.0), Point(near, 49.95, 0.9, 11.0),
                             Point(far_cfg, 49.5, 0.7, 5.0)]
        space = self.space

        class _Dist:
            def between(self, a, b):
                pair = frozenset((a.key, b.key))
                if pair == frozenset((centre.key, near.key)):
                    return 0.01, "callset"
                if centre.key in pair or near.key in pair:
                    return 0.11, "callset"
                return 0.0, "callset"

            def far(self, a, b, min_callset, min_settings):
                return self.between(a, b)[0] >= min_callset

        ctx.dist = _Dist()
        seats = asyncio.run(step.pick())
        self.assertEqual(seats[0].key, centre.key)
        self.assertEqual(seats[1].key, far_cfg.key)
        from search.piles import enabled as piles_enabled
        self.assertEqual(step.picked[1]["group"], "pile 2" if piles_enabled() else "far")

    def test_an_elite_child_differs_in_at_least_four_settings(self):
        ctx, cl = self.climber()
        step = StepA(ctx)
        elite = cl.centre
        mate_p = dict(elite.params)
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            mate_p[name] = ctx.ladders.shift(name, elite.params[name], 4) or ctx.ladders.shift(name, elite.params[name], -4)
        mate = self.space._assemble(mate_p)
        child = step.elite_child(elite, mate)
        self.assertIsNotNone(child)
        diffs = sum(child.params.get(n) != elite.params.get(n) for n in ctx.space.names)
        self.assertGreaterEqual(diffs, 4)

    def test_the_elite_window_probes_one_to_three_settings(self):
        from search.piles import Pile
        ctx, cl = self.climber()
        step = StepA(ctx)
        ctx.ga_opening = True
        ctx.best = lambda: (cl.centre, cl.score)
        ctx._piles = [Pile(cl.centre, cl.score)]
        ctx._pile_seen = {cl.centre.key}
        said = []
        ctx.say = lambda msg, *args: said.append(msg % args if args else msg)
        children = asyncio.run(step.breed_piles(1, 12, why="GA"))
        self.assertEqual(len(children), 12)
        for child in children:
            diffs = sum(child.params.get(n) != cl.centre.params.get(n) for n in ctx.space.names)
            self.assertGreaterEqual(diffs, 1)
            self.assertLessEqual(diffs, 3)
        self.assertIn("elite probes", said[-1])

    def test_an_elite_probe_close_to_one_in_flight_is_still_sent(self):
        from search.piles import Pile, enabled as piles_enabled
        if not piles_enabled():
            self.skipTest("elite probes of the best run while piles are on")
        ctx, cl = self.climber()
        step = StepA(ctx)
        ctx.step_a = step
        ctx.ga_opening = True
        ctx.best = lambda: (cl.centre, cl.score)
        ctx._piles = [Pile(cl.centre, cl.score)]
        ctx._pile_seen = {cl.centre.key}
        near = step.elite_probe(cl.centre, cl.centre)
        self.assertIsNotNone(near)
        ctx.ga_hills = [near]
        ctx.ga_live = [(near, None)]
        asyncio.run(submit_ga(ctx, 30))
        sent = [c for batch in ctx.sent for c in batch]
        self.assertGreaterEqual(len(sent), 20)
        close = 0
        for child in sent:
            diffs = sum(child.params.get(n) != cl.centre.params.get(n) for n in ctx.space.names)
            self.assertGreaterEqual(diffs, 1)
            self.assertLessEqual(diffs, 3)
            near_diffs = sum(child.params.get(n) != near.params.get(n) for n in ctx.space.names)
            if near_diffs < 4:
                close += 1
        self.assertGreaterEqual(close, 10)

    def test_a_mutation_moves_one_rung(self):
        ctx, cl = self.climber()
        step = StepA(ctx)
        parent = cl.centre
        from search.separate import _LINKED
        moved = 0
        for _ in range(100):
            child = step.child(parent, parent)
            self.assertIsNotNone(child)
            for n in _LINKED:
                self.assertEqual(child.params.get(n), parent.params.get(n))
            for n in ctx.space.names:
                if child.params.get(n) == parent.params.get(n):
                    continue
                moved += 1
                one = {ctx.ladders.shift(n, parent.params[n], 1), ctx.ladders.shift(n, parent.params[n], -1)}
                self.assertIn(child.params[n], one)
        self.assertGreater(moved, 0)

    def test_a_child_keeps_each_verified_group_with_one_parent(self):
        from search.separate import LINKED_GROUPS
        ctx, cl = self.climber()
        ctx.rng.seed(1)
        step = StepA(ctx)
        left = cl.centre
        right_p = dict(left.params)
        for name in ctx.space.names:
            shifted = ctx.ladders.shift(name, left.params[name], 3)
            if shifted is not None and shifted != left.params[name]:
                right_p[name] = shifted
        if right_p.get("pcr_indel_model") == left.params.get("pcr_indel_model"):
            right_p["pcr_indel_model"] = "NONE" if left.params.get("pcr_indel_model") != "NONE" else "AGGRESSIVE"
        right = self.space._assemble(right_p)
        saw_left = {i: False for i in range(len(LINKED_GROUPS))}
        saw_right = {i: False for i in range(len(LINKED_GROUPS))}
        for _ in range(40):
            kid = step.child(left, right)
            self.assertIsNotNone(kid)
            elite = step.elite_child(left, right)
            self.assertIsNotNone(elite)
            for born in (kid, elite):
                for i, group in enumerate(LINKED_GROUPS):
                    from_left = all(born.params.get(name) == left.params.get(name) for name in group)
                    from_right = all(born.params.get(name) == right.params.get(name) for name in group)
                    self.assertTrue(from_left or from_right)
                    if born is kid:
                        saw_left[i] = saw_left[i] or from_left
                        saw_right[i] = saw_right[i] or from_right
        self.assertTrue(all(saw_left.values()))
        self.assertTrue(all(saw_right.values()))

    def test_step_b_mutation_rises_from_0_3_to_0_7(self):
        from search.candidates import mutation_share, step_b_mutation
        ctx, _ = self.climber()
        ctx.b_span = (0.0, 100.0)
        ctx.phase = "A: candidates"
        ctx.now = 50.0
        self.assertIsNone(step_b_mutation(ctx))
        ctx.phase = "B: climbing"
        ctx.now = 0.0
        self.assertAlmostEqual(step_b_mutation(ctx), 0.3)
        ctx.now = 50.0
        self.assertAlmostEqual(step_b_mutation(ctx), 0.5)
        ctx.now = 100.0
        self.assertAlmostEqual(step_b_mutation(ctx), 0.7)
        ctx.now = 140.0
        self.assertAlmostEqual(step_b_mutation(ctx), 0.7)
        ctx.phase = "C: GA"
        ctx.step_c_gen = None
        self.assertAlmostEqual(mutation_share(ctx), 0.3)
        ctx.step_c_gen = 1
        self.assertAlmostEqual(mutation_share(ctx), 0.3)
        ctx.step_c_gen = 2
        self.assertAlmostEqual(mutation_share(ctx), 0.6)

    def test_a_plateau_probes_once_and_then_breeds_children_far_from_the_climbers(self):
        ctx, cl = self.climber()
        centre = cl.centre
        twin = self.space._assemble({**centre.params,
                                     "max_alternate_alleles": centre.params["max_alternate_alleles"] + 1,
                                     "max_num_haplotypes_in_population": 30})
        base = 50.0
        cl.score = base
        ctx.climbers = [cl]
        ctx.best = lambda: (centre, base)
        ctx.scored = lambda: [(centre, base), (twin, base)]
        probes = next_probes(ctx, centre, base)
        self.assertTrue(probes)
        for p in probes:
            diffs = [n for n in ctx.space.names if p.params.get(n) != centre.params.get(n)]
            self.assertEqual(diffs, [n for n in diffs if n not in
                             ("max_alternate_alleles", "max_num_haplotypes_in_population")])
            self.assertEqual(len(diffs), 1)
            one = {ctx.ladders.shift(diffs[0], centre.params[diffs[0]], 1),
                   ctx.ladders.shift(diffs[0], centre.params[diffs[0]], -1)}
            self.assertIn(p.params[diffs[0]], one)
        self.assertEqual(asyncio.run(submit_ga(ctx, 10)), [])          # probes are queued; nothing breeds until a breeder exists
        self.assertTrue(cl.pending)
        held = len(cl.pending)
        near = self.space._assemble({**centre.params, "min_pruning": centre.params["min_pruning"] + 1})
        moved = dict(centre.params)
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            moved[name] = ctx.ladders.shift(name, centre.params[name], 4) or ctx.ladders.shift(name, centre.params[name], -4)
        far_child = self.space._assemble(moved)

        class _Dist:
            def between(self, a, b):
                return setting_distance(self_space, a, b), "settings"

            def far(self, a, b, min_callset, min_settings):
                return setting_distance(self_space, a, b) >= min_settings

        self_space = self.space
        ctx.dist = _Dist()

        class _Breed:
            async def breed(self, gen, n, why="", accept=None, elite=False):
                out = []
                for c in (near, far_child):
                    ok = True if c is not None and accept is None else (accept(c) if c is not None else False)
                    if hasattr(ok, "__await__"):
                        ok = await ok
                    if ok:
                        out.append(c)
                return out[:n]

        ctx.step_a = _Breed()
        futs = asyncio.run(submit_ga(ctx, 10))
        self.assertEqual(len(cl.pending), held)            # the probes are not queued again
        self.assertEqual(len(futs), 1)
        self.assertEqual(ctx.sent[-1][0].key, far_child.key)
        ctx._probe_base, ctx._probe_done = base, True
        again = asyncio.run(submit_ga(ctx, 10))
        self.assertEqual(len(cl.pending), held)
        self.assertEqual(len(again), 1)                            # that GA run is back: the hill is not a seat
        self.assertEqual(ctx.sent[-1][0].key, far_child.key)

    def test_a_ga_child_is_not_sent_on_the_same_hill_as_another_ga_child(self):
        ctx, cl = self.climber()
        ctx.climbers = [cl]
        ctx.best = lambda: (cl.centre, cl.score)
        ctx.scored = lambda: [(cl.centre, cl.score)]
        first = dict(cl.centre.params)
        second = dict(cl.centre.params)
        names = ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                 "heterozygosity", "min_base_quality_score", "min_mapping_quality_score")
        for i, name in enumerate(names):
            first[name] = ctx.ladders.shift(name, first[name], 4) or first[name]
            second[name] = ctx.ladders.shift(name, second[name], 4) or second[name]
        for name in ("min_pruning", "heterozygosity"):
            second[name] = ctx.ladders.shift(name, second[name], -4) or second[name]
        a = self.space._assemble(first)
        b = self.space._assemble(second)
        self.assertIsNotNone(a)
        self.assertIsNotNone(b)
        diffs = sum(a.params.get(n) != b.params.get(n) for n in ctx.space.names)
        self.assertLess(diffs, 4)

        class _Breed:
            async def breed(self, gen, n, why="", accept=None, elite=False):
                out = []
                for c in (a, b):
                    ok = True if c is not None and accept is None else (accept(c) if c is not None else False)
                    if hasattr(ok, "__await__"):
                        ok = await ok
                    if ok:
                        out.append(c)
                return out[:n]

        ctx.step_a = _Breed()
        futs = asyncio.run(submit_ga(ctx, 10))
        self.assertEqual(len(futs), 1)
        self.assertEqual(ctx.sent[-1][0].key, a.key)
        self.assertTrue(held_by_ga(ctx, b))

    def test_an_empty_seat_takes_the_best_scored_config_that_is_far(self):
        ctx, cl = self.climber()
        centre = cl.centre
        ctx.climbers = [cl]
        ctx.best = lambda: (centre, 50.0)
        near = self.space._assemble({**centre.params, "min_pruning": centre.params["min_pruning"] + 1})
        moved = dict(centre.params)
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            moved[name] = ctx.ladders.shift(name, centre.params[name], 4) or ctx.ladders.shift(name, centre.params[name], -4)
        far_cfg = self.space._assemble(moved)
        low = self.space._assemble({**moved, "max_alternate_alleles": 12})

        space = self.space

        class _Dist:
            def between(self, a, b):
                return setting_distance(space, a, b), "settings"

            def far(self, a, b, min_callset, min_settings):
                return setting_distance(space, a, b) >= min_settings

        ctx.dist = _Dist()
        self.assertIsNotNone(far_cfg)
        ctx.scored = lambda: [(centre, 50.0), (near, 49.9), (far_cfg, 49.4), (low, 40.0)]
        got = far_replacement(ctx)
        self.assertIsNotNone(got)
        self.assertEqual(got[0].key, far_cfg.key)
        self.assertEqual(got[2], "far (scored)")
        ctx.scored = lambda: [(centre, 50.0), (near, 49.9)]
        self.assertIsNone(far_replacement(ctx))
        scouts = far_scouts(ctx, 2)
        self.assertEqual(len(scouts), 2)
        from common.tuning import tune
        ms = float(tune("candidates.min_setting_distance"))
        for c in scouts:
            self.assertGreaterEqual(setting_distance(self.space, c, centre), ms)
            self.assertGreater(sum(c.params.get(n) != centre.params.get(n) for n in ctx.space.names), 2)

    def test_a_small_loss_steps_past_the_rung_and_a_loss_of_0_10_closes_the_side(self):
        ctx, cl = self.climber()
        centre = cl.centre
        name = "assembly_region_padding"
        i0 = int(ctx.ladders.index(name, centre.params[name]))
        twin = self.space._assemble({**centre.params, "max_alternate_alleles": 8, "max_num_haplotypes_in_population": 30})
        down = single(ctx, centre, name, i0 - 1)
        small = single(ctx, centre, name, i0 + 1)
        base = 50.0
        ctx.scored = lambda: [(centre, base), (twin, base), (down, base), (small, base - 0.02)]
        probe = next(p for p in next_probes(ctx, centre, base) if p.params[name] != centre.params[name])
        self.assertEqual(int(ctx.ladders.index(name, probe.params[name])), i0 + 2)
        closed = single(ctx, centre, name, i0 + 1)
        ctx.scored = lambda: [(centre, base), (twin, base), (down, base), (closed, base - 0.10)]
        again = [p for p in next_probes(ctx, centre, base) if p.params[name] != centre.params[name]]
        self.assertTrue(again)
        self.assertLess(int(round(ctx.ladders.index(name, again[0].params[name]))), i0)

    def test_a_gain_keeps_other_flat_settings_and_two_flat_grids_release_a_trailing_seat(self):
        ctx, low = self.climber()
        high = Climber(ctx, 1, low.centre, 90.0)
        ctx.climbers = [low, high]
        low.score = 10.0
        low.note_gain(["min_base_quality_score"], 0.0)
        low.note_gain(["min_pruning"], 0.2)                          # the assembly gain does not forget the read failure
        self.assertIn("min_base_quality_score", low.flat)
        self.assertNotIn("min_pruning", low.flat)
        self.assertFalse(low.trailing())
        low.note_gain(["min_pruning"], 0.0)
        low.note_gain(["min_base_quality_score"], 0.0)               # the last two grids are both flat
        self.assertTrue(low.trailing())
        self.assertGreaterEqual(low.gap_to_best(), 0.10)

    def test_a_returned_run_above_the_centre_takes_it(self):
        ctx, cl = self.climber()
        higher = self.space._assemble({**cl.centre.params, "min_pruning": cl.centre.params["min_pruning"] + 1})
        cl.owed[higher.key] = higher
        base = cl.score
        ctx.f = lambda c, higher=higher, base=base: base + 1.0 if c.key == higher.key else base
        self.assertTrue(cl.claim_higher())
        self.assertEqual(cl.centre.key, higher.key)
        self.assertEqual(cl.score, base + 1.0)

    def test_a_returned_run_inside_the_pile_becomes_the_centre_when_the_gap_fails(self):
        from search.piles import Pile, enabled as piles_enabled
        if not piles_enabled():
            self.skipTest("pile fingerprints are off on this tag")
        ctx, cl = self.climber()
        cl.origin = "pile 1"
        other = self.space._assemble({**cl.centre.params, "assembly_region_padding": cl.centre.params["assembly_region_padding"] + 20})
        higher = self.space._assemble({**cl.centre.params, "min_pruning": cl.centre.params["min_pruning"] + 1})
        ctx._piles = [Pile(cl.centre, cl.score), Pile(other, cl.score - 1)]
        fingerprint = ctx._piles[0].config.key

        class _Dist:
            def between(self, a, b):
                if other.key in (a.key, b.key) and higher.key in (a.key, b.key):
                    return 0.04, "callset"
                if other.key in (a.key, b.key):
                    return 0.20, "callset"
                return 0.0, "callset"

        ctx.dist = _Dist()
        cl.owed[higher.key] = higher
        base = cl.score
        ctx.f = lambda c, higher=higher, base=base: base + 1.0 if c.key == higher.key else base
        self.assertTrue(cl.claim_higher())
        self.assertEqual(cl.centre.key, higher.key)
        self.assertEqual(cl.score, base + 1.0)
        self.assertEqual(len(ctx._piles), 1)
        self.assertEqual(ctx._piles[0].config.key, higher.key)
        self.assertIn(fingerprint, ctx._piles[0].scores)
        self.assertIn(other.key, ctx._piles[0].scores)

    def test_a_live_grid_keeps_the_centre_and_writes_the_setting_itself(self):
        from search.climb import adopt_reference
        from search.piles import Pile, enabled as piles_enabled, note
        if not piles_enabled():
            self.skipTest("pile fingerprints are off on this tag")
        ctx, cl = self.climber()
        cl.origin = "pile 1"
        cl.batch_live = True
        cl.kind = "read"
        cl.group = ["min_pruning"]
        cl.wave_before = cl.score
        cl.wave_centre = cl.centre
        higher = self.space._assemble({**cl.centre.params, "min_pruning": cl.centre.params["min_pruning"] + 1})
        pile = Pile(cl.centre, cl.score)
        ctx._piles = [pile]
        ctx.climbers = [cl]
        note(ctx, [pile], higher, cl.score + 1.0, cl.score + 1.0)
        self.assertEqual(pile.config.key, cl.centre.key)
        self.assertIn(higher.key, pile.scores)
        self.assertFalse(adopt_reference(ctx, [cl]))
        self.assertEqual(cl.centre.key, pile.config.key)
        self.assertFalse(cl.claim_higher())
        cl.wave_results = [(higher, cl.score + 1.0)]
        cl.line_keys = {higher.key}
        cl._finish_grid_wave()
        self.assertEqual(cl.centre.key, higher.key)
        self.assertEqual(ctx._piles[0].config.key, higher.key)

    def test_step_b_breeds_a_hill_outside_the_peak_interval(self):
        """Generation 1 stays inside the interval. Step B's free cores breed from a hill outside it."""
        from common.types import FULL
        ctx, cl = self.climber()
        centre = cl.centre
        moved = dict(centre.params)
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            moved[name] = ctx.ladders.shift(name, centre.params[name], 4) or ctx.ladders.shift(name, centre.params[name], -4)
        far_cfg = self.space._assemble(moved)
        self.assertGreaterEqual(setting_distance(self.space, centre, far_cfg), 0.05)
        step = StepA(ctx)
        step._range = {n: [centre.params.get(n)] for n in ctx.space.names if ctx.ladders.movable(n)}
        self.assertFalse(step._inside(far_cfg))

        class _Row:
            def __init__(self, score, core, fp):
                self.score, self.core, self.region_fp, self.gate_passed = score, core, fp, True

        class _Client:
            latest = {(centre.key, FULL): _Row(50.0, 0.9, 10), (far_cfg.key, FULL): _Row(49.0, 0.8, 20)}
            configs = {centre.key: centre, far_cfg.key: far_cfg}

        class _Dist:
            def between(self, a, b):
                return setting_distance(self_space, a, b), "settings"

            def far(self, a, b, min_callset, min_settings):
                return setting_distance(self_space, a, b) >= min_settings

        self_space = self.space
        ctx.client = _Client()
        ctx.dist = _Dist()
        ctx.step_a = step
        ctx.best = lambda: (centre, 50.0)
        ctx.phase = "B"
        self.assertEqual(len(step.pool()), 1)
        self.assertEqual(asyncio.run(step.breed(1, 8)), [])
        asyncio.run(submit_ga(ctx, 8))
        sent = [c for batch in ctx.sent for c in batch]
        self.assertTrue(sent)
        self.assertTrue(getattr(step, "_bound_pool", True))
        self.assertTrue(any(not step._inside(c) for c in sent))


class StepCTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.ctx = FakeCtx(self.space, lambda c: 0.0)
        self.rungs = self.ctx.ladders.rungs[MDBL]

    def test_gaps_skipped_by_growth_and_coarse_rungs_beyond(self):
        # the climb went 3 -> 4 -> 6 -> 9 and 9 lost: the peak is 6, 5 and 9 were tried, 7 and 8 never
        peak = self.space._assemble({MDBL: self.rungs[6]})
        known = {n: {} for n in self.space.names if self.ctx.ladders.movable(n)}
        known[MDBL] = {5: 98.0, 9: 97.0}
        want = lines.plan(self.ctx, peak, known, coarse=3)[MDBL]
        self.assertIn(7, want)
        self.assertIn(8, want)
        self.assertNotIn(5, want)
        self.assertNotIn(9, want)
        last = len(self.rungs) - 1
        self.assertEqual([i for i in want if i > 9], lines.spread(10, last, 3))           # beyond 9, coarse
        self.assertEqual([i for i in want if i < 5], lines.spread(0, 4, 3))               # beyond 5, coarse

    def test_a_second_high_point_behind_a_valley(self):
        peak = self.space._assemble({MDBL: self.rungs[3]})
        known = {MDBL: {2: 98.9, 4: 98.9, 6: 95.0, 9: 95.0, 12: 98.85, 14: 98.95, 16: 98.85}}
        highs = far.second_highs(self.ctx, peak, 99.0, known)
        self.assertEqual([(n, i, lo, hi) for _, n, i, lo, hi in highs], [(MDBL, 14, 13, 15)])
        known[MDBL][14] = 98.5                                     # more than far.second_margin below the peak
        self.assertEqual(far.second_highs(self.ctx, peak, 99.0, known), [])

    def landscape(self, far_alone, bonus):
        """MDBL: a peak at rung 3 (99), a valley, a second hill around rung 14; on that hill each min_pruning rung up
        from the start adds `bonus`."""
        start_mp = self.space._assemble({}).params[MP]
        rungs = self.rungs

        def f(c):
            i = rungs.index(c.params[MDBL])
            if i <= 4:
                base = 99.0 - 0.1 * abs(i - 3)
            elif i < 12:
                base = 95.0
            else:
                base = far_alone - 0.1 * abs(i - 14) + bonus * max(0, c.params[MP] - start_mp)
            return base
        return f

    def run_c(self, f):
        ctx = FakeCtx(self.space, f)
        start = self.space._assemble({MDBL: self.rungs[3]})
        cl = Climber(ctx, 0, start, f(start))
        cl.state = "peak"
        ctx.climbers.append(cl)
        asyncio.run(far.run(ctx, [cl], until=1e9))
        return ctx, cl

    def test_a_better_point_beyond_the_valley_moves_the_peak(self):
        ctx, cl = self.run_c(self.landscape(100.5, 0.0))
        self.assertEqual(cl.centre.params[MDBL], self.rungs[14])
        self.assertEqual(cl.score, 100.5)

    def test_a_lower_second_high_wins_after_its_short_climb(self):
        ctx, cl = self.run_c(self.landscape(98.95, 0.3))           # 14 alone is below the peak; its hill is higher
        self.assertEqual(cl.centre.params[MDBL], self.rungs[14])
        self.assertGreater(cl.score, 99.2)
        self.assertTrue(any(d.origin.startswith("step C landing") for d in ctx.climbers))


class VirtualRoundTest(unittest.TestCase):
    def test_a_whole_round_climbs_and_picks(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=1)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        searcher = engine.searcher(r, 300)
        seen = []

        async def wrapped(client):
            seen.append(client)
            await searcher(client)

        out = run_round(ROUND, Stub(), ex, FakeScorer(runner), wrapped, DryRunSubmitter(), r,
                        LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=72 * 60)
        eng = seen[0].notes["engine"]
        from search.engine_tag import VERSION as TAG
        self.assertEqual(eng["version"], TAG)
        # The fine-tune runs until the pick. A simulated job can jump the clock past the round's stop
        # while that climb is in flight, so the last snapshot is still the fine-tune.
        self.assertIn(eng["phase"], ("picked", "C: fine-tune"))
        first = eng["basins"][0]
        self.assertGreaterEqual(first["best"], first["opened"])
        self.assertTrue(first["best_key"] and first["first_key"])
        self.assertGreaterEqual(round(out.best.score, 4), first["best"])
        labels = [w["label"] for w in eng["waves"]]
        self.assertNotIn("generation 1", labels)
        self.assertNotIn("generation 2", labels)
        self.assertNotIn("generation 3", labels)


class PileTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())

    def ctx(self, dist):
        ctx = FakeCtx(self.space, lambda c: 0.0)
        ctx.dist = dist
        return ctx

    def test_a_higher_score_at_0_05_replaces_the_fingerprint(self):
        from search.piles import refresh
        ctx = self.ctx(_Map({}))
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 4})
        ctx.dist = _Map({(a.key, b.key): 0.05})
        piles = refresh(ctx, [(a, 80.0), (b, 81.0)])
        self.assertEqual(len(piles), 1)
        self.assertEqual(piles[0].config.key, b.key)
        self.assertEqual(piles[0].score, 81.0)

    def test_a_higher_score_within_0_10_replaces_the_reference(self):
        from search.piles import refresh
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 6})
        ctx = self.ctx(_Map({(a.key, b.key): 0.08}))
        piles = refresh(ctx, [(a, 80.0), (b, 81.0)])
        self.assertEqual(len(piles), 1)
        self.assertEqual(piles[0].config.key, b.key)
        self.assertEqual(piles[0].score, 81.0)

    def test_a_lower_score_within_0_10_joins_without_replacing(self):
        from search.piles import refresh
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 6})
        ctx = self.ctx(_Map({(a.key, b.key): 0.07}))
        piles = refresh(ctx, [(a, 80.0), (b, 79.5)])
        self.assertEqual(len(piles), 1)
        self.assertEqual(piles[0].config.key, a.key)
        self.assertIn(b.key, piles[0].scores)

    def test_a_fingerprint_farther_than_0_10_opens_a_second_pile(self):
        from search.piles import refresh
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 8})
        ctx = self.ctx(_Map({(a.key, b.key): 0.11}))
        piles = refresh(ctx, [(a, 80.0), (b, 79.0)])
        self.assertEqual([p.config.key for p in piles], [a.key, b.key])

    def test_the_climber_develops_from_the_highest_within_0_10(self):
        from search.climb import Climber, adopt_reference
        from search.piles import Pile
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 6})
        ctx = self.ctx(_Map({(a.key, b.key): 0.08}))
        ctx.ga_opening = True
        ctx.T = lambda seconds: seconds
        ctx.ga_moved = 0.0
        ctx.now = 0.0
        pile = Pile(b, 81.0)
        pile.members.append(a)
        pile.scores[a.key] = 80.0
        ctx._piles = [pile]
        ctx._pile_seen = {a.key, b.key}
        cl = Climber(ctx, 0, a, 80.0, origin="pile 1")
        ctx.climbers.append(cl)
        self.assertTrue(adopt_reference(ctx, [cl]))
        self.assertEqual(cl.centre.key, b.key)
        self.assertEqual(cl.score, 81.0)

    def test_a_near_copy_does_not_take_a_second_parent_slot(self):
        from search.piles import Pile, mating_pool
        ctx = FakeCtx(self.space, lambda c: 0.0)
        a = self.space._assemble({"min_pruning": 2, "assembly_region_padding": 50})
        b = self.space._assemble({"min_pruning": 3, "assembly_region_padding": 50})
        pile = Pile(a, 80.0)
        pile.members.append(b)
        pile.scores[b.key] = 79.0
        self.assertEqual([c.key for c in mating_pool(ctx, pile)], [a.key])

    def test_a_pile_reports_its_three_highest_scores(self):
        from search.piles import Pile, pile_rows
        ctx = FakeCtx(self.space, lambda c: 0.0)
        pile = None
        for i, score in enumerate((80.0, 83.0, 81.0, 82.0)):
            c = self.space._assemble({"min_pruning": 2 + i})
            if pile is None:
                pile = Pile(c, score)
            else:
                pile.members.append(c)
                pile.scores[c.key] = score
                if score > pile.score:
                    pile.config, pile.score = c, score
        ctx._piles = [pile]
        rows = pile_rows(ctx)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["pile"], 1)
        self.assertEqual([item["score"] for item in rows[0]["top"]], [83.0, 82.0, 81.0])
        self.assertEqual(len(rows[0]["top"]), 3)

    def test_the_mating_pool_keeps_32_parents(self):
        from search.piles import POOL, Pile, mating_pool
        ctx = FakeCtx(self.space, lambda c: 0.0)
        pile = None
        for i in range(POOL + 1):
            c = self.space._assemble({
                "assembly_region_padding": 50 + i,
                "min_assembly_region_size": 10 + i,
                "phred_scaled_global_read_mismapping_rate": 10 + i,
                "base_quality_score_threshold": i,
                "max_reads_per_alignment_start": 25 + i,
            })
            if pile is None:
                pile = Pile(c, 100.0 - i)
            else:
                pile.members.append(c)
                pile.scores[c.key] = 100.0 - i
        chosen = mating_pool(ctx, pile)
        self.assertEqual(len(chosen), POOL)
        self.assertEqual(chosen[0].key, pile.config.key)
        self.assertNotIn(pile.members[-1].key, [c.key for c in chosen])

    def test_a_climber_takes_a_centre_inside_its_pile_when_the_gap_fails(self):
        from search.climb import Climber, adopt_reference
        from search.piles import Pile, centre_allowed
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 8})
        step = self.space._assemble({"min_pruning": 5})
        outside = self.space._assemble({"min_pruning": 9})
        ctx = self.ctx(_Map({
            (a.key, b.key): 0.20, (a.key, step.key): 0.04, (b.key, step.key): 0.08,
            (a.key, outside.key): 0.11, (b.key, outside.key): 0.04,
        }))
        ctx.ga_opening = True
        ctx.T = lambda seconds: seconds
        ctx.ga_moved = 0.0
        ctx.now = 0.0
        pile = Pile(a, 80.0)
        pile.members.append(step)
        pile.scores[step.key] = 81.0
        ctx._piles = [pile, Pile(b, 79.0)]
        ctx._pile_seen = {a.key, b.key, step.key}
        cl = Climber(ctx, 0, a, 80.0, origin="pile 1")
        ctx.climbers.append(cl)
        self.assertTrue(centre_allowed(cl, step))
        self.assertFalse(centre_allowed(cl, outside))
        self.assertTrue(adopt_reference(ctx, [cl]))
        self.assertEqual(cl.centre.key, step.key)
        self.assertEqual(cl.score, 81.0)
        self.assertEqual(len(ctx._piles), 1)
        self.assertEqual(ctx._piles[0].config.key, step.key)
        self.assertIn(b.key, ctx._piles[0].scores)

    def test_a_new_peak_moves_the_fingerprint_and_merges_the_overlapped_pile(self):
        from search.climb import Climber
        from search.piles import Pile, note
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 8})
        step = self.space._assemble({"min_pruning": 5})
        ctx = self.ctx(_Map({
            (a.key, b.key): 0.20, (a.key, step.key): 0.04, (b.key, step.key): 0.08,
        }))
        ctx._piles = [Pile(b, 79.0), Pile(a, 80.0)]
        cl = Climber(ctx, 0, a, 80.0, origin="pile 2")
        ctx.climbers.append(cl)
        note(ctx, ctx._piles, step, 81.0, 81.0)
        self.assertEqual(len(ctx._piles), 1)
        self.assertEqual(ctx._piles[0].config.key, step.key)
        self.assertEqual(ctx._piles[0].score, 81.0)
        self.assertIn(a.key, ctx._piles[0].scores)
        self.assertIn(b.key, ctx._piles[0].scores)
        self.assertEqual(cl.origin, "pile 1")

    def test_a_peaked_fingerprint_is_not_seated_again(self):
        from search.piles import Pile, unseated
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 8})
        ctx = self.ctx(_Map({(a.key, b.key): 0.12}))
        piles = [Pile(a, 85.7893), Pile(b, 85.3558)]
        done = Climber(ctx, 0, a, 85.7893, origin="pile 1")
        done.state = "peak"
        ctx.climbers = [done]
        self.assertEqual([origin for _, _, origin in unseated(ctx, piles)], ["pile 2"])

    def test_the_same_fingerprint_stays_finished_when_its_score_ticks_up(self):
        from search.piles import Pile, unseated
        a = self.space._assemble({"min_pruning": 2})
        ctx = self.ctx(_Map({}))
        pile = Pile(a, 85.7900)
        done = Climber(ctx, 0, a, 85.7893, origin="pile 1")
        done.state = "peak"
        ctx.climbers = [done]
        self.assertEqual(unseated(ctx, [pile]), [])

    def test_a_fingerprint_similar_to_a_closed_one_is_not_seated(self):
        from search.piles import Pile, history_rows, unseated
        a = self.space._assemble({"min_pruning": 2})
        higher = self.space._assemble({"min_pruning": 4})
        ctx = self.ctx(_Map({(a.key, higher.key): 0.04}))
        pile = Pile(higher, 86.0)
        done = Climber(ctx, 0, a, 85.7893, origin="pile 1")
        done.state = "peak"
        ctx.climbers = [done]
        self.assertEqual(unseated(ctx, [pile]), [])
        self.assertTrue(history_rows(ctx)[0]["closed"])

    def test_six_piles_keep_two_candidates_for_a_fast_switch(self):
        from search.climb import Climber
        from search.piles import note, unseated

        class _Far:
            def between(self, a, b):
                return (0.0 if a.key == b.key else 0.20), "callset"

        ctx = self.ctx(_Far())
        ctx.climbers = []
        piles = []
        configs = []
        pads = ctx.ladders.rungs["assembly_region_padding"]
        for i in range(6):
            config = self.space._assemble({"assembly_region_padding": pads[i * 3]})
            configs.append(config)
            note(ctx, piles, config, 80.0 - 0.1 * i, 80.0)
        self.assertEqual([pile.config.key for pile in piles], [config.key for config in configs])
        extra = self.space._assemble({"assembly_region_padding": pads[-1]})
        note(ctx, piles, extra, 79.4, 80.0)
        self.assertEqual(len(piles), 6)
        for i, config in enumerate(configs[:4]):
            ctx.climbers.append(Climber(ctx, i, config, 80.0 - 0.1 * i, origin=f"pile {i + 1}"))
        waiting = unseated(ctx, piles)
        self.assertEqual([item[2] for item in waiting], ["pile 5", "pile 6"])
        ctx.climbers[0].state = "peak"
        self.assertEqual(unseated(ctx, piles)[0][2], "pile 5")
        nxt = self.space._assemble({"assembly_region_padding": pads[-2]})
        note(ctx, piles, nxt, 79.4, 80.0)
        self.assertEqual(piles[0].config.key, nxt.key)
        self.assertEqual([pile.config.key for pile in piles[1:]], [config.key for config in configs[1:]])

    def test_a_far_config_above_the_lowest_pile_takes_that_climber(self):
        from search.climb import Climber, switch_to_candidate
        from search.piles import note, unseated

        class _Far:
            def between(self, a, b):
                return (0.0 if a.key == b.key else 0.20), "callset"

        ctx = self.ctx(_Far())
        ctx.ga_opening = True
        pads = ctx.ladders.rungs["assembly_region_padding"]
        scores = [80.0, 80.4, 80.5, 80.7, 79.8, 79.7]
        piles = []
        configs = []
        for i, score in enumerate(scores):
            config = self.space._assemble({"assembly_region_padding": pads[i * 3]})
            configs.append(config)
            note(ctx, piles, config, score, 80.7)
        ctx._piles = piles
        ctx._pile_seen = {pile.config.key for pile in piles}
        ctx.climbers = [Climber(ctx, i, configs[i], scores[i], origin=f"pile {i + 1}") for i in range(4)]
        below = self.space._assemble({"assembly_region_padding": pads[-1]})
        note(ctx, piles, below, 79.5, 80.7)
        self.assertEqual([pile.config.key for pile in piles], [config.key for config in configs])
        far = self.space._assemble({"assembly_region_padding": pads[-2]})
        note(ctx, piles, far, 84.0, 84.0)
        self.assertEqual(piles[5].config.key, far.key)
        self.assertTrue(switch_to_candidate(ctx, ctx.climbers))
        self.assertEqual(ctx.climbers[0].centre.key, far.key)
        self.assertEqual(ctx.climbers[0].origin, "pile 6")
        self.assertEqual(ctx.climbers[1].origin, "pile 2")
        self.assertIn("pile 1", [item[2] for item in unseated(ctx, piles)])

    def test_a_higher_candidate_takes_the_lowest_climber(self):
        from search.climb import Climber, switch_to_candidate
        from search.piles import Pile
        low_c = self.space._assemble({"min_pruning": 2})
        high_c = self.space._assemble({"min_pruning": 8})
        cand = self.space._assemble({"min_dangling_branch_length": 8})
        below = self.space._assemble({"min_dangling_branch_length": 4})
        ctx = self.ctx(_Map({}))
        ctx.ga_opening = True
        ctx._piles = [Pile(low_c, 80.0), Pile(high_c, 90.0), Pile(below, 79.0)]
        ctx._pile_seen = {low_c.key, high_c.key, below.key}
        low = Climber(ctx, 0, low_c, 80.0, origin="pile 1")
        high = Climber(ctx, 1, high_c, 90.0, origin="pile 2")
        ctx.climbers = [low, high]
        self.assertFalse(switch_to_candidate(ctx, [low, high]))
        self.assertEqual(low.centre.key, low_c.key)
        ctx._piles[2] = Pile(cand, 85.0)
        ctx._pile_seen.add(cand.key)
        self.assertTrue(switch_to_candidate(ctx, [low, high]))
        self.assertEqual(low.centre.key, cand.key)
        self.assertEqual(low.score, 85.0)
        self.assertEqual(low.origin, "pile 3")
        self.assertEqual(high.centre.key, high_c.key)
        self.assertEqual(high.origin, "pile 2")

    def test_a_candidate_does_not_take_a_centre_another_climber_holds(self):
        from search.climb import Climber, switch_to_candidate
        from search.piles import Pile
        held = self.space._assemble({"min_pruning": 2})
        low_c = self.space._assemble({"min_pruning": 8})
        own = self.space._assemble({"min_dangling_branch_length": 8})
        ctx = self.ctx(_Map({}))
        ctx.ga_opening = True
        pile = Pile(own, 84.0)
        pile.members.append(held)
        pile.scores[held.key] = 90.0
        ctx._piles = [Pile(held, 85.0), Pile(low_c, 80.0), pile]
        ctx._pile_seen = {held.key, low_c.key, own.key}
        holder = Climber(ctx, 0, held, 85.0, origin="pile 1")
        low = Climber(ctx, 1, low_c, 80.0, origin="pile 2")
        ctx.climbers = [holder, low]
        self.assertTrue(switch_to_candidate(ctx, [holder, low]))
        self.assertEqual(low.centre.key, own.key)
        self.assertEqual(low.origin, "pile 3")
        self.assertEqual(holder.centre.key, held.key)
        self.assertEqual(holder.state, "climbing")
        low.centre, low.score, low.origin = low_c, 80.0, "pile 2"
        pile.scores[own.key] = 79.0
        pile.score = 79.0
        self.assertFalse(switch_to_candidate(ctx, [holder, low]))
        self.assertEqual(low.centre.key, low_c.key)

    def test_a_candidate_does_not_take_a_centre_a_peaked_climber_holds(self):
        from search.climb import Climber, switch_to_candidate
        from search.piles import Pile
        held = self.space._assemble({"min_pruning": 2})
        low_c = self.space._assemble({"min_pruning": 8})
        own = self.space._assemble({"min_dangling_branch_length": 8})
        ctx = self.ctx(_Map({(held.key, own.key): 0.20, (held.key, low_c.key): 0.20, (own.key, low_c.key): 0.20}))
        ctx.ga_opening = True
        pile = Pile(own, 84.0)
        pile.members.append(held)
        pile.scores[held.key] = 90.0
        ctx._piles = [Pile(held, 85.0), Pile(low_c, 80.0), pile]
        ctx._pile_seen = {held.key, low_c.key, own.key}
        holder = Climber(ctx, 0, held, 85.0, origin="pile 1")
        holder.state = "peak"
        low = Climber(ctx, 1, low_c, 80.0, origin="pile 2")
        ctx.climbers = [holder, low]
        self.assertTrue(switch_to_candidate(ctx, [low]))
        self.assertEqual(low.centre.key, own.key)
        self.assertEqual(low.origin, "pile 3")
        self.assertEqual(holder.state, "peak")
        self.assertEqual(holder.centre.key, held.key)
        low.centre, low.score, low.origin = low_c, 80.0, "pile 2"
        pile.scores[own.key] = 79.0
        pile.score = 79.0
        self.assertFalse(switch_to_candidate(ctx, [low]))
        self.assertEqual(low.centre.key, low_c.key)

    def test_fingerprint_history_keeps_a_pile_that_moved(self):
        from search.piles import history_rows, note
        a = self.space._assemble({"min_pruning": 2})
        b = self.space._assemble({"min_pruning": 4})
        ctx = self.ctx(_Map({(a.key, b.key): 0.04}))
        logged = []
        ctx.say = lambda msg, *args: logged.append(msg % args if args else msg)
        piles = []
        note(ctx, piles, a, 80.0, 80.0)
        note(ctx, piles, b, 81.0, 81.0)
        rows = history_rows(ctx)
        self.assertEqual([row["key"] for row in rows], [a.key[:8], b.key[:8]])
        self.assertIn(a.key[:8], logged[-1])
        self.assertIn(b.key[:8], logged[-1])

    def test_a_config_near_a_remembered_fingerprint_does_not_open_a_pile(self):
        from search.piles import note
        a = self.space._assemble({"min_pruning": 2})
        moved = self.space._assemble({"min_pruning": 4})
        again = self.space._assemble({"min_pruning": 8})
        ctx = self.ctx(_Map({
            (a.key, moved.key): 0.04, (a.key, again.key): 0.08, (moved.key, again.key): 0.12}))
        piles = []
        note(ctx, piles, a, 80.0, 81.0)
        note(ctx, piles, moved, 81.0, 81.0)
        note(ctx, piles, again, 79.0, 81.0)
        self.assertEqual([p.config.key for p in piles], [moved.key])


class _Map:
    def __init__(self, pairs):
        self.pairs = {}
        for (a, b), d in pairs.items():
            self.pairs[frozenset((a, b))] = d

    def between(self, a, b):
        if a.key == b.key:
            return 0.0, "callset"
        return self.pairs[frozenset((a.key, b.key))], "callset"

    def far(self, a, b, min_callset, min_settings):
        return self.between(a, b)[0] >= min_callset


if __name__ == "__main__":
    unittest.main()
