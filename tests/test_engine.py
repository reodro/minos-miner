"""search engine v9: ladders, callset distance, Step A's groups, Step B's grid and growth, Step C's children, and a
whole virtual round."""
import asyncio
import gzip
import os
import random
import tempfile
import unittest

from common.types import AnswerSet, Config, Regime, Round, Variant
from evaluator import FakeRunner, FakeScorer, SimExecutor
from evaluator.sketch import callset_sketch
from main import DryRunSubmitter, LoopSettings, run_round
from search import engine, far
from search.candidates import Point, pareto_layers
from search.climb import Climber
from search.distance import callset_distance, setting_distance
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges

ROUND = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")


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

    def propose(self, batch, tag, batches):
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


def pair_ends(cl):
    """Riders the centre's run loses because a cheap setting sits at the end of its ladder (or a pair collapses)."""
    lad, c = cl.ctx.ladders, cl.start
    singles = sum(1 for n in cl.cheap_names for k in (1, -1) if lad.shift(n, c.params[n], k) is not None)
    full = 2 * len(cl.cheap_names)
    pairs = sum(1 for i, n in enumerate(cl.cheap_names) for m in cl.cheap_names[i + 1:]
                for a in (1, -1) for b in (1, -1)
                if lad.shift(n, c.params[n], a) is not None and lad.shift(m, c.params[m], b) is not None)
    full_pairs = 4 * len(cl.cheap_names) * (len(cl.cheap_names) - 1) // 2
    return (full - singles) + (full_pairs - pairs)


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

    def test_a_read_pass_grid_is_every_combination_of_five_settings(self):
        ctx, cl = self.climber()
        kind, group = cl.next_group()
        self.assertEqual((kind, len(group)), ("read", 5))
        expected = 1
        for n in group:
            expected *= len(cl.values(n, cl.start.params[n]))
        asyncio.run(cl.grid(kind, group, until=1e9))
        heads = [b[0] for b in ctx.sent]
        self.assertEqual(len({h.key for h in heads}), len(heads))
        self.assertEqual(len(heads), expected)                        # every combination, the centre included
        cheap = len(cl.cheap_names)
        self.assertEqual(len(ctx.sent[0]) - 1, 2 * cheap + 4 * cheap * (cheap - 1) // 2 - pair_ends(cl))   # the centre's run: singles and pairs
        self.assertTrue(all(len(b) - 1 <= 2 * cheap for b in ctx.sent[1:]))                     # the others: singles

    def test_the_climb_reaches_the_peak_and_grows_its_step(self):
        ctx, cl = self.climber()
        asyncio.run(cl.run(until=1e9))
        self.assertEqual(cl.score, 100.0)
        self.assertEqual({n: cl.centre.params[n] for n in self.target}, self.target)
        self.assertEqual(cl.state, "peak")



MDBL, MP = "min_dangling_branch_length", "min_pruning"


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
        want = far.plan(self.ctx, peak, known, coarse=3)[MDBL]
        self.assertIn(7, want)
        self.assertIn(8, want)
        self.assertNotIn(5, want)
        self.assertNotIn(9, want)
        last = len(self.rungs) - 1
        self.assertEqual([i for i in want if i > 9], far.spread(10, last, 3))           # beyond 9, coarse
        self.assertEqual([i for i in want if i < 5], far.spread(0, 4, 3))               # beyond 5, coarse

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
        self.assertEqual(eng["version"], "v9")
        self.assertEqual(eng["phase"], "picked")
        first = eng["basins"][0]
        self.assertGreater(first["best"], first["opened"])                  # Step B climbed
        self.assertTrue(first["best_key"] and first["first_key"])
        self.assertGreaterEqual(round(out.best.score, 4), first["best"])
        self.assertGreater(len(eng["waves"]), 3)                            # wave 0 and three generations


if __name__ == "__main__":
    unittest.main()
