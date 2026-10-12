"""v11: hills (grouping, identity across refits, funding, convergence) and a whole virtual round."""
import random
import unittest

from common.tuning import DEFAULTS, TUNING
from search.space import Space
from search.v9.steps import Ladders
from search.v10.forest import Forest
from search.v11.engine import Hill, Hills, converged, endgame_hills, valley_verdict
from tests.helpers import ranges
from tests.test_v10 import ROUND, Stub

MP, PAD, MBQ, MMQ, APT = ("min_pruning", "assembly_region_padding", "min_base_quality_score",
                          "min_mapping_quality_score", "active_probability_threshold")


class HillsTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.lad = Ladders(self.space)
        self.forest = Forest(self.space, self.lad, trees=10, jobs=1)      # for its list of searched settings

    def far_apart(self, base, k):
        """A config k settings away from `base`, each moved to the other end of its ladder."""
        p = dict(base.params)
        names = [n for n in self.forest.names if len(self.lad.rungs[n]) > 4][:k]
        for n in names:
            r = self.lad.rungs[n]
            p[n] = r[-1] if self.lad.index(n, p[n]) < len(r) / 2 else r[0]
        return self.space._assemble(p)

    def test_near_configs_share_a_hill_and_far_ones_found_another(self):
        a = self.space._assemble({})
        a2 = self.space._assemble({**a.params, PAD: a.params[PAD] + 2})       # one small move: same hill
        b = self.far_apart(a, 6)                                               # six settings far: another hill
        h = Hills(self.space)
        h.update([(a, 90.0, 0.9, 100), (a2, 89.9, 0.9, 100), (b, 89.5, 0.9, 100)], {}, now=60.0, gain=0.02)
        self.assertEqual(len(h.hills), 2)
        sizes = sorted(len(x.members) for x in h.hills.values())
        self.assertEqual(sizes, [1, 2])

    def test_hills_are_grouped_by_calls_not_settings(self):
        import base64, struct
        from search.v9.distance import WEIGHTS

        def sketch(vals):
            raw = struct.pack(f">{len(vals)}I", *vals)
            return {"k": 128, "h": {c: base64.b64encode(raw).decode() for c in WEIGHTS}}
        a = self.space._assemble({})
        b = self.far_apart(a, 6)                       # far in settings ...
        c = self.space._assemble({**a.params, PAD: a.params[PAD] + 1})    # ... and one rung away
        same, other = sketch(list(range(1, 129))), sketch(list(range(1000, 1128)))
        h = Hills(self.space)
        h.update([(a, 90.0, 0.9, 100), (b, 89.9, 0.9, 100), (c, 89.8, 0.9, 100)],
                 {a.key: same, b.key: same, c.key: other}, now=60.0, gain=0.02)
        groups = sorted(sorted(m[0].key for m in x.members) for x in h.hills.values())
        self.assertEqual(groups, sorted([sorted([a.key, b.key]), [c.key]]))   # same calls: one hill; other calls: two

    def test_a_hill_keeps_its_id_and_its_gain_clock(self):
        a = self.space._assemble({})
        b = self.far_apart(a, 6)
        h = Hills(self.space)
        h.update([(a, 90.0, 0.9, 100), (b, 89.5, 0.9, 100)], {}, now=60.0, gain=0.02)
        ids = {x.founder.key: x.id for x in h.hills.values()}
        a3 = self.space._assemble({**a.params, PAD: a.params[PAD] + 1})
        h.update([(a3, 90.1, 0.9, 100), (a, 90.0, 0.9, 100), (b, 89.5, 0.9, 100)], {}, now=600.0, gain=0.02)
        by = {x.id: x for x in h.hills.values()}
        self.assertEqual(by[ids[a.key]].founder.key, a3.key)       # the same hill, a new best
        self.assertEqual(by[ids[a.key]].last_gain, 600.0)          # +0.1 is a gain
        self.assertEqual(by[ids[b.key]].last_gain, 60.0)           # b did not move

    def test_funding_drops_a_stale_low_hill(self):
        h = Hills(self.space)
        c = self.space._assemble({})
        h.hills = {0: Hill(0, c, 90.0, 0, 90.0, 900), 1: Hill(1, c, 89.0, 0, 89.0, 0), 2: Hill(2, c, 89.8, 0, 89.8, 0)}
        got = [x.id for x in h.fund(now=1000.0, count=3, drop=0.5, patience=600.0)]
        self.assertEqual(got, [0, 2])                              # hill 1: 1.0 below and quiet for 1000 s

    def test_converged_needs_every_funded_hill_quiet(self):
        c = self.space._assemble({})
        quiet, busy = Hill(0, c, 90, 0, 90, 0), Hill(1, c, 89, 0, 89, 1500)
        self.assertFalse(converged([quiet, busy], now=1800, after=1200, window=600))
        self.assertTrue(converged([quiet], now=1800, after=1200, window=600))
        self.assertFalse(converged([quiet], now=1000, after=1200, window=600))


class ValleyTest(unittest.TestCase):
    def test_the_verdict(self):
        self.assertEqual(valley_verdict(82.6, 82.4, [82.35], 0.1, 5), "same")         # within 0.1 of the lower
        self.assertEqual(valley_verdict(82.6, 82.4, [82.0] * 5, 0.1, 5), "valley")
        self.assertEqual(valley_verdict(82.6, 82.4, [82.0] * 4, 0.1, 5), "unknown")    # not enough to call it

    def sketched(self):
        import base64, struct
        from search.v9.distance import WEIGHTS

        def sketch(vals):
            raw = struct.pack(f">{len(vals)}I", *sorted(vals))
            return {"k": 128, "h": {c: base64.b64encode(raw).decode() for c in WEIGHTS}}
        space = Space(ranges())
        a = space._assemble({})
        b = space._assemble({**a.params, PAD: a.params[PAD] + 40})
        mids = [space._assemble({**a.params, PAD: a.params[PAD] + 10 + i}) for i in range(5)]
        A, B = list(range(1, 129)), list(range(1001, 1129))
        M = list(range(1, 65)) + list(range(1001, 1065))           # half of each: 0.5 from both, they are 1.0 apart
        raw = {a.key: sketch(A), b.key: sketch(B), **{m.key: sketch(M) for m in mids}}
        return space, a, b, mids, raw

    def test_a_high_config_between_two_hills_merges_them(self):
        space, a, b, mids, raw = self.sketched()
        h = Hills(space)
        h.update([(a, 90.0, 0.9, 100), (b, 89.9, 0.9, 100), (mids[0], 89.85, 0.9, 100)], raw, now=60.0, gain=0.02)
        homes = [x for x in h.hills.values() if any(m[0].key == a.key for m in x.members)]
        self.assertTrue(any(m[0].key == b.key for m in homes[0].members))     # a and b: one hill
        self.assertTrue(h.events and "no valley" in h.events[0])

    def test_low_configs_between_two_hills_keep_them_apart(self):
        space, a, b, mids, raw = self.sketched()
        h = Hills(space)
        pts = [(a, 90.0, 0.9, 100), (b, 89.9, 0.9, 100)] + [(m, 88.0, 0.9, 100) for m in mids]
        h.update(pts, raw, now=60.0, gain=0.02)
        tops = sorted(h.hills.values(), key=lambda x: -x.best)[:2]
        self.assertEqual({tops[0].founder.key, tops[1].founder.key}, {a.key, b.key})   # still two hills
        self.assertTrue(any("valley, kept apart" in e for e in h.events))


class EndgameTest(unittest.TestCase):
    def test_only_hills_within_the_window_of_the_leader_are_climbed(self):
        c = Space(ranges())._assemble({})
        hs = [Hill(0, c, 82.40, 0, 0, 0), Hill(1, c, 82.65, 0, 0, 0), Hill(2, c, 82.36, 0, 0, 0), Hill(3, c, 82.10, 0, 0, 0)]
        self.assertEqual([h.id for h in endgame_hills(hs, 4, 0.3)], [1, 0, 2])
        self.assertEqual([h.id for h in endgame_hills(hs, 2, 0.3)], [1, 0])
        self.assertEqual(endgame_hills([], 4, 0.3), [])


class VirtualRoundTest(unittest.TestCase):
    def test_a_whole_round_funds_several_hills(self):
        from evaluator import FakeRunner, FakeScorer, SimExecutor
        from main import DryRunSubmitter, LoopSettings, run_round
        from search.v11.engine import searcher
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=3)
        ex = SimExecutor(runner, 64, duration=lambda job: 150.0)
        seen = []
        inner = searcher(r, 64, seed=1)

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
        self.assertEqual(eng["version"], "v11")
        rf = eng["rf"]
        self.assertGreater(rf["refits"], 3)
        self.assertGreaterEqual(len(rf["hills"]), 2)                    # more than one hill funded to the end
        self.assertGreater(sum(h["picks"] > 0 for h in rf["hills"]), 1)  # and more than one actually worked
        self.assertIn("cross", rf["picks"])
        self.assertTrue(eng["basins"])                                  # the endgame climbed
        self.assertIsNotNone(out.best)


if __name__ == "__main__":
    unittest.main()
