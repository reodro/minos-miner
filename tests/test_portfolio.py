"""search.portfolio: the opening's seeds - profile first, spaced apart, capped; and what a reviewed run adds."""
import json
import tempfile
import unittest
from pathlib import Path

from common.types import Config
from search import portfolio
from search.space import Space
from tests.helpers import ranges


class PortfolioTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges(), "gatk")
        self.d = self.space._assemble({}).params
        self.tmp = tempfile.TemporaryDirectory()
        self.ours = Path(self.tmp.name) / "ours.jsonl"

    def tearDown(self):
        self.tmp.cleanup()

    def member(self, profile, rank=1, src="live:20260910T1000000", **change):
        return {"src": src, "profile": profile, "rank": rank, "config": {**self.d, **change}}

    def test_the_legacy_portfolio_loads_and_fits_the_space(self):
        m = portfolio.members(ours=self.ours)
        self.assertGreater(len(m), 700)
        seeds = portfolio.pick(self.space, "B", 50, pool=m)
        self.assertEqual(len(seeds), 50)
        self.assertEqual(len({c.key for c in seeds}), 50)

    FAR1 = dict(min_mapping_quality_score=60, max_reads_per_alignment_start=300, min_dangling_branch_length=20)
    FAR2 = dict(min_mapping_quality_score=0, max_reads_per_alignment_start=25, min_pruning=10)

    def test_profile_first_spacing_and_cap(self):
        pool = [self.member("A", min_pruning=10, heterozygosity=0.01, **{k: v for k, v in self.FAR1.items()}),
                self.member("B", **self.FAR1),
                self.member("B", heterozygosity=0.00101, **self.FAR1),        # a near-copy of the one before
                self.member("B", **self.FAR2)]
        got = portfolio.pick(self.space, "B", 10, pool=pool)
        self.assertEqual(len(got), 3)                                     # the near-copy is skipped
        self.assertEqual([c.params["min_mapping_quality_score"] for c in got[:2]], [60, 0])   # B first
        self.assertEqual(got[-1].params["heterozygosity"], 0.01)          # A fills up last
        self.assertEqual(len(portfolio.pick(self.space, "B", 1, pool=pool)), 1)
        self.assertEqual(portfolio.pick(self.space, "B", 0, pool=pool), [])

    def test_a_reviewed_run_adds_spaced_near_best_configs(self):
        configs = {}
        for change, score in [(self.FAR1, 88.70), ({**self.FAR1, "min_pruning": 3}, 88.69), (self.FAR2, 88.60),
                              ({"min_pruning": 9, "heterozygosity": 0.01}, 87.00)]:
            c = Config("gatk", {**self.d, **change})
            configs[c.key] = {"config": c.to_dict(), "local": score}
        n = portfolio.add_run(self.space, configs, "20261002T092900Z", "B", path=self.ours)
        rows = [json.loads(l) for l in self.ours.read_text().splitlines()]
        self.assertEqual(n, len(rows))
        self.assertEqual([r["local"] for r in rows], [88.70, 88.60])      # the best first; 88.69 is a near-copy
        self.assertNotIn(87.0, [r["local"] for r in rows])                 # outside the 0.3 window
        self.assertTrue(all(r["src"] == "minos:20261002T092900Z" and r["profile"] == "B" for r in rows))
        none = Path(self.ours).parent / "no_kings.json"
        self.assertEqual(portfolio.members(ours=self.ours, kings=none)[0]["src"], "minos:20261002T092900Z")  # ours first
        kings = Path(self.ours).parent / "kings.json"
        kings.write_text(json.dumps({"members": [{"src": "legacy-king:20260915T10490000", "config": {}}]}))
        got = portfolio.members(ours=self.ours, kings=kings)
        self.assertEqual([got[0]["src"], got[1]["src"]], ["legacy-king:20260915T10490000", "minos:20261002T092900Z"])
        self.assertEqual(portfolio.add_run(self.space, configs, "20261002T092900Z", "B", path=self.ours), 0)
        self.assertTrue(all(r["group"] == "score" for r in rows))

    def test_a_reviewed_run_also_keeps_high_core_and_fewest_fp(self):
        configs = {}
        specs = [
            (self.FAR1, 88.70, 0.80, 1.0),
            (self.FAR2, 80.00, 0.95, 1.1),                                          # below the score window, highest core
            ({"min_pruning": 9, "heterozygosity": 0.01, "assembly_region_padding": 400}, 79.00, 0.945, 0.2),
            ({"min_pruning": 1, "heterozygosity": 0.0001, "assembly_region_padding": 10}, 78.00, 0.70, 0.1),
        ]
        for change, score, core, fp in specs:
            c = Config("gatk", {**self.d, **change})
            configs[c.key] = {"config": c.to_dict(), "local": score, "core": core, "fp_per_target": fp}
        portfolio.add_run(self.space, configs, "20261002T092900Z", "B", path=self.ours)
        rows = [json.loads(l) for l in self.ours.read_text().splitlines()]
        self.assertEqual({r["group"] for r in rows}, {"score", "core", "fp"})
        self.assertEqual(next(r["local"] for r in rows if r["group"] == "core"), 80.00)
        # fewest false positives among the cores within CORE_BAND of the best: not the 0.70-core config with fewer
        self.assertEqual([r["local"] for r in rows if r["group"] == "fp"], [79.00])
        seeds = portfolio.pick(self.space, "B", 3, pool=rows)
        self.assertEqual(len(seeds), 3)

    def test_a_replay_never_seeds_itself(self):
        rows = [self.member("B", src="minos:20261002T092900Z", min_pruning=9),
                self.member("B", src="live:20261002T0929000", **self.FAR1),
                self.member("B", src="minos:20261003T101500Z", **self.FAR2)]
        seeds = portfolio.pick(self.space, "B", 5, pool=rows, exclude_round="20261002T092900Z")
        self.assertEqual(len(seeds), 1)
        self.assertEqual(seeds[0].params["min_pruning"], self.FAR2["min_pruning"])

    def test_no_bam_no_profile(self):
        self.assertIsNone(portfolio.bam_profile("/nonexistent.bam", "chr20", 1, 1000))


if __name__ == "__main__":
    unittest.main()
