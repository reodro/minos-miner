import random
import unittest

from common.scoring import finalize
from common.types import CLASSES, ClassCounts, Config, ScoreResult
from search.score2 import PASSED, Second, gate_slope, second
from utils.scoring import plausibility_gate

TRUTH = {"titv_truth_snp": 2.08, "hethom_truth_snp": 1.89}


def result(key: str, counts: dict, region_fp: int, n: int = 500, ratios: dict = None) -> ScoreResult:
    r = ScoreResult("2026-01-01T00:00:00+00:00", key, 1, {k: ClassCounts(**v) for k, v in counts.items()},
                    region_fp, n, dict(ratios if ratios is not None else {**TRUTH, "titv_query_snp": 2.0,
                                                                          "hethom_query_snp": 1.8}))
    return finalize(r)


def counts(rng: random.Random) -> dict:
    return {k: {"tp": rng.randint(50, 900), "fp": rng.randint(0, 200), "fn": rng.randint(0, 400)} for k in CLASSES}


def config(i: int, **extra) -> Config:
    return Config("gatk", {"min_base_quality_score": 10 + i % 30, "heterozygosity": 0.001 * (1 + i), **extra})


class SecondScoreTest(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(7)

    def scored(self, n: int):
        out = []
        for i in range(n):
            c = config(i)
            r = result(c.key, counts(self.rng), self.rng.randint(100, 1500))
            out.append((r, second(r, c)))
        return out

    def test_keeps_the_strict_order_of_the_base_score(self):
        rows = self.scored(300)
        for ra, a in rows:
            for rb, b in rows:
                if ra.score > rb.score:
                    self.assertGreater(a, b)

    def test_rises_and_falls_with_the_base_score_along_a_path(self):
        rows = self.scored(400)
        for (ra, a), (rb, b) in zip(rows, rows[1:]):
            if rb.score != ra.score:
                self.assertEqual(rb.score > ra.score, b > a)

    def test_breaks_every_tie_of_the_base_score(self):
        same = counts(self.rng)
        rows = []
        for i in range(50):
            c = config(i, max_alternate_alleles=3 + i)
            rows.append(second(result(c.key, same, 900), c))
        self.assertEqual(len({r.s for r in rows}), 1)
        self.assertEqual(len(set(rows)), len(rows))
        self.assertEqual(len({r.informative for r in rows}), 1)

    def test_the_same_config_has_the_same_second_score(self):
        c = config(1)
        r = result(c.key, counts(self.rng), 700)
        self.assertEqual(second(r, c), second(r, Config("gatk", dict(c.params))))

    def test_a_gated_config_sits_below_every_passing_one_and_slopes_toward_the_gate(self):
        c = counts(self.rng)
        ok = second(result("ok", c, 900), config(1))
        near = second(result("near", c, 900, ratios={**TRUTH, "titv_query_snp": 2.08 + 0.7,
                                                     "hethom_query_snp": 1.8}), config(2))
        far = second(result("far", c, 900, ratios={**TRUTH, "titv_query_snp": 2.08 + 1.5,
                                                   "hethom_query_snp": 1.8}), config(3))
        empty = second(result("empty", c, 900, ratios={**TRUTH, "titv_query_snp": 0.0,
                                                       "hethom_query_snp": 0.0}), config(4))
        self.assertEqual((near.s, far.s, empty.s), (0.0, 0.0, 0.0))
        self.assertEqual(ok.gate, PASSED)
        self.assertGreater(ok, near)
        self.assertGreater(near, far)
        self.assertGreater(far, empty)

    def test_equal_violations_are_ordered_by_the_ungated_score(self):
        bad = {**TRUTH, "titv_query_snp": 3.0, "hethom_query_snp": 1.8}
        c = counts(self.rng)
        fewer_fp = second(result("a", c, 300, ratios=bad), config(1))
        more_fp = second(result("b", c, 1200, ratios=bad), config(2))
        self.assertEqual(fewer_fp.gate[:2], more_fp.gate[:2])
        self.assertGreater(fewer_fp, more_fp)

    def test_the_gate_slope_reads_the_ratios_like_the_validator(self):
        rng = random.Random(3)
        for _ in range(500):
            ratios = {"titv_truth_snp": rng.choice([0.0, rng.uniform(1.5, 2.5)]),
                      "titv_query_snp": rng.choice([0.0, rng.uniform(0.5, 3.5)]),
                      "hethom_truth_snp": rng.choice([0.0, rng.uniform(1.0, 2.5)]),
                      "hethom_query_snp": rng.choice([0.0, rng.uniform(0.0, 4.5)])}
            r = result("r", counts(rng), 800, ratios=ratios)
            slope = gate_slope(r)
            self.assertEqual(plausibility_gate(ratios), r.gate_passed)
            if r.gate_passed:
                self.assertEqual(slope, PASSED)
            else:
                self.assertTrue(slope[0] < 0 or slope[1] < 0, (ratios, slope))

    def test_no_score_has_no_second_score(self):
        r = result("x", counts(self.rng), 800)
        r.score = None
        self.assertIsNone(second(r, config(1)))
        self.assertIsNone(second(None, config(1)))

    def test_informative_drops_only_the_code(self):
        r = result("x", counts(self.rng), 800)
        s = second(r, config(1))
        self.assertIsInstance(s, Second)
        self.assertEqual(s.informative, (s.s, s.gate))


if __name__ == "__main__":
    unittest.main()
