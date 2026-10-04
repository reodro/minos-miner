import math
import unittest

from common.scoring import score_from_parts, v2_score
from common.types import CLASSES, ClassCounts, Config, Round, ScoreResult
from tests.helpers import ranges


class TypesTest(unittest.TestCase):
    def test_config_key_ignores_param_order(self):
        a = Config("gatk", {"min_pruning": 3, "pcr_indel_model": "NONE"})
        b = Config("gatk", {"pcr_indel_model": "NONE", "min_pruning": 3})
        self.assertEqual(a.key, b.key)
        self.assertNotEqual(a.key, Config("gatk", {"min_pruning": 4, "pcr_indel_model": "NONE"}).key)

    def test_submission_shape(self):
        self.assertEqual(Config("gatk", {"min_pruning": 3}).to_submission(),
                         {"tool": "gatk", "gatk_options": {"min_pruning": 3}})

    def test_round_tag_and_region(self):
        r = Round("2026-09-28T08:17:00+00:00", "chr14", 45912001, 53254000, "", "")
        self.assertEqual(r.tag, "20260928T081700Z")
        self.assertEqual(r.region, "chr14:45912001-53254000")
        self.assertEqual(Round.parse_region("chr14:45,912,001-53,254,000"), ("chr14", 45912001, 53254000))

    def test_score_result_roundtrip(self):
        r = ScoreResult("r", "k", 2, {k: ClassCounts(1, 2, 3) for k in CLASSES}, 10, 600, core=0.5, score=50.0)
        self.assertEqual(ScoreResult.from_dict(r.to_dict()), r)


class ParamsTest(unittest.TestCase):
    def setUp(self):
        self.r = ranges()

    def test_defaults_are_valid(self):
        for tool in self.r.tools:
            self.assertEqual(self.r.errors(self.r.defaults(tool)), [], tool)

    def test_rejects_out_of_range_unknown_and_wrong_type(self):
        self.assertTrue(self.r.errors(Config("gatk", {"min_pruning": 11})))
        self.assertTrue(self.r.errors(Config("gatk", {"no_such_param": 1})))
        self.assertTrue(self.r.errors(Config("gatk", {"min_pruning": 2.5})))
        self.assertTrue(self.r.errors(Config("gatk", {"recover_all_dangling_branches": 1})))
        self.assertTrue(self.r.errors(Config("gatk", {"dont_use_soft_clipped_bases": True})))  # platform allows only False
        self.assertTrue(self.r.errors(Config("freebayes", {})))

    def test_cross_constraints_gatk_rejects(self):
        d = self.r.defaults("gatk").params
        self.assertTrue(self.r.errors(Config("gatk", {**d, "min_assembly_region_size": 215, "max_assembly_region_size": 162})))
        self.assertTrue(self.r.errors(Config("gatk", {**d, "base_quality_score_threshold": 5})))
        self.assertEqual(self.r.errors(Config("gatk", {**d, "base_quality_score_threshold": 6})), [])


class ScoringTest(unittest.TestCase):
    """Checks the wrapper around the upstream formula, not the formula itself."""

    def counts(self, tp=100, fp=0, fn=0):
        return {k: ClassCounts(tp, fp, fn) for k in CLASSES}

    def test_perfect_callset_scores_100(self):
        core, score, gate = v2_score(self.counts(), region_fp=0, n_targets=600)
        self.assertAlmostEqual(core, 1.0)
        self.assertAlmostEqual(score, 100.0)
        self.assertTrue(gate)

    def test_region_fp_enters_germline_term(self):
        _, score, _ = v2_score(self.counts(), region_fp=600, n_targets=600)
        self.assertAlmostEqual(score, 100 * (0.70 + 0.30 * math.exp(-1 / 8)))

    def test_score_from_parts_uses_the_upstream_weights(self):
        _, score, _ = v2_score(self.counts(), region_fp=600, n_targets=600)
        self.assertAlmostEqual(score_from_parts(1.0, 1.0), score)

    def test_failed_gate_scores_zero(self):
        ratios = {"titv_truth_snp": 2.0, "titv_query_snp": 1.0, "hethom_truth_snp": 1.5, "hethom_query_snp": 1.5}
        _, score, gate = v2_score(self.counts(), 0, 600, ratios)
        self.assertFalse(gate)
        self.assertEqual(score, 0.0)


if __name__ == "__main__":
    unittest.main()
