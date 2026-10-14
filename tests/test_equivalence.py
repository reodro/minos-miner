import unittest

from common.types import Config
from evaluator.equivalence import BQST, CONT, MBQ, MMQ, MRPAS, RoundFacts, _class_low, canonical

QUALS = frozenset({2, 14, 15, 16, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40})   # round 20261004T092900Z
FACTS = RoundFacts(QUALS, frozenset(range(61)), 22)


class EquivalenceTest(unittest.TestCase):
    """Two settings share a class exactly when no value GATK compares them with lies between them."""

    def test_min_base_quality_classes_follow_the_closed_interval_rule(self):
        cq = FACTS.compared_quals                       # the BAM's qualities plus 0 and 20 (overlapping mates)
        self.assertEqual(cq[:6], [0, 2, 14, 15, 16, 20])
        for a in range(60):
            for b in range(a, 60):
                same = a == b or not any(a <= v <= b for v in cq)
                self.assertEqual(same, _class_low(a, cq, True) == _class_low(b, cq, True), (a, b))

    def test_quality_threshold_classes_follow_the_half_open_rule_after_the_mapq_cap(self):
        for mmq in (0, 20, 45):
            capped = sorted({min(q, m) for q in FACTS.compared_quals for m in FACTS.mapqs if m >= mmq})
            for a in range(6, 60):
                for b in range(a, 60):
                    same = a == b or not any(a <= v < b for v in capped)
                    self.assertEqual(same, max(_class_low(a, capped, False), 6) == max(_class_low(b, capped, False), 6))

    def test_canonical(self):
        c = Config("gatk", {MRPAS: 143, MBQ: 12, BQST: 13, MMQ: 20, "min_pruning": 2})
        self.assertEqual(canonical(c, FACTS).params, {MRPAS: 22, MBQ: 3, BQST: 6, MMQ: 20, "min_pruning": 2})
        self.assertIs(canonical(c, None), c)                                       # no facts: unchanged
        kept = Config("gatk", {MRPAS: 20, MBQ: 14, BQST: 4, MMQ: 20})              # below the max / a compared value /
        self.assertIs(canonical(kept, FACTS), kept)                                # a threshold GATK rejects
        self.assertEqual(canonical(Config("gatk", {MRPAS: 300}), FACTS).key,
                         canonical(Config("gatk", {MRPAS: 25}), FACTS).key)        # no start has more than 22 reads

    def test_contamination_fractions_that_remove_no_reads_are_one_run(self):
        # Depths of 20261010T044200Z: max 71, so every fraction below 1/71 removes 0 reads.
        depths = frozenset(i for i in range(1, 72) if i not in (64, 65))
        facts = RoundFacts(QUALS, frozenset(range(61)), 22, depths=depths)
        zero = canonical(Config("gatk", {CONT: 0.01}), facts)
        self.assertEqual(zero.params[CONT], 0.0)
        self.assertEqual(canonical(Config("gatk", {CONT: 0.0125}), facts).params[CONT], 0.0)
        kept = canonical(Config("gatk", {CONT: 0.01875}), facts)          # removes a read at depth 54; not the zero class
        self.assertEqual(kept.params[CONT], 0.01875)
        self.assertNotEqual(canonical(Config("gatk", {CONT: 0.045}), facts).params[CONT],
                            canonical(Config("gatk", {CONT: 0.0475}), facts).params[CONT])
        bare = canonical(Config("gatk", {CONT: 0.01}), FACTS)              # no depths recorded: the fraction stays
        self.assertEqual(bare.params[CONT], 0.01)

    def test_a_result_that_ran_as_another_config_is_also_found_under_its_own_name(self):
        import shutil
        import tempfile
        from pathlib import Path
        from evaluator.local import _link_as
        d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d, True)
        (d / "canon_full.vcf.gz").write_text("v")
        (d / "canon_full.vcf.gz.tbi").write_text("i")
        _link_as(str(d / "canon_full.vcf.gz"), "canon", "asked")
        self.assertEqual((d / "asked_full.vcf.gz").read_text(), "v")       # main.review finds <key>_full.vcf.gz
        self.assertEqual((d / "asked_full.vcf.gz.tbi").read_text(), "i")
        _link_as("fake://canon", "canon", "asked")                          # not a file: nothing to do


if __name__ == "__main__":
    unittest.main()
