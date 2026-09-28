import gzip
import os
import unittest

from answer_set import Provider
from answer_set.fingerprint import Fingerprint, best_match
from answer_set.germline import alleles, keep, parse_query
from common.settings import SETTINGS
from common.types import Regime, Round, Variant

ROUND_0817 = os.path.join(SETTINGS.rounds_dir, "20260928T081700Z")


class GermlineRuleTest(unittest.TestCase):
    def test_keep(self):
        self.assertTrue(keep("PASS", ["A"], gq=20, dp=10))
        self.assertFalse(keep("PASS", ["A"], gq=19, dp=10))
        self.assertFalse(keep("PASS", ["A"], gq=20, dp=9))
        self.assertFalse(keep("VQSRTrancheSNP99.80to100.00", ["A"], gq=99, dp=30))
        self.assertFalse(keep("PASS", ["A", "*"], gq=99, dp=30))

    def test_alleles(self):
        self.assertEqual(alleles("chr1", 5, "T", ["C"], "1/1"), [Variant("chr1", 5, "T", "C", "1/1")])
        self.assertEqual(alleles("chr1", 5, "TA", ["T", "TAA"], "1|2"),
                         [Variant("chr1", 5, "TA", "T", "0/1"), Variant("chr1", 5, "TA", "TAA", "0/1")])
        self.assertEqual(alleles("chr1", 5, "T", ["C"], "0/0"), [])

    def test_parse_query(self):
        lines = ["100\tA\tG\tPASS\t0/1\t30\t99", "200\tA\tG\tPASS\t0/1\t5\t99", "300\tA\tG,*\tPASS\t1/1\t30\t99"]
        self.assertEqual(parse_query("chr1", lines), [Variant("chr1", 100, "A", "G", "0/1")])


class FingerprintTest(unittest.TestCase):
    def test_best_match(self):
        sites = {1: ("A", "G", [0, 1, 2]), 2: ("C", "T", [1, 1, 0]), 3: ("G", "A", [2, 0, -1])}
        fp = best_match(["s0", "s1", "s2"], sites, {1: 1, 2: 1, 3: 0})
        self.assertEqual((fp.donor, fp.score), ("s1", 1.0))
        self.assertTrue(Fingerprint("x", 0.99, "y", 0.7, 50).decisive)
        self.assertFalse(Fingerprint("x", 0.69, "y", 0.68, 50).decisive)


def _norm(pos, ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


@unittest.skipUnless(os.environ.get("RUN_SLOW") and os.path.exists(os.path.join(ROUND_0817, "input.bam")),
                     "set RUN_SLOW=1 (needs round 08:17 and 1000G data, ~1 min)")
class ProviderEndToEndTest(unittest.TestCase):
    def test_1000g_round_germline_is_exact(self):
        r = Round("2026-09-28T08:17:00+00:00", "chr14", 45912001, 53254000,
                  bam=os.path.join(ROUND_0817, "input.bam"), reference="")
        a = Provider().build(r)
        self.assertEqual((a.regime, a.donor), (Regime.KG, "NA18878"))
        truth = {}
        for line in gzip.open(os.path.join(ROUND_0817, "truth.vcf.gz"), "rt"):
            if line.startswith("#") or "SYNTHETIC" in line.split("\t")[7]:
                continue
            p = line.rstrip("\n").split("\t")
            for v in alleles(p[0], int(p[1]), p[3], p[4].split(","), p[9].split(":")[0]):
                truth[_norm(v.pos, v.ref, v.alt)] = v.gt
        ours = {_norm(v.pos, v.ref, v.alt): v.gt for v in a.germline}
        self.assertEqual(ours, truth)


if __name__ == "__main__":
    unittest.main()
