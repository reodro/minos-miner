import gzip
import os
import unittest

from answer_set import Provider
from answer_set.fingerprint import Fingerprint, best_match
from answer_set.germline import canonical_gt, keep, parse_query
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

    def test_parse_query_keeps_whole_records(self):
        lines = ["100\tA\tG\tPASS\t0/1\t30\t99",           # kept
                 "200\tA\tG\tPASS\t0/1\t5\t99",            # DP < 10
                 "300\tA\tG,*\tPASS\t1/1\t30\t99",          # spanning-deletion allele
                 "400\tTA\tTAA,T,TAAA\tPASS\t3|1\t30\t99",   # multi-allelic: one record, every ALT, unphased GT
                 "500\tC\tT\tPASS\t0/0\t30\t99"]            # not the donor's variant
        self.assertEqual(parse_query("chr1", lines), [Variant("chr1", 100, "A", "G", "0/1"),
                                                      Variant("chr1", 400, "TA", "TAA,T,TAAA", "1/3")])
        self.assertEqual(canonical_gt("2|1"), "1/2")


class FingerprintTest(unittest.TestCase):
    def test_best_match(self):
        sites = {1: ("A", "G", [0, 1, 2]), 2: ("C", "T", [1, 1, 0]), 3: ("G", "A", [2, 0, -1])}
        fp = best_match(["s0", "s1", "s2"], sites, {1: 1, 2: 1, 3: 0})
        self.assertEqual((fp.donor, fp.score), ("s1", 1.0))
        self.assertTrue(Fingerprint("x", 0.99, "y", 0.7, 50).decisive)
        self.assertFalse(Fingerprint("x", 0.69, "y", 0.68, 50).decisive)


@unittest.skipUnless(os.environ.get("RUN_SLOW") and os.path.exists(os.path.join(ROUND_0817, "input.bam")),
                     "set RUN_SLOW=1 (needs round 08:17 and 1000G data, ~1 min)")
class ProviderEndToEndTest(unittest.TestCase):
    def test_1000g_round_germline_is_exact(self):
        r = Round("2026-09-28T08:17:00+00:00", "chr14", 45912001, 53254000,
                  bam=os.path.join(ROUND_0817, "input.bam"), reference="")
        a = Provider().build(r)
        self.assertEqual((a.regime, a.donor), (Regime.KG, "NA18878"))
        truth = set()
        for line in gzip.open(os.path.join(ROUND_0817, "truth.vcf.gz"), "rt"):
            if line.startswith("#") or "SYNTHETIC" in line.split("\t")[7]:
                continue
            p = line.rstrip("\n").split("\t")
            truth.add((int(p[1]), p[3], p[4], canonical_gt(p[9].split(":")[0])))
        self.assertEqual({(v.pos, v.ref, v.alt, v.gt) for v in a.germline}, truth)       # whole records, as written


if __name__ == "__main__":
    unittest.main()


class SyntheticClassifyTest(unittest.TestCase):
    """Parent-based decisions: absent in the parent -> injected (p=1), present -> real donor variant (dropped)."""

    class FakeParent:
        def __init__(self, counts):
            self.c = counts

        def counts(self, pos, ref, alt):
            return self.c.get(pos, (0, 0))

    def test_rules(self):
        from answer_set.synthetic import Candidate, classify
        from common.types import Round, Variant
        r = Round("2026-09-28T08:17:00+00:00", "chr1", 1, 1000, "", "")
        cands = [Candidate(10, "A", "G", 10, 20),      # parent ref only -> injected
                 Candidate(20, "C", "T", 10, 20),      # parent shows alt -> real donor variant
                 Candidate(30, "G", "A", 10, 20),      # parent too thin -> uncertain
                 Candidate(40, "T", "C", 1, 30),       # too weak in the round
                 Candidate(50, "A", "AT", 19, 20),     # germline allele
                 Candidate(60, "C", "A", 19, 20)]      # injected, hom
        parent = self.FakeParent({10: (0, 25), 20: (12, 25), 30: (0, 3), 60: (0, 30)})
        germ = [Variant("chr1", 50, "A", "AT,ATT", "1/2")]
        out, stats = classify(r, cands, parent, germ)
        got = {s.variant.pos: (s.p, s.variant.gt) for s in out}
        self.assertEqual(got, {10: (1.0, "0/1"), 30: (0.5, "0/1"), 60: (1.0, "1/1")})
        self.assertEqual((stats["donor_real"], stats["weak"], stats["germline"]), (1, 1, 1))
        out, _ = classify(r, cands[:1], None, [])                      # no parent at all -> uncertain
        self.assertEqual([s.p for s in out], [0.5])

    def test_germline_keys_match_bcftools_norm(self):
        from answer_set.synthetic import germline_keys
        from common.types import Variant
        self.assertEqual(germline_keys([Variant("chr1", 100, "TCA", "TCACA,T", "1/2")]),
                         {(100, "T", "TCA"), (100, "TCA", "T")})

    def test_canonical_gt(self):
        from answer_set.giab import canonical_gt
        self.assertEqual([canonical_gt(g) for g in ("1|0", "0/1", "2|1", "1/1")], ["0/1", "0/1", "1/2", "1/1"])
