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


class K2OutputTest(unittest.TestCase):
    """K2Provider reads the builder's files as written; scoring is pointed at those exact files."""

    def _write(self, d, name, lines):
        with gzip.open(os.path.join(d, name), "wt") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
            f.writelines(line + "\n" for line in lines)

    def test_read_output(self):
        import json
        import tempfile
        from pathlib import Path
        from answer_set.k2.provider import read_output
        with tempfile.TemporaryDirectory() as d:
            self._write(d, "truth.vcf.gz", ["chr21\t100\t.\tA\tG\t60\tPASS\t.\tGT\t0/1",
                                            "chr21\t200\t.\tCA\tC\t60\tPASS\tSYNTHETIC\tGT\t1/1",
                                            "chr21\t300\t.\tTAA\tT,TA\t60\tPASS\tSYNTHETIC\tGT\t1/2"])
            self._write(d, "mutations.vcf.gz", [f"chr21\t{p}\t.\tA\tG\t60\tPASS\tSYNTHETIC\tGT\t0/1" for p in (1, 2, 3)])
            with gzip.open(os.path.join(d, "scores.tsv.gz"), "wt") as f:
                f.write("pos\tref\talt\tp\n200\tCA\tC\t0.9100\n300\tTAA\tT\t0.8000\n300\tTA\tT\t0.7000\n")
            Path(d, "report.json").write_text(json.dumps({"mode": "GIAB:HG004", "model": "giab track g006", "threshold": 0.5}))
            r = Round("2026-09-30T03:29:00+00:00", "chr21", 1, 1000, bam="/x.bam", reference="/ref.fa")
            a = read_output(r, Path(d))
            self.assertEqual((a.regime, a.donor, a.n_targets), (Regime.GIAB, "HG004", 3))
            self.assertEqual(a.germline, [Variant("chr21", 100, "A", "G", "0/1")])
            self.assertEqual([(c.variant.pos, c.p) for c in a.synthetic], [(200, 0.91), (300, 0.7)])   # 1/2: weaker allele
            self.assertEqual(a.notes["truth_vcf"], os.path.join(d, "truth.vcf.gz"))
            self.assertEqual((a.notes["region"], a.notes["reference"]), ("chr21:1-1000", "/ref.fa"))
            r2 = Round(r.round_id, "chr21", 1, 1000, bam="/x.bam", reference="/ref.fa", n_targets=611)
            self.assertEqual(read_output(r2, Path(d)).n_targets, 611)                                # platform N wins

    def test_regime_of(self):
        from answer_set.k2.provider import regime_of
        self.assertEqual(regime_of("1000G:HG03559"), (Regime.KG, "HG03559"))
        self.assertEqual(regime_of("1000G:NA18878 WITHOUT DRAGEN (download failed, reads only)"), (Regime.KG, "NA18878"))
        self.assertEqual(regime_of("unrecognised (reads only)"), (Regime.UNKNOWN, None))


class BuilderFactoryTest(unittest.TestCase):
    def test_names(self):
        from answer_set import BUILDER_NAMES, K2Provider, builder
        self.assertIsInstance(builder(), type(builder(BUILDER_NAMES[0])))           # None -> the default
        self.assertIsInstance(builder("provider"), Provider)
        self.assertIsInstance(builder("k2"), K2Provider)
        with self.assertRaises(ValueError):
            builder("nope")

    def test_env_override(self):
        from unittest import mock
        from answer_set import K2Provider, builder
        from answer_set.builders import SearchableOracle
        with mock.patch.dict(os.environ, {"ANSWER_SET_BUILDER": "oracle"}):
            self.assertIsInstance(builder(), SearchableOracle)
            self.assertIsInstance(builder("k2"), K2Provider)                          # an explicit name still wins
        with mock.patch.dict(os.environ, {"ANSWER_SET_BUILDER": ""}):
            self.assertIsInstance(builder(), K2Provider)                              # empty: the default

    @unittest.skipUnless(os.path.exists(os.path.join(SETTINGS.rounds_dir, "20261001T164100Z", "truth.vcf.gz")),
                         "needs released round 2026-10-01 16:41")
    def test_oracle_is_searchable(self):
        from answer_set.builders import SearchableOracle
        from backtest.data import load
        a = SearchableOracle().build(load("20261001T164100Z"))
        self.assertEqual(a.regime, Regime.GIAB)                                       # main searches only if not UNKNOWN
        self.assertTrue(a.germline and a.synthetic)


class K2RobustnessTest(unittest.TestCase):
    """A weaker build must never be silent: retries, input manifest, warnings."""

    def test_retry(self):
        from answer_set.k2.util import retry
        calls, waits = [], []

        def flaky():
            calls.append(1)
            if len(calls) < 3:
                raise OSError("connection reset")
            return "ok"
        self.assertEqual(retry(flaky, 3, 10, sleep=waits.append), ("ok", ["OSError: connection reset"] * 2))
        self.assertEqual(waits, [10, 20])
        with self.assertRaises(RuntimeError) as e:
            retry(lambda: 1 / 0, 2, 0, sleep=lambda s: None)
        self.assertIn("failed after 2 attempts: ZeroDivisionError", str(e.exception))

    def test_manifest(self):
        import tempfile
        from pathlib import Path
        from answer_set.k2.inputs import check, record
        with tempfile.TemporaryDirectory() as d:
            root, man = Path(d), Path(d, "inputs.json")
            a, b = Path(d, "a.vcf.gz"), Path(d, "b.bed")
            a.write_text("AAAA"); b.write_text("BB")
            self.assertEqual(check([a], man, root), [f"no input manifest {man} (python -m answer_set.k2.inputs record)"])
            self.assertEqual(record([a, b], man, root), 2)
            self.assertEqual(check([a, b], man, root, full=True), [])
            a.write_text("AAAAA")                                                     # other size
            os.utime(b, ns=(1, 1))                                                    # touched
            self.assertEqual(check([a, b, Path(d, "c")], man, root),
                             ["changed: a.vcf.gz", "modified since recorded: b.bed", "missing: c"])
            b.write_text("XX"); os.utime(b, ns=(1, 1))                                 # same size, new content
            self.assertEqual(check([b], man, root, full=True), ["changed: b.bed"])

    def test_warnings(self):
        from answer_set.k2.provider import warnings_of
        self.assertEqual(warnings_of({"inputs": "ok", "readdiff": {"edited": 5}}), [])
        w = warnings_of({"inputs": ["changed: giab/HG002.vcf.gz"], "readdiff": {"error": "failed after 3 attempts"}})
        self.assertEqual(len(w), 2)
        self.assertIn("HG002", w[0]); self.assertIn("read comparison failed", w[1])
        self.assertIn("needed 2 tries", warnings_of({"readdiff": {"failed_tries": ["OSError: x"]}})[0])


class K2AssembleTest(unittest.TestCase):
    """Model scores -> truth records: donor background copied, threshold, <= 2 alleles per position, BED rule."""

    def test_giab_round(self):
        from answer_set.k2.assemble import assemble_records
        donor = {(100, "A", "G"): "het", (150, "C", "CT"): "hom"}
        keys = [(100, "A", "G"), (200, "T", "C"), (300, "G", "A"), (300, "G", "GA"), (300, "G", "T"),
                (400, "A", "C"), (500, "C", "G")]
        p_any = [0.9, 0.95, 0.9, 0.8, 0.7, 0.3, 0.99]
        p_hom = [0.0, 0.9, 0.9, 0.1, 0.1, 0.0, 0.0]
        recs = assemble_records(keys, p_any, p_hom, 0.5, True, donor, allowed=lambda k: k[0] < 500)
        got = {r[:5] for r in recs}
        self.assertEqual(got, {(100, "A", "G", "het", False), (150, "C", "CT", "hom", False),    # donor, as published
                               (200, "T", "C", "hom", True),                                     # model call, P(hom) 0.9
                               (300, "G", "A", "het", True), (300, "G", "GA", "het", True)})     # top 2 of 3 -> 1/2
        # 400: below threshold; 500: outside the donor's confident region

    def test_kg_round(self):
        from answer_set.k2.assemble import assemble_records
        recs = assemble_records([(100, "A", "G")], [0.9], [0.2], 0.5, False, {})
        self.assertEqual([r[:5] for r in recs], [(100, "A", "G", "het", False)])        # 1000G: nothing SYNTHETIC


ROUND_1641 = os.path.join(SETTINGS.rounds_dir, "20261001T164100Z")


@unittest.skipUnless(os.environ.get("RUN_SLOW") and os.path.exists(os.path.join(ROUND_1641, "input.bam")),
                     "set RUN_SLOW=1 (needs round 2026-10-01 16:41, the answer_set/k2 data, Docker, internet; ~1 min)")
class K2EndToEndTest(unittest.TestCase):
    def test_giab_round_matches_released_truth(self):
        import tempfile
        from answer_set import K2Provider
        from answer_set.k2.util import read_alleles, ref_path
        from backtest.data import load
        r = load("20261001T164100Z", reference=str(ref_path("chr18")))
        with tempfile.TemporaryDirectory() as d:
            a = K2Provider(work_dir=d).build(r)
            ours = {k for k, _, _ in read_alleles(a.notes["truth_vcf"])}
            real = {k: "SYNTHETIC" in i for k, _, i in read_alleles(os.path.join(r.released_dir, "truth.vcf.gz"),
                                                                     keep_info=True)}
        self.assertEqual((a.regime, a.donor), (Regime.GIAB, "HG002"))
        self.assertNotIn("error", a.notes["report"]["readdiff"])                   # Google originals reached
        self.assertEqual(a.notes["warnings"], [])                                  # inputs as recorded, nothing degraded
        self.assertEqual({k for k, syn in real.items() if not syn} - ours, set())  # donor background copied exactly
        m = len(ours & real.keys())
        p, rc = m / len(ours), m / len(real)
        self.assertGreaterEqual(2 * p * rc / (p + rc), 0.998)                      # measured 0.9991


if __name__ == "__main__":
    unittest.main()
