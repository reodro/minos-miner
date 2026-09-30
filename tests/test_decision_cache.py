import gzip
import os
import shutil
import tempfile
import unittest

from common.settings import SETTINGS
from common.types import FULL, AnswerSet, CallSet, Regime
from evaluator.decision_cache import (DecisionCacheScorer, Record, in_spans, lines_in_window, make_clusters,
                                      read_records, window_pad)

SLOW_TAG = "20260928T081700Z"
BENCH = "/minos/minos-miner/test/out/bench"


def rec(pos, ref="A", alt="G", gt="0/1", flt="PASS"):
    return Record(pos, pos + len(ref) - 1, ref, alt, gt, flt)


class ClusterTest(unittest.TestCase):
    def test_gap_below_threshold_merges_at_threshold_splits(self):
        cl = make_clusters([(100, 100), (1099, 1099), (2099, 2099)], [], gap=1000)
        self.assertEqual([(c.start, c.end) for c in cl], [(100, 1099), (2099, 2099)])

    def test_span_end_uses_ref_length(self):
        # a 50 bp deletion at 100 ends at 149; 1148 is 999 bp after it -> same cluster
        cl = make_clusters([(1148, 1148)], [rec(100, ref="A" * 50, alt="A")], gap=1000)
        self.assertEqual([(c.start, c.end) for c in cl], [(100, 1148)])

    def test_query_records_assigned_to_their_cluster(self):
        q = [rec(10), rec(20, gt="1/1"), rec(5000, alt="T")]
        cl = make_clusters([(15, 15), (9000, 9000)], q, gap=1000)
        self.assertEqual(len(cl), 3)
        self.assertEqual(cl[0].query, ((10, "A", "G", "0/1", "PASS"), (20, "A", "G", "1/1", "PASS")))
        self.assertEqual(cl[1].query, ((5000, "A", "T", "0/1", "PASS"),))
        self.assertEqual(cl[2].query, ())          # truth-only cluster

    def test_key_depends_on_query_content(self):
        truth = [(100, 100)]
        k = lambda q: make_clusters(truth, q)[0].key  # noqa: E731
        self.assertEqual(k([rec(120)]), k([rec(120)]))
        self.assertNotEqual(k([rec(120)]), k([rec(120, gt="1/1")]))
        self.assertNotEqual(k([rec(120)]), k([rec(120, flt="LowQual")]))
        self.assertNotEqual(k([rec(120)]), k([rec(120, alt="T")]))
        self.assertNotEqual(k([rec(120)]), k([]))

    def test_untouched_clusters_keep_their_key(self):
        truth = [(100, 100), (10_000, 10_000)]
        a = make_clusters(truth, [rec(110)])
        b = make_clusters(truth, [rec(110), rec(10_050, alt="C")])
        self.assertEqual(a[0].key, b[0].key)
        self.assertNotEqual(a[1].key, b[1].key)

    def test_windows_of_distinct_clusters_do_not_overlap(self):
        for gap in (300, 1000, 1001):
            pad = window_pad(gap)
            cl = make_clusters([(100, 100), (100 + gap, 100 + gap)], [], gap=gap)
            self.assertEqual(len(cl), 2)
            self.assertLess(cl[0].end + pad, cl[1].start - pad)

    def test_in_spans_and_lines_in_window(self):
        spans = [(10, 20), (50, 60)]
        self.assertEqual([p for p in (9, 10, 20, 21, 50, 61) if in_spans(p, spans)], [10, 20, 50])
        lines = [(5, "a"), (10, "b"), (10, "c"), (15, "d"), (30, "e")]
        self.assertEqual(lines_in_window(lines, 10, 15), ("b", "c", "d"))
        self.assertEqual(lines_in_window(lines, 16, 29), ())


class ReadRecordsTest(unittest.TestCase):
    def test_region_filter_sort_and_phasing(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d)
        path = os.path.join(d, "q.vcf.gz")
        with gzip.open(path, "wt") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
            f.write("chr1\t200\t.\tAC\tA\t50\tPASS\t.\tGT:DP\t0|1:10\n")
            f.write("chr1\t100\t.\tA\tG\t50\tLowQual\t.\tGT\t1/1\n")
            f.write("chr1\t5000\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\n")
            f.write("chr2\t150\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\n")
        self.assertEqual(read_records(path, "chr1", 100, 1000),
                         [Record(100, 100, "A", "G", "1/1", "LowQual"), Record(200, 201, "AC", "A", "0/1", "PASS")])


class ScorerContractTest(unittest.TestCase):
    def test_only_full_fidelity(self):
        with self.assertRaises(NotImplementedError):
            DecisionCacheScorer(work_dir=tempfile.mkdtemp()).score(
                CallSet("r", "k", "/x.vcf.gz", "", 0, fidelity="chr1:1-100"), AnswerSet("r", 1, Regime.GIAB, None, [], [], 1))

    def test_needs_truth_files(self):
        with self.assertRaises(NotImplementedError):
            DecisionCacheScorer(work_dir=tempfile.mkdtemp()).score(
                CallSet("r", "k", "/x.vcf.gz", "", 0, fidelity=FULL), AnswerSet("r", 1, Regime.GIAB, None, [], [], 1))

    def test_one_cache_per_answer_set_version(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d)
        truth = os.path.join(d, "t.vcf.gz")
        with gzip.open(truth, "wt") as f:
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n")
        notes = {"truth_vcf": truth, "mutations_vcf": truth, "region": "chr1:1-100", "reference": "/r/chr1.fa"}
        s = DecisionCacheScorer(work_dir=d)
        v1 = s._cache_for(AnswerSet("r", 1, Regime.GIAB, None, [], [], 1, notes=notes))
        self.assertIs(v1, s._cache_for(AnswerSet("r", 1, Regime.GIAB, None, [], [], 1, notes=notes)))
        self.assertIsNot(v1, s._cache_for(AnswerSet("r", 2, Regime.GIAB, None, [], [], 1, notes=notes)))
        self.assertEqual(v1.truth_spans, [(5, 5)])
        self.assertEqual(v1.sdf, "/r/chr1.sdf")


@unittest.skipUnless(os.environ.get("RUN_SLOW") and os.path.exists(os.path.join(SETTINGS.rounds_dir, SLOW_TAG))
                     and os.path.isdir(BENCH),
                     f"set RUN_SLOW=1 (needs round {SLOW_TAG}, bench callsets, docker hap.py; ~2 min)")
class DecisionCacheExactnessTest(unittest.TestCase):
    def test_equals_happy_scorer(self):
        from answer_set import OracleProvider
        from backtest.data import load
        from evaluator import HappyScorer
        r = load(SLOW_TAG, reference=os.path.join(SETTINGS.data_dir, "reference/chr14/chr14.fa"))
        answer = OracleProvider().build(r)
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d)
        names = [f"baseline/{SLOW_TAG}/{n}" for n in ("stock", "a", "b")] + \
                [f"mutants/{SLOW_TAG}/rand{i:02d}" for i in (0, 1)]
        scorer = DecisionCacheScorer(work_dir=os.path.join(d, "dcache"))
        for n in names + names[:1]:      # the repeat must be served fully from the cache
            src = os.path.join(BENCH, n + ".vcf.gz")
            # HappyScorer writes hap.py's files next to the query: give it a private copy
            copy = os.path.join(d, n.replace("/", "_") + ".vcf.gz")
            shutil.copy(src, copy)
            shutil.copy(src + ".tbi", copy + ".tbi")
            with self.subTest(callset=n):
                self.assertEqual(scorer.score(CallSet(r.round_id, n, src, "", 0), answer),
                                 HappyScorer().score(CallSet(r.round_id, n, copy, "", 0), answer))
        st = scorer.stats
        self.assertEqual(st["callsets"], len(names) + 1)
        self.assertGreaterEqual(st["cached_callsets"], 1)
        self.assertLessEqual(st["happy_runs"], len(names))


if __name__ == "__main__":
    unittest.main()
