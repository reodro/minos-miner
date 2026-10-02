import gzip
import os
import unittest

from common.settings import SETTINGS
from common.types import Regime, Round, Variant


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
            retry(lambda: (_ for _ in ()).throw(OSError("down")), 2, 0, sleep=lambda s: None)
        self.assertIn("failed after 2 attempts: OSError: down", str(e.exception))

    def test_retry_bug_is_not_retried(self):
        from answer_set.k2.util import retry
        calls = []

        def bug():
            calls.append(1)
            return None < 1                                   # the TypeError readdiff hit on an unmapped read
        with self.assertRaises(RuntimeError) as e:
            retry(bug, 3, 10, sleep=lambda s: self.fail("waited before a retry of a bug"))
        self.assertEqual(len(calls), 1)
        self.assertIn("a bug, not retried", str(e.exception))

    def test_retry_budget(self):
        from answer_set.k2.util import retry
        now, waits = [0.0], []

        def slow():
            now[0] += 100                                     # every try takes 100 s and fails
            raise OSError("timeout")
        with self.assertRaises(RuntimeError) as e:
            retry(slow, 5, 10, sleep=lambda s: (waits.append(s), now.__setitem__(0, now[0] + s)), budget_s=250,
                  clock=lambda: now[0])
        self.assertEqual(waits, [10, 20])                     # tries end at 100, 210, 330; a 30 s wait: past 250
        self.assertIn("out of time after 3 attempts", str(e.exception))

    def test_bcftools_timeout_removes_container(self):
        import subprocess
        from unittest import mock
        from answer_set.k2 import util
        calls = []

        def run(args, **kw):
            calls.append(args)
            if args[:2] == ["docker", "run"]:
                raise subprocess.TimeoutExpired(args, kw["timeout"])
            return subprocess.CompletedProcess(args, 0)
        with mock.patch.object(util.subprocess, "run", run):
            with self.assertRaises(RuntimeError) as e:
                util.bcftools("bcftools view x", [], timeout_s=5)
        self.assertIn("timed out after 5 s", str(e.exception))
        name = calls[0][calls[0].index("--name") + 1]
        self.assertEqual(calls[1], ["docker", "rm", "-f", name])

    def test_original_reads_fetch_is_bounded(self):
        """A stuck download fails the fetch within its budget instead of blocking the round."""
        import subprocess
        import tempfile
        from pathlib import Path
        from unittest import mock
        from answer_set.k2 import giab_origin
        with tempfile.TemporaryDirectory() as d:
            Path(d, "HG001.novaseq.pcr-free.50x.dedup.grch38.bam.bai").write_bytes(b"x")
            now = [0.0]

            def run(args, **kw):
                now[0] += kw["timeout"]                       # every piece hangs until killed
                raise subprocess.TimeoutExpired(args, kw["timeout"])
            with mock.patch.object(giab_origin, "CACHE", Path(d)), \
                    mock.patch.object(giab_origin.subprocess, "run", run), \
                    mock.patch.object(giab_origin.time, "monotonic", lambda: now[0]), \
                    mock.patch.object(giab_origin.time, "sleep", lambda s: now.__setitem__(0, now[0] + s)), \
                    mock.patch.object(giab_origin.cf, "ThreadPoolExecutor",
                                      lambda n: __import__("concurrent.futures").futures.ThreadPoolExecutor(1)):
                with self.assertRaises(RuntimeError):
                    giab_origin.fetch("HG001", "chr1", 1000, 2000, chunks=2, budget_s=400)
            self.assertLessEqual(now[0], 400 + giab_origin.PIECE_TIMEOUT_S)
            self.assertEqual([x.name for x in Path(d).iterdir()], ["HG001.novaseq.pcr-free.50x.dedup.grch38.bam.bai"])

    def test_readdiff_skips_unmapped_reads(self):
        """An unmapped round read (placed at its mate, no reference_end) that is an exact piece of an original read
        broke the whole read comparison (TypeError sorting read ends) in 1 of 28 rounds."""
        import random
        import tempfile
        import pysam
        from answer_set.k2 import readdiff
        rnd = random.Random(7)
        seq = "".join(rnd.choice("ACGT") for _ in range(151))
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}

        def bam(path, reads):
            with pysam.AlignmentFile(path, "wb", header=header) as o:
                for name, s, start, unmapped in reads:
                    r = pysam.AlignedSegment()
                    r.query_name, r.query_sequence, r.reference_id, r.reference_start = name, s, 0, start
                    r.query_qualities = pysam.qualitystring_to_array("I" * len(s))
                    if unmapped:
                        r.flag = 4
                    else:
                        r.cigarstring, r.mapping_quality = f"{len(s)}M", 60
                    o.write(r)
            pysam.index(path)
        with tempfile.TemporaryDirectory() as d:
            bam(f"{d}/orig.bam", [("o1", seq, 5000, False)])
            bam(f"{d}/round.bam", [("r1", seq[10:138], 5010, False), ("r2", seq[20:148], 5020, True)])
            found, stats = readdiff.recover(f"{d}/round.bam", f"{d}/orig.bam", "chr1", 4000, 6000, workers=1)
        self.assertEqual(found, {})
        self.assertEqual(stats, {"unedited": 1})

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
