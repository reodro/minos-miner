"""answer_set.k2.nest: N without a platform N - from a raised pool, else past draws of the same type that fit the pool."""
import gzip
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pysam

from answer_set.k2 import nest

HEADER = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
BINNED = [2, 3, 4, 5, 6, 10, 20, 30]
FULL = list(range(2, 40))
MB = "chr1:1-1000000"                     # 1 Mb: pool per Mb = pool


def _bam(path, quals):
    with pysam.AlignmentFile(path, "wb", header=HEADER) as f:
        for i in range(40):
            r = pysam.AlignedSegment()
            r.query_name, r.reference_id, r.reference_start = f"r{i}", 0, 100 + i
            r.query_sequence, r.cigarstring, r.mapping_quality = "A" * 50, "50M", 60
            r.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + quals[(i + j) % len(quals)]) for j in range(50)))
            f.write(r)


def _vcf(path, n, info="."):
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i in range(n):
            f.write(f"chr1\t{i + 1}\t.\tA\tC\t60\tPASS\t{info}\n")


def _round(root, tag, n, quals, pool):
    d = os.path.join(root, tag)
    os.makedirs(d, exist_ok=True)
    _bam(os.path.join(d, "input.bam"), quals)
    _vcf(os.path.join(d, "mutations.vcf.gz"), n)
    _vcf(os.path.join(d, "truth.vcf.gz"), pool, "SYNTHETIC")
    Path(d, "round.json").write_text(json.dumps({"region": MB}))


def _tag(i):
    h, m = divmod(i * 72, 60)
    return f"2026101{7 + h // 24}T{h % 24:02d}{m:02d}00Z"


class NestTest(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.cache = mock.patch.object(nest, "ROUND_CACHE", Path(self.root) / "rounds.json")
        self.cache.start()
        # on 2026-10-17/18, every 72 min: 12 binned rounds N 300..410 and 12 full ones N 500..610, base pool 610/Mb;
        # then 6 full rounds with raised pools: N 700..750, pool = 8 N / 0.97
        for i in range(24):
            binned = i % 2 == 0
            _round(self.root, _tag(i), (300 if binned else 500) + 10 * (i // 2), BINNED if binned else FULL, 610)
        for i in range(24, 30):
            n = 700 + 10 * (i - 24)
            _round(self.root, _tag(i), n, FULL, round(8 * n / 0.97))
        self.binned = os.path.join(self.root, _tag(0), "input.bam")
        self.full = os.path.join(self.root, _tag(1), "input.bam")
        self.later = "20261020T000000Z"            # everything revealed

    def tearDown(self):
        self.cache.stop()

    def test_bam_type(self):
        self.assertEqual(nest.bam_type(self.binned), "binned")
        self.assertEqual(nest.bam_type(self.full), "full")
        self.assertIsNone(nest.bam_type(os.path.join(self.root, "missing.bam")))

    def test_raised_pool_gives_n_from_the_pool(self):
        # 8000/Mb is far above the base: N = q x pool / 8, q = 0.97 from the past raised rounds
        n, how = nest.estimate(self.later, self.full, MB, round(8000 * nest.POOL_SCALE), self.root)
        self.assertEqual(how["rule"], "raised pool")
        self.assertAlmostEqual(n, 0.97 * 1000, delta=2)

    def test_base_pool_takes_the_draws_that_fit(self):
        # base pool 610 -> cap 610 / 8 / 0.95 = 80: none of the past draws fit -> type median, capped
        n, how = nest.estimate(self.later, self.binned, MB, round(610 * nest.POOL_SCALE), self.root)
        self.assertEqual((n, how["rule"]), (80, "type median, capped"))
        # 4800 on 8 Mb (600/Mb, base): cap = 4800 / 0.9976 / 8 / 0.95 = 633
        n, how = nest.estimate(self.later, self.full, "chr1:1-8000000", 4800, self.root)
        self.assertEqual((how["rule"], how["cap"]), ("draws that fit the pool", 633))
        self.assertEqual(n, 555)                                       # full draws that fit: 500..610, not 700..750

    def test_type_is_kept_apart(self):
        n, _ = nest.estimate(self.later, self.binned, "chr1:1-8000000", 4800, self.root)
        self.assertEqual(n, 355)                                       # binned draws 300..410

    def test_a_short_pool_count_sets_no_cap(self):
        n, how = nest.estimate(self.later, self.full, MB, 300, self.root)   # 300/Mb: the build lost flags
        self.assertEqual((how["per_mb"], how["rule"], how["cap"]), (None, "draws that fit the pool", None))
        self.assertEqual(n, 585)                                       # all 18 full draws: 500..610 and 700..750

    def test_unrevealed_rounds_are_not_used(self):
        n, _ = nest.estimate("20261017T120000Z", self.binned, MB, 300, self.root)   # 4 binned rounds out by then
        self.assertIsNone(n)


class ProviderTest(unittest.TestCase):
    """read_output uses nest's N unless the platform gave one; an unrecognised round gives no pool (no cap)."""

    def _out(self, d, mode="GIAB:HG004", report_extra=None):
        import json
        for name, lines in (("truth.vcf.gz", ["chr21\t200\t.\tCA\tC\t60\tPASS\tSYNTHETIC\tGT\t1/1"]),
                            ("mutations.vcf.gz", ["chr21\t1\t.\tA\tG\t60\tPASS\tSYNTHETIC\tGT\t0/1"])):
            with gzip.open(os.path.join(d, name), "wt") as f:
                f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
                f.write("\n".join(lines) + "\n")
        with gzip.open(os.path.join(d, "scores.tsv.gz"), "wt") as f:
            f.write("pos\tref\talt\tp\n200\tCA\tC\t0.9\n")
        Path(d, "report.json").write_text(json.dumps({"mode": mode, **(report_extra or {})}))

    def test_estimate_used_unless_the_platform_gave_n(self):
        from answer_set.k2.provider import read_output
        from common.types import Round
        with tempfile.TemporaryDirectory() as d:
            self._out(d)
            seen = []
            fake = lambda tag, bam, region, pool, rounds_dir=None: (seen.append(pool) or (512, {"type": "full"}))
            with mock.patch.object(nest, "estimate", side_effect=fake):
                r = Round("2026-10-19T09:30:00+00:00", "chr21", 1, 1000, bam="/x.bam", reference="/r.fa",
                          n_targets=480, n_source="recent rounds")
                a = read_output(r, Path(d))
                self.assertEqual((a.n_targets, a.notes["n_targets_source"], seen), (512, "type + pool cap", [1]))
                r = Round(r.round_id, "chr21", 1, 1000, bam="/x.bam", reference="/r.fa", n_targets=611, n_source="platform")
                self.assertEqual(read_output(r, Path(d)).n_targets, 611)
                self.assertEqual(len(seen), 1)                                   # not even asked
                self._out(d, mode="unrecognised (reads only)")                   # no inserted flags: pool 0 -> no cap
                read_output(Round(r.round_id, "chr21", 1, 1000, bam="/x.bam", reference="/r.fa"), Path(d))
                self.assertEqual(seen[-1], 0)
            with mock.patch.object(nest, "estimate", return_value=(None, {})):   # too few rounds: caller's N stays
                r = Round(r.round_id, "chr21", 1, 1000, bam="/x.bam", reference="/r.fa", n_targets=480, n_source="recent rounds")
                self.assertEqual(read_output(r, Path(d)).n_targets, 480)


if __name__ == "__main__":
    unittest.main()
