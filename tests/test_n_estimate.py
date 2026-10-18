"""main.n_estimate: N for a round without a published N is the median of recent revealed rounds' target counts."""
import gzip
import os
import tempfile
import unittest
import unittest.mock

from common.tuning import TUNING
from main.n_estimate import recent_n


def _round(root, tag, n):
    d = os.path.join(root, tag)
    os.makedirs(d, exist_ok=True)
    with gzip.open(os.path.join(d, "mutations.vcf.gz"), "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")
        for i in range(n):
            f.write(f"chr1\t{i + 1}\t.\tA\tC\n")


class RecentNTest(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()
        # rounds every 72 minutes; counts 300, 310, ... so the median is easy to read
        for i, tag in enumerate(["20261018T010000Z", "20261018T021200Z", "20261018T032400Z", "20261018T043600Z",
                                 "20261018T054800Z", "20261018T070000Z", "20261018T081200Z"]):
            _round(self.root, tag, 300 + 10 * i)

    def test_median_of_rounds_revealed_before_the_round_started(self):
        # 20261018T090000Z: rounds starting at or before 05:30 are revealed (3.5 h) -> 01:00 .. 04:36 (300..330)
        self.assertEqual(recent_n("20261018T090000Z", self.root, k=10), 315)

    def test_only_the_last_k(self):
        self.assertEqual(recent_n("20261018T120000Z", self.root, k=3), 350)   # 05:48, 07:00, 08:12 -> 340, 350, 360

    def test_too_few_revealed_rounds_give_no_estimate(self):
        self.assertIsNone(recent_n("20261018T050000Z", self.root, k=10))     # only 01:00 is out by then

    def test_zero_turns_it_off(self):
        with unittest.mock.patch.dict(TUNING.values, {"main.n_recent_rounds": 0}):
            self.assertIsNone(recent_n("20261018T120000Z", self.root))

    def test_the_setting_sets_k(self):
        with unittest.mock.patch.dict(TUNING.values, {"main.n_recent_rounds": 3}):
            self.assertEqual(recent_n("20261018T120000Z", self.root), 350)


if __name__ == "__main__":
    unittest.main()
