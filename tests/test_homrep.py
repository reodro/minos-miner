"""answer_set.k2.homrep: repeat-aware het/hom for 1000G inserted indels."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pysam

from answer_set.k2 import homrep

# 1-based: position 30 is the anchor before a 14 bp CA repeat (CACACACACACACA)
SEQ = "GATTACAGGCTTACAGATTACCGTAGGCAT" + "CACACACACACACA" + "GGATTACAGCTAGCTAGGCTTACGATCGATCCGATAGCTAGGCATGTTGACA"
ANCHOR = 30


class Always:
    """stand-in model: always hom"""
    def predict_proba(self, X):
        return np.array([[0.0, 1.0]] * len(X))


class HomRepTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        (d / "ref.fa").write_text(">chrT\n" + SEQ + "\n")
        pysam.faidx(str(d / "ref.fa"))
        self.fa = pysam.FastaFile(str(d / "ref.fa"))
        p = d / "r.bam"
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": len(SEQ)}]}
        with pysam.AlignmentFile(str(p), "wb", header=hdr) as out:
            for i in range(6):
                a = pysam.AlignedSegment(out.header)
                a.query_name, a.reference_id, a.reference_start, a.mapping_quality = f"r{i}", 0, 0, 60
                a.query_sequence, a.cigartuples = SEQ, [(0, len(SEQ))]
                a.query_qualities = pysam.qualitystring_to_array("I" * len(SEQ))
                out.write(a)
        pysam.index(str(p))
        self.bam = pysam.AlignmentFile(str(p))

    def tearDown(self):
        self.fa.close(); self.bam.close(); self.tmp.cleanup()

    def test_repeat_at(self):
        self.assertEqual(homrep.repeat_at(self.fa, "chrT", (ANCHOR, "TCA", "T")), (2, 14, 7))

    def test_new_feats_only_for_long_repeat_indels(self):
        k2 = (ANCHOR, "TCA", "T")                         # 2 bp deletion in the CA repeat
        f = homrep.new_feats(self.bam, self.fa, "chrT", k2)
        self.assertEqual(len(f), len(homrep.NEW))
        self.assertEqual(f[:5], [2, 14, 7, 1.0, 1.0])
        self.assertEqual(f[5:7], [0, 6])                  # no read carries it, all 6 show the reference
        self.assertIsNone(homrep.new_feats(self.bam, self.fa, "chrT", (ANCHOR, "TC", "T")))     # 1 bp: k029 keeps it
        self.assertIsNone(homrep.new_feats(self.bam, self.fa, "chrT", (5, "ACA", "A")))         # not in a repeat

    def test_apply_only_called_non_person_candidates(self):
        keys = [(ANCHOR, "TCA", "T")] * 3
        X = np.zeros((3, 2))
        p_hom = np.array([0.1, 0.1, 0.1])
        person = np.array([np.nan, np.nan, 1.0])          # the third is the person's GATK genotype
        p_any = np.array([0.9, 0.1, 0.9])                 # the second is not called
        out, n = homrep.apply(p_hom, keys, X, person, p_any, 0.3, self.bam, self.fa, "chrT", Always())
        self.assertEqual(n, 1)
        self.assertEqual(out.tolist(), [1.0, 0.1, 0.1])


if __name__ == "__main__":
    unittest.main()
