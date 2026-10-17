"""answer_set.k2.readguard: later steps' inserted indels the reads do not carry are dropped."""
import tempfile
import unittest
from pathlib import Path

import pysam

from answer_set.k2 import readguard

SEQ = "GATTACAGGCTTACAGATTACCGTAGGCATCGATCGGATTACAGCTAGCTAGGCTTACGATCGATCCGATAGCTAGGCATG"   # 82 bp
POS = 41                      # 1-based anchor of the test insertion A>AT (SEQ[40] == "A")


class ReadGuardTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        (d / "ref.fa").write_text(">chrT\n" + SEQ + "\n")
        pysam.faidx(str(d / "ref.fa"))
        self.fa = pysam.FastaFile(str(d / "ref.fa"))
        self.d = d

    def tearDown(self):
        self.fa.close()
        self.tmp.cleanup()

    def bam(self, n_ref, n_ins):
        """reads over the whole reference: n_ref plain, n_ins with a T inserted after POS"""
        p = self.d / f"r{n_ref}_{n_ins}.bam"
        hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": len(SEQ)}]}
        with pysam.AlignmentFile(str(p), "wb", header=hdr) as out:
            for i in range(n_ref + n_ins):
                a = pysam.AlignedSegment(out.header)
                a.query_name, a.reference_id, a.reference_start, a.mapping_quality = f"r{i}", 0, 0, 60
                if i < n_ref:
                    a.query_sequence, a.cigartuples = SEQ, [(0, len(SEQ))]
                else:
                    a.query_sequence = SEQ[:POS] + "T" + SEQ[POS:]
                    a.cigartuples = [(0, POS), (1, 1), (0, len(SEQ) - POS)]
                a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
                out.write(a)
        pysam.index(str(p))
        return pysam.AlignmentFile(str(p))

    def test_counts(self):
        key = (POS, SEQ[POS - 1], SEQ[POS - 1] + "T")
        self.assertEqual(readguard.counts(self.bam(7, 3), self.fa, "chrT", key, {key: "het"}), (3, 7))

    def test_drops_unsupported_later_record_only(self):
        key = (POS, SEQ[POS - 1], SEQ[POS - 1] + "T")
        rec = key + ("het", True, 1.0)
        out, n = readguard.apply([rec], self.bam(10, 0), self.fa, "chrT", {key: 0.1}, 0.5)
        self.assertEqual((out, n), ([], 1))                       # a later step's record no read carries
        out, n = readguard.apply([rec], self.bam(10, 0), self.fa, "chrT", {key: 0.9}, 0.5)
        self.assertEqual((out, n), ([rec], 0))                    # the model's own call is never touched
        out, n = readguard.apply([rec], self.bam(5, 5), self.fa, "chrT", {}, 0.5)
        self.assertEqual((out, n), ([rec], 0))                    # carried by half the reads: kept
        person = key + ("het", False, 2.0)
        self.assertEqual(readguard.apply([person], self.bam(10, 0), self.fa, "chrT", {}, 0.5), ([person], 0))


if __name__ == "__main__":
    unittest.main()
