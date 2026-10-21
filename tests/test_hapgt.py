"""answer_set.k2.hapgt: reads judged against local haplotypes (an indel with or without its neighbours)."""
import os
import tempfile
import unittest

import pysam

from answer_set.k2 import hapgt

REF = "ACGTTGCAAGCTTCGATCGGATCCATGCAAGTCCGATTAGCCATGGCATTCAGGTACCTTGAACGTAGCTAGCATCGGATACGTTCAGTCAGTAGGCATCA" * 3


def _read(name, seq, start, cigar):
    r = pysam.AlignedSegment()
    r.query_name, r.reference_id, r.reference_start = name, 0, start
    r.query_sequence, r.cigarstring, r.mapping_quality = seq, cigar, 60
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return r


class HapCountsTest(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        fa = os.path.join(self.d, "ref.fa")
        with open(fa, "w") as f:
            f.write(">chr1\n" + REF + "\n")
        pysam.faidx(fa)
        self.fa = pysam.FastaFile(fa)
        # event: delete REF[100] (1-based 101) after anchor 100; neighbour: 2-bp deletion after anchor 110
        self.k = (100, REF[99:101], REF[99])
        self.nb = (110, REF[109:112], REF[109])
        both = REF[60:100] + REF[101:110] + REF[112:161]                       # the deletion and the neighbour
        alone = REF[60:100] + REF[101:160]                                     # the deletion only
        plain = REF[60:160]
        bam = os.path.join(self.d, "r.bam")
        with pysam.AlignmentFile(bam, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                                   "SQ": [{"SN": "chr1", "LN": len(REF)}]}) as f:
            rs = [_read(f"b{i}", both, 60, "40M1D9M2D49M") for i in range(6)]
            rs += [_read(f"a{i}", alone, 60, "40M1D59M") for i in range(6)]
            for r in sorted(rs, key=lambda r: r.reference_start):
                f.write(r)
        pysam.index(bam)
        self.bam = pysam.AlignmentFile(bam)
        self.plain = plain

    def test_neighbour_makes_every_read_a_carrier(self):
        alone = hapgt.hap_counts(self.bam, self.fa, "chr1", self.k, 0, [])
        self.assertEqual((alone["carry"], alone["other"]), (6, 6))           # judged alone: looks het (6 of 12)
        with_nb = hapgt.hap_counts(self.bam, self.fa, "chr1", self.k, 0, [self.nb])
        self.assertEqual((with_nb["carry"], with_nb["without"], with_nb["other"]), (12, 0, 0))   # really hom

    def test_window_neighbours(self):
        nb, pos = hapgt.neighbours([], [(105, "A", "C"), (300, "A", "C")], {(95, "G", "T"): "het"})
        self.assertEqual([x[1] for x in hapgt.window_neigh(nb, pos, self.k, 0)], ["person", "called"])

    def test_edit_rejects_overlaps(self):
        self.assertIsNone(hapgt._edit("ACGTACGT", 1, [(2, "CG", "C"), (3, "G", "T")]))
        self.assertEqual(hapgt._edit("ACGTACGT", 1, [(2, "C", "T")]), "ATGTACGT")


if __name__ == "__main__":
    unittest.main()
