"""answer_set.k2.calldna: which candidates the DNA call model re-scores."""
import unittest
from unittest import mock

import numpy as np

from answer_set.k2 import calldna


class _Const:
    """stand-in model: probability 0.77 for every row, and remembers the rows it saw"""
    def predict_proba(self, X):
        self.rows = X
        return np.column_stack([np.full(len(X), 0.23), np.full(len(X), 0.77)])


class CallDnaTest(unittest.TestCase):
    def test_only_non_person_indels_are_rescored(self):
        keys = [(100, "A", "AT"), (200, "C", "G"), (300, "GT", "G"), (400, "T", "TA")]
        p = np.array([0.1, 0.2, 0.3, 0.4])
        X = np.arange(8, dtype=float).reshape(4, 2)
        person = np.array([0.0, 0.0, np.nan, 1.0])          # the last indel is one of the person's GATK calls
        m = _Const()
        with mock.patch.object(calldna, "feats", lambda bam, fa, chrom, k: [1.0, 2.0]):
            out, n = calldna.apply(p, keys, X, person, None, None, "chr1", m)
        self.assertEqual(n, 2)                              # the two non-person indels
        np.testing.assert_allclose(out, [0.77, 0.2, 0.77, 0.4])
        np.testing.assert_allclose(m.rows, [[0, 1, 1, 2], [4, 5, 1, 2]])   # own columns, then the DNA measurements
        np.testing.assert_allclose(p, [0.1, 0.2, 0.3, 0.4])                # input untouched

    def test_no_indel_candidates(self):
        out, n = calldna.apply(np.array([0.5]), [(1, "A", "C")], np.zeros((1, 2)), np.array([0.0]), None, None, "chr1",
                               _Const())
        self.assertEqual(n, 0)
        np.testing.assert_allclose(out, [0.5])


if __name__ == "__main__":
    unittest.main()
