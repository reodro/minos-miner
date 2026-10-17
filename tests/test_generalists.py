"""Cross-round generalists as wave-0 seeds (search/portfolio.generalists)."""
import json
import os
import tempfile
import unittest
from pathlib import Path

from search import portfolio
from search.space import Space
from tests.helpers import ranges


class GeneralistsTest(unittest.TestCase):
    def test_best_first_snapped_and_capped(self):
        space = Space(ranges())
        rows = [{"params": {"min_mapping_quality_score": 40, "assembly_region_padding": 50}, "fitness": 0.2},
                {"params": {"min_mapping_quality_score": 41, "assembly_region_padding": 50}, "fitness": 0.1},
                {"params": {"min_mapping_quality_score": 40, "assembly_region_padding": 50}, "fitness": 0.05},
                {"params": {"min_mapping_quality_score": 42}, "fitness": 0.0}]
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "generalists.json"
            p.write_text(json.dumps(rows))
            got = portfolio.generalists(space, 2, p)
            self.assertEqual(len(got), 2)
            self.assertEqual(got[0].params["min_mapping_quality_score"], 40)
            self.assertEqual(got[1].params["min_mapping_quality_score"], 41)       # the duplicate is skipped
            self.assertEqual(portfolio.generalists(space, 0, p), [])

    def test_no_file_no_seeds(self):
        self.assertEqual(portfolio.generalists(Space(ranges()), 10, Path("/nonexistent/g.json")), [])


if __name__ == "__main__":
    unittest.main()
