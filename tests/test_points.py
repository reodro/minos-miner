"""dashboard/points.py: every scored config of a round, projected to 3D for the island view (2026-10-19)."""
import json
import sqlite3
import tempfile
import unittest
import unittest.mock
from pathlib import Path

from dashboard import points
from tests.helpers import ranges
from search.space import Space

MMQ = "min_mapping_quality_score"


class PointsTest(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        space = Space(ranges())
        db = sqlite3.connect(self.dir / "R1.sqlite")
        db.execute("CREATE TABLE results (round_id TEXT, config_key TEXT, config TEXT, fidelity TEXT, "
                   "answer_set_version INTEGER, score REAL, score_sd REAL, gate_passed INTEGER, result TEXT, "
                   "runtime_s REAL, tag TEXT, created REAL)")
        for i, (v, tag) in enumerate([(20, "A"), (30, "B0"), (40, "B101"), (50, "G"), (55, "R"), (58, "D")]):
            c = space._assemble({MMQ: v})
            db.execute("INSERT INTO results VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                       ("R1", c.key, json.dumps(c.to_dict()), "full", 1, 80 + i, None, 1, "{}", 1.0, tag, i))
        db.commit()
        db.close()
        points._cache.clear()
        self.patch = unittest.mock.patch.object(points, "_space", space)
        self.patch.start()

    def tearDown(self):
        self.patch.stop()

    def test_every_scored_config_with_its_source_island_and_climber(self):
        out = points.points("R1", self.dir)
        self.assertEqual(out["n"], 6)
        self.assertEqual(len(out["xyz"]), 18)
        self.assertEqual(len(out["basis"]), 3)
        by = {s: (i, c) for s, i, c in zip(out["score"], out["island"], out["climber"])}
        self.assertEqual(by[81.0], (0, 0))
        self.assertEqual(by[82.0], (1, 101))
        self.assertEqual(sorted(out["sources"][k] for k in out["src"]), ["climber", "climber", "ga", "released", "rf", "step_a"])

    def test_a_missing_round_is_empty(self):
        self.assertEqual(points.points("nope", self.dir)["n"], 0)


if __name__ == "__main__":
    unittest.main()
