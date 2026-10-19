"""The dashboard's climber parts and island distances (search/climb.py part_of / row, search/islands.py), 2026-10-19."""
import unittest

from search.climb import Climber
from search.islands import island_distances
from search.space import Space
from tests.helpers import ranges
from tests.test_islands import JumpCtx

MMQ = "min_mapping_quality_score"


class PartsTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.ctx = JumpCtx(self.space, lambda c: 80.0)
        self.cl = Climber(self.ctx, 0, self.space._assemble({MMQ: 30}), 80.0)

    def at(self, v):
        return self.space._assemble({MMQ: v})

    def test_each_move_is_credited_to_the_part_that_found_it(self):
        self.ctx.ga_keys, self.ctx.rf_keys = {self.at(33).key, self.at(34).key}, {self.at(34).key}
        self.cl.move(self.at(31), 80.1, "read grid")
        self.cl.move(self.at(32), 80.3, "jump")
        self.cl.move(self.at(33), 80.4, "ga")
        self.cl.move(self.at(34), 80.6, "ga")             # an RF-UCB config counts as a GA child, credited to RF
        P = self.cl.row()["parts"]
        self.assertEqual([P[k]["moves"] for k in ("grid", "jump", "ga", "rf")], [1, 1, 1, 1])
        self.assertAlmostEqual(P["jump"]["gain"], 0.2)
        self.assertAlmostEqual(P["rf"]["gain"], 0.2)
        self.assertEqual([t[2] for t in self.cl.trail], ["grid", "jump", "ga", "rf"])
        self.assertEqual([e[3] for e in self.cl.history[-4:]], ["grid", "jump", "ga", "rf"])   # the timeline colours

    def test_the_row_carries_the_centre_encoding_and_runs_in_flight(self):
        row = self.cl.row()
        self.assertEqual(len(row["vec"]), len(self.space.names))
        self.assertEqual(row["parts"]["grid"], {"moves": 0, "gain": 0.0, "last_s": None, "sent": 0, "inflight": 0})

    def test_a_lower_move_adds_no_gain(self):
        self.cl.move(self.at(40), 79.0, "kick")
        self.assertEqual(self.cl.row()["parts"]["jump"]["gain"], 0.0)


class IslandDistanceTest(unittest.TestCase):
    def test_distance_and_settings_that_differ(self):
        rows = [{"island": 0, "best_vec": [0.0, 0.5, 1.0]}, {"island": 1, "best_vec": [0.0, 0.1, 0.7]},
                {"island": 2, "best_vec": None}]
        out = island_distances(rows)
        self.assertEqual(len(out), 1)
        self.assertEqual((out[0]["a"], out[0]["b"], out[0]["settings"]), (0, 1, 2))
        self.assertAlmostEqual(out[0]["distance"], 0.5)


if __name__ == "__main__":
    unittest.main()
