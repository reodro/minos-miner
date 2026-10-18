"""dashboard.field.refresh_rank: a review's rank row placed again on the field as published now (2026-10-18)."""
import unittest

from dashboard.field import refresh_rank


class FakeCache:
    def __init__(self, scores):
        self.rows = scores

    def scores(self, round_id):
        return self.rows.get(round_id, [])


class RefreshRankTest(unittest.TestCase):
    def test_a_live_row_ranked_before_its_field_existed_gets_its_rank_once_the_field_is_out(self):
        row = {"round_id": "r1", "true_score": 90.97, "rank": 1, "field": 0}
        field = [91.12, 91.0, 90.97, 90.9, 90.8]                       # our own submission is one of the rows
        out = refresh_rank(row, FakeCache({"r1": field}), replay=False)
        self.assertEqual((out["rank"], out["field"]), (3, 5))
        self.assertAlmostEqual(out["gap_to_best"], 0.15, places=4)

    def test_no_published_field_means_no_rank_not_rank_one(self):
        out = refresh_rank({"round_id": "r2", "true_score": 90.0, "rank": 1, "field": 0}, FakeCache({}), replay=False)
        self.assertIsNone(out["rank"])
        self.assertTrue(out["filling"])

    def test_a_replay_joins_the_field(self):
        out = refresh_rank({"round_id": "r3", "true_score": 80.0, "rank": 2, "field": 3},
                           FakeCache({"r3": [81.0, 79.0]}), replay=True)
        self.assertEqual((out["rank"], out["field"]), (2, 3))

    def test_a_row_without_a_true_score_is_left_alone(self):
        row = {"round_id": "r4", "rank": None}
        self.assertIs(refresh_rank(row, FakeCache({"r4": [1.0]}), replay=False), row)


if __name__ == "__main__":
    unittest.main()
