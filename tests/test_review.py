import json
import tempfile
import unittest
from pathlib import Path

from common.types import Variant
from main.review import answer_set_accuracy, history, pick_configs, summarize


def v(pos, gt="0/1", alt="G"):
    return Variant("chr20", pos, "A", alt, gt)


class AccuracyTest(unittest.TestCase):
    def test_germline_exact_needs_genotype_but_sites_do_not(self):
        ours = [v(1), v(2, "1/1"), v(3)]
        truth = [v(1, "1|0"), v(2, "0/1"), v(4)]                 # phased GT counts as the same genotype
        a = answer_set_accuracy(ours, [], truth, [], 10, 20)
        self.assertEqual(a["germline_exact"]["match"], 1)
        self.assertAlmostEqual(a["germline_exact"]["precision"], 1 / 3, places=3)
        self.assertEqual(a["germline_sites"]["match"], 2)
        self.assertEqual(a["n_targets"]["ratio"], 0.5)

    def test_synthetic_all_and_confident(self):
        ours = [(v(10), 1.0), (v(11), 0.5), (v(12), 0.99)]
        truth = [v(10), v(11), v(13)]
        a = answer_set_accuracy([], ours, [], truth, None, None)
        self.assertEqual(a["synthetic"]["match"], 2)
        self.assertAlmostEqual(a["synthetic"]["recall"], 2 / 3, places=3)
        self.assertEqual(a["synthetic_confident"]["ours"], 2)
        self.assertEqual(a["synthetic_confident"]["match"], 1)
        self.assertIsNone(a["n_targets"]["ratio"])
        self.assertIsNone(a["germline_exact"]["precision"])      # nothing of ours: undefined, not 0


class PickAndSummaryTest(unittest.TestCase):
    def test_pick_final_baseline_then_top(self):
        configs = {k: {"local": s} for k, s in [("a", 90), ("b", 80), ("c", 85), ("base", 70), ("x", None), ("z", 0.0)]}
        self.assertEqual(pick_configs(configs, "b", "base", 2), ["b", "base", "a", "c"])
        self.assertEqual(pick_configs(configs, "a", "gone", 1), ["a"])           # unknown keys are skipped
        self.assertNotIn("z", pick_configs(configs, "a", "base", 10))            # scored 0: not worth a true score

    def test_summary(self):
        local = {"f": 90.0, "b": 80.0, "t1": 88.0, "t2": 85.0}
        true = {"f": 87.0, "b": 84.0, "t1": 88.0, "t2": 86.0}
        s = summarize(local, true, "f", "b")
        self.assertEqual(s["final_true"], 87.0)
        self.assertEqual(s["best_true_config"], "t1")
        self.assertEqual(s["pick_gap"], 1.0)
        self.assertEqual(s["gain_over_baseline"], 3.0)
        self.assertEqual(s["local_minus_true"], 3.0)
        self.assertEqual(s["n"], 4)
        self.assertAlmostEqual(s["spearman"], 0.8, places=4)

    def test_summary_without_truth(self):
        s = summarize({"f": 90.0}, {}, "f", None)
        self.assertIsNone(s["final_true"])
        self.assertIsNone(s["pick_gap"])
        self.assertIsNone(s["spearman"])


class HistoryTest(unittest.TestCase):
    def test_reviewed_and_ranked_runs(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "rounds.jsonl"
            rows = [{"type": "round", "run_id": "T1-1", "round_id": "r1"},
                    {"type": "round", "round_id": "r0"},                     # old format: no run id, never reviewed
                    {"type": "review", "run_id": "T1-1"},
                    {"type": "round", "run_id": "T2-2", "round_id": "r2"},
                    {"type": "rank", "run_id": "T1-1"}]
            p.write_text("\n".join(json.dumps(r) for r in rows) + "\nnot json\n")
            rounds, reviewed, ranked = history(p)
            self.assertEqual([r["run_id"] for r in rounds], ["T1-1", "T2-2"])
            self.assertEqual(reviewed, {"T1-1"})
            self.assertEqual(ranked, {"T1-1"})


if __name__ == "__main__":
    unittest.main()
