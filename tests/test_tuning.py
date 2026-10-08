"""common.tuning: defaults, file overrides, rejection with reasons, round settings waiting for the next round."""
import json
import os
import tempfile
import unittest
from pathlib import Path

from common.tuning import DEFAULTS, KNOBS, TUNING, Tuning, describe, read_overrides, tune, validate, write_overrides


class TuningTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.path = Path(self.dir.name) / "tuning.json"
        self.t = Tuning()
        self._n = 0

    def tearDown(self):
        self.dir.cleanup()

    def write(self, data):
        self._n += 1
        self.path.write_text(data if isinstance(data, str) else json.dumps(data))
        os.utime(self.path, ns=(10 ** 18 + self._n, 10 ** 18 + self._n))   # a new stamp even within one tick
        self.t._maybe_reload()                                                 # main notices within a second

    def test_defaults_without_a_file(self):
        self.assertIsNone(TUNING.path)                       # tests never read the repo's tuning.json
        self.assertEqual(tune("scout.share"), 0.0)
        self.assertEqual(self.t.get("climb.shake_k"), [2, 3, 5, 8])
        self.t.use_file(self.path)                           # no file yet: still the defaults
        self.assertEqual(self.t.values, DEFAULTS)

    def test_a_live_override_applies_and_its_removal_restores_the_default(self):
        self.t.use_file(self.path)
        self.write({"scout.share": 0.2, "merges.mutated": 10})
        self.assertEqual((self.t.get("scout.share"), self.t.get("merges.mutated")), (0.2, 10))
        self.assertEqual([c["name"] for c in self.t.changes], ["merges.mutated", "scout.share"])
        self.write({})
        self.assertEqual(self.t.get("scout.share"), 0.0)

    def test_bad_values_are_rejected_and_the_value_in_force_stays(self):
        self.t.use_file(self.path)
        self.write({"scout.share": 0.3})
        self.write({"scout.share": 2, "riders.global": 1.5, "climb.shake_k": [2, "x"], "nope.x": 1})
        self.assertEqual(self.t.get("scout.share"), 0.3)    # previous value, not the default
        self.assertEqual(self.t.get("riders.global"), 2)
        why = {r["name"]: r["reason"] for r in self.t.rejected}
        self.assertIn("outside", why["scout.share"])
        self.assertIn("whole number", why["riders.global"])
        self.assertIn("not a number", why["climb.shake_k"])
        self.assertEqual(why["nope.x"], "unknown setting")

    def test_broken_json_keeps_everything(self):
        self.t.use_file(self.path)
        self.write({"scout.share": 0.25})
        self.write("{not json")
        self.assertEqual(self.t.get("scout.share"), 0.25)
        self.assertIn("not valid JSON", self.t.snapshot()["error"])

    def test_checks_across_settings(self):
        ok, rejected = validate({"schedule.reduce_to_1_start_min": 30.0, "evaluator.small_batch": 2.5,
                                 "evaluator.big_batch": 2.0})
        self.assertEqual(ok, {})
        self.assertEqual({r["name"] for r in rejected},
                         {"schedule.reduce_to_1_start_min", "evaluator.small_batch", "evaluator.big_batch"})
        ok, rejected = validate({"schedule.reduce_to_3_start_min": 30, "schedule.reduce_to_3_by_min": 36})
        self.assertEqual((ok, rejected), ({"schedule.reduce_to_3_start_min": 30.0,
                                           "schedule.reduce_to_3_by_min": 36.0}, []))

    def test_round_settings_wait_for_the_next_round(self):
        self.t.use_file(self.path)
        self.t.round_started()
        self.write({"domains.max": 4, "scout.share": 0.1})
        self.assertEqual(self.t.get("scout.share"), 0.1)     # live: now
        self.assertEqual(self.t.get("domains.max"), 24)      # round: waits
        self.assertEqual(self.t.snapshot()["pending"], {"domains.max": 4})
        self.t.round_ended()
        self.assertEqual(self.t.get("domains.max"), 4)
        self.assertEqual(self.t.snapshot()["pending"], {})
        self.t.round_started()
        self.write({"domains.max": 5})
        self.t.round_started()                                # the next round's start applies it
        self.assertEqual(self.t.get("domains.max"), 5)

    def test_write_and_read_back(self):
        write_overrides(self.path, {"pick.band": 0.1})
        self.assertEqual(read_overrides(self.path), ({"pick.band": 0.1}, None))
        self.assertEqual(sorted(os.listdir(self.dir.name)), ["tuning.json"])   # no temp file left behind

    def test_every_default_passes_its_own_check(self):
        for k in KNOBS.values():
            self.assertTrue(k.check(k.default)[0], k.name)
        self.assertEqual(validate(dict(DEFAULTS)), (DEFAULTS, []))
        self.assertEqual(len(describe()), len(KNOBS))


if __name__ == "__main__":
    unittest.main()
