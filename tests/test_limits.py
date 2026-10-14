"""The validator's per-run time limit (minos_subnet PR #40): a config that ran longer here is flagged, not dropped."""
import unittest
from types import SimpleNamespace

from common.limits import VALIDATOR_RUN_S, over_validator_limit
from search.client import EvalClient


class ValidatorLimitTest(unittest.TestCase):
    def test_the_limit(self):
        self.assertEqual(VALIDATOR_RUN_S, 900)
        self.assertTrue(over_validator_limit(901))
        self.assertFalse(over_validator_limit(900))
        self.assertFalse(over_validator_limit(None))
        self.assertFalse(over_validator_limit("x"))

    def _client(self):
        stub = SimpleNamespace(round=SimpleNamespace(round_id="r"), configs={"k": object()}, _pending={},
                               callsets={}, latest={}, on_result=[], _note_best=lambda r: None,
                               store=SimpleNamespace(add=lambda *a: None, add_sketch=lambda *a: None))
        return stub

    def test_a_slow_run_is_flagged_and_still_recorded(self):
        stub = self._client()
        cs = SimpleNamespace(round_id="r", config_key="k", job_id="j", ok=True, error=None, fidelity="full",
                             runtime_s=1200.0, sketch=None)
        result = SimpleNamespace(score=80.0)
        with self.assertLogs("search.client", level="WARNING") as logs:
            EvalClient.complete(stub, cs, result)
        self.assertIn("k", stub.slow)
        self.assertIn(("k", "full"), stub.latest)
        self.assertTrue(any("over the validator" in l for l in logs.output))

    def test_a_normal_run_is_not_flagged(self):
        stub = self._client()
        cs = SimpleNamespace(round_id="r", config_key="k", job_id="j", ok=True, error=None, fidelity="full",
                             runtime_s=180.0, sketch=None)
        EvalClient.complete(stub, cs, SimpleNamespace(score=80.0))
        self.assertNotIn("k", getattr(stub, "slow", set()))


if __name__ == "__main__":
    unittest.main()
