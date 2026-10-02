import unittest
from dataclasses import replace

from common.types import AnswerSet, CallSet, Regime, ScoreResult
from evaluator.fleet.node import Node


class _CountingScorer:
    def __init__(self):
        self.calls = []

    def score(self, cs, aset):
        self.calls.append(cs.config_key)
        return ScoreResult(cs.round_id, cs.config_key, aset.version, {}, 0, 1, score=50.0)


class _NoRunner:
    def run(self, *a, **k):
        raise AssertionError("not used")


class NodeScoringUnitTest(unittest.TestCase):
    def setUp(self):
        self.scorer = _CountingScorer()
        self.node = Node(2, _NoRunner(), lambda i: self.scorer, score_threads=1)
        self.node._scorers.append(self.scorer)
        self.sent = []
        self.node.send = self.sent.append
        aset = AnswerSet("r1", 1, Regime.UNKNOWN, None, [], [], 10)
        self.node._rounds["r1"] = (None, aset)

    def _score(self, unit):
        self.node._to_score.put(unit)
        import threading, time
        t = threading.Thread(target=self.node._scorer, args=(0,), daemon=True)
        t.start()
        deadline = time.time() + 5
        while self.node._to_score.qsize() and time.time() < deadline:
            time.sleep(0.01)
        time.sleep(0.2)
        self.node._stop.set()
        t.join(timeout=3)

    def test_unit_scores_each_config_once_and_reuses_identical_records(self):
        cs = lambda key, digest: CallSet("r1", key, f"/x/{key}.vcf.gz", digest, 0.0)   # noqa: E731
        self._score([cs("a", "d1"), cs("b", "d2"), cs("c", "d1")])   # c has a's records
        self.assertEqual(sorted(self.scorer.calls), ["a", "b"])
        got = {m["callset"]["config_key"]: m["callset"]["score"]["config_key"] for m in self.sent}
        self.assertEqual(got, {"a": "a", "b": "b", "c": "c"})        # c's result is its own (key carried over)

    def test_records_seen_earlier_in_the_round_are_not_scored_again(self):
        first = CallSet("r1", "a", "/x/a.vcf.gz", "d1", 0.0)
        self.node._remember(first, self.node._rounds["r1"][1], self.scorer.score(first, self.node._rounds["r1"][1]))
        self.scorer.calls.clear()
        self._score([replace(first, config_key="z", job_id="j9")])
        self.assertEqual(self.scorer.calls, [])
        self.assertEqual(self.sent[0]["callset"]["score"]["config_key"], "z")


if __name__ == "__main__":
    unittest.main()
