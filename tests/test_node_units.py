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
        self.node._rounds["r1"] = None
        self.node._asets["r1"] = aset

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
        self.node._remember(first, self.node._asets["r1"], self.scorer.score(first, self.node._asets["r1"]))
        self.scorer.calls.clear()
        self._score([replace(first, config_key="z", job_id="j9")])
        self.assertEqual(self.scorer.calls, [])
        self.assertEqual(self.sent[0]["callset"]["score"]["config_key"], "z")

    def test_units_made_before_the_answer_set_are_held_and_scored_when_it_arrives(self):
        import tempfile
        from pathlib import Path
        self.node._rounds["r2"] = None                          # GATK can run: the round is prepared
        unit = [CallSet("r2", "a", "/x/a.vcf.gz", "d1", 0.0, job_id="j1")]
        self._score(unit)
        self.assertEqual(self.scorer.calls, [])                 # nothing to score it against yet
        self.assertEqual(self.sent, [])
        self.assertEqual(self.node._held, {"r2": [unit]})
        d = Path(tempfile.mkdtemp())
        for f in ("truth.vcf.gz", "mutations.vcf.gz"):
            (d / f).write_text("x")
        self.node.set_answer_set({"round_id": "r2", "answer_set": {"version": 1, "n_targets": 5, "notes": {
            "truth_vcf": str(d / "truth.vcf.gz"), "mutations_vcf": str(d / "mutations.vcf.gz")}}})
        self.assertEqual(self.node._held, {})
        self.assertIn("truth_id", self.node._asets["r2"].notes)
        self.node._stop.clear()
        self._score(self.node._to_score.get(timeout=1))         # the released unit, scored like any other
        self.assertEqual(self.scorer.calls, ["a"])
        self.assertEqual(self.sent[0]["callset"]["job_id"], "j1")

    def test_aborted_jobs_are_never_reported(self):
        self.node.abort_jobs(["j7"])
        self._score([CallSet("r1", "a", "/x/a.vcf.gz", "d1", 0.0, job_id="j7"),
                     CallSet("r1", "b", "/x/b.vcf.gz", "d2", 0.0, job_id="j8")])
        self.assertEqual(self.scorer.calls, ["b"])
        self.assertEqual([m["callset"]["job_id"] for m in self.sent], ["j8"])


class NodeForgetTest(unittest.TestCase):
    def test_a_node_on_mains_machine_never_deletes_the_round_dir(self):
        """A node sharing main's data dir (host "local": the bench) forgets a round without deleting its directory,
        which is main's stored round; a node with its own copy deletes it."""
        import os
        import tempfile
        from pathlib import Path
        from common.settings import SETTINGS
        from common.types import Round
        base = Path(SETTINGS.cache_dir) / "test_rounds"                 # under our own folders, where deleting is allowed
        base.mkdir(parents=True, exist_ok=True)
        for shared, kept in ((True, True), (False, False)):
            d = Path(tempfile.mkdtemp(dir=base))
            (d / "input.bam").write_text("x")
            node = Node(1, _NoRunner(), lambda i: _CountingScorer(), score_threads=1)
            node._shared_data = shared
            node._rounds["r1"] = Round("r1", "chr20", 1, 10, str(d / "input.bam"), "")
            node.forget_rounds(set())
            self.assertEqual(os.path.isdir(d), kept, f"shared={shared}")
            if kept:
                import shutil
                shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
