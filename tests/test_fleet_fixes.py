"""Fleet fixes of 2026-10-17: rounds kept as a setting and never pruned while they have work; the `keep` list a node
receives is current when queued; a run the node never reports is reclaimed; copies are capped; a node never drops a
finished callset silently; the agent's watchdog; the live TRUE scorer only scores whole callsets."""
import json
import tempfile
import threading
import time
import unittest
import unittest.mock
from pathlib import Path

from common.tuning import TUNING
from common.types import CallSet, Config, EvalJob, Round
from evaluator.fleet.coordinator import FleetExecutor, _Session
from tests.helpers import ranges
from tests.test_fleet import Planner


def _round(i: int) -> Round:
    return Round(f"2026-01-0{i}T00:00:00+00:00", "chr20", 1, 1000, f"/r{i}/input.bam", "")


def _ex() -> FleetExecutor:
    ex = FleetExecutor(Planner(), inventory=Path("/nonexistent"))
    ex.refresh_hosts = lambda: None                              # no inventory, no ssh
    return ex


class Sent(list):
    """A session's outgoing messages, captured instead of written to ssh."""


def _session(workers=4, rounds=(), solo_s=100.0) -> _Session:
    s = _Session(host=None, workers=workers, ready=True, rounds_ready=set(rounds), solo_s=solo_s)
    s.proc = object()                                           # send() needs a process; nothing is written
    sent = Sent()
    s.send = lambda msg, sent=sent: sent.append(msg) or True     # type: ignore[method-assign]
    s.sent = sent                                                # type: ignore[attr-defined]
    return s


class KeepRoundsTest(unittest.TestCase):
    def test_the_oldest_idle_round_is_pruned_past_the_setting(self):
        ex = _ex()
        with unittest.mock.patch.dict(TUNING.values, {"evaluator.keep_rounds": 2}):
            for i in (1, 2, 3):
                ex.publish_round(_round(i))
        self.assertEqual(list(ex._rounds), [_round(2).round_id, _round(3).round_id])

    def test_a_round_with_queued_or_running_jobs_is_never_pruned(self):
        ex = _ex()
        d = ranges().defaults("gatk").params
        with unittest.mock.patch.dict(TUNING.values, {"evaluator.keep_rounds": 2}):
            ex.publish_round(_round(1))
            ex.submit([EvalJob("q1", _round(1), Config("gatk", dict(d)))])   # queued: no node has the round
            ex.publish_round(_round(2))
            ex.publish_round(_round(3))
            self.assertIn(_round(1).round_id, ex._rounds)            # busy: kept beyond the setting
            self.assertNotIn(_round(2).round_id, ex._rounds)         # the oldest idle one went instead
            ex.end_round(_round(1).round_id)                         # its work is over
            ex.publish_round(_round(4))
        self.assertNotIn(_round(1).round_id, ex._rounds)

    def test_more_rounds_for_multi_round_tuning(self):
        ex = _ex()
        with unittest.mock.patch.dict(TUNING.values, {"evaluator.keep_rounds": 9}):
            for i in range(1, 10):
                ex.publish_round(_round(i))
        self.assertEqual(len(ex._rounds), 9)

    def test_a_reopened_round_takes_results_again(self):
        ex = _ex()
        ex.publish_round(_round(1))
        ex.end_round(_round(1).round_id)
        self.assertIn(_round(1).round_id, ex._ended)
        ex.publish_round(_round(1))
        self.assertNotIn(_round(1).round_id, ex._ended)


class KeepListTest(unittest.TestCase):
    def test_the_keep_list_is_read_when_the_message_is_queued(self):
        """A slow BAM push for round 1 finishes after rounds 2 and 3 were published: round 1's message must still
        name them, or the node deletes them while their jobs run (2026-10-17)."""
        import evaluator.fleet.coordinator as C
        ex = _ex()
        s = _session()
        ex._sessions = {"n": s}
        release = threading.Event()

        def round_message(h, r, bam_url, bai_url):
            if r.round_id == _round(1).round_id:
                release.wait(5)                                     # the slow push
            return {"type": "round", "round": {"round_id": r.round_id}}

        with unittest.mock.patch.object(C, "round_message", round_message):
            ex.publish_round(_round(1))                             # its thread blocks in the push
            time.sleep(0.2)
            ex.publish_round(_round(2))
            ex.publish_round(_round(3))
            time.sleep(0.3)
            release.set()
            time.sleep(0.5)
        last = [m for m in s.sent if m.get("type") == "round"][-1]
        self.assertEqual(last["round"]["round_id"], _round(1).round_id)    # arrived last ...
        self.assertEqual(set(last["keep"]), {_round(i).round_id for i in (1, 2, 3)})   # ... and keeps 2 and 3


class LostRunTest(unittest.TestCase):
    def test_a_run_the_node_never_reports_is_queued_again_and_dropped_there(self):
        import evaluator.fleet.coordinator as C
        ex = _ex()
        r = _round(1)
        ex._rounds[r.round_id] = C._Published(r)
        d = ranges().defaults("gatk").params
        a, b = _session(rounds={r.round_id}, solo_s=100.0), _session(rounds={r.round_id}, solo_s=100.0)
        a.workers, b.workers = 1, 0                                  # the job goes to a; b joins later
        ex._sessions = {"a": a, "b": b}
        ex._now = lambda: 1000.0
        ex.submit([EvalJob("L", r, Config("gatk", dict(d)))])
        self.assertIn("L", a.jobs)
        ex._now = lambda: 1000.0 + 500.0                           # long, but under the 10-minute floor
        ex.poll(0)
        self.assertIn("L", a.jobs)
        b.workers = 1
        ex._now = lambda: 1000.0 + 700.0                           # > 6 x 100 s and > 600 s: lost
        with unittest.mock.patch.object(C, "round_message",
                                        lambda h, rr, u, v: {"type": "round", "round": {"round_id": rr.round_id}}):
            ex.poll(0)
            time.sleep(0.3)
        self.assertNotIn("L", a.jobs)
        self.assertIn({"type": "abort_jobs", "job_ids": ["L"]}, a.sent)
        self.assertIn("L", b.jobs)                                   # queued again; the other node runs it
        self.assertIn(r.round_id, [m["round"]["round_id"] for m in a.sent if m.get("type") == "round"])  # re-shipped
        cs = CallSet(r.round_id, Config("gatk", dict(d)).key, "/x.vcf.gz", "d", 1.0, job_id="L")
        ex._result(b, cs)
        self.assertEqual([c.job_id for c in ex.poll(0)], ["L"])


class CopyCapTest(unittest.TestCase):
    def _overdue_fleet(self, n_batches: int, workers: int):
        ex = _ex()
        r = _round(1)
        d = ranges().defaults("gatk").params
        slow = _session(workers=n_batches, rounds={r.round_id}, solo_s=100.0)
        ex._sessions = {"slow": slow}
        ex._now = lambda: 1000.0
        ex.submit([EvalJob(f"L{i}", r, Config("gatk", {**d, "min_pruning": 2 + i % 8, "max_alternate_alleles": 1 + i}),
                           batch="small") for i in range(n_batches)])
        fast = _session(workers=workers, rounds={r.round_id}, solo_s=50.0)
        ex._sessions["fast"] = fast
        ex._now = lambda: 1000.0 + 2.0 * 100.0                     # every batch is overdue
        return ex

    def test_copies_hold_at_most_the_set_share_of_the_workers(self):
        ex = self._overdue_fleet(10, 90)                             # 100 workers in all
        with unittest.mock.patch.dict(TUNING.values, {"evaluator.max_copy_share": 0.05}):
            ex._speculate()
        self.assertEqual(ex.speculated, 5)

    def test_share_zero_turns_copies_off(self):
        ex = self._overdue_fleet(4, 10)
        with unittest.mock.patch.dict(TUNING.values, {"evaluator.max_copy_share": 0.0}):
            ex._speculate()
        self.assertEqual(ex.speculated, 0)


class NodeTest(unittest.TestCase):
    def _node(self):
        from evaluator.fake import FakeRunner, FakeScorer
        from evaluator.fleet.node import Node
        from evaluator.fleet.inventory import REPO
        from common.params import ParamRanges
        runner = FakeRunner(ParamRanges.from_file(str(REPO / "tests" / "fixtures" / "parameter_ranges.json")))
        node = Node(2, runner, lambda i: FakeScorer(runner), score_threads=1)
        out = []
        node.send = out.append                                       # type: ignore[method-assign]
        return node, out

    def test_a_finished_callset_of_a_forgotten_round_is_reported_as_a_failure(self):
        from evaluator.fleet import codec
        node, out = self._node()
        node._scorers.append(node.scorer_factory(0))
        t = threading.Thread(target=node._scorer, args=(0,), daemon=True)
        t.start()
        cs = CallSet("2026-01-01T00:00:00+00:00", "k", "/x.vcf.gz", "d", 1.0, job_id="J1")
        node._to_score.put([cs])                                     # its round is not (any more) on this node
        time.sleep(1.5)
        node._stop.set()
        got = [codec.callset_from(m["callset"]) for m in out if m.get("type") == "result"]
        self.assertEqual([(c.job_id, c.ok) for c in got], [("J1", False)])
        self.assertIn("forgotten", got[0].error)

    def test_units_held_for_an_answer_set_are_reported_when_their_round_is_forgotten(self):
        from evaluator.fleet import codec
        node, out = self._node()
        rid = "2026-01-01T00:00:00+00:00"
        node._held[rid] = [[CallSet(rid, "k", "/x.vcf.gz", "d", 1.0, job_id="J2")]]
        node.forget_rounds(set())
        got = [codec.callset_from(m["callset"]) for m in out if m.get("type") == "result"]
        self.assertEqual([(c.job_id, c.ok) for c in got], [("J2", False)])

    def test_the_watchdog_waits_for_main_s_first_ping(self):
        import evaluator.fleet.node as N
        node, _ = self._node()
        node._heard = time.time() - 10 * N.WATCHDOG_S
        self.assertFalse(node._silent())                             # an older main never pings
        node._pinged = True
        self.assertTrue(node._silent())
        node._heard = time.time()
        self.assertFalse(node._silent())

    def test_an_agent_that_died_has_its_gatk_jobs_reaped_and_a_live_one_does_not(self):
        import os
        import evaluator.fleet.node as N
        with tempfile.TemporaryDirectory() as d, unittest.mock.patch.object(N, "AGENTS", d), \
                unittest.mock.patch.object(N, "reap") as reap:
            Path(d, "minos-node-dead01.pid").write_text("999999999")
            Path(d, "minos-node-live01.pid").write_text(str(os.getpid()))
            N.register("minos-node-new001")
            reap.assert_called_once_with("minos-node-dead01")
            self.assertEqual(sorted(os.listdir(d)), ["minos-node-live01.pid", "minos-node-new001.pid"])
            N.unregister("minos-node-new001")
            self.assertEqual(os.listdir(d), ["minos-node-live01.pid"])


class TrueLiveTest(unittest.TestCase):
    def test_only_a_whole_bgzf_file_counts(self):
        import gzip
        from dashboard.true_live import BGZF_EOF, _complete
        with tempfile.TemporaryDirectory() as d:
            whole = Path(d, "a.vcf.gz")
            whole.write_bytes(gzip.compress(b"##fileformat=VCFv4.2\n") + BGZF_EOF)
            cut = Path(d, "b.vcf.gz")
            cut.write_bytes(whole.read_bytes()[:-5])
            self.assertTrue(_complete(whole))
            self.assertFalse(_complete(cut))
            self.assertFalse(_complete(Path(d, "missing.vcf.gz")))


if __name__ == "__main__":
    unittest.main()
