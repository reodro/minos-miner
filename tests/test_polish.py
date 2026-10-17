"""The test-only polish mode (search/polish.py): climbs from given centres, moves at once, names the pick."""
import asyncio
import json
import os
import tempfile
import unittest
import unittest.mock

from common.types import FULL, ScoreResult
from search import polish
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges

MMQ = "min_mapping_quality_score"


class FakeClient:
    """Scores a config one tick after it is sent: higher the closer min_mapping_quality_score is to 50."""

    def __init__(self, space):
        self.space, self.latest, self.configs, self.pick_key = space, {}, {}, None
        self.left = 400.0
        self.round = type("R", (), {"tag": "20990101T000000Z"})()

    def time_left(self):
        self.left -= 2.0
        return self.left

    def idle_slots(self):
        return 50

    def f(self, c):
        return 80.0 - abs(int(c.params.get(MMQ)) - 50) * 0.01

    def propose_many(self, configs, *, tag="", batches=None, priority=0):
        loop = asyncio.get_event_loop()
        out = []
        for c in configs:
            fut = loop.create_future()
            self.configs[c.key] = c
            r = ScoreResult.__new__(ScoreResult)
            r.score, r.gate_passed, r.fidelity = self.f(c), True, FULL
            loop.call_later(0.01, lambda c=c, r=r, fut=fut: (self.latest.__setitem__((c.key, FULL), r), fut.set_result(None)))
            out.append(fut)
        return out


class PolishTest(unittest.TestCase):
    def test_it_climbs_from_the_centre_and_names_the_best_as_the_pick(self):
        space = Space(ranges())
        centre = space._assemble({MMQ: 40})
        spec = {"target": 79.99, "centres": [{"key": centre.key, "config": centre.to_dict(), "score": 79.90}],
                "known": {centre.key: 79.90}}
        client = FakeClient(space)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "spec.json")
            json.dump(spec, open(path, "w"))
            with unittest.mock.patch.object(polish, "TICK_S", 0.02), unittest.mock.patch.object(polish, "STOP_S", 20.0), \
                    unittest.mock.patch("search.polish.os.path.dirname", side_effect=lambda p: d):
                asyncio.run(polish.run(client, space, path))
        best = max(r.score for r in client.latest.values())
        self.assertGreater(best, 79.95)                     # several moves up the hill, not one
        self.assertEqual(client.latest[(client.pick_key, FULL)].score, best)

    def test_a_lower_line_climbs_on_its_own_while_the_best_is_elsewhere(self):
        space = Space(ranges())
        top = space._assemble({MMQ: 50, "assembly_region_padding": 60})
        low = space._assemble({MMQ: 30, "assembly_region_padding": 40})
        spec = {"target": 81.0, "known": {top.key: 80.0, low.key: 79.8},
                "centres": [{"key": top.key, "config": top.to_dict(), "score": 80.0, "label": "best"},
                            {"key": low.key, "config": low.to_dict(), "score": 79.8, "label": "indel 8bp"}]}
        client = FakeClient(space)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "spec.json")
            json.dump(spec, open(path, "w"))
            with unittest.mock.patch.object(polish, "TICK_S", 0.02), unittest.mock.patch.object(polish, "STOP_S", 20.0), \
                    unittest.mock.patch("search.polish.os.path.dirname", side_effect=lambda p: d):
                asyncio.run(polish.run(client, space, path))
            out = json.load(open([os.path.join(d, "test", "runs", f) for f in os.listdir(os.path.join(d, "test", "runs"))][0]))
        lines = {ln["label"]: ln for ln in out["lines"]}
        self.assertGreater(lines["indel 8bp"]["best"], 79.85)      # it climbed although it never led the round
        self.assertEqual(lines["best"]["best"], 80.0)

    def test_a_frozen_guarded_line_keeps_its_setting_and_its_class(self):
        from common.types import ClassCounts
        space = Space(ranges())
        c = space._assemble({MMQ: 40, "min_assembly_region_size": 28})
        spec = {"known": {}, "centres": [{"key": c.key, "config": c.to_dict(), "score": 80.0, "label": "8bp",
                                          "freeze": {"min_assembly_region_size": 28}, "keep_class": "indel_8bp",
                                          "keep_min_f1": 0.8}]}
        p = polish.Polish(space, spec)
        p.plan_centre(c, 0)
        sent = [x for q in p.queue for x in q[1]]
        self.assertTrue(sent)
        self.assertTrue(all(x.params.get("min_assembly_region_size") == 28 for x in sent))
        good = type("R", (), {"counts": {"indel_8bp": ClassCounts(tp=80, fp=10, fn=10)}})()   # F1 0.889
        bad = type("R", (), {"counts": {"indel_8bp": ClassCounts(tp=70, fp=20, fn=20)}})()    # F1 0.778
        self.assertTrue(p.keeps(0, good))
        self.assertFalse(p.keeps(0, bad))

    def test_the_cheap_sweep_covers_the_calling_threshold_ladder(self):
        space = Space(ranges())
        p = polish.Polish(space, {"centres": [], "known": {}})
        centre = space._assemble({MMQ: 40})
        n = "standard_min_confidence_threshold_for_calling"
        vals = p.cheap_values(centre, n)
        self.assertGreaterEqual(len(vals), len(Ladders(space).rungs[n]) // 3)


if __name__ == "__main__":
    unittest.main()
