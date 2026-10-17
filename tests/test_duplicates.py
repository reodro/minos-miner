"""The duplicate guard (search/duplicates.py): park predicted repeats of a callset, never drop them on a prediction."""
import asyncio
import random
import unittest
import unittest.mock

from common.tuning import TUNING
from common.types import FULL, ScoreResult
from search import engine
from search.climb import Climber
from search.duplicates import PARKED, DuplicateGuard
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges
from tests.test_islands import FakeBase, JumpCtx

H = "max_num_haplotypes_in_population"
MMQ = "min_mapping_quality_score"
PAD = "assembly_region_padding"


class GuardTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        lad = Ladders(self.space)
        self.h = lad.rungs[H][2:4]                  # two neighbouring rungs: the change under test
        self.mmq = lad.rungs[MMQ][:12]
        self.g = DuplicateGuard(self.space.names, random.Random(0))
        self.knobs = unittest.mock.patch.dict(TUNING.values, {"dups.min_ties": 5, "dups.audit": 0.0, "dups.near": 0.2})
        self.knobs.start()

    def tearDown(self):
        self.knobs.stop()

    def cfg(self, h, mmq, pad=50):
        return self.space._assemble({H: h, MMQ: mmq, PAD: pad})

    def teach(self, n, score=80.0):
        """n contexts where the H change tied."""
        for i in range(n):
            for h in self.h:
                c = self.cfg(h, self.mmq[i])
                self.g.note_sent(c)
                self.g.note_score(c, score, 82.0)

    def test_a_change_is_a_rule_only_after_enough_ties_and_no_difference(self):
        self.teach(4)
        base = self.cfg(self.h[0], self.mmq[10])
        self.g.note_sent(base)
        self.assertFalse(self.g.decide(self.cfg(self.h[1], self.mmq[10]), "B0", 82.0))       # 4 ties: no rule
        self.teach(5)
        self.assertTrue(self.g.decide(self.cfg(self.h[1], self.mmq[10]), "B0", 82.0))
        self.assertEqual(self.g.still_parked(), 1)

    def test_a_different_score_breaks_the_rule_and_sends_what_it_parked(self):
        self.teach(5)
        base = self.cfg(self.h[0], self.mmq[10])
        self.g.note_sent(base)
        parked = self.cfg(self.h[1], self.mmq[10])
        self.assertTrue(self.g.decide(parked, "B0", 82.0))
        a, b = self.cfg(self.h[0], self.mmq[11]), self.cfg(self.h[1], self.mmq[11])
        self.g.note_score(a, 79.0, 82.0)
        out = self.g.note_score(b, 79.5, 82.0)                       # the same change, a different score
        self.assertEqual([c.key for c, _ in out], [parked.key])
        self.assertEqual(self.g.still_parked(), 0)
        self.assertFalse(self.g.decide(self.cfg(self.h[1], self.mmq[9], pad=60), "B0", 82.0))  # broken for good

    def test_a_twin_near_the_best_sends_the_config_now_or_when_it_returns(self):
        self.teach(5)
        scored = self.cfg(self.h[0], self.mmq[10])
        self.g.note_sent(scored)
        self.g.note_score(scored, 81.9, 82.0)                         # within 0.2 of the best
        self.assertFalse(self.g.decide(self.cfg(self.h[1], self.mmq[10]), "B0", 82.0))
        running = self.cfg(self.h[0], self.mmq[11])
        self.g.note_sent(running)
        waiting = self.cfg(self.h[1], self.mmq[11])
        self.assertTrue(self.g.decide(waiting, "B0", 82.0))           # twin still running: parked on it
        out = self.g.note_score(running, 82.1, 82.1)                  # the twin comes back as the new best
        self.assertEqual([c.key for c, _ in out], [waiting.key])

    def test_a_twin_far_below_the_best_keeps_the_config_parked(self):
        self.teach(5)
        running = self.cfg(self.h[0], self.mmq[11])
        self.g.note_sent(running)
        waiting = self.cfg(self.h[1], self.mmq[11])
        self.assertTrue(self.g.decide(waiting, "B0", 82.0))
        self.assertEqual(self.g.note_score(running, 80.0, 82.0), [])
        self.assertEqual(self.g.still_parked(), 1)

    def test_audits_send_a_share_of_the_predicted_repeats(self):
        self.teach(5)
        self.g.note_sent(self.cfg(self.h[0], self.mmq[10]))
        with unittest.mock.patch.dict(TUNING.values, {"dups.audit": 1.0}):
            self.assertFalse(self.g.decide(self.cfg(self.h[1], self.mmq[10]), "B0", 82.0))
        self.assertEqual(self.g.n_audits, 1)

    def test_a_hash_collision_never_makes_a_twin(self):
        self.teach(5)
        self.g.note_sent(self.cfg(self.h[0], self.mmq[10]))
        with unittest.mock.patch.object(DuplicateGuard, "_proj", staticmethod(lambda v, d: 0)):
            g = DuplicateGuard(self.space.names, random.Random(0))
            g.rules = self.g.rules
            g.note_sent(self.cfg(self.h[0], self.mmq[10], pad=60))    # same bucket, not one change away
            self.assertFalse(g.decide(self.cfg(self.h[1], self.mmq[10]), "B0", 82.0))

    def test_switched_off_nothing_is_parked(self):
        self.teach(5)
        self.g.note_sent(self.cfg(self.h[0], self.mmq[10]))
        with unittest.mock.patch.dict(TUNING.values, {"dups.enabled": 0}):
            self.assertFalse(self.g.decide(self.cfg(self.h[1], self.mmq[10]), "B0", 82.0))


class ContextTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        self.knobs = unittest.mock.patch.dict(TUNING.values, {"dups.min_ties": 1, "dups.audit": 0.0, "dups.near": 0.2})
        self.knobs.start()
        lad = Ladders(self.space)
        self.h = lad.rungs[H][2:4]
        self.mmq = lad.rungs[MMQ][:6]

    def tearDown(self):
        self.knobs.stop()
        self.loop.close()

    def cfg(self, h, mmq):
        return self.space._assemble({H: h, MMQ: mmq, PAD: 50})

    def ctx_with_rule(self):
        base = FakeBase()
        ctx = engine.Context(base, self.space, random.Random(0), None)
        ctx._top = 90.0
        for h in self.h:                                               # one tie establishes the rule (min_ties 1)
            c = self.cfg(h, self.mmq[0])
            ctx.propose([c], "G", ["big"])
            base.finish(c.key, 80.0)
        ctx._top = 90.0
        return base, ctx

    def test_a_parked_config_is_done_at_once_and_not_sent(self):
        base, ctx = self.ctx_with_rule()
        twin = self.cfg(self.h[0], self.mmq[1])
        ctx.propose([twin], "B0", ["small"])
        rep = self.cfg(self.h[1], self.mmq[1])
        rider = self.space._assemble({**rep.params, "indel_heterozygosity": 0.0001})
        fs = ctx.propose([rep, rider], "B0", ["small", "big"])
        self.assertTrue(fs[0].done())
        self.assertEqual(fs[0].result().error, PARKED)
        self.assertNotIn(rep.key, base.pending)
        self.assertIn(rider.key, base.pending)                         # not predicted: sent
        self.assertIn(rep.key, ctx.proposed)                           # never proposed again in a loop

    def test_nothing_is_parked_while_the_pick_is_certified(self):
        base, ctx = self.ctx_with_rule()
        ctx.propose([self.cfg(self.h[0], self.mmq[1])], "B0", ["small"])
        ctx.phase = "C: certify"
        rep = self.cfg(self.h[1], self.mmq[1])
        ctx.propose([rep], "C", ["small"])
        self.assertIn(rep.key, base.pending)

    def test_a_released_config_is_sent_and_can_move_a_climber(self):
        base, ctx = self.ctx_with_rule()
        twin = self.cfg(self.h[0], self.mmq[1])
        ctx.propose([twin], "B0", ["small"])
        rep = self.cfg(self.h[1], self.mmq[1])
        ctx.propose([rep], "B0", ["small"])
        self.assertNotIn(rep.key, base.pending)
        base.finish(twin.key, 89.9)                                    # the twin is near the best
        ctx.send_released()
        self.assertIn(rep.key, base.pending)
        self.assertIn(rep.key, ctx.ga_keys)


class ClimberParkedTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()

    def test_a_parked_run_is_neither_a_crash_nor_a_loss(self):
        from search.client import Outcome
        failed = {}
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.note_failure = lambda source, reason: failed.__setitem__(source, failed.get(source, 0) + 1)
        centre = self.space._assemble({MMQ: 40, PAD: 50})
        cl = Climber(ctx, 0, centre, 80.0)
        c = self.space._assemble({MMQ: 40, PAD: 40})
        fut = self.loop.create_future()
        fut.set_result(Outcome(c, None, PARKED, "B0"))
        cl._note_sent([fut], [c], "jump")
        cl.leap_meta[c.key] = ((PAD, 2, -1), 80.0, centre, c)
        cl.waiting[c.key] = c
        self.assertEqual(cl.drop_failed(), 0)
        self.assertEqual(failed, {})
        self.assertNotIn(c.key, cl.waiting)
        self.assertEqual(cl.crash_edge, {})
