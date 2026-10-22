"""Centre jumps (search/leaps.py) and islands (search/islands.py)."""
import asyncio
import collections
import random
import unittest
import unittest.mock

from common.tuning import TUNING
from common.types import FULL, ScoreResult
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search import engine, islands, leaps
from search.candidates import OFFICIAL_SEED, StepA
from search.climb import Climber
from search.space import ASSEMBLY_PARAMS, ACTIVITY_PARAMS, FULL_RUN_PARAMS, Space
from search.steps import Ladders
from tests.helpers import ranges
from tests.test_engine import ROUND, Stub

MMQ = "min_mapping_quality_score"


class JumpCtx:
    """A climber's context whose configs score only once they are sent: f(config) decides the score."""

    def __init__(self, space, f, history=()):
        self.space, self.f = space, f
        self.ladders = Ladders(space)
        self.climbers = []
        self.now = 0.0
        self.cap = 40
        self.done = {}
        for c, s in history:
            self.done[c.key] = (c, s)

    def over(self):
        return False

    def clock(self):
        return self.now

    def elapsed(self):
        return self.now

    def say(self, *a):
        pass

    def owner(self, c, me):
        return None

    def score(self, c):
        hit = self.done.get(c.key)
        return None if hit is None else hit[1]

    def scored(self):
        return list(self.done.values())

    def fresh(self, c):
        return c is not None and c.key not in self.done

    def propose(self, batch, tag, batches, priority=0):
        loop = asyncio.get_event_loop()
        out = []
        for c in batch:
            self.done[c.key] = (c, self.f(c))
            fut = loop.create_future()
            fut.set_result(None)
            out.append(fut)
        return out


def _rungs(ctx, name, a, b):
    return abs(ctx.ladders.index(name, b) - ctx.ladders.index(name, a))


class LeapTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.centre = self.space._assemble({MMQ: 40})
        self.names = [n for n in FULL_RUN_PARAMS + ASSEMBLY_PARAMS + ACTIVITY_PARAMS if n in self.space.names]

    def test_a_leap_moves_one_expensive_setting_two_or_more_rungs_nearest_first(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        found = leaps.candidates(ctx, self.centre, self.names, lambda c: False)
        self.assertTrue(found)
        gaps = []
        for name, k, c in found:
            moved = [n for n in self.space.names if c.params.get(n) != self.centre.params.get(n)]
            self.assertEqual(moved, [name])
            self.assertIn(name, self.names)
            self.assertGreaterEqual(abs(k), 2)
            gaps.append(abs(k))
        self.assertEqual(gaps, sorted(gaps))
        self.assertIn(13, gaps)

    def test_a_scored_or_refused_leap_is_not_planned(self):
        down2 = self.space._assemble({**self.centre.params, MMQ: 38})
        ctx = JumpCtx(self.space, lambda c: 80.0, history=[(down2, 79.0)])
        up2 = self.space._assemble({**self.centre.params, MMQ: 42})
        keys = {c.key for _, _, c in leaps.candidates(ctx, self.centre, self.names, lambda c: c.key == up2.key)}
        self.assertNotIn(down2.key, keys)
        self.assertNotIn(up2.key, keys)
        self.assertIn(self.space._assemble({**self.centre.params, MMQ: 35}).key, keys)

    def test_expected_value_waits_for_twenty_scores(self):
        hist = [(self.space._assemble({MMQ: 20 + i}), 80.0 + i) for i in range(10)]
        self.assertEqual(leaps.expected(JumpCtx(self.space, lambda c: 0.0, hist), [MMQ]), {})

    def test_half_the_leaps_go_by_expected_value(self):
        """The good scores all sit at mapping quality 27: the second leap is the one that lands there."""
        rng = random.Random(1)
        hist = []
        for i in range(60):
            good = i < 12
            p = {MMQ: 27 if good else rng.choice([38, 39, 41, 42, 50, 55]), "min_pruning": 2 + i % 3,
                 "max_num_haplotypes_in_population": 20 + i}                  # 60 distinct configs
            c = self.space._assemble(p)
            if c is not None:
                hist.append((c, 82.0 if good else 70.0 + rng.random()))
        ctx = JumpCtx(self.space, lambda c: 0.0, hist)
        centre = self.space._assemble({MMQ: 40, "min_pruning": 5})          # no leap lands on a scored config
        out = leaps.order(ctx, centre, self.names, lambda c: False)
        nearest = leaps.candidates(ctx, centre, self.names, lambda c: False)
        self.assertEqual(out[0].key, nearest[0][2].key)
        self.assertEqual(out[1].params[MMQ], 27)
        self.assertEqual(len(out), len(nearest))


class ClimberJumpTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()                 # asyncio.run in earlier tests leaves no current loop
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()

    def certify(self, f, centre):
        ctx = JumpCtx(self.space, f)
        ctx.done[centre.key] = (centre, f(centre))
        cl = Climber(ctx, 0, centre, f(centre))
        cl.climb_until = 1e9
        ctx.climbers.append(cl)
        for _ in range(20):
            if cl.state != "climbing":
                break
            if not cl.batch_live and not cl.begin_jumps(1e9):
                break
            while cl.pending:
                cl.submit_one()
            self.assertTrue(cl.gate_done())
            cl.close_wave()
        return cl

    def test_a_centre_jump_that_gains_moves_the_climber_before_any_peak(self):
        """17:06 run: 58 and 60 lose, 46 wins. Every one-rung neighbour loses here; 9-21 rungs down on MMQ gains."""
        centre = self.space._assemble({MMQ: 59})

        def f(c):
            return 82.6953 if 38 <= c.params.get(MMQ) <= 50 else (82.6168 if c.key == centre.key else 82.60)

        cl = self.certify(f, centre)
        self.assertTrue(38 <= cl.centre.params.get(MMQ) <= 50)
        self.assertAlmostEqual(cl.score, 82.6953)

    def test_the_peak_is_named_once_every_centre_jump_has_lost(self):
        centre = self.space._assemble({MMQ: 59})
        with unittest.mock.patch.dict(TUNING.values, {"climb.after_peak": 0}):       # on, it keeps jumping
            cl = self.certify(lambda c: 82.6168 if c.key == centre.key else 82.60, centre)
        self.assertEqual(cl.state, "peak")
        self.assertEqual(cl.centre.key, centre.key)
        down = sorted(59 - c.params.get(MMQ) for c, _ in cl.ctx.done.values()
                      if [n for n in self.space.names if c.params.get(n) != centre.params.get(n)] == [MMQ]
                      and c.params.get(MMQ) < 59)
        classes = {leaps.size_class(k) for k in down if k > 1}
        self.assertEqual(classes, {0, 1, 2, 3})                         # one jump in every distance class
        two = [c for c, _ in cl.ctx.done.values()
               if len([n for n in self.space.names if c.params.get(n) != centre.params.get(n)]) == 2]
        self.assertTrue(two)                                             # then two settings at once

    def test_no_centre_jump_when_the_share_is_zero(self):
        centre = self.space._assemble({MMQ: 59})

        def f(c):
            return 82.6953 if c.params.get(MMQ) == 46 else (82.6168 if c.key == centre.key else 82.60)

        with unittest.mock.patch.dict(TUNING.values, {"climb.leap_share": 0.0}):
            cl = self.certify(f, centre)
        self.assertEqual(cl.state, "peak")
        self.assertEqual(cl.centre.key, centre.key)

    def test_island_climbers_are_numbered_apart(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.climber_base = 100
        cl = Climber(ctx, 2, self.space._assemble({}), 80.0)
        self.assertEqual(cl.id, 102)


class BanditTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.names = [n for n in FULL_RUN_PARAMS + ASSEMBLY_PARAMS + ACTIVITY_PARAMS if n in self.space.names]

    def test_the_priors_favour_far_jumps_on_rugged_settings(self):
        arms = leaps.Arms()
        self.assertGreater(arms.mean((MMQ, 2)), arms.mean((MMQ, 0)))
        self.assertGreater(arms.mean((MMQ, 2)), arms.mean(("min_pruning", 2)))

    def test_an_arm_learns_from_its_own_jumps(self):
        arms = leaps.Arms()
        before = arms.mean(("min_pruning", 1))
        for _ in range(5):
            arms.learn(("min_pruning", 1), 0.05, 0.005)
        self.assertGreater(arms.mean(("min_pruning", 1)), 3 * before)
        for _ in range(30):
            arms.learn((MMQ, 2), 0.0, 0.005)
        self.assertLess(arms.mean((MMQ, 2)), leaps.Arms().mean((MMQ, 2)) / 2)

    def test_picks_go_mostly_to_rugged_settings_beyond_three_rungs(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.rng = random.Random(3)
        centre = self.space._assemble({MMQ: 40, "max_num_haplotypes_in_population": 60})
        got = leaps.pick(ctx, centre, self.names, lambda c: False, 40, ctx.rng)
        self.assertEqual(len(got), 40)
        self.assertEqual(len({c.key for _, _, _, c in got}), 40)
        rugged = [g for g in got if g[0] in leaps.RUGGED]
        share = len([n for n in self.names if n in leaps.RUGGED]) / len(self.names)
        self.assertGreater(len(rugged), 1.5 * share * len(got))         # well over their share of the settings
        self.assertLess(sum(1 for g in got if g[1] == 0), 10)
        for name, ci, k, c in got:
            self.assertEqual(leaps.size_class(k), ci)
            self.assertEqual([n for n in self.space.names if c.params.get(n) != centre.params.get(n)], [name])

    def test_the_certificate_takes_one_jump_per_class_and_direction(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        centre = self.space._assemble({MMQ: 40})
        got = leaps.coverage(ctx, centre, self.names, lambda c: False, random.Random(0))
        cells = [(name, ci, k > 0) for name, ci, k, _ in got]
        self.assertEqual(len(cells), len(set(cells)))
        self.assertIn((MMQ, 3, False), cells)
        self.assertIn((MMQ, 2, True), cells)

    def test_two_setting_jumps_move_exactly_two_settings(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        centre = self.space._assemble({MMQ: 40, "max_num_haplotypes_in_population": 60})
        got = leaps.pairs(ctx, centre, self.names, lambda c: False, random.Random(0))
        self.assertTrue(got)
        for name, ci, k, c in got:
            self.assertEqual(len([n for n in self.space.names if c.params.get(n) != centre.params.get(n)]), 2)

    def test_each_direction_learns_on_its_own(self):
        """Mapping quality won 7 times in two runs, always down."""
        arms = leaps.Arms()
        down, up = (MMQ, 2, -1), (MMQ, 2, 1)
        start = arms.mean(up)
        for _ in range(4):
            arms.learn(down, 0.03, 0.005)
            arms.learn(up, 0.0, 0.005)
        self.assertGreater(arms.mean(down), 2 * arms.mean(up))
        self.assertLess(arms.mean(up), start)

    def test_picks_follow_the_direction_that_won(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.rng = random.Random(5)
        arms = leaps.arms_of(ctx)
        for ci in range(len(leaps.CLASSES)):
            for _ in range(6):
                arms.learn((MMQ, ci, -1), 0.04, 0.005)
                arms.learn((MMQ, ci, 1), 0.0, 0.005)
        centre = self.space._assemble({MMQ: 40})
        got = [g for g in leaps.pick(ctx, centre, self.names, lambda c: False, 30, ctx.rng) if g[0] == MMQ]
        self.assertTrue(got)
        self.assertGreater(sum(1 for g in got if g[2] < 0), 3 * max(1, sum(1 for g in got if g[2] > 0)))

    def test_a_winning_jump_is_tried_again_then_twice_as_far(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        centre = self.space._assemble({MMQ: 40})
        got = leaps.momentum(ctx, centre, MMQ, -1, 3, lambda c: False)
        self.assertEqual([(g[0], g[2]) for g in got], [(MMQ, -3), (MMQ, -6)])
        self.assertEqual([g[3].params[MMQ] for g in got], [37, 34])

    def test_a_climber_sends_its_winning_move_first_in_the_next_wave(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.rng = random.Random(1)
        old = self.space._assemble({MMQ: 43})
        new = self.space._assemble({MMQ: 40})
        cl = Climber(ctx, 0, new, 80.0)
        ctx.climbers.append(cl)
        cl.note_momentum(old, new)
        self.assertEqual(cl.momentum, (MMQ, -1, 3))
        planned = cl.plan_leaps(new, 10)
        self.assertEqual([c.params[MMQ] for c in planned[:2]], [37, 34])
        self.assertIsNone(cl.momentum)

    def test_the_islands_share_their_jump_arms(self):
        arms = leaps.Arms()
        seats = {}
        a, b = islands.Island(0, 2, seats, arms), islands.Island(1, 2, seats, arms)
        self.assertIs(a.arms, b.arms)

    def test_three_jump_wins_double_the_share_and_two_misses_halve_it(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        cl = Climber(ctx, 0, self.space._assemble({}), 80.0)
        start = cl.leap_share
        for _ in range(3):
            cl._wave_leaps, cl._wave_leap_wins = 5, 1
            cl.adapt_leap_share()
        from search.climb import LEAP_SHARE_MAX
        self.assertAlmostEqual(cl.leap_share, min(LEAP_SHARE_MAX, 2 * start))
        for _ in range(2):
            cl._wave_leaps, cl._wave_leap_wins = 5, 0
            cl.adapt_leap_share()
        self.assertAlmostEqual(cl.leap_share, max(0.1, min(LEAP_SHARE_MAX, 2 * start) / 2))


class FakeBase:
    """The shared client: a config is pending until finish(); a scored one comes back at once."""

    def __init__(self, idle=0):
        self.round = ROUND
        self.ranges = self.store = None
        self.executor = object()
        self.budget = self.deadline_ts = None
        self.answer_set = None
        self.clock = lambda: self.now
        self.now = 0.0
        self.spent = 0.0
        self.latest, self.callsets, self.configs = {}, {}, {}
        self.on_result, self.on_rescore = [], []
        self.notes, self.pick_key = {}, None
        self.pending = {}
        self.idle = idle

    exhausted = False

    def time_left(self):
        return 1e9

    def idle_slots(self):
        return self.idle

    def propose_many(self, configs, *, fidelity=FULL, priority=0, tag="", batch="big", batches=None):
        loop = asyncio.get_event_loop()
        out = []
        for c in configs:
            fut = loop.create_future()
            hit = self.latest.get((c.key, fidelity))
            if hit is not None:
                fut.set_result(type("O", (), {"config": c, "result": hit})())
            else:
                self.configs[c.key] = c
                self.pending.setdefault(c.key, []).append(fut)
                self.spent += 1.0
            out.append(fut)
        return out

    def finish(self, key, score):
        c = self.configs[key]
        r = ScoreResult.__new__(ScoreResult)
        r.score, r.gate_passed, r.fidelity = score, True, FULL
        self.latest[(key, FULL)] = r
        outcome = type("O", (), {"config": c, "result": r})()
        for fut in self.pending.pop(key, []):
            fut.set_result(outcome)
        for cb in self.on_result:
            cb(outcome)


class IslandTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()

    def views(self, base, cap=810, k=2):
        out = []
        for i, q in enumerate(islands.quotas(cap, k)):
            out.append(islands.IslandClient(base, i, q, out))
        return out

    def config(self, mmq):
        return self.space._assemble({MMQ: mmq})

    def test_one_island_per_324_workers(self):
        self.assertEqual([islands.count(n) for n in (100, 324, 647, 648, 810, 1056)], [1, 1, 1, 2, 2, 3])
        self.assertEqual(islands.quotas(810, 2), [405, 405])
        self.assertEqual(islands.quotas(1056, 3), [352, 352, 352])

    def test_an_island_sees_only_the_configs_it_proposed(self):
        base = FakeBase()
        a, b = self.views(base)
        x, y = self.config(40), self.config(45)
        a.propose_many([x])
        b.propose_many([y])
        base.finish(x.key, 81.0)
        base.finish(y.key, 82.0)
        self.assertEqual(list(a.latest), [(x.key, FULL)])
        self.assertEqual(list(b.latest), [(y.key, FULL)])
        self.assertEqual(set(a.configs), {x.key})
        heard = []
        b.on_result.append(lambda o: heard.append(o.config.key))
        fut = b.propose_many([x])[0]                     # another island's score comes back at once
        self.assertTrue(fut.done())
        self.assertIn((x.key, FULL), b.latest)
        self.assertEqual(heard, [x.key])
        self.assertEqual(base.spent, 2.0)

    def test_a_config_still_running_elsewhere_is_joined(self):
        base = FakeBase()
        a, b = self.views(base)
        x = self.config(40)
        a.propose_many([x])
        fb = b.propose_many([x])[0]
        self.assertFalse(fb.done())
        self.assertEqual(b.in_flight, 1)
        base.finish(x.key, 81.0)
        self.assertTrue(fb.done())
        self.assertIn((x.key, FULL), b.latest)
        self.assertEqual(b.in_flight, 0)

    def test_idle_workers_leave_room_for_an_island_under_its_share(self):
        base = FakeBase()
        a, b = self.views(base)
        a.propose_many([self.config(1 + (i % 60)) for i in range(405)])
        b.propose_many([self.space._assemble({MMQ: 1 + (i % 60), "min_pruning": 2 + i // 60}) for i in range(300)])
        base.idle = 105
        self.assertEqual(a.idle_slots(), 0)
        self.assertEqual(b.idle_slots(), 105)

    def test_a_quiet_island_leaves_its_share_to_the_others(self):
        base = FakeBase()
        a, b = self.views(base)
        a.propose_many([self.config(1 + (i % 60)) for i in range(405)])
        base.idle = 405
        self.assertEqual(a.idle_slots(), 0)
        base.now = islands.ACTIVE_S + 1
        a.last_propose = base.now
        self.assertEqual(a.idle_slots(), 405)

    def test_the_pick_is_the_best_island_pick(self):
        base = FakeBase()
        a, b = self.views(base)
        x, y = self.config(40), self.config(45)
        a.propose_many([x])
        b.propose_many([y])
        base.finish(x.key, 81.0)
        base.finish(y.key, 82.0)
        a.pick_key = x.key
        self.assertEqual(base.pick_key, y.key)                       # island b's best counts before it picks
        b.pick_key = y.key
        self.assertEqual(base.pick_key, y.key)

    def test_the_notes_show_every_island(self):
        base = FakeBase()
        a, b = self.views(base)
        a.notes["engine"] = {"best": 82.1, "phase": "B", "basins": [{"id": 0, "best": 82.1}], "events": [[2.0, "a"]]}
        b.notes["engine"] = {"best": 82.3, "phase": "C", "basins": [{"id": 100, "best": 82.3}], "events": [[1.0, "b"]]}
        islands.merge_notes(base, [a, b], 810)
        eng = base.notes["engine"]
        self.assertEqual(eng["phase"], "C")
        self.assertEqual(eng["capacity"], 810)
        self.assertEqual([r["island"] for r in eng["islands"]], [0, 1])
        self.assertEqual([(r["id"], r["island"]) for r in eng["basins"]], [(0, 0), (100, 1)])
        self.assertEqual([e[1] for e in eng["events"]], ["b", "a"])

    def test_a_seat_close_to_another_islands_seat_is_passed_over(self):
        seats = {}
        first, second = islands.Island(0, 2, seats), islands.Island(1, 2, seats)

        class _Ctx:
            space = self.space

            def say(self, *a):
                pass

        ctx = _Ctx()
        s0 = self.space._assemble({MMQ: 40})
        near = self.space._assemble({MMQ: 41, "min_pruning": 3})
        far = self.space._assemble({MMQ: 20, "min_pruning": 5, "min_dangling_branch_length": 10,
                                    "assembly_region_padding": 300, "heterozygosity": 0.005,
                                    "max_num_haplotypes_in_population": 50, "min_assembly_region_size": 30,
                                    "max_assembly_region_size": 200, "active_probability_threshold": 0.01})
        self.assertEqual(first.apart(ctx, [s0]), [s0])
        self.assertEqual(second.apart(ctx, [far, near]), [far])          # a later start close to island 0 goes
        self.assertEqual(second.apart(ctx, [near, far]), [near, far])    # the island's best start always stays

    def test_portfolio_seeds_are_dealt_by_rank(self):
        ranked = [self.space._assemble({MMQ: 10 + i, "min_pruning": 5}) for i in range(12)]

        def seeds_for(index):
            class _Ctx:
                pass

            ctx = _Ctx()
            ctx.space, ctx.cap, ctx.profile = self.space, 10, None
            ctx.say = lambda *a: None
            ctx.island = islands.Island(index, 2, {})
            ctx.client = type("Client", (), {"round": type("Rnd", (), {"tag": "20261012T115400Z"})()})()
            with unittest.mock.patch("search.candidates.portfolio.pick", lambda space, profile, n, **kw: ranked[:n]):
                return StepA(ctx).seeds()

        zero, one = seeds_for(0), seeds_for(1)
        official = self.space._assemble(dict(OFFICIAL_SEED)).key
        self.assertIn(official, {c.key for c in zero})
        self.assertNotIn(official, {c.key for c in one})
        self.assertEqual([c.params[MMQ] for c in one], [11, 13, 15, 17])
        self.assertEqual([c.params[MMQ] for c in zero if c.key in {r.key for r in ranked}], [10, 12, 14, 16])


class IslandRoundTest(unittest.TestCase):
    def test_two_islands_climb_and_the_best_island_is_picked(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 600, duration=lambda job: 150.0)
        seen = []

        async def wrapped(client):
            seen.append(client)
            await engine.searcher(r, 600)(client)

        with unittest.mock.patch.dict(TUNING.values, {"islands.workers": 300}):
            out = run_round(ROUND, Stub(), ex, FakeScorer(runner), wrapped, DryRunSubmitter(), r,
                            LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                            deadline_ts=72 * 60)
        client = seen[0]
        eng = client.notes["engine"]
        self.assertEqual([row["island"] for row in eng["islands"]], [0, 1])
        self.assertEqual({row["island"] for row in eng["basins"]}, {0, 1})
        self.assertTrue(all(row["id"] >= 100 for row in eng["basins"] if row["island"] == 1))
        self.assertIsNotNone(client.pick_key)
        best_island = max((row["best"] for row in eng["islands"] if row["best"] is not None))
        round_best = max(r.score for r in client.latest.values() if r.score is not None and r.gate_passed is not False)
        picked = client.latest[(client.pick_key, FULL)].score
        self.assertAlmostEqual(picked, round_best, places=6)            # the pick is the highest local score
        self.assertGreaterEqual(picked, best_island - 1e-4)                  # the notes round to 4 decimals


if __name__ == "__main__":
    unittest.main()


class CrashCtx(JumpCtx):
    """JumpCtx whose runs return Outcomes; f(config) None is a GATK crash (no score)."""

    def __init__(self, space, f):
        super().__init__(space, f)
        self.failed = {}

    def note_failure(self, source, reason):
        self.failed[source] = self.failed.get(source, 0) + 1

    def fresh(self, c):
        return c is not None and c.key not in self.done and c.key not in getattr(self, "sent_keys", set())

    def propose(self, batch, tag, batches, priority=0):
        from search.client import Outcome
        loop = asyncio.get_event_loop()
        self.sent_keys = getattr(self, "sent_keys", set())
        out = []
        for c in batch:
            self.sent_keys.add(c.key)
            s = self.f(c)
            fut = loop.create_future()
            if s is None:
                fut.set_result(Outcome(c, None, "rc=3 java.lang.IllegalStateException: Graph must have ref source"))
            else:
                self.done[c.key] = (c, s)
                r = ScoreResult.__new__(ScoreResult)
                r.score, r.gate_passed, r.fidelity = s, True, FULL
                fut.set_result(Outcome(c, r))
            out.append(fut)
        return out


class CrashTest(unittest.TestCase):
    """044200Z (2026-10-16): crashes doubled with jumps, because a crashed jump never reached the arms."""

    def setUp(self):
        self.space = Space(ranges())
        self.centre = self.space._assemble({MMQ: 40, "assembly_region_padding": 50})
        self.pad = "assembly_region_padding"
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()

    def climber(self, f):
        ctx = CrashCtx(self.space, f)
        ctx.rng = random.Random(2)
        ctx.done[self.centre.key] = (self.centre, 80.0)
        cl = Climber(ctx, 0, self.centre, 80.0)
        cl.climb_until = 1e9
        ctx.climbers.append(cl)
        return ctx, cl

    def run_certificate(self, ctx, cl):
        for _ in range(10):
            if cl.state != "climbing" or (not cl.batch_live and not cl.begin_jumps(1e9)):
                break
            while cl.pending:
                cl.submit_one()
            cl.close_wave()

    def test_a_crashed_jump_is_a_loss_and_leaves_the_waiting_list(self):
        crash = lambda c: None if c.params.get(self.pad) < 45 else 79.0
        ctx, cl = self.climber(crash)
        self.run_certificate(ctx, cl)
        down = [(a, st) for a, st in leaps.arms_of(ctx).stats.items() if a[0] == self.pad and a[2] == -1]
        self.assertTrue(down)
        self.assertTrue(all(st[1] >= 1 and st[0] == 0 for _, st in down))     # tried, never won
        self.assertFalse([k for k, c in cl.waiting.items() if c.params.get(self.pad) < 45])
        self.assertGreater(ctx.failed.get("jump", 0), 0)
        self.assertIn((self.pad, -1), cl.crash_edge)

    def test_no_jump_goes_past_a_crash_from_the_same_centre(self):
        crash = lambda c: None if c.params.get(self.pad) < 45 else 79.0
        ctx, cl = self.climber(crash)
        self.run_certificate(ctx, cl)
        edge = cl.crash_edge[(self.pad, -1)]
        lower = self.space._assemble({**self.centre.params, self.pad: 20})
        self.assertTrue(cl.beyond_crash(self.centre, lower))
        upper = self.space._assemble({**self.centre.params, self.pad: 60})
        self.assertFalse(cl.beyond_crash(self.centre, upper))
        far_centre = self.space._assemble({**self.centre.params, self.pad: 120})
        self.assertFalse(cl.beyond_crash(far_centre, self.space._assemble({**far_centre.params, self.pad: 100})))
        self.assertLess(edge[1], edge[0])

    def test_the_context_counts_failed_runs_by_source(self):
        from search.client import Outcome
        base = FakeBase()
        ctx = engine.Context(base, self.space, random.Random(0), None)
        ctx._learn_score(Outcome(self.centre, None, "rc=3 java.lang.IllegalStateException: Graph | batch crashed", "G"))
        ctx._learn_score(Outcome(self.centre, None, "rc=3 java.lang.IllegalStateException: Graph", "A"))
        ctx._learn_score(Outcome(self.centre, None, "rc=3 x", "B0"))                # climbers count their own
        self.assertEqual(ctx.failures, {"GA": 1, "step A": 1})
        self.assertEqual(list(ctx.fail_reasons.values()), [2])


class StepCFillTest(unittest.TestCase):
    """044200Z (2026-10-16): Step C kept 6-34% of Step B's cores busy; the fine-tune's idle workers now take jumps."""

    def setUp(self):
        self.space = Space(ranges())

    def test_idle_workers_take_jumps_then_multi_setting_moves(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.rng = random.Random(4)
        ctx.pick_at = lambda: 1e9
        centre = self.space._assemble({MMQ: 40, "max_num_haplotypes_in_population": 60})
        cl = Climber(ctx, 0, centre, 80.0, origin="step C")
        ctx.climbers.append(cl)
        n = cl.fill_step_c(3000)                                         # more than single jumps can cover
        self.assertGreater(n, 100)
        singles = [c for c in cl.leaps if len([x for x in self.space.names if c.params.get(x) != centre.params.get(x)]) == 1]
        self.assertEqual(cl.leaps[:len(singles)], singles)              # single jumps first
        self.assertEqual(len({c.key for c in cl.leaps}), len(cl.leaps))
        multi = [c for c in cl.leaps if len([x for x in self.space.names if c.params.get(x) != centre.params.get(x)]) >= 2]
        self.assertTrue(multi)
        self.assertTrue(all(c.key in cl.leap_meta for c in cl.leaps))

    def test_a_blocked_fill_is_bounded_and_not_retried(self):
        # 2026-10-22 03:09 live: the crash guard blocked nearly every candidate; each fill made 20 tries per idle worker
        # (~11,000) and was called again at once - the event loop starved for the rest of the round
        from search import leaps
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.rng = random.Random(4)
        ctx.pick_at = lambda: 1e9
        calls = [0]

        def blocked(c):
            calls[0] += 1
            return False
        ctx.fresh = blocked
        cl = Climber(ctx, 0, self.space._assemble({MMQ: 40}), 80.0, origin="step C")
        ctx.climbers.append(cl)
        self.assertEqual(cl.fill_step_c(600), 0)
        # single jumps (bounded by the space, ~2,000 here) + at most MUTATION_TRIES mutation tries; before: + 20 x 600
        self.assertLessEqual(calls[0], 2000 + leaps.MUTATION_TRIES)
        before = calls[0]
        self.assertEqual(cl.fill_step_c(600), 0)                        # nothing new: not tried again
        self.assertEqual(calls[0], before)

    def test_no_fill_once_a_run_could_not_return_before_the_pick(self):
        from search.climb import FILL_CUTOFF_S
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.now = 1000.0
        ctx.pick_at = lambda: 1000.0 + FILL_CUTOFF_S - 1
        cl = Climber(ctx, 0, self.space._assemble({MMQ: 40}), 80.0, origin="step C")
        self.assertEqual(cl.fill_step_c(100), 0)


class WaveSpeedTest(unittest.TestCase):
    """044200Z islands (2026-10-16): grid waves took 4.0 minutes (2.8 before islands); heads queued behind GA runs."""

    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()

    def climber_with_wave(self, back, heads=10):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        cl = Climber(ctx, 0, self.space._assemble({MMQ: 40}), 80.0)
        cl.batch_live = True
        futs = []
        for i in range(heads):
            f = self.loop.create_future()
            if i < back:
                f.set_result(None)
            futs.append(f)
        cl.gate_futs = futs
        return cl

    def test_the_ga_leaves_room_for_the_next_heads_once_a_wave_is_half_back(self):
        from search.climb import head_reserve
        self.assertEqual(head_reserve([self.climber_with_wave(back=2)]), 0)
        self.assertEqual(head_reserve([self.climber_with_wave(back=5)]), 10)
        self.assertEqual(head_reserve([self.climber_with_wave(back=40, heads=60)]), 24)
        with unittest.mock.patch.dict(TUNING.values, {"climb.head_reserve": 0}):
            self.assertEqual(head_reserve([self.climber_with_wave(back=5)]), 0)

    def test_the_ga_opening_ends_when_a_child_beats_the_lowest_seat(self):
        from search.climb import ga_beat_lowest
        ctx = JumpCtx(self.space, lambda c: 80.0)
        low = Climber(ctx, 1, self.space._assemble({MMQ: 30}), 82.05, origin="far")
        high = Climber(ctx, 0, self.space._assemble({MMQ: 40}), 82.22)
        child = self.space._assemble({MMQ: 35})
        ctx.done[child.key] = (child, 82.10)
        ctx.ga_keys = {child.key}
        ctx.best = lambda: (high.centre, 82.22)
        self.assertTrue(ga_beat_lowest(ctx, [low, high]))
        ctx.done[child.key] = (child, 82.00)
        self.assertFalse(ga_beat_lowest(ctx, [low, high]))


class RelinkTest(unittest.TestCase):
    """044200Z islands (2026-10-16): the better peak was eight settings from a peak no single move could leave."""

    def setUp(self):
        self.space = Space(ranges())
        self.peak = self.space._assemble({MMQ: 49, "max_assembly_region_size": 102, "max_num_haplotypes_in_population": 8,
                                          "min_pruning": 2, "phred_scaled_global_read_mismapping_rate": 47})
        self.better = self.space._assemble({MMQ: 47, "max_assembly_region_size": 105, "max_num_haplotypes_in_population": 54,
                                            "min_pruning": 1, "phred_scaled_global_read_mismapping_rate": 36})

    def moved(self, c):
        return [n for n in self.space.names if c.params.get(n) != self.peak.params.get(n)]

    @unittest.skip("built on validator-edge values outside the published ranges, removed 2026-10-18 (search/space.py VALIDATOR_EDGES)")

    def test_subsets_of_the_differences_toward_a_better_config(self):
        ctx = JumpCtx(self.space, lambda c: 80.0, history=[(self.better, 82.6664)])
        got = leaps.relink(ctx, self.peak, 82.6298, lambda c: False, random.Random(0))
        diff = set(self.moved(self.better))
        self.assertTrue(got)
        sizes = collections.Counter(len(self.moved(c)) for _, _, _, c in got)
        self.assertEqual(set(sizes), {2, 3, 4})
        for _, _, _, c in got:
            self.assertTrue(set(self.moved(c)) <= diff)
            self.assertTrue(all(c.params[n] == self.better.params[n] for n in self.moved(c)))

    @unittest.skip("built on validator-edge values outside the published ranges, removed 2026-10-18 (search/space.py VALIDATOR_EDGES)")

    def test_a_target_far_below_or_too_close_is_not_used(self):
        near = self.space._assemble({**self.peak.params, MMQ: 48})
        ctx = JumpCtx(self.space, lambda c: 80.0, history=[(self.better, 82.40), (near, 82.70)])
        self.assertEqual(leaps.relink(ctx, self.peak, 82.6298, lambda c: False, random.Random(0)), [])

    @unittest.skip("built on validator-edge values outside the published ranges, removed 2026-10-18 (search/space.py VALIDATOR_EDGES)")

    def test_the_other_islands_best_is_a_target(self):
        ctx = JumpCtx(self.space, lambda c: 80.0)
        got = leaps.relink(ctx, self.peak, 82.6298, lambda c: False, random.Random(0),
                           elsewhere=[(self.better, 82.6664)])
        self.assertTrue(got)


class FollowUpTest(unittest.TestCase):
    """The 2026-10-16 follow-ups: late seats, mid-wave claims, risky heads, the adoption switch, carried arms."""

    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()

    def test_no_far_seat_in_the_last_twenty_minutes_of_step_b(self):
        from search.climb import LATE_SEAT_MIN, late_for_seat
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.now = 1000.0
        self.assertTrue(late_for_seat(ctx, 1000.0 + 60 * LATE_SEAT_MIN - 1))
        self.assertFalse(late_for_seat(ctx, 1000.0 + 60 * LATE_SEAT_MIN + 1))

    def test_a_run_from_an_earlier_wave_takes_the_centre_during_a_wave(self):
        centre = self.space._assemble({MMQ: 40})
        better = self.space._assemble({MMQ: 41, "min_pruning": 3})
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.done[centre.key] = (centre, 80.0)
        cl = Climber(ctx, 0, centre, 80.0)
        cl.climb_until = 1e9
        cl._watch(better)                                   # sent in wave 0
        cl._wave_id = 1
        cl.batch_live = True                                # wave 1 is live
        same_wave = self.space._assemble({MMQ: 42})
        cl._watch(same_wave)
        ctx.done[same_wave.key] = (same_wave, 81.0)
        ctx.done[better.key] = (better, 80.5)
        self.assertTrue(cl.claim_higher())
        self.assertEqual(cl.centre.key, better.key)         # the same-wave 81.0 waits for its wave to close

    def test_a_head_next_to_a_crash_goes_without_riders(self):
        centre = self.space._assemble({MMQ: 40, "assembly_region_padding": 50})
        cl = Climber(JumpCtx(self.space, lambda c: 80.0), 0, centre, 80.0)
        lad = cl.ctx.ladders
        i = lad.index("assembly_region_padding", 50)
        cl.crash_edge[("assembly_region_padding", -1)] = (i, i - 2)
        head = self.space._assemble({**centre.params, "assembly_region_padding": lad.shift("assembly_region_padding", 50, -1)})
        self.assertFalse(cl.beyond_crash(centre, head))
        self.assertTrue(cl.beyond_crash(centre, head, margin=1))

    def test_the_own_hill_adoption_can_be_switched_off(self):
        from search.climb import adopt_ga_best
        ctx = JumpCtx(self.space, lambda c: 80.0)
        low = Climber(ctx, 0, self.space._assemble({"min_pruning": 2, "assembly_region_padding": 50}), 10.0)
        high_start = self.space._assemble({"min_pruning": 6, "assembly_region_padding": 55})
        high = Climber(ctx, 1, high_start, 90.0)
        ga = self.space._assemble({**high_start.params, "min_pruning": 7})
        ctx.climbers = [low, high]
        ctx.ga_keys = {ga.key}
        ctx.done[ga.key] = (ga, 100.0)
        ctx.best = lambda: (ga, 100.0)
        with unittest.mock.patch.dict(TUNING.values, {"climb.own_hill_adopt": 0}):
            self.assertFalse(adopt_ga_best(ctx, [low, high]))
        self.assertEqual(high.centre.key, high_start.key)
        self.assertTrue(adopt_ga_best(ctx, [low, high]))
        self.assertEqual(high.centre.key, ga.key)

    def test_arms_carry_to_the_next_round_halved_and_capped(self):
        import tempfile, os
        arms = leaps.Arms()
        for _ in range(4):
            arms.learn((MMQ, 2, -1), 0.03, 0.005)
        for _ in range(60):
            arms.learn((MMQ, 2, 1), 0.0, 0.005)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "search", "jump_arms.json")
            leaps.save_arms(arms, path)
            nxt = leaps.load_arms(path)
            self.assertEqual(nxt.carried[(MMQ, 2, -1)][:2], [4.0, 4.0])
            self.assertLessEqual(nxt.carried[(MMQ, 2, 1)][1], leaps.CARRY_MAX_TRIES)
            self.assertGreater(nxt.mean((MMQ, 2, -1)), leaps.Arms().mean((MMQ, 2, -1)))
            leaps.save_arms(nxt, path)                      # nothing new learnt: the carried counts halve
            self.assertEqual(leaps.load_arms(path).carried[(MMQ, 2, -1)][:2], [2.0, 2.0])
        self.assertEqual(leaps.load_arms("/nonexistent/arms.json").carried, {})


class GuardAndStallTest(unittest.TestCase):
    """The v4 fleet follow-ups (2026-10-16): crash guard, stall after two flat waves, one take per wave, the logs."""

    PAD = "assembly_region_padding"

    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)

    def tearDown(self):
        self.loop.close()

    def cfg(self, pad, mmq=40):
        return self.space._assemble({MMQ: mmq, self.PAD: pad})

    def taught(self, crashes=60, ok=300):
        from search.crashguard import CrashGuard
        g = CrashGuard(self.space, Ladders(self.space))
        rng = random.Random(1)
        lad = g.ladders.rungs[self.PAD]
        low = [v for v in lad if v < 25] or lad[:2]
        high = [v for v in lad if v >= 40]
        for i in range(crashes):
            g.note(self.cfg(rng.choice(low), 30 + i % 20).params, True)
        for i in range(ok):
            g.note(self.cfg(rng.choice(high), 30 + i % 20).params, False)
        return g, low, high

    def test_the_guard_learns_the_crash_region_and_says_nothing_before_enough_crashes(self):
        from search.crashguard import CrashGuard, MIN_CRASHES
        g, low, high = self.taught(crashes=MIN_CRASHES - 1)
        self.assertIsNone(g.model)
        self.assertFalse(g.blocks(self.cfg(low[0])))
        g, low, high = self.taught()
        self.assertIsNotNone(g.model)
        self.assertTrue(g.blocks(self.cfg(low[0], 33)))
        self.assertTrue(g.risky(self.cfg(low[0], 33)))
        self.assertFalse(g.blocks(self.cfg(high[0], 33)))
        self.assertFalse(g.risky(self.cfg(high[0], 33)))
        self.assertEqual(len(g.skipped), 1)
        g.blocks(self.cfg(low[0], 33))
        self.assertEqual(len(g.skipped), 1)                          # a config is counted once
        self.assertIn("1 configs not sent", g.summary())

    def test_audits_send_some_withheld_configs_and_a_wrong_rule_stops_blocking(self):
        from search import crashguard as cg
        g, low, high = self.taught()
        lad = g.ladders.rungs[self.PAD]
        risky = [self.cfg(v, m) for v in low for m in range(10, 60)]
        sent = [c for c in risky if not g.blocks(c)]
        self.assertTrue(sent)                                        # about a tenth are audited
        self.assertLess(len(sent), len(risky) * 0.3)
        self.assertGreater(g.n_audits, 0)
        for c in sent:                                               # the audited configs all come back fine
            g.note(c.params, False, c.key)
        self.assertGreaterEqual(g.n_off, 1)
        withheld = [c for c in risky if c.key in g.skipped]
        free = [c for c in withheld if g.p_crash(c) < cg.SKIP_P or g._leaf(g.features(c.params)) in g.off]
        self.assertTrue(free)               # the rule is off, or the next fit learnt the audits: no longer blocked
        self.assertFalse(any(g.blocks(c) for c in free))
        self.assertIn("rules stopped blocking", g.summary())

    def test_audited_crashes_keep_the_rule(self):
        g, low, high = self.taught()
        risky = [self.cfg(v, m) for v in low for m in range(10, 60)]
        sent = [c for c in risky if not g.blocks(c)]
        for c in sent:
            g.note(c.params, True, c.key)
        self.assertEqual(g.n_off, 0)

    def test_the_tree_walk_matches_predict_proba(self):
        g, low, high = self.taught()
        rng = random.Random(3)
        lad = g.ladders.rungs[self.PAD]
        for _ in range(50):
            c = self.cfg(rng.choice(lad), rng.choice([20, 30, 40, 50]))
            x = g.features(c.params)
            self.assertAlmostEqual(g._leaf_rate(x), float(g.model.predict_proba([x])[0][1]), places=9)

    def test_fresh_refuses_a_config_the_guard_predicts_to_crash(self):
        from search.client import Outcome
        ctx = engine.Context(FakeBase(), self.space, random.Random(0), None)
        g, low, high = self.taught()
        ctx.guard = g
        self.assertFalse(ctx.fresh(self.cfg(low[0], 33)))
        self.assertTrue(ctx.fresh(self.cfg(high[0], 33)))
        fresh_ctx = engine.Context(FakeBase(), self.space, random.Random(0), None)
        r = ScoreResult.__new__(ScoreResult)
        r.score, r.gate_passed, r.fidelity = 80.0, True, FULL
        fresh_ctx._learn_score(Outcome(self.cfg(low[0]), None, "rc=3 java.lang.IllegalStateException: Graph", "G"))
        fresh_ctx._learn_score(Outcome(self.cfg(low[0]), None, "run cancelled: round over", "G"))
        fresh_ctx._learn_score(Outcome(self.cfg(high[0]), r, None, "G"))
        self.assertEqual((len(fresh_ctx.guard.crashed), len(fresh_ctx.guard.ok)), (1, 1))   # a cancel is not a crash

    def test_a_head_the_guard_finds_risky_goes_without_riders(self):
        g, low, high = self.taught()
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.guard = g
        centre = self.cfg(high[0])
        cl = Climber(ctx, 0, centre, 80.0)
        sent = []
        ctx.propose = lambda batch, tag, batches, priority=0: (sent.append(list(batch)),
                                                              [self.loop.create_future() for _ in batch])[1]
        cl.pending = [self.cfg(low[0], 33)]
        cl._submit_head()
        self.assertEqual(len(sent[0]), 1)

    def test_two_flat_waves_start_the_certificate(self):
        cl = Climber(JumpCtx(self.space, lambda c: 80.0), 0, self.cfg(50), 80.0)
        cl.recent = [0.02, 0.0]
        self.assertFalse(cl.flat_waves())
        cl.recent = [0.02, 0.0, 0.001]
        self.assertTrue(cl.flat_waves())
        with unittest.mock.patch.dict(TUNING.values, {"climb.stall_waves": 0}):
            self.assertFalse(cl.flat_waves())

    def climber_with_returns(self):
        centre = self.space._assemble({MMQ: 40})
        ctx = JumpCtx(self.space, lambda c: 80.0)
        said = []
        ctx.say = lambda fmt, *a: said.append(fmt % a)
        ctx.done[centre.key] = (centre, 80.0)
        cl = Climber(ctx, 0, centre, 80.0)
        cl.climb_until = 1e9
        return ctx, cl, said

    def test_returned_runs_from_earlier_waves_can_be_taken_one_after_another(self):
        ctx, cl, said = self.climber_with_returns()
        a = self.space._assemble({MMQ: 41, "min_pruning": 3})
        b = self.space._assemble({MMQ: 39, "min_pruning": 4})
        cl._watch(a)
        cl._watch(b)                                        # both sent in wave 0
        cl._wave_id = 1
        cl.batch_live = True
        ctx.done[a.key] = (a, 80.5)
        self.assertTrue(cl.claim_higher())
        self.assertEqual(cl.centre.key, a.key)
        ctx.done[b.key] = (b, 81.0)
        self.assertTrue(cl.claim_higher())                  # no longer held until the take's wave closes
        self.assertEqual(cl.centre.key, b.key)

    def test_a_returned_centre_jump_is_logged_as_one(self):
        ctx, cl, said = self.climber_with_returns()
        j = self.space._assemble({MMQ: 20})
        cl._watch(j)
        cl.leap_sent.add(j.key)
        ctx.done[j.key] = (j, 80.5)
        self.assertTrue(cl.claim_higher())
        self.assertTrue(any(m.startswith("centre jump c0 (returned)") for m in said), said)

    def test_it_stays_there_is_silent_when_the_holder_is_already_as_high(self):
        from search.climb import adopt_ga_best
        for ga_score, loud in ((50.0, False), (95.0, True)):
            ctx = JumpCtx(self.space, lambda c: 80.0)
            said = []
            ctx.say = lambda fmt, *a: said.append(fmt % a)
            low = Climber(ctx, 0, self.space._assemble({"min_pruning": 2, self.PAD: 50}), 10.0)
            high_start = self.space._assemble({"min_pruning": 6, self.PAD: 55})
            high = Climber(ctx, 1, high_start, 90.0)
            ga = self.space._assemble({**high_start.params, "min_pruning": 7})
            ctx.climbers = [low, high]
            ctx.ga_keys = {ga.key}
            ctx.done[ga.key] = (ga, ga_score)
            ctx.best = lambda: (high_start, 90.0)
            with unittest.mock.patch.dict(TUNING.values, {"climb.own_hill_adopt": 0}):
                adopt_ga_best(ctx, [low, high])
            self.assertEqual(any("it stays there" in m for m in said), loud, (ga_score, said))

    def test_the_head_reserve_is_counted_in_worker_seconds(self):
        from search.climb import note_reserve
        ctx = JumpCtx(self.space, lambda c: 80.0)
        note_reserve(ctx, 10)
        ctx.now = 6.0
        note_reserve(ctx, 0)
        ctx.now = 20.0
        note_reserve(ctx, 4)
        self.assertAlmostEqual(ctx.reserve_ws, 60.0)


class MoveAtOnceTest(unittest.TestCase):
    """2026-10-16: a better result waited 13-41 minutes per run before its climber moved (half for the rest of its
    wave, a third for a cube on the old centre). A run that clears climb.accept_min now moves the centre at once."""

    def setUp(self):
        self.space = Space(ranges())
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        self.lad = Ladders(self.space)

    def tearDown(self):
        self.loop.close()

    def fut(self, done):
        f = self.loop.create_future()
        if done:
            f.set_result(None)
        return f

    def live_compass(self, group, results, pending=()):
        """A climber with a live compass wave: `results` are back (config, score), `pending` still running."""
        centre = self.space._assemble({MMQ: 40, "assembly_region_padding": 50})
        ctx = JumpCtx(self.space, lambda c: 80.0)
        ctx.done[centre.key] = (centre, 80.0)
        cl = Climber(ctx, 0, centre, 80.0)
        cl.climb_until = 1e9
        cl.jumps_enabled = False
        cl._wave_id = 1
        cl.batch_live, cl.phase, cl.kind, cl.group = True, "compass", "read", list(group)
        cl.wave_before, cl.wave_centre, cl.wave_results = 80.0, centre, []
        cl._wave_t0, cl._wave_n = 0.0, 0
        cl.gate_futs = [self.fut(True) for _ in results] + [self.fut(False) for _ in pending]
        for c, s in results:
            ctx.done[c.key] = (c, s)
            cl._watch(c)
        for c in pending:
            cl._watch(c)
        return ctx, cl, centre

    def step(self, centre, name, k):
        return self.space._assemble({**centre.params, name: self.lad.shift(name, centre.params[name], k)})

    def centre(self):
        return self.space._assemble({MMQ: 40, "assembly_region_padding": 50})

    def test_a_run_that_clears_the_bar_moves_the_centre_before_the_wave_is_back(self):
        c0 = self.centre()
        up, slow = self.step(c0, MMQ, 1), self.step(c0, MMQ, -1)
        ctx, cl, centre = self.live_compass([MMQ], [(up, 80.5)], pending=[slow])
        self.assertFalse(cl.gate_done())
        self.assertTrue(cl.early_gain())
        cl.close_wave()
        self.assertEqual(cl.centre.key, up.key)
        self.assertIn(slow.key, cl.owed)                    # the slow run still counts when it returns
        cl._wave_id, cl.batch_live = 2, True
        ctx.done[slow.key] = (slow, 81.0)
        self.assertTrue(cl.claim_higher())
        self.assertEqual(cl.centre.key, slow.key)

    def test_a_gain_below_the_bar_waits_for_the_wave(self):
        c0 = self.centre()
        up, slow = self.step(c0, MMQ, 1), self.step(c0, MMQ, -1)
        ctx, cl, centre = self.live_compass([MMQ], [(up, 80.003)], pending=[slow])
        self.assertFalse(cl.early_gain())
        with unittest.mock.patch.dict(TUNING.values, {"climb.move_at_once": 0}):
            ctx.done[up.key] = (up, 80.5)
            self.assertFalse(cl.early_gain())

    def test_no_cube_on_the_old_centre_once_a_compass_point_clears_the_bar(self):
        c0 = self.centre()
        pad = "assembly_region_padding"
        a, b = self.step(c0, MMQ, 1), self.step(c0, pad, 1)
        ctx, cl, centre = self.live_compass([MMQ, pad], [(a, 80.5), (b, 80.002)])
        cl.close_wave()
        self.assertNotEqual(cl.phase, "confirm")
        self.assertEqual(cl.centre.key, a.key)

    def test_the_cube_still_runs_when_no_single_point_clears_the_bar(self):
        c0 = self.centre()
        pad = "assembly_region_padding"
        a, b = self.step(c0, MMQ, 1), self.step(c0, pad, 1)
        ctx, cl, centre = self.live_compass([MMQ, pad], [(a, 80.003), (b, 80.002)])
        cl.close_wave()
        self.assertEqual(cl.phase, "confirm")
        self.assertEqual(cl.centre.key, centre.key)

    def test_a_grid_wave_closes_at_the_quorum_but_a_certificate_waits_for_every_neighbour(self):
        ctx, cl, centre = self.live_compass([MMQ], [])
        cl.gate_futs = [self.fut(True) for _ in range(4)] + [self.fut(False)]
        self.assertTrue(cl.gate_done())                     # 4 of 5 = climb.quorum 0.8
        cl.kind = "jumps"
        self.assertFalse(cl.gate_done())

    def test_step_c_closes_its_wave_on_a_higher_round_best(self):
        ctx, cl, centre = self.live_compass([MMQ], [], pending=[self.step(self.centre(), MMQ, 1)])
        cl.origin = "step C"
        other = self.space._assemble({MMQ: 30, "assembly_region_padding": 70})
        ctx.best = lambda: (other, 80.4)
        self.assertTrue(cl.early_gain())
        ctx.best = lambda: (other, 80.002)
        self.assertFalse(cl.early_gain())


class GaOpeningOffTest(unittest.TestCase):
    """climb.ga_open_min 0 (the default since 2026-10-18): Step B starts with the grids, no whole-fleet GA first."""

    def test_no_ga_opening_when_it_is_off(self):
        import asyncio
        from search.climb import open_with_ga
        space = Space(ranges())
        ctx = JumpCtx(space, lambda c: 80.0)
        cl = Climber(ctx, 0, space._assemble({MMQ: 40}), 82.0)
        with unittest.mock.patch("search.climb.free_workers", return_value=10):
            asyncio.run(open_with_ga(ctx, [cl], until=1e9))
        self.assertFalse(getattr(ctx, "ga_opening", False))
