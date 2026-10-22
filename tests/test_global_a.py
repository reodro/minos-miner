"""Island mode (search/global_a.py): the band, and a whole simulated round with islands.global_a on."""
import unittest
import unittest.mock

from common.tuning import TUNING
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search import engine, global_a
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges
from tests.test_engine import ROUND, Stub


class BandTest(unittest.TestCase):
    def test_band_fixes_pcr_and_keeps_numbers_near_the_centre(self):
        space = Space(ranges())
        lad = Ladders(space)
        centre = dict(space._assemble({}).params, pcr_indel_model="AGGRESSIVE", min_pruning=3,
                      phred_scaled_global_read_mismapping_rate=40)
        b = global_a.banded(space, centre)
        self.assertEqual(b.specs["pcr_indel_model"].allowed, ("AGGRESSIVE",))
        self.assertNotIn("pcr_indel_model", b.names)
        blad = Ladders(b)
        for n in ("min_pruning", "phred_scaled_global_read_mismapping_rate"):
            i = lad.index(n, centre[n])
            for v in blad.rungs[n]:
                self.assertLessEqual(abs(lad.index(n, v) - i), global_a.BAND, (n, v))
        self.assertEqual(space.specs["pcr_indel_model"].allowed, Space(ranges()).specs["pcr_indel_model"].allowed)
        self.assertEqual(b.names[:0], [])
        self.assertTrue(set(b.names) <= set(space.names))

    def test_local_peaks_are_apart_and_highest_first(self):
        ys = [1, 5, 4, 3, 6, None, 2, 9, 1]
        self.assertEqual(global_a.local_peaks(ys), [7, 4, 1])


class StartsTest(unittest.TestCase):
    def test_starts_scale_with_the_fleet(self):
        self.assertEqual(global_a.n_starts(300), 8)
        self.assertEqual(global_a.n_starts(971), 24)
        self.assertEqual(global_a.n_starts(5000), global_a.MAX_STARTS)


class HeadsTest(unittest.TestCase):
    def test_heads_have_bands_that_do_not_overlap(self):
        from types import SimpleNamespace
        space = Space(ranges())
        ga = global_a.GlobalA.__new__(global_a.GlobalA)
        ga.lad = Ladders(space)
        ga.spans = {}                                   # no hills known: the rung-gap fallback
        base = space._assemble({}).params
        mk = lambda **kw: SimpleNamespace(params=dict(base, **kw), key=str(sorted(kw.items())))
        a = mk(pcr_indel_model="NONE", min_pruning=2, min_base_quality_score=20)
        b = mk(pcr_indel_model="NONE", min_pruning=2, min_base_quality_score=22)    # 2 rungs from a: same band
        c = mk(pcr_indel_model="NONE", min_pruning=2, min_base_quality_score=30)    # 10 rungs from a
        d = mk(pcr_indel_model="AGGRESSIVE", min_pruning=2, min_base_quality_score=20)
        heads = ga.heads([(a, 3.0), (b, 2.9), (c, 2.8), (d, 2.7)], 3)
        self.assertEqual([h.key for h, _ in heads], [a.key, c.key, d.key])
        self.assertTrue(ga.in_band(a, b))
        self.assertFalse(ga.in_band(a, c))
        # with hills: min_base_quality_score 18-24 and 28-34 are two hills; 20 and 22 share one, 30 is the other
        ga.spans = {"min_base_quality_score": [(20, 18, 24), (30, 28, 34)]}
        heads = ga.heads([(a, 3.0), (b, 2.9), (c, 2.8), (d, 2.7)], 3)
        self.assertEqual([h.key for h, _ in heads], [a.key, c.key, d.key])   # b shares a's hill; d: another pcr
        self.assertEqual(ga.bands_for(c), {"min_base_quality_score": (28, 34)})
        self.assertTrue(ga.in_band(a, b))
        self.assertFalse(ga.in_band(a, c))

    def test_hills_merge_plateau_bumps_and_drop_low_hills(self):
        apt = [86.06, 86.28, 86.28, 86.17, 86.15, 85.98, 86.02, 86.22, 86.20, 77.92, 77.44, 76.9, 75.5, 75.9]
        self.assertEqual(global_a.hills(apt), [(1, 0, 8)])      # one plateau; the hill behind the cliff is too low
        self.assertEqual(global_a.hills([1, 5, 2, 0, 5.2, 4]), [(4, 4, 4), (1, 1, 1)])


class CertifyOutsideBandTest(unittest.TestCase):
    def test_certify_tries_one_rung_outside_the_band_and_climbing_does_not(self):
        from types import SimpleNamespace
        from search.climb import one_rung_neighbors
        space = Space(ranges())
        head = space._assemble(dict(space._assemble({}).params, min_base_quality_score=20, pcr_indel_model="NONE"))
        band = global_a.banded(space, head.params, {"min_base_quality_score": (19, 20)})   # hill 19-20, +-2: 18-22
        self.assertEqual(sorted(band.specs["min_base_quality_score"].allowed or
                                range(int(band.specs["min_base_quality_score"].min),
                                      int(band.specs["min_base_quality_score"].max) + 1)), [18, 19, 20, 21, 22])
        centre = space._assemble(dict(head.params, min_base_quality_score=22))         # climbed to the band's edge
        ctx = SimpleNamespace(space=band, ladders=Ladders(band), full_space=space, certifying=False)
        mbq = lambda cs: sorted({c.params["min_base_quality_score"] for c in cs})
        self.assertEqual(mbq(one_rung_neighbors(ctx, centre)), [21, 22])        # climbing stays in the band (22: moves
                                                                                # of other settings)
        ctx.certifying = True
        self.assertEqual(mbq(one_rung_neighbors(ctx, centre)), [21, 22, 23])    # certify: one rung outside
        self.assertEqual({c.params["pcr_indel_model"] for c in one_rung_neighbors(ctx, centre)}, {"NONE"})


class IslandModeRoundTest(unittest.TestCase):
    def test_fill_uses_the_free_workers_after_a1(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a": 1, "islands.workers": 150,
                                                      "islands.global_a_fill": 1}):
            with self.assertLogs("search.global_a", level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        sent = {l.split("global A ")[1].split(":")[0]: int(l.split(": ")[1].split(" sent")[0])
                for l in logs.output if "global A A" in l and " sent" in l}
        self.assertGreaterEqual(sent.get("A3", 0), 250)     # a wave after A1 fills most of the 300 workers
        self.assertIsNotNone(out.best)

    def test_quorum_round(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 100.0 + (hash(job) % 7) * 20.0)   # uneven runs
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a": 1, "islands.workers": 150,
                                                      "islands.global_a_waves": 1, "islands.global_a_quorum": 0.8}):
            with self.assertLogs("search.global_a", level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        a1 = next(l for l in logs.output if "global A A1" in l)
        sent, scored = int(a1.split(": ")[1].split(" sent")[0]), int(a1.split(", ")[1].split(" scored")[0])
        self.assertGreaterEqual(scored, int(0.8 * sent) - 5)
        self.assertIsNotNone(out.best)

    def test_one_wave_round(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a": 1, "islands.workers": 150,
                                                      "islands.global_a_waves": 1, "islands.global_a_fill": 1}):
            with self.assertLogs("search.global_a", level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertIn("global A A1", text)
        self.assertNotIn("global A A3", text)
        self.assertNotIn("global A A4", text)
        self.assertIn("segments", text)
        self.assertIsNotNone(out.best)

    def test_additive_climb_round(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a": 1, "islands.workers": 150,
                                                      "islands.global_a_waves": 1, "climb.additive": 1,
                                                      "climb.group_size": 9}):
            with self.assertLogs("search.engine", level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertIn("additive jump on", text)              # winners combined, not the 3^m cube
        self.assertNotIn(": cube on", text)
        self.assertIsNotNone(out.best)

    def test_probe_extra_seats_borrow_round(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        knobs = {"islands.global_a": 1, "islands.workers": 150, "islands.global_a_waves": 1,
                 "islands.global_a_probe": 150.0, "islands.global_a_extra": 1, "islands.global_seats": 4,
                 "candidates.workers_per_candidate": 40, "climb.borrow": 1, "climb.additive": 1, "climb.group_size": 9,
                 "islands.seat_apart": 3, "climb.min_apart": 2}
        with unittest.mock.patch.dict(TUNING.values, knobs):
            with self.assertLogs(level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertIn("global A probe:", text)                 # starts ran alone first
        seats = [l for l in logs.output if "step A: global (island mode)" in l]
        self.assertTrue(seats and all("segment(s)" in l for l in seats))
        self.assertTrue(any(int(l.split("island mode): ")[1].split(" segment")[0]) >= 2 for l in seats))
        self.assertIsNotNone(out.best)

    def test_min_dist_round(self):
        import math
        from search import climb
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        knobs = {"islands.global_a": 1, "islands.workers": 150, "islands.global_a_waves": 1, "islands.global_seats": 4,
                 "candidates.workers_per_candidate": 40, "climb.borrow": 1, "climb.additive": 1, "climb.group_size": 9,
                 "islands.seat_apart": 3, "climb.min_apart": 2, "climb.min_dist": 0.5, "climb.rf_ucb_share": 0.5}
        gaps, real = [], climb.Climber.move

        def move(cl, c, s, why):                   # after each move: the mover's distance to every other climber
            real(cl, c, s, why)
            if why.startswith("C ") or cl.state != "climbing" or cl.origin == "step C":
                return
            v = climb.encoded(cl.ctx, cl.centre)
            gaps.extend(math.dist(v, climb.encoded(cl.ctx, o.centre)) for o in cl.ctx.climbers
                        if o is not cl and o.state == "climbing" and o.centre is not None)
        with unittest.mock.patch.dict(TUNING.values, knobs), unittest.mock.patch.object(climb.Climber, "move", move):
            with self.assertLogs(level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        self.assertTrue(gaps)                                      # climbers moved with others climbing
        self.assertGreaterEqual(min(gaps), 0.5)
        self.assertIsNotNone(out.best)

    # The live tuning.json of 2026-10-22 (islands.workers scaled to the simulated fleet). 02:34 live: step B of island 1
    # crashed in the GA hand-off (a renamed local shadowed the near-climber list) - no test ran these knobs together.
    LIVE = {"climb.leap_share": 0.5, "climb.grid_fill": 0.5, "climb.rf_ucb_share": 0.5, "climb.all_steps": 1,
            "climb.grid_cores": 0.4, "climb.jump_cores": 0.3, "islands.by_climber": 1, "climb.heads_first": 1,
            "climb.flat_shift": 0.5, "islands.global_a": 1, "islands.global_k": 2, "islands.global_a_fill": 1,
            "islands.global_a_waves": 1, "climb.group_size": 9, "islands.global_a_quorum": 0.8, "climb.additive": 1,
            "islands.global_a_probe": 150.0, "islands.global_a_extra": 1, "islands.global_seats": 4,
            "candidates.workers_per_candidate": 40, "climb.borrow": 1, "islands.seat_apart": 3, "climb.min_apart": 2,
            "evaluator.split_k": 1, "far.step_c_top": 2, "climb.min_dist": 0.5, "islands.workers": 150}

    def test_multi_setting_moves_round(self):
        # offline test T1: each grid candidate moves 3-5 settings of its group (a smaller group: all of them)
        from search import climb
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        sizes, real = [], climb.Climber.multi_heads

        def spy(cl, centre, group, n):
            out = real(cl, centre, group, n)
            names = [x for x in group if any(o for o in cl.compass_offsets(x))]
            sizes.extend((sum(c.params.get(x) != centre.params.get(x) for x in cl.all_names), len(names)) for c in out)
            return out
        knobs = dict(self.LIVE, **{"climb.move_k_min": 3, "climb.move_k_max": 5})
        with unittest.mock.patch.dict(TUNING.values, knobs), unittest.mock.patch.object(climb.Climber, "multi_heads", spy):
            with self.assertLogs(level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertNotIn("step B failed", text)
        self.assertNotIn("Traceback", text)
        self.assertTrue(sizes)                                          # the variant built the grid waves
        self.assertTrue(all(min(3, m) <= k <= 5 for k, m in sizes), [x for x in sizes if not min(3, x[1]) <= x[0] <= 5][:10])
        self.assertIsNotNone(out.best)

    def test_t1_independent_climbers_round(self):
        # offline test T1: one island per (seat, variant) climber, no GA/RF hand-offs; only the B islands build
        # multi-setting grid candidates
        import json, os, tempfile
        from search import climb
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        base = r.defaults("gatk").params
        other = Space(r).random(__import__("random").Random(3)).params
        seats = [{"params": p, "move_k": mk} for p in (base, other) for mk in (None, None, [3, 5])]
        f = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump({"seats": seats}, f); f.close()
        calls, real = [], climb.Climber.multi_heads

        def spy(cl, centre, group, n):
            calls.append(getattr(getattr(cl.ctx, "island", None), "index", None))
            return real(cl, centre, group, n)
        knobs = dict(self.LIVE, **{"islands.global_a": 0, "climb.borrow": 0, "climb.min_dist": 0.0, "climb.min_apart": 0,
                                   "far.step_c_top": 0, "climb.rf_ucb_share": 0.0})
        env = {"SEARCH_SEATS": f.name, "MINOS_T1_ISOLATE": "1"}
        try:
            with unittest.mock.patch.dict(TUNING.values, knobs), unittest.mock.patch.dict(os.environ, env), \
                    unittest.mock.patch.object(climb.Climber, "multi_heads", spy):
                with self.assertLogs(level="INFO") as logs:
                    out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(),
                                    r, LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                    deadline_ts=72 * 60)
        finally:
            os.unlink(f.name)
        text = "\n".join(logs.output)
        self.assertIn("T1 test: 6 independent climbers", text)
        self.assertNotIn("takes the GA best", text)
        self.assertNotIn("step B failed", text)
        self.assertNotIn("Traceback", text)
        self.assertTrue(calls)
        self.assertEqual(set(calls) - {2, 5}, set())                    # only the B islands (2 and 5)
        self.assertIsNotNone(out.best)

    def _a1(self, knobs, duration, crashes=lambda job: False, workers=300):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, workers, duration=duration, crashes=crashes)
        with unittest.mock.patch.dict(TUNING.values, dict(self.LIVE, **knobs)):
            with self.assertLogs(level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        return out, "\n".join(logs.output)

    def test_idle_fix_1_starts_in_flight_replace_a_crashed_start(self):
        crash = lambda job: (job.config.params.get("pcr_indel_model") == "HOSTILE"
                             and job.config.params.get("min_pruning") == 3)
        out, text = self._a1({"islands.global_a_inline": 1}, lambda job: 150.0, crash)
        self.assertNotIn("global A probe:", text)                      # no probe wait
        line = next(l for l in text.splitlines() if "global A A1 (starts in flight)" in l)
        self.assertIn("HOSTILE/3", line)                                # the crashed start was replaced
        self.assertNotIn(" 0 start(s) crashed", line)
        self.assertIsNotNone(out.best)

    def test_idle_fix_2_an_emptying_fleet_ends_the_sweep(self):
        # live 04:58: busy 993 -> 735 -> 317 -> 114 from t+8 to t+14 min while A1 waited for its 80% quorum
        from types import SimpleNamespace as NS
        stub = lambda idle: NS(client=NS(idle_slots=lambda: idle, capacity=lambda: 993))
        end = lambda knob, idle, back, sent: global_a.GlobalA.emptying(stub(idle), back, sent)
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a_tail": 0.5}):
            self.assertTrue(end(0.5, 876, 4000, 5300))      # 75% back, 876 of 993 idle (t+14): the wave ends
            self.assertFalse(end(0.5, 258, 4000, 5300))     # 75% back but 735 still busy (t+10): keep waiting
            self.assertFalse(end(0.5, 900, 3000, 5300))     # idle, but only 57% back: keep waiting
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a_tail": 0.0}):
            self.assertFalse(end(0.0, 990, 5000, 5300))     # off

    def test_generalized_step_a_side_by_side(self):
        # offline test G: both Step A versions in one round, each on half the workers; G covers every searched setting
        import os
        from search import global_g
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        sent = {}
        real = global_g.GeneralizedA.run

        async def spy(g):
            out = await real(g)
            sent.update({n: {json.dumps(c.params.get(n)) for _, c in g.plan()} for n in g.values})
            sent["_cells"] = g.report["cells"]
            return out
        import json
        knobs = dict(self.LIVE, **{"islands.global_a_share": 1.0})
        with unittest.mock.patch.dict(TUNING.values, knobs), unittest.mock.patch.dict(os.environ, {"MINOS_STEPA_AB": "20"}), \
                unittest.mock.patch.object(global_g.GeneralizedA, "run", spy):
            with self.assertLogs(level="INFO") as logs:
                run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                          LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")], deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertIn("Step A test:", text)
        self.assertIn("global A A1", text)                           # today's Step A ran
        self.assertIn("generalized A:", text)                        # and the generalized one
        self.assertNotIn("Traceback", text)
        self.assertGreaterEqual(sent["_cells"], 30)                  # ~36 sub-domains (+ anchors)
        many = [n for n in sent if n != "_cells" and len(sent[n]) >= 3]
        self.assertGreaterEqual(len(many), 12)                       # most settings take several values

    def test_territories_split_two_hills(self):
        # nearest-better clustering: two separate hills come out as two territories, each rooted at its top
        from common.types import Config
        from search import territory
        values = {"x": list(range(21)), "y": list(range(21))}
        pts = []
        for cx, cy, top in ((3, 3, 90.0), (17, 17, 89.0)):
            for dx in (-2, -1, 0, 1, 2):
                for dy in (-2, -1, 0, 1, 2):
                    pts.append((Config("gatk", {"x": cx + dx, "y": cy + dy}), top - 0.3 * (abs(dx) + abs(dy))))
        t = territory.territories(pts, values)
        self.assertEqual(len(t), 2)
        self.assertEqual((t[0]["root"].params, t[1]["root"].params), ({"x": 3, "y": 3}, {"x": 17, "y": 17}))
        self.assertEqual((t[0]["size"], t[1]["size"]), (25, 25))
        self.assertEqual(territory.band(values, t[1]["members"], ["x"]), {"x": (13, 20)})

    def test_generalized_step_a_round(self):
        # islands.global_g: scan -> territories -> race -> islands, under the live settings
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        with unittest.mock.patch.dict(TUNING.values, dict(self.LIVE, **{"islands.global_g": 1, "islands.global_a_tail": 0.5})):
            with self.assertLogs(level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertIn("generalized A:", text)
        self.assertIn("engine: generalized Step A:", text)
        self.assertNotIn("Traceback", text)
        self.assertNotIn("step B failed", text)
        line = next(l for l in logs.output if "engine: generalized Step A:" in l)
        if "race of 1 " not in line:                                 # more than global_k territories: a race ran
            self.assertIn("race over", text)
            self.assertIn("lost the race", text)
        self.assertIsNotNone(out.best)

    def test_build_b_rounds(self):
        # Build B1 (scan -> today's value-based heads) and B2 (scan -> best complete configs as heads), live settings
        r = ranges()
        for heads in (0, 1):
            runner = FakeRunner(r, seed=0, batching=True, hills=2)
            ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
            knobs = dict(self.LIVE, **{"islands.global_a_scan": 1, "islands.global_a_heads": heads})
            with unittest.mock.patch.dict(TUNING.values, knobs):
                with self.assertLogs(level="INFO") as logs:
                    out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(),
                                    r, LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                    deadline_ts=72 * 60)
            text = "\n".join(logs.output)
            self.assertIn("generalized A:", text, f"heads {heads}")
            self.assertNotIn("global A A1:", text)                    # the one-setting sweep did not run
            self.assertIn("island mode:", text)
            self.assertNotIn("Traceback", text)
            self.assertNotIn("step B failed", text)
            self.assertIsNotNone(out.best)

    def test_live_knobs_round(self):
        r = ranges()
        for seed in (0, 1, 2):
            runner = FakeRunner(r, seed=seed, batching=True, hills=2)
            ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
            with unittest.mock.patch.dict(TUNING.values, self.LIVE):
                with self.assertLogs(level="INFO") as logs:
                    out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(),
                                    r, LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                    deadline_ts=72 * 60)
            text = "\n".join(logs.output)
            self.assertNotIn("step B failed", text, f"seed {seed}")
            self.assertNotIn("Traceback", text, f"seed {seed}")
            self.assertIsNotNone(out.best)

    def test_step_c_top2_round(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        knobs = {"islands.global_a": 1, "islands.workers": 150, "islands.global_a_waves": 1, "far.step_c_top": 2,
                 "islands.global_seats": 3, "candidates.workers_per_candidate": 40}
        with unittest.mock.patch.dict(TUNING.values, knobs):
            with self.assertLogs("search.engine", level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        tops = [l for l in logs.output if "the round's top 2 climbers" in l]
        self.assertEqual(len(tops), 2)                                   # both islands took part
        self.assertEqual(len({l.split("climbers ")[1].split(";")[0] for l in tops}), 1)   # and agree on the top 2
        self.assertEqual(sum(1 for l in logs.output if "step C: fine-tune" in l and "opened at" in l), 2)
        self.assertIsNotNone(out.best)

    def test_short_tail_round(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a": 1, "islands.workers": 150,
                                                      "islands.global_a_short": 1}):
            with self.assertLogs("search.global_a", level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), engine.searcher(r, 300), DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertIn("A4_1_pm3+1", text)
        self.assertNotIn("A4_3", text)
        self.assertIsNotNone(out.best)

    def test_a_whole_round_in_island_mode(self):
        r = ranges()
        runner = FakeRunner(r, seed=0, batching=True, hills=2)
        ex = SimExecutor(runner, 300, duration=lambda job: 150.0)
        seen = []
        with unittest.mock.patch.dict(TUNING.values, {"islands.global_a": 1, "islands.workers": 150}):
            searcher = engine.searcher(r, 300)

            async def wrapped(client):
                seen.append(client)
                await searcher(client)

            with self.assertLogs("search.global_a", level="INFO") as logs:
                out = run_round(ROUND, Stub(), ex, FakeScorer(runner), wrapped, DryRunSubmitter(), r,
                                LoopSettings(stop_before_deadline_s=60), baselines=[r.defaults("gatk")],
                                deadline_ts=72 * 60)
        text = "\n".join(logs.output)
        self.assertIn("global A A1", text)
        self.assertIn("segments", text)
        self.assertIsNotNone(out.best)
        notes = seen[0].notes
        self.assertTrue(notes, "the round left no engine notes")


if __name__ == "__main__":
    unittest.main()
