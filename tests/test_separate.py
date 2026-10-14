"""v9.29 separate mountains: the line, the view, the dummy, and the parent cap."""
import asyncio
import random
import unittest
from unittest.mock import patch

from search.piles import Pile, absorb_parents
from search.separate import (
    LINKED_GROUPS, _cover, _cross, _far_from, aim_compass, best_grid_point, claim_fingerprint, climb_probes,
    compass_axes, compass_points, cube_points, dummy_config, fill_free, gate, hill_children, line_point,
    valley_status, winning_faces,
)
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges


class _Ctx:
    def __init__(self):
        self.space = Space(ranges())
        self.ladders = Ladders(self.space)
        self.rng = random.Random(0)

    def fresh(self, config):
        return config is not None


def _cfg(ctx, **updates):
    params = dict(ctx.space._assemble({}).params)
    params.update(updates)
    return ctx.space._assemble(params)


class _Dist:
    def __init__(self):
        self.table = {}

    def link(self, a, b, dist):
        self.table[(a.key, b.key)] = dist
        self.table[(b.key, a.key)] = dist

    def between(self, a, b):
        if a.key == b.key:
            return 0.0, "callset"
        return self.table.get((a.key, b.key), 0.9), "callset"


class _Board:
    """Scores, callset distances, and the small futures a valley probe sends."""

    def __init__(self):
        self.space = Space(ranges())
        self.ladders = Ladders(self.space)
        self.rng = random.Random(0)
        self.dist = _Dist()
        self.scores = {}
        self.configs = {}
        self._piles = []
        self._pile_seen = set()
        self.climbers = []
        self.sent = []

    def add(self, config, score):
        self.scores[config.key] = float(score)
        self.configs[config.key] = config
        self._pile_seen.add(config.key)
        return config

    def score(self, config):
        if config is None:
            return None
        return self.scores.get(config.key)

    def scored(self):
        return [(self.configs[key], score) for key, score in self.scores.items()]

    def fresh(self, config):
        return config is not None

    def say(self, *_args):
        return None

    def clock(self):
        return 0.0

    def over(self):
        return False

    def T(self, seconds):
        return seconds

    def propose(self, configs, _tag, _groups, priority=2):
        self.sent.extend(configs)
        loop = asyncio.get_running_loop()
        done = []
        for _config in configs:
            fut = loop.create_future()
            fut.set_result(None)
            done.append(fut)
        return done


class SeparateTest(unittest.TestCase):
    def setUp(self):
        self.ctx = _Ctx()

    def test_the_midpoint_of_an_ordered_setting_sits_between_the_centres(self):
        higher = _cfg(self.ctx, min_pruning=2)
        lower = _cfg(self.ctx, min_pruning=10)
        mid = line_point(self.ctx, higher, lower, 0.5)
        self.assertEqual(mid.params["min_pruning"], 6)

    def test_a_categorical_setting_keeps_the_higher_centre_through_the_halfway_point(self):
        name = "pcr_indel_model"
        allowed = self.ctx.space.specs[name].allowed
        higher = _cfg(self.ctx, **{name: allowed[0]})
        lower = _cfg(self.ctx, **{name: allowed[-1]})
        self.assertEqual(line_point(self.ctx, higher, lower, 0.5).params[name], allowed[0])
        self.assertEqual(line_point(self.ctx, higher, lower, 0.6).params[name], allowed[-1])

    def test_the_compass_tries_one_rung_either_side(self):
        centre = _cfg(self.ctx, min_pruning=3, assembly_region_padding=10)
        points = compass_points(self.ctx, centre, ["min_pruning", "assembly_region_padding"])
        pruning = sorted(p.params["min_pruning"] for p in points if p.params["assembly_region_padding"] == 10)
        padding = sorted(p.params["assembly_region_padding"] for p in points if p.params["min_pruning"] == 3)
        self.assertEqual(pruning, [2, 4])
        self.assertEqual(padding, [9, 11])
        self.assertTrue(all(
            sum(p.params[n] != centre.params[n] for n in self.ctx.space.names) == 1 for p in points))

    def test_two_live_faces_send_only_their_joint_config(self):
        centre = _cfg(self.ctx, min_pruning=3, assembly_region_padding=10)
        board = _Board()
        board.add(centre, 80.0)
        up = board.add(_cfg(board, min_pruning=4, assembly_region_padding=10), 81.0)
        side = board.add(_cfg(board, min_pruning=3, assembly_region_padding=11), 80.5)
        flat = board.add(_cfg(board, min_pruning=2, assembly_region_padding=10), 79.0)
        names = ["min_pruning", "assembly_region_padding"]
        faces = winning_faces(board, centre, names, [up, side, flat])
        self.assertEqual(faces, {"min_pruning": 4, "assembly_region_padding": 11})
        joint = cube_points(self.ctx, centre, faces)
        self.assertEqual(len(joint), 1)
        self.assertEqual(joint[0].params["min_pruning"], 4)
        self.assertEqual(joint[0].params["assembly_region_padding"], 11)

    def test_the_next_compass_takes_the_highest_score_on_this_mountain(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        higher = board.add(_cfg(board, min_pruning=4), 90.0)
        board.dist.link(centre, higher, 0.05)
        board.best = lambda: (higher, 90.0)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        asyncio.run(aim_compass(board, cl))
        self.assertEqual(cl.centre.key, higher.key)
        self.assertEqual(cl.score, 90.0)

    def test_a_far_fingerprint_does_not_become_the_compass_before_the_valley(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        other = board.add(_cfg(board, min_pruning=8), 95.0)
        board.dist.link(centre, other, 0.50)
        board.best = lambda: (other, 95.0)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        asyncio.run(aim_compass(board, cl))
        self.assertEqual(cl.centre.key, centre.key)

    def test_the_last_twelve_minutes_reach_one_rung_further(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        board.pick_at = lambda: 1000.0
        centre = _cfg(board, min_pruning=3, assembly_region_padding=10)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        cl.compass_only = True
        cl.line_keys = set()
        early = [p.params["assembly_region_padding"] for p in cl.compass_heads(centre, ["assembly_region_padding"])]
        self.assertEqual(sorted(early), [9, 11])
        board.clock = lambda: 400.0
        late = [p.params["assembly_region_padding"] for p in cl.compass_heads(centre, ["assembly_region_padding"])]
        self.assertEqual(sorted(late), [8, 9, 11, 12])
        outer = _cfg(board, min_pruning=5, assembly_region_padding=10)
        inner = _cfg(board, min_pruning=4, assembly_region_padding=10)
        self.assertEqual(cl.live_axes(centre, ["min_pruning"], [(outer, 81.0)], 80.0), [])
        self.assertEqual(
            cl.live_axes(centre, ["min_pruning"], [(inner, 80.06), (outer, 81.0)], 80.0),
            ["min_pruning"])
        side = _cfg(board, min_pruning=3, assembly_region_padding=11)
        cube = cl.confirm_heads(centre, ["min_pruning", "assembly_region_padding"], set())
        self.assertTrue(cube)
        self.assertTrue(all(c.params["min_pruning"] in (2, 3, 4) for c in cube))
        self.assertNotIn(side.key, {c.key for c in cube})

    def test_the_fine_tune_keeps_a_far_fingerprint_off_the_compass(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        near = board.add(_cfg(board, min_pruning=4), 82.0)
        far = board.add(_cfg(board, min_pruning=8), 99.0)
        board.dist.link(centre, near, 0.05)
        board.dist.link(centre, far, 0.50)
        board.dist.link(near, far, 0.50)
        board.best = lambda: (near, 82.0)
        cl = Climber(board, 0, centre, 80.0, origin="step C")
        cl.compass_only = True
        asyncio.run(aim_compass(board, cl, valley=False))
        self.assertEqual(cl.centre.key, near.key)
        board.best = lambda: (far, 99.0)
        asyncio.run(aim_compass(board, cl, valley=False))
        self.assertEqual(cl.centre.key, near.key)

    def test_the_compass_tries_the_axis_that_moved_the_score(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3, assembly_region_padding=10), 80.0)
        pruning = board.add(_cfg(board, min_pruning=4, assembly_region_padding=10), 80.2)
        padding = board.add(_cfg(board, min_pruning=3, assembly_region_padding=11), 83.0)
        board.dist.link(centre, pruning, 0.02)
        board.dist.link(centre, padding, 0.02)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        cl.compass_only = True
        axes = compass_axes(board, cl)
        self.assertEqual(axes[0], "assembly_region_padding")
        self.assertIn("min_pruning", axes)
        self.assertLessEqual(len(axes), 4)
        cl.flat.add("assembly_region_padding")
        self.assertNotIn("assembly_region_padding", compass_axes(board, cl))

    def test_a_tied_neighbour_does_not_take_a_compass_slot(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3, assembly_region_padding=10), 80.0)
        tied = board.add(_cfg(board, min_pruning=4, assembly_region_padding=10), 80.0)
        moved = board.add(_cfg(board, min_pruning=3, assembly_region_padding=11), 83.0)
        board.dist.link(centre, tied, 0.02)
        board.dist.link(centre, moved, 0.02)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        cl.compass_only = True
        axes = compass_axes(board, cl)
        self.assertEqual(len(axes), 4)
        self.assertEqual(axes[0], "assembly_region_padding")
        self.assertNotIn("min_pruning", axes)

    def test_a_centre_move_reopens_a_missed_axis(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        nxt = _cfg(board, min_pruning=4)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        cl.compass_only = True
        cl.flat.update(cl.read_names + cl.assembly_names)
        cl.move(nxt, 81.0, "grid")
        self.assertEqual(cl.flat, set())
        cl.note_gain(["min_pruning"], 0.0)
        self.assertIn("min_pruning", cl.flat)
        self.assertNotIn("assembly_region_padding", cl.flat)
        self.assertEqual(cl.state, "climbing")

    def test_a_finished_centre_does_not_end_step_b(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        cl.compass_only = True
        cl.jumps_enabled = False
        cl.flat.update(cl.read_names + cl.assembly_names)
        self.assertFalse(cl.begin_grid(1e9))
        self.assertEqual(cl.state, "climbing")
        cl.wave_before = cl.score
        cl.kind, cl.group = "read", ["min_pruning"]
        cl.wave_results = []
        cl._wave_t0 = 0.0
        cl.gate_futs = []
        cl._finish_grid_wave()
        self.assertEqual(cl.state, "climbing")

    def test_a_face_inside_0_05_is_not_live(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = _cfg(board, min_pruning=3, assembly_region_padding=10)
        noise = _cfg(board, min_pruning=4, assembly_region_padding=10)
        live = _cfg(board, min_pruning=3, assembly_region_padding=11)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        cl.compass_only = True
        cl.line_keys = set()
        group = ["min_pruning", "assembly_region_padding"]
        self.assertEqual(
            cl.live_axes(centre, group, [(noise, 80.04), (live, 80.2)], 80.0),
            ["assembly_region_padding"])

    def test_a_hill_child_stays_0_10_from_the_climber(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        parent = board.add(_cfg(board, min_pruning=8), 70.0)
        board.dist.link(centre, parent, 0.50)
        board._far_bases = [parent]
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        close = _cfg(board, min_pruning=4)
        board.dist.link(centre, close, 0.05)
        self.assertFalse(_far_from(board, close, [centre]))
        kids = hill_children(board, cl, 3)
        self.assertTrue(kids)
        names = [n for n in board.space.names if board.ladders.movable(n)]
        for kid in kids:
            diffs = sum(kid.params.get(n) != centre.params.get(n) for n in names)
            self.assertGreaterEqual(diffs, 4)
            self.assertGreaterEqual(board.dist.between(kid, centre)[0], 0.10)

    def test_the_hill_lane_crosses_two_piles(self):
        from search.climb import Climber
        from search.piles import Pile
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        left_p = dict(centre.params)
        right_p = dict(centre.params)
        moved = []
        for name in ("min_pruning", "min_dangling_branch_length", "assembly_region_padding",
                     "heterozygosity", "active_probability_threshold", "min_base_quality_score"):
            shifted = board.ladders.shift(name, centre.params[name], 4) or board.ladders.shift(name, centre.params[name], -4)
            if shifted is None:
                continue
            left_p[name] = shifted
            other = board.ladders.shift(name, centre.params[name], -4) or board.ladders.shift(name, centre.params[name], 6)
            right_p[name] = other or shifted
            moved.append(name)
        left = board.add(board.space._assemble(left_p), 70.0)
        right = board.add(board.space._assemble(right_p), 71.0)
        board.dist.link(centre, left, 0.40)
        board.dist.link(centre, right, 0.40)
        board.dist.link(left, right, 0.40)
        board._piles = [Pile(left, 70.0), Pile(right, 71.0)]
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        kids = hill_children(board, cl, 4)
        self.assertTrue(kids)
        differed = [name for name in moved if left.params.get(name) != right.params.get(name)]
        self.assertGreaterEqual(len(differed), 4)
        for kid in kids:
            from_left = sum(kid.params.get(name) == left.params.get(name) for name in differed)
            from_right = sum(kid.params.get(name) == right.params.get(name) for name in differed)
            self.assertGreater(from_left, 0)
            self.assertGreater(from_right, 0)
            self.assertGreaterEqual(board.dist.between(kid, centre)[0], 0.10)

    def test_a_cross_keeps_a_verified_group_with_one_parent(self):
        board = _Board()
        board.rng = random.Random(1)
        left = _cfg(board)
        right_p = dict(left.params)
        for name in board.space.names:
            shifted = board.ladders.shift(name, left.params[name], 3)
            if shifted is not None and shifted != left.params[name]:
                right_p[name] = shifted
        if right_p.get("pcr_indel_model") == left.params.get("pcr_indel_model"):
            right_p["pcr_indel_model"] = "NONE" if left.params.get("pcr_indel_model") != "NONE" else "AGGRESSIVE"
        right = board.space._assemble(right_p)
        saw_left = {i: False for i in range(len(LINKED_GROUPS))}
        saw_right = {i: False for i in range(len(LINKED_GROUPS))}
        for _ in range(40):
            kid = _cross(board, left, right)
            self.assertIsNotNone(kid)
            for i, group in enumerate(LINKED_GROUPS):
                from_left = all(kid.params.get(name) == left.params.get(name) for name in group)
                from_right = all(kid.params.get(name) == right.params.get(name) for name in group)
                self.assertTrue(from_left or from_right)
                saw_left[i] = saw_left[i] or from_left
                saw_right[i] = saw_right[i] or from_right
        self.assertTrue(all(saw_left.values()))
        self.assertTrue(all(saw_right.values()))

    def test_one_pile_does_not_shift_a_linked_knob(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board), 80.0)
        parent_p = dict(centre.params)
        for name in ("min_base_quality_score", "pair_hmm_gap_continuation_penalty", "min_pruning"):
            shifted = board.ladders.shift(name, centre.params[name], 4)
            if shifted is not None:
                parent_p[name] = shifted
        parent = board.add(board.space._assemble(parent_p), 70.0)
        board.dist.link(centre, parent, 0.50)
        board._far_bases = [parent]
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        kids = hill_children(board, cl, 6)
        self.assertTrue(kids)
        for kid in kids:
            self.assertEqual(kid.params.get("min_base_quality_score"), parent.params.get("min_base_quality_score"))
            self.assertEqual(kid.params.get("pair_hmm_gap_continuation_penalty"),
                             parent.params.get("pair_hmm_gap_continuation_penalty"))
            self.assertGreaterEqual(board.dist.between(kid, centre)[0], 0.10)

    def test_the_population_mates_and_may_land_on_the_climber(self):
        from search.separate import mating_room, population_children
        board = _Board()
        left = board.add(_cfg(board, min_base_quality_score=12, pair_hmm_gap_continuation_penalty=8,
                               min_pruning=2), 80.0)
        right = board.add(_cfg(board, min_base_quality_score=20, pair_hmm_gap_continuation_penalty=14,
                                min_pruning=6), 79.0)
        board.dist.link(left, right, 0.02)
        kids = population_children(board, 12)
        self.assertTrue(kids)
        group = ("min_base_quality_score", "pair_hmm_gap_continuation_penalty")
        for kid in kids:
            from_left = all(kid.params.get(name) == left.params.get(name) for name in group)
            from_right = all(kid.params.get(name) == right.params.get(name) for name in group)
            self.assertTrue(from_left or from_right)
            self.assertFalse(from_left and from_right)
        self.assertEqual(mating_room(1000), 984)
        self.assertEqual(mating_room(4), 4)

    def test_step_b_breeds_inside_the_piles_and_looks_for_a_new_fingerprint(self):
        from search.piles import Pile
        from search.separate import pile_ga
        board = _Board()
        left = board.add(_cfg(board, min_base_quality_score=12, pair_hmm_gap_continuation_penalty=8,
                               active_probability_threshold=0.002, min_mapping_quality_score=20,
                               min_pruning=2), 80.0)
        right = board.add(_cfg(board, min_base_quality_score=20, pair_hmm_gap_continuation_penalty=14,
                                active_probability_threshold=0.004, min_mapping_quality_score=40,
                                min_pruning=6), 79.0)
        board.dist.link(left, right, 0.02)
        pile = Pile(left, 80.0)
        pile.members.append(right)
        pile.scores[right.key] = 79.0
        board._piles = [pile]
        cl = type("C", (), {"centre": left})()
        kids = pile_ga(board, cl, 8)
        self.assertTrue(kids)
        group = ("min_base_quality_score", "pair_hmm_gap_continuation_penalty",
                 "active_probability_threshold", "min_mapping_quality_score")
        parents = (left, right)
        for kid in kids:
            matched = [parent for parent in parents
                       if all(kid.params.get(name) == parent.params.get(name) for name in group)]
            self.assertEqual(len(matched), 1)
        far = [kid for kid in kids if board.dist.between(kid, left)[0] >= 0.10]
        self.assertTrue(far)

    def test_a_near_population_score_becomes_the_next_centre(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        near = board.add(_cfg(board, min_pruning=4), 81.0)
        far = board.add(_cfg(board, min_pruning=8), 95.0)
        board.dist.link(centre, near, 0.04)
        board.dist.link(centre, far, 0.50)
        board.dist.link(near, far, 0.50)
        board.best = lambda: (far, 95.0)
        cl = Climber(board, 0, centre, 80.0, origin="pile 1")
        cl.compass_only = True
        asyncio.run(aim_compass(board, cl))
        self.assertEqual(cl.centre.key, near.key)

    def test_step_a_presents_a_score_range_and_updates_it_in_a_small_batch(self):
        from search.separate import candidate_range, refine_range
        board = _Board()
        board.add(_cfg(board, min_pruning=2), 80.0)
        board.add(_cfg(board, min_pruning=4), 79.6)
        board.add(_cfg(board, min_pruning=8), 79.0)
        rows = candidate_range(board)
        self.assertEqual([score for _config, score in rows], [80.0, 79.6])
        board._step_a_open = True
        board._dummies_sent = True
        board._step_a_rest = [_cfg(board, min_pruning=3, assembly_region_padding=pad) for pad in range(20, 60)]
        take = refine_range(board, 80)
        self.assertTrue(take)
        self.assertLessEqual(len(take), 8)
        self.assertLess(len(board._step_a_rest), 40)
        self.assertEqual(board.candidate_range[0][1], 80.0)

        class _Pending:
            def done(self):
                return False

        board._range_futs = [_Pending()]
        self.assertEqual(refine_range(board, 80), [])

    def test_the_fine_tune_does_not_send_the_hill_lane(self):
        from search.climb import Climber
        board = _Board()
        board.elapsed = lambda: 0.0
        board.owner = lambda config, me=None: None
        centre = board.add(_cfg(board, min_pruning=3), 80.0)
        board._far_bases = [_cfg(board, min_pruning=8)]
        cl = Climber(board, 0, centre, 80.0, origin="step C")
        cl.compass_only = True
        cl.tried_axes = set(cl.read_names[1:] + cl.assembly_names)
        sent = []
        board.propose = lambda configs, tag, batches, priority=0: sent.extend(configs)
        asyncio.run(fill_free(board, cl, 4))
        self.assertFalse(getattr(board, "_hill_live", None))
        self.assertTrue(sent)
        self.assertTrue(all(
            sum(c.params.get(n) != centre.params.get(n) for n in cl.read_names + cl.assembly_names) <= 3
            for c in sent))
        self.assertTrue(climb_probes(board, cl, 1))

    def test_absorbing_parents_raises_the_cap_from_32_to_64(self):
        higher = Pile(_cfg(self.ctx, min_pruning=3), 80.0)
        lower = Pile(_cfg(self.ctx, min_pruning=8), 79.0)
        absorb_parents(self.ctx, higher, lower)
        self.assertEqual(higher.pool_cap, 64)
        self.assertEqual(higher.absorbed, 1)
        self.assertIn(lower.config.key, higher.scores)

    def test_a_scored_dip_on_the_path_is_a_valley(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        sample = board.add(_cfg(board, min_pruning=6), 69.0)
        board.dist.link(higher, lower, 0.40)
        board.dist.link(higher, sample, 0.20)
        board.dist.link(lower, sample, 0.20)
        decision, probe = valley_status(board, higher, lower)
        self.assertEqual(decision, "dip")
        self.assertIsNone(probe)

    def test_a_covered_middle_with_no_dip_is_the_same_slope(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        sample = board.add(_cfg(board, min_pruning=6), 80.5)
        board.dist.link(higher, lower, 0.40)
        board.dist.link(higher, sample, 0.20)
        board.dist.link(lower, sample, 0.20)
        decision, probe = valley_status(board, higher, lower)
        self.assertEqual(decision, "slope")
        self.assertIsNone(probe)

    def test_an_empty_middle_asks_for_the_midpoint(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        board.dist.link(higher, lower, 0.40)
        decision, probe = valley_status(board, higher, lower)
        self.assertEqual(decision, "need")
        self.assertEqual(probe.params["min_pruning"], 6)

    def test_a_midpoint_inside_the_band_asks_for_the_quarter(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        mid = board.add(line_point(board, higher, lower, 0.5), 79.95)
        board.dist.link(higher, lower, 0.40)
        board.dist.link(higher, mid, 0.20)
        board.dist.link(lower, mid, 0.20)
        decision, probe = valley_status(board, higher, lower)
        self.assertEqual(decision, "need")
        self.assertEqual(probe.key, line_point(board, higher, lower, 0.25).key)

    def test_a_quarter_still_inside_the_band_is_the_same_slope(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        mid = board.add(line_point(board, higher, lower, 0.5), 79.95)
        quarter = board.add(line_point(board, higher, lower, 0.25), 79.96)
        board.dist.link(higher, lower, 0.40)
        board.dist.link(higher, mid, 0.20)
        board.dist.link(lower, mid, 0.20)
        board.dist.link(higher, quarter, 0.10)
        board.dist.link(lower, quarter, 0.30)
        decision, probe = valley_status(board, higher, lower)
        self.assertEqual(decision, "slope")
        self.assertIsNone(probe)

    def test_a_midpoint_off_the_callset_path_stays_undecided(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        mid = board.add(line_point(board, higher, lower, 0.5), 70.0)
        board.dist.link(higher, lower, 0.40)
        board.dist.link(higher, mid, 0.30)
        board.dist.link(lower, mid, 0.30)
        decision, probe = valley_status(board, higher, lower)
        self.assertEqual(decision, "undecided")
        self.assertIsNone(probe)

    def test_the_centre_takes_the_finished_grid_and_not_the_child(self):
        board = _Board()
        board.add(_cfg(board, min_pruning=2), 80.0)
        child = board.add(_cfg(board, min_pruning=8), 90.0)
        rung = board.add(_cfg(board, min_pruning=4), 81.0)
        missed = board.add(_cfg(board, min_pruning=6), 79.0)
        self.assertGreater(board.score(child), 80.0)
        chosen = best_grid_point(board, [rung, missed], 80.0)
        self.assertEqual(chosen[0].key, rung.key)
        self.assertIsNone(best_grid_point(board, [missed], 80.0))
        self.assertIsNone(best_grid_point(board, [], 80.0))

    def test_callsets_within_0_10_meet(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        board.dist.link(higher, lower, 0.08)
        decision, probe = valley_status(board, higher, lower)
        self.assertEqual(decision, "meet")
        self.assertIsNone(probe)

    def test_a_downward_bowl_snaps_the_dummy_to_the_crest(self):
        ctx = self.ctx
        def spaced(name):
            rungs = ctx.ladders.rungs[name]
            if ctx.space.specs[name].type != "int" or len(rungs) < 9:
                return False
            units = [ctx.space._to_unit(ctx.space.specs[name], rung) for rung in rungs]
            return max(units) - min(units) > 0.5 and min(b - a for a, b in zip(units, units[1:])) > 0.05

        ints = [name for name in ctx.space.names if spaced(name)]
        along, across = ints[0], ints[1]
        along_rungs = ctx.ladders.rungs[along]
        across_rungs = ctx.ladders.rungs[across]
        sampled = along_rungs[::2][:5]
        crest_rung = along_rungs[along_rungs.index(sampled[1]) + 1]
        held = across_rungs[len(across_rungs) // 2]
        crest_u = ctx.space._to_unit(ctx.space.specs[along], crest_rung)
        held_u = ctx.space._to_unit(ctx.space.specs[across], held)

        def score_of(config):
            vec = ctx.space.encode(config)
            ia, ib = ctx.space.names.index(along), ctx.space.names.index(across)
            return 80.0 - 400.0 * (vec[ia] - crest_u) ** 2 - 400.0 * (vec[ib] - held_u) ** 2

        samples = []
        peak = None
        for value in sampled:
            for side in across_rungs[:: max(1, len(across_rungs) // 4)][:4]:
                config = _cfg(ctx, **{along: value, across: side})
                samples.append((ctx.space.encode(config), score_of(config)))
                if peak is None or samples[-1][1] > peak[0]:
                    peak = (samples[-1][1], config)
        dummy = dummy_config(ctx, peak[1], samples)
        self.assertIsNotNone(dummy)
        self.assertEqual(dummy.params[along], crest_rung)
        self.assertEqual(dummy.params[across], held)

    def test_an_upward_bowl_sends_no_dummy(self):
        ctx = self.ctx
        ints = [name for name in ctx.space.names
                if ctx.space.specs[name].type == "int" and len(ctx.ladders.rungs[name]) >= 5]
        along = ints[0]
        centre = _cfg(ctx, **{along: ctx.ladders.rungs[along][len(ctx.ladders.rungs[along]) // 2]})
        samples = [(ctx.space.encode(centre), 10.0)]
        for rung in ctx.ladders.rungs[along][::2][:6]:
            config = _cfg(ctx, **{along: rung})
            vec = ctx.space.encode(config)
            ia = ctx.space.names.index(along)
            samples.append((vec, 10.0 + 200.0 * (vec[ia] - ctx.space.encode(centre)[ia]) ** 2))
        self.assertIsNone(dummy_config(ctx, centre, samples))


class CoverTest(unittest.IsolatedAsyncioTestCase):
    async def test_free_workers_breed_while_a_small_batch_is_out(self):
        gate = asyncio.get_running_loop().create_future()
        bred = []

        class Client:
            def idle_slots(self):
                return 0 if bred else 900

        class Ctx:
            client = Client()

            def clock(self):
                return 0

            def over(self):
                return False

            def publish(self):
                return None

        async def fake_fill(ctx, cl, n):
            bred.append(n)
            gate.set_result(None)

        with patch("search.separate.fill_free", fake_fill):
            await _cover(Ctx(), [gate], 1e9)
        self.assertEqual(bred, [900])


class _Climber:
    def __init__(self, centre, score, ident):
        self.centre = centre
        self.score = score
        self.id = ident
        self.state = "climbing"
        self.origin = "step A"
        self.why = None

    def move(self, config, score, why):
        self.centre = config
        self.score = score
        self.why = why


class ValleySeatTest(unittest.IsolatedAsyncioTestCase):
    async def test_one_climber_sits_on_the_highest_fingerprint(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        sample = board.add(_cfg(board, min_pruning=8), 69.0)
        board.dist.link(higher, lower, 0.40)
        board.dist.link(higher, sample, 0.20)
        board.dist.link(lower, sample, 0.20)
        board._piles = [Pile(higher, 90.0), Pile(lower, 80.0)]
        seated = await gate(board, [higher, lower])
        self.assertEqual([config.key for config in seated], [higher.key])
        self.assertEqual(board.sent, [])

    async def test_gate_drops_the_lower_fingerprint_when_the_path_does_not_dip(self):
        board = _Board()
        higher = board.add(_cfg(board, min_pruning=2), 90.0)
        lower = board.add(_cfg(board, min_pruning=10), 80.0)
        sample = board.add(_cfg(board, min_pruning=8), 80.4)
        board.dist.link(higher, lower, 0.40)
        board.dist.link(higher, sample, 0.20)
        board.dist.link(lower, sample, 0.20)
        board._piles = [Pile(higher, 90.0), Pile(lower, 80.0)]
        seated = await gate(board, [higher, lower])
        self.assertEqual([config.key for config in seated], [higher.key])
        self.assertEqual(board.sent, [])

    async def test_a_fingerprint_below_the_lowest_climber_is_not_probed(self):
        board = _Board()
        high = board.add(_cfg(board, min_pruning=2), 85.0)
        low = board.add(_cfg(board, min_pruning=4), 70.0)
        waiting = board.add(_cfg(board, min_pruning=10), 60.0)
        board.dist.link(high, waiting, 0.40)
        board.dist.link(low, waiting, 0.40)
        board._piles = [Pile(high, 85.0), Pile(low, 70.0), Pile(waiting, 60.0)]
        climbers = [_Climber(high, 85.0, 0), _Climber(low, 70.0, 1)]
        board.climbers = climbers
        moved = await claim_fingerprint(board, climbers)
        self.assertFalse(moved)
        self.assertEqual(climbers[1].centre.key, low.key)
        self.assertEqual(board.sent, [])

    async def test_a_dip_moves_the_lowest_climber_onto_the_new_fingerprint(self):
        board = _Board()
        high = board.add(_cfg(board, min_pruning=2), 85.0)
        low = board.add(_cfg(board, min_pruning=4), 70.0)
        found = board.add(_cfg(board, min_pruning=10), 88.0)
        sample = board.add(_cfg(board, min_pruning=8), 74.0)
        board.dist.link(found, high, 0.40)
        board.dist.link(found, low, 0.55)
        board.dist.link(found, sample, 0.20)
        board.dist.link(high, sample, 0.20)
        board._piles = [Pile(high, 85.0), Pile(low, 70.0), Pile(found, 88.0)]
        climbers = [_Climber(high, 85.0, 0), _Climber(low, 70.0, 1)]
        board.climbers = climbers
        moved = await claim_fingerprint(board, climbers)
        self.assertTrue(moved)
        self.assertEqual(climbers[1].centre.key, found.key)
        self.assertEqual(climbers[1].origin, "pile 3")
        self.assertEqual(climbers[0].centre.key, high.key)
        self.assertEqual(board.sent, [])
