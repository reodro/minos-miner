import json
import random
import unittest

from common.store import ResultStore
from common.types import ClassCounts, Config
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search import pattern, random_stream
from search.pattern import (bam_max_quality, block_moves, block_swaps, confidence_steps, core_steps,
                            cross_parents, expensive_distance, pattern_regions, region_centers, worst_class)
from search.space import Space
from tests.helpers import ranges
from tests.test_parts import ROUND, StubProvider


class PatternTest(unittest.TestCase):
    def run_one(self, strat, seed, max_evals=200):
        r = ranges()
        runner = FakeRunner(r, seed=seed)
        out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner), strat, DryRunSubmitter(),
                        r, LoopSettings(max_evals=max_evals), baselines=[])
        return out, runner

    def test_beats_random_and_respects_budget(self):
        r = ranges()
        best = {"random": 0.0, "pattern": 0.0}
        for seed in range(3):
            for name, strat in (("random", random_stream(r, in_flight=8, seed=seed)),
                                ("pattern", pattern(r, in_flight=8, seed=seed, max_bq=None))):
                out, _ = self.run_one(strat, seed)
                self.assertLessEqual(out.evaluated + out.failed, 200 + 8 * 7)   # in-flight families may finish after
                best[name] += out.best.score
        self.assertGreater(best["pattern"], best["random"])

    def test_cap_keeps_base_quality_cutoffs_scorable(self):
        out, runner = self.run_one(pattern(ranges(), in_flight=8, seed=1, max_bq=30), seed=1, max_evals=120)
        self.assertTrue(runner.configs)
        for c in runner.configs.values():
            self.assertLessEqual(c.params["min_base_quality_score"], 30)
            self.assertLessEqual(c.params["base_quality_score_threshold"], 30)

    def test_crossover_mixes_two_scored_configs_and_obeys_the_cap(self):
        space = Space(ranges())
        low = space._assemble({"min_pruning": 2, "assembly_region_padding": 50, "min_base_quality_score": 10})
        high = space._assemble({"min_pruning": 8, "assembly_region_padding": 400, "min_base_quality_score": 40})
        kids = cross_parents(space, random.Random(0), [low, high], 80, cap_at=30)
        self.assertTrue(kids)
        mixed = False
        for c in kids:
            self.assertLessEqual(c.params["min_base_quality_score"], 30)
            self.assertIn(c.params["min_pruning"], (2, 8))
            self.assertIn(c.params["assembly_region_padding"], (50, 400))
            pair = (c.params["min_pruning"], c.params["assembly_region_padding"])
            mixed = mixed or pair in ((2, 400), (8, 50))
        self.assertTrue(mixed)

    def test_block_move_changes_one_stage_and_obeys_the_cap(self):
        space = Space(ranges())
        base = space._assemble({})
        kids = block_moves(space, random.Random(0), base, 12, cap_at=30)
        stages = (
            {"min_base_quality_score", "min_mapping_quality_score", "max_reads_per_alignment_start"},
            {"active_probability_threshold", "min_assembly_region_size", "max_assembly_region_size",
             "assembly_region_padding"},
            {"min_pruning", "min_dangling_branch_length", "recover_all_dangling_branches",
             "max_num_haplotypes_in_population"},
        )
        self.assertTrue(kids)
        for c in kids:
            self.assertLessEqual(c.params["min_base_quality_score"], 30)
            changed = {n for n in space.expensive
                       if not space._same(space.specs[n], c.params[n], base.params[n])}
            self.assertGreaterEqual(len(changed), 2)
            self.assertEqual(sum(changed <= stage for stage in stages), 1)

    def test_block_swap_copies_one_stage(self):
        space = Space(ranges())
        a = space._assemble({"min_pruning": 2, "min_dangling_branch_length": 4,
                             "max_num_haplotypes_in_population": 128, "assembly_region_padding": 50})
        b = space._assemble({"min_pruning": 8, "min_dangling_branch_length": 20,
                             "max_num_haplotypes_in_population": 16, "assembly_region_padding": 400})
        kids = block_swaps(space, random.Random(1), [a, b], 40, cap_at=30)
        swapped = [c for c in kids if c.params["min_pruning"] == 8 and c.params["assembly_region_padding"] == 50]
        self.assertTrue(swapped)
        for c in swapped:
            self.assertEqual(c.params["min_dangling_branch_length"], 20)
            self.assertEqual(c.params["max_num_haplotypes_in_population"], 16)
            self.assertEqual(c.params["min_base_quality_score"], a.params["min_base_quality_score"])

    def test_core_step_copies_one_indel_setting_and_keeps_the_reads(self):
        space = Space(ranges())
        high = space._assemble({"min_base_quality_score": 16, "min_mapping_quality_score": 46,
                                "max_reads_per_alignment_start": 116, "min_pruning": 2,
                                "min_dangling_branch_length": 4, "assembly_region_padding": 100,
                                "max_num_haplotypes_in_population": 128,
                                "pair_hmm_gap_continuation_penalty": 10})
        low = space._assemble({"min_base_quality_score": 30, "min_mapping_quality_score": 60,
                               "max_reads_per_alignment_start": 25, "min_pruning": 4,
                               "min_dangling_branch_length": 20, "assembly_region_padding": 50,
                               "max_num_haplotypes_in_population": 64,
                               "pair_hmm_gap_continuation_penalty": 7})
        kids = core_steps(space, high, low, cap_at=30)
        self.assertTrue(kids)
        changed_names = []
        for c in kids:
            self.assertEqual(c.params["min_base_quality_score"], 16)
            self.assertEqual(c.params["min_mapping_quality_score"], 46)
            self.assertEqual(c.params["max_reads_per_alignment_start"], 116)
            changed = [n for n in space.names if c.params[n] != high.params[n]]
            self.assertEqual(len(changed), 1)
            changed_names.append(changed[0])
        self.assertIn("min_dangling_branch_length", changed_names)
        self.assertIn("min_pruning", changed_names)
        self.assertIn("assembly_region_padding", changed_names)
        self.assertIn("pair_hmm_gap_continuation_penalty", changed_names)
        self.assertNotIn("min_base_quality_score", changed_names)
        self.assertEqual(core_steps(space, high, high, cap_at=30), [])

    def test_overcalled_long_indel_raises_the_calling_threshold_only(self):
        space = Space(ranges())
        base = space._assemble({"standard_min_confidence_threshold_for_calling": 80.0,
                                "min_base_quality_score": 14, "pcr_indel_model": "NONE",
                                "min_mapping_quality_score": 54, "max_reads_per_alignment_start": 50})
        counts = {
            "indel_1bp": ClassCounts(1572, 159, 96),
            "indel_2_3bp": ClassCounts(699, 148, 57),
            "indel_4_7bp": ClassCounts(455, 87, 36),
            "indel_8bp": ClassCounts(184, 116, 34),
            "snp_het": ClassCounts(9938, 339, 249),
            "snp_hom": ClassCounts(4966, 6, 40),
        }
        hole = worst_class(counts)
        self.assertEqual(hole[1], "indel_8bp")
        kids = confidence_steps(space, base, counts, cap_at=30)
        self.assertTrue(kids)
        for c in kids:
            self.assertGreater(c.params["standard_min_confidence_threshold_for_calling"], 80.0)
            self.assertEqual(c.params["min_base_quality_score"], 14)
            self.assertEqual(c.params["pcr_indel_model"], "NONE")
            self.assertEqual(c.params["min_mapping_quality_score"], 54)
            changed = [n for n in space.names if c.params[n] != base.params[n]]
            self.assertEqual(changed, ["standard_min_confidence_threshold_for_calling"])
        # A missed class, precision above recall, walks the threshold down.
        missed = dict(counts)
        missed["indel_8bp"] = ClassCounts(100, 10, 120)
        down = confidence_steps(space, base, missed, cap_at=30)
        self.assertTrue(down)
        self.assertTrue(all(c.params["standard_min_confidence_threshold_for_calling"] < 80.0 for c in down))
        # A SNP hole does not spend the reserved steps.
        snp = {name: ClassCounts(1000, 0, 0) for name in counts}
        snp["snp_het"] = ClassCounts(1000, 400, 400)
        self.assertEqual(confidence_steps(space, base, snp, cap_at=30), [])

    def test_region_centers_drop_a_close_config_and_a_weak_one(self):
        space = Space(ranges())
        best = space._assemble({})
        far = space._assemble({"min_pruning": 8, "min_base_quality_score": 30,
                               "max_reads_per_alignment_start": 25})
        close = space._assemble({"min_dangling_branch_length": 6})
        weak = space._assemble({"min_pruning": 10, "min_mapping_quality_score": 60,
                                "max_reads_per_alignment_start": 25})
        centers = region_centers(space, [(90.0, best), (89.0, far), (88.9, close), (80.0, weak)],
                                 4, 0.35, band=1.5)
        self.assertEqual([c.key for c in centers], [best.key, far.key])
        self.assertGreaterEqual(expensive_distance(space, best, far), 0.35)

    def test_each_region_proposes_only_inside_its_cell(self):
        r = ranges()
        store = ResultStore()
        runner = FakeRunner(r, seed=0)
        strat = pattern_regions(r, in_flight=8, seed=0, max_bq=None, n_regions=2, separation=0.2,
                                scout_min=8, band=50.0)
        out = run_round(ROUND, StubProvider(), SimExecutor(runner, 8), FakeScorer(runner),
                        strat, DryRunSubmitter(), r, LoopSettings(max_evals=160), baselines=[], store=store)
        self.assertIsNotNone(out.best)
        space = Space(r)
        centers = strat.centers
        self.assertGreaterEqual(len(centers), 2)
        lanes = {}
        for tag, raw, _score in store.db.execute("SELECT tag, config, score FROM results"):
            if tag.startswith("region"):
                lanes.setdefault(tag, []).append(Config.from_dict(json.loads(raw)))
        self.assertIn("region0", lanes)
        self.assertIn("region1", lanes)

        def owner(config):
            distances = [expensive_distance(space, config, center) for center in centers]
            return min(range(len(distances)), key=lambda i: (distances[i], i))

        for tag, configs in lanes.items():
            ident = int(tag.replace("region", ""))
            for config in configs:
                self.assertEqual(owner(config), ident)

    def test_auto_cap_without_a_bam_searches_uncapped(self):
        self.assertIsNone(bam_max_quality("/nonexistent.bam", "chr20:1-1000"))
        out, _ = self.run_one(pattern(ranges(), in_flight=8, seed=2), seed=2, max_evals=60)
        self.assertIsNotNone(out.best)


if __name__ == "__main__":
    unittest.main()
