import unittest

from common.types import Config
from evaluator import FakeRunner, FakeScorer, SimExecutor
from main import DryRunSubmitter, LoopSettings, run_round
from search.basin import (OFFICIAL_SEED, cheap_family, signature, structural_names)
from search.basin import basin_race
from search.space import GENOTYPING_PARAMS, PINNED, Space
from tests.helpers import ranges
from tests.test_parts import ROUND, StubProvider


def space() -> Space:
    pinned = dict(PINNED["gatk"])
    pinned["pcr_indel_model"] = "NONE"
    return Space(ranges(), "gatk", pinned=pinned)


class BasinTest(unittest.TestCase):
    def test_signature_ignores_genotyping_and_splits_a_real_structural_step(self):
        sp = space()
        base = sp._assemble(dict(OFFICIAL_SEED))
        self.assertIsNotNone(base)
        nudged = sp._assemble({**base.params, "standard_min_confidence_threshold_for_calling": 30.0,
                               "indel_heterozygosity": 0.0005, "contamination_fraction_to_filter": 0.0,
                               "max_alternate_alleles": 2})
        self.assertEqual(signature(sp, base), signature(sp, nudged))
        far = sp._assemble({**base.params, "min_pruning": 10, "assembly_region_padding": 400,
                            "max_reads_per_alignment_start": 300})
        self.assertNotEqual(signature(sp, base), signature(sp, far))
        self.assertNotIn("standard_min_confidence_threshold_for_calling", structural_names(sp))

    def test_family_changes_genotyping_only_and_keeps_pcr_none(self):
        sp = space()
        base = sp._assemble(dict(OFFICIAL_SEED))
        kids = cheap_family(sp, base, limit=5)
        self.assertGreaterEqual(len(kids), 4)
        self.assertLessEqual(len(kids), 5)
        for kid in kids:
            self.assertEqual(kid.params["pcr_indel_model"], "NONE")
            changed = [name for name in sp.names if kid.params.get(name) != base.params.get(name)]
            self.assertTrue(changed)
            self.assertTrue(set(changed) <= set(GENOTYPING_PARAMS))
            self.assertEqual(signature(sp, kid), signature(sp, base))

    def test_a_short_race_submits_and_stays_inside_the_budget(self):
        r = ranges()
        runner = FakeRunner(r, seed=0)
        sub = DryRunSubmitter()
        out = run_round(ROUND, StubProvider(), SimExecutor(runner, 4), FakeScorer(runner),
                        basin_race(r, in_flight=4), sub, r,
                        LoopSettings(max_evals=40), baselines=[])
        self.assertGreater(out.evaluated, 0)
        self.assertLessEqual(out.evaluated + out.failed, 40 + 4 * 6)
        self.assertTrue(sub.sent)
        # The submitted key is one the strategy actually ran.
        self.assertIn(sub.sent[-1][1], runner.configs)

    def test_basin_race_builds_a_strategy(self):
        strat = basin_race(ranges(), in_flight=8)
        self.assertEqual(strat.__name__, "run")
