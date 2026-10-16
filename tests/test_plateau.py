"""search/plateau.py: per-setting intervals from one-setting slices, the second score, and its use in the plateau
probe (centre of the flat; the rung toward the interval first, the opposite rung after)."""
import unittest

from search import plateau
from search.climb import next_probes
from search.space import Space
from tests.helpers import ranges
from tests.test_engine import FakeCtx

MP = "min_pruning"                        # ladder 2..10: index = value - 2
MDBL = "min_dangling_branch_length"
HET = "heterozygosity"


class PlateauTest(unittest.TestCase):
    def setUp(self):
        self.space = Space(ranges())
        self.base = self.space._assemble({MP: 5})
        # a slice of MP on another background (heterozygosity changed): the score peaks at MP 8 (index 6)
        bg = self.space._assemble({HET: 0.005})
        self.slice = [self.space._assemble({**bg.params, MP: v}) for v in range(2, 11)]
        self.tie = self.space._assemble({**self.base.params, MDBL: self.base.params[MDBL] - 1})
        slice_keys = {c.key for c in self.slice}

        def f(c):
            if c.key in slice_keys:
                return 90.0 - 0.2 * abs(c.params[MP] - 8)
            if c.key in (self.base.key, self.tie.key):
                return 95.0                                   # the flat
            return 80.0
        self.ctx = FakeCtx(self.space, f)
        self.ctx.sent = [self.slice, [self.base, self.tie]]

    def test_the_interval_is_the_peak_of_the_kept_slices(self):
        iv = plateau.intervals(self.ctx)
        self.assertEqual(iv[MP], {6})
        self.assertNotIn(MDBL, iv)                            # no slice of MDBL: silent, no direction invented

    def test_a_flat_slice_is_not_kept(self):
        flat = [self.space._assemble({**self.slice[0].params, MP: v}) for v in range(2, 11)]
        ctx = FakeCtx(self.space, lambda c: 90.0 + 0.01 * c.params[MP])   # spread 0.08 <= 0.10
        ctx.sent = [flat]
        self.assertNotIn(MP, plateau.intervals(ctx))

    def test_the_second_score_counts_rungs_outside_the_interval(self):
        n = len(self.ctx.ladders.rungs[MP])
        self.assertEqual(plateau.second_score(self.ctx, self.space._assemble({MP: 8})), 0.0)
        self.assertAlmostEqual(plateau.second_score(self.ctx, self.space._assemble({MP: 4})), -4 / n)

    def test_the_flat_centres_on_the_tied_config_nearest_the_intervals(self):
        near = self.space._assemble({**self.base.params, MP: 7})
        ctx = FakeCtx(self.space, lambda c: 95.0 if c.key in (self.base.key, near.key) else
                      90.0 - 0.2 * abs(c.params[MP] - 8))
        ctx.sent = [self.slice, [self.base, near]]
        self.assertEqual(plateau.centre_of_flat(ctx, 95.0, self.base).key, near.key)

    def test_the_probe_sends_the_rung_toward_the_interval_first_and_the_opposite_rung_after(self):
        out = next_probes(self.ctx, self.base, 95.0)
        mp = [(i, c.params[MP]) for i, c in enumerate(out) if c.params[MP] != 5]
        self.assertEqual(sorted(v for _, v in mp), [4, 6])
        toward = next(i for i, v in mp if v == 6)
        away = next(i for i, v in mp if v == 4)
        self.assertLess(toward, away)
        self.assertEqual(away, len(out) - 1)                  # after every first rung of every locked setting

    def test_switched_off_it_is_the_old_probe(self):
        from unittest import mock
        with mock.patch("search.plateau.enabled", return_value=False):
            out = next_probes(self.ctx, self.base, 95.0)
            self.assertEqual(plateau.centre_of_flat(self.ctx, 95.0, self.base).key, self.base.key)
        self.assertEqual(len([c for c in out if c.params[MP] != 5]), 1)    # one nearest rung, no direction


if __name__ == "__main__":
    unittest.main()
