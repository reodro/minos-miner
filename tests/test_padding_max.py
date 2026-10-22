"""space.padding_max (2026-10-19): the search's ceiling for assembly_region_padding, inside the published range."""
import unittest
import unittest.mock

from common.tuning import TUNING
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges

PAD = "assembly_region_padding"


class PaddingMaxTest(unittest.TestCase):
    def space(self, v):
        with unittest.mock.patch.dict(TUNING.values, {"space.padding_max": v}):
            return Space(ranges())

    def test_off_keeps_the_published_maximum(self):
        sp = self.space(0)
        self.assertEqual(Ladders(sp).rungs[PAD][-1], ranges().tools["gatk"][PAD].max)

    def test_the_ladder_stops_at_the_ceiling_and_keeps_the_floor(self):
        rungs = Ladders(self.space(60)).rungs[PAD]
        self.assertEqual((rungs[0], rungs[-1]), (5, 60))

    def test_a_seed_above_the_ceiling_is_clamped(self):
        c = self.space(60)._assemble({PAD: 200})
        self.assertEqual(c.params[PAD], 60)


class FloorAndAptTest(unittest.TestCase):
    def space(self, **v):
        with unittest.mock.patch.dict(TUNING.values, v):
            return Space(ranges())

    def test_padding_floor_raises_the_ladder_and_seeds(self):
        sp = self.space(**{"space.padding_min": 12, "space.padding_max": 60})
        rungs = Ladders(sp).rungs[PAD]
        self.assertEqual((rungs[0], rungs[-1], len(rungs)), (12, 60, 49))
        self.assertEqual(sp._assemble({PAD: 8}).params[PAD], 12)

    def test_apt_ceiling(self):
        apt = "active_probability_threshold"
        sp = self.space(**{"space.apt_max": 0.023})
        self.assertLessEqual(Ladders(sp).rungs[apt][-1], 0.023)
        self.assertEqual(sp._assemble({apt: 0.04}).params[apt], 0.023)

    def test_off_changes_nothing(self):
        sp = self.space(**{"space.padding_min": 0, "space.apt_max": 0.0})
        self.assertEqual(Ladders(sp).rungs[PAD][0], 5)
        self.assertEqual(Ladders(sp).rungs["active_probability_threshold"][-1], ranges().tools["gatk"]["active_probability_threshold"].max)


if __name__ == "__main__":
    unittest.main()
