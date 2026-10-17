"""The GA model filter (search/model_filter.py): keeps the children a per-class model of the round ranks highest."""
import random
import unittest
import unittest.mock

from common.types import CLASSES, ClassCounts, ScoreResult
from search import model_filter as mf
from search.space import Space
from search.steps import Ladders
from tests.helpers import ranges

MMQ = "min_mapping_quality_score"


def result(c):
    q = 1.0 - abs(int(c.params.get(MMQ)) - 45) / 60.0          # best at min_mapping_quality_score 45
    counts = {k: ClassCounts(tp=int(800 * q) + 100, fp=int(60 * (1 - q)), fn=900 - int(800 * q)) for k in CLASSES}
    r = ScoreResult("r", c.key, 1, counts, region_fp=int(400 * (1 - q)) + 50, n_targets=500)
    from common.scoring import finalize
    return finalize(r)


class FakeCtx:
    def __init__(self, space):
        self.space, self.ladders = space, Ladders(space)
        self.said = []

    def say(self, *a):
        self.said.append(a)


class ModelFilterTest(unittest.TestCase):
    def test_keeps_the_children_the_model_ranks_best(self):
        space = Space(ranges())
        rng = random.Random(0)
        scored = []
        for _ in range(300):
            c = space._assemble({MMQ: rng.randint(0, 60), "assembly_region_padding": rng.randint(10, 200)})
            scored.append((c, result(c)))
        ctx = FakeCtx(space)
        m = mf.model_of(ctx)
        with unittest.mock.patch.object(mf, "FIRST_FIT", 200):
            self.assertTrue(m.maybe_fit(scored))
        kids = [space._assemble({MMQ: v, "assembly_region_padding": 80}) for v in (5, 20, 44, 46, 58, 30)]
        mf.scored_results = lambda ctx: scored                    # the context's scored configs
        got = mf.keep_best(ctx, kids, 2)
        self.assertEqual(sorted(int(c.params[MMQ]) for c in got), [44, 46])

    def test_no_model_yet_keeps_the_first_n(self):
        space = Space(ranges())
        ctx = FakeCtx(space)
        mf.scored_results = lambda ctx: []
        kids = [space._assemble({MMQ: v}) for v in (5, 20, 44)]
        self.assertEqual(mf.keep_best(ctx, kids, 2), kids[:2])


if __name__ == "__main__":
    unittest.main()
