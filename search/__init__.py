"""Search: the engine that proposes configs through EvalClient (engine.py: asynchronous basin climbing). Owner: search.

Contract with main: main asks `searcher(ranges, workers)` for the round's search - one async function of the round's
EvalClient - builds the EvalClient (round, executor, scorer, answer set, budget, deadline, clock) and runs it until the
deadline. The search only proposes configs and reads results - it never submits to the platform (main keeps the best
config submitted). What the search is (one algorithm, or several combined) is decided here alone: changing it never
needs a change in main.
"""
from search.base import RandomSampler, Strategy, random_stream  # noqa: F401
from search.client import EvalClient, Outcome  # noqa: F401
from search.basin import basin_race  # noqa: F401
from search.pattern import pattern, pattern_regions  # noqa: F401
from search.regions import gradient_regions  # noqa: F401
from search.rf_ucb import rf_ucb  # noqa: F401
from search.tpe import tpe  # noqa: F401


def searcher(ranges, workers: int) -> Strategy:
    """The search main runs every round, in replays and live alike. Owner: search - changing it never touches main.
    MINOS_ENGINE selects the engine: "v3" (engine3: isolated iterated local search with racing, the default) or "v2"
    (engine: asynchronous basin climbing), so the two can be compared on the same round."""
    import os
    if os.environ.get("MINOS_ENGINE", "v3") == "v2":
        from search.engine import searcher as engine
        return engine(ranges, workers)
    from search.engine3 import searcher as engine3
    return engine3(ranges, workers)
