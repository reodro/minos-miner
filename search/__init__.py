"""Search: the engine main runs proposes configs through EvalClient: search/v10 (RF_UCB), search/v11 (several hills) or search/v9. Owner: search.

Contract with main: main asks `searcher(ranges, workers)` for the round's search - one async function of the round's
EvalClient - builds the EvalClient (round, executor, scorer, answer set, budget, deadline, clock) and runs it until the
deadline. The search only proposes configs and reads results; main submits the best config, or the one the search
names in `client.pick_key`. Engine history: search/README.md and git.
"""
from search.base import RandomSampler, Strategy, random_stream  # noqa: F401
from search.client import EvalClient, Outcome  # noqa: F401


def searcher(ranges, workers: int) -> Strategy:
    """The search main runs every round, in replays and live alike: MINOS_ENGINE=v9, v10 (RF_UCB, the default) v11 (RF_UCB over several hills) or v12 (trust regions)."""
    import os
    which = os.environ.get("MINOS_ENGINE", "v10")
    if which == "v9":
        from search.v9.engine import searcher as engine
    elif which == "v11":
        from search.v11.engine import searcher as engine
    elif which == "v12":
        from search.v12.engine import searcher as engine
    else:
        from search.v10.engine import searcher as engine
    return engine(ranges, workers)
