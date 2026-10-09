"""Search: the engine main runs proposes configs through EvalClient: search/v10 (RF_UCB) or search/v9. Owner: search.

Contract with main: main asks `searcher(ranges, workers)` for the round's search - one async function of the round's
EvalClient - builds the EvalClient (round, executor, scorer, answer set, budget, deadline, clock) and runs it until the
deadline. The search only proposes configs and reads results; main submits the best config, or the one the search
names in `client.pick_key`. Engine history: search/README.md and git.
"""
from search.base import RandomSampler, Strategy, random_stream  # noqa: F401
from search.client import EvalClient, Outcome  # noqa: F401


def searcher(ranges, workers: int) -> Strategy:
    """The search main runs every round, in replays and live alike: MINOS_ENGINE=v9 or v10 (RF_UCB, the default)."""
    import os
    if os.environ.get("MINOS_ENGINE", "v10") == "v9":
        from search.v9.engine import searcher as engine
    else:
        from search.v10.engine import searcher as engine
    return engine(ranges, workers)
