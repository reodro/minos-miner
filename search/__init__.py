"""Search: strategies that propose configs through EvalClient. Owner: search.

Contract with main: main asks `strategies(name, ranges, workers)` for the lineup, builds the EvalClient (round,
executor, scorer, answer set, budget, deadline, clock) and runs the strategies as asyncio tasks until the deadline; a
strategy only proposes configs and reads results - it never submits to the platform (main keeps the best config
submitted). The lineup lives in search/lineup.py: adding or retuning strategies never needs a change in main.
"""
from search.base import RandomSampler, Strategy, random_stream  # noqa: F401
from search.client import EvalClient, Outcome  # noqa: F401
from search.lineup import STRATEGY_NAMES, strategies  # noqa: F401
from search.rf_ucb import rf_ucb  # noqa: F401
