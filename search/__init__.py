"""Search: strategies that propose configs through EvalClient. Owner: search.

Contract with main: main builds the EvalClient (round, executor, scorer, answer set, budget) and runs strategies
as asyncio tasks until the deadline; a strategy only proposes configs and reads results - it never submits to the
platform (main keeps the best config submitted).
"""
from search.base import RandomSampler, Strategy, random_stream  # noqa: F401
from search.client import EvalClient, Outcome  # noqa: F401
from search.rf_ucb import rf_ucb  # noqa: F401
