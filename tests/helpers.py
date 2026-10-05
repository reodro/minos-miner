import asyncio
import os

from common.params import ParamRanges

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


def ranges() -> ParamRanges:
    return ParamRanges.from_file(os.path.join(FIXTURES, "parameter_ranges.json"))


def together(*searches):
    """Several search functions run side by side as one searcher (main runs one search; tests of search behaviour
    sometimes combine a few)."""
    async def run(client):
        await asyncio.gather(*(s(client) for s in searches))
    return run
