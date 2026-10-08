"""Search strategies. Owner: search.

A strategy is an async function `strategy(client)` that runs until the client is exhausted or the round ends
(the round loop cancels it). It can stream, wait on barriers, or run sequential chains - see search/client.py.
RF_UCB, QGA, ladder, restarts, local polish ... are each one strategy; several run at once.
"""
from __future__ import annotations

import asyncio
import random
from typing import Awaitable, Callable, List

from common.params import ParamRanges
from common.types import Config
from search.v422.space import PINNED
from search.v422.client import EvalClient

Strategy = Callable[[EvalClient], Awaitable[None]]


class RandomSampler:
    def __init__(self, ranges: ParamRanges, tool: str = "gatk", seed: int = 0):
        self.ranges, self.specs, self.tool, self.rng = ranges, ranges.tools[tool], tool, random.Random(seed)

    def _value(self, s):
        if s.allowed:
            return self.rng.choice(s.allowed)
        if s.type == "bool":
            return self.rng.random() < 0.5
        if s.type == "int":
            return self.rng.randint(int(s.min), int(s.max))
        return round(self.rng.uniform(float(s.min), float(s.max)), 6)

    def sample(self, n: int) -> List[Config]:
        """n configs that pass ParamRanges.errors, including the cross-parameter constraints GATK enforces."""
        out = []
        while len(out) < n:
            c = Config(self.tool, {k: PINNED.get(self.tool, {}).get(k, self._value(s)) for k, s in self.specs.items()})
            if not self.ranges.errors(c):
                out.append(c)
        return out


def random_stream(ranges: ParamRanges, in_flight: int = 8, seed: int = 0, tag: str = "random") -> Strategy:
    """Baseline every real strategy must beat: keeps `in_flight` uniform samples running, refilling as each ends."""
    sampler = RandomSampler(ranges, seed=seed)

    async def run(client: EvalClient) -> None:
        running = set()
        while True:
            while len(running) < in_flight and not client.exhausted:
                running.add(client.propose(sampler.sample(1)[0], tag=tag))
            if not running:
                return
            _, running = await asyncio.wait(running, return_when=asyncio.FIRST_COMPLETED)

    return run
