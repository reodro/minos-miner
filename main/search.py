"""Search strategies. Owner: main. RF_UCB, QGA, ladder, restarts etc. implement Searcher."""
from __future__ import annotations

import random
from typing import Iterable, List, Protocol

from common.params import ParamRanges
from common.types import Config, ScoreResult


class Searcher(Protocol):
    def propose(self, n: int) -> List[Config]:
        """Next n configs to evaluate. Must be valid per ParamRanges."""

    def observe(self, results: Iterable[ScoreResult]) -> None:
        """New scores for configs this searcher proposed (or seeded)."""

    def rescored(self, results: Iterable[ScoreResult]) -> None:
        """The answer set changed: these replace every earlier score."""


class RandomSearcher:
    """Uniform sampling inside the published ranges. The baseline every real searcher must beat."""

    def __init__(self, ranges: ParamRanges, tool: str = "gatk", seed: int = 0):
        self.specs = ranges.tools[tool]
        self.tool = tool
        self.rng = random.Random(seed)

    def _sample(self, s):
        if s.allowed:
            return self.rng.choice(s.allowed)
        if s.type == "bool":
            return self.rng.random() < 0.5
        if s.type == "int":
            return self.rng.randint(int(s.min), int(s.max))
        return round(self.rng.uniform(float(s.min), float(s.max)), 6)

    def propose(self, n: int) -> List[Config]:
        return [Config(self.tool, {k: self._sample(s) for k, s in self.specs.items()}) for _ in range(n)]

    def observe(self, results: Iterable[ScoreResult]) -> None:
        pass

    def rescored(self, results: Iterable[ScoreResult]) -> None:
        pass
