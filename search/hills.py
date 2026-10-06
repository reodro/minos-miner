"""Which hill a config is on: one rule for the whole engine. Owner: search.

Two configs are on different hills when they differ on MIN_FAR or more live settings by more than 10% of the
setting's range (enums and booleans: any difference). Ownership, adoption, the scout's exclusion zone, merging and
founding all ask this module, so they can never disagree about where a hill ends.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

from common.types import Config
from search.axes import far_axes
from search.space import Space

MIN_FAR = 4


class Hills:
    def __init__(self, space: Space, names: Callable[[], Sequence[str]]):
        """`names`: the live settings right now (parked dead settings do not separate hills)."""
        self.space, self.names = space, names

    def distance(self, a: Config, b: Config) -> int:
        return far_axes(self.space, self.names(), a, b)

    def same(self, a: Config, b: Config) -> bool:
        return self.distance(a, b) < MIN_FAR

    def nearest(self, c: Config, centres: Dict[int, Config]) -> Tuple[Optional[int], int]:
        best: Tuple[Optional[int], int] = (None, 10 ** 9)
        for i, z in centres.items():
            d = self.distance(c, z)
            if d < best[1]:
                best = (i, d)
        return best

    def owner(self, c: Config, centres: Dict[int, Config]) -> Optional[int]:
        """The hill whose centre `c` is on (the nearest one, when several), or None when it is on no known hill."""
        i, d = self.nearest(c, centres)
        return i if i is not None and d < MIN_FAR else None

    def cluster(self, ranked: Sequence[Tuple[Config, float]], k: int, band: float,
                min_far: int = MIN_FAR) -> List[Tuple[Config, float]]:
        """The best config of up to `k` distinct hills among `ranked` (best first), none more than `band` below the top:
        a config starts a hill when it is at least `min_far` settings from every hill already chosen."""
        out: List[Tuple[Config, float]] = []
        if not ranked:
            return out
        top = ranked[0][1]
        for c, s in ranked:
            if len(out) >= k or s < top - band:
                break
            if all(self.distance(c, o) >= min_far for o, _ in out):
                out.append((c, s))
        return out
