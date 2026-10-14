"""A value that has already tied an exact score is the same callset and is not sent again.

A setting value becomes tied only when two scored configs share that exact score and that value is the one
setting still unexplained between them. A later config is refused when every value it changes is one of those
tied values on the config it is compared with. An untried value of the same setting is still sent: on
20261009T070600Z treating the whole ladder as tied refused climbs of 0.05 to 0.11.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Set, Tuple


class TiedCallsets:
    """Exemplars of scored configs, and the values that tied each one."""

    def __init__(self, names: Sequence[str]):
        self.names = list(names)
        self.seen: Set[str] = set()
        self._kept: List[list] = []
        self._by_score: Dict[float, List[int]] = defaultdict(list)
        self._buckets: Dict[frozenset, Dict[tuple, List[int]]] = defaultdict(lambda: defaultdict(list))
        self._where: List[Optional[Tuple[frozenset, tuple]]] = []
        self._active: List[frozenset] = []

    def note(self, key: str, params: dict, score: float) -> None:
        """Record one scored config. A one-setting tie, after values already tied, marks that value."""
        if key in self.seen:
            return
        self.seen.add(key)
        v = tuple(params.get(n) for n in self.names)
        i = len(self._kept)
        self._kept.append([v, score, {}])
        self._by_score[score].append(i)
        self._where.append(None)
        moved: List[int] = []
        for j in self._by_score[score][:-1]:
            u, _, extra_j = self._kept[j]
            extra_i = self._kept[i][2]
            unexplained = [d for d in range(len(self.names))
                           if v[d] != u[d]
                           and v[d] not in extra_j.get(d, ())
                           and u[d] not in extra_i.get(d, ())]
            if len(unexplained) == 1:
                bit = unexplained[0]
                extra_j.setdefault(bit, set()).add(v[bit])
                extra_i.setdefault(bit, set()).add(u[bit])
                moved.append(j)
        if self._kept[i][2] or moved:
            self._put(i)
            for j in moved:
                self._put(j)

    def holds(self, params: dict) -> bool:
        """True when `params` differs from a scored config only by values that already tied it."""
        if not self._active:
            return False
        v = tuple(params.get(n) for n in self.names)
        n = len(self.names)
        for key in self._active:
            proj = tuple(v[i] for i in range(n) if i not in key)
            for i in self._buckets[key].get(proj, ()):
                extra = self._kept[i][2]
                exemplar = self._kept[i][0]
                diffs = [d for d in range(n) if v[d] != exemplar[d]]
                if diffs and all(v[d] in extra.get(d, ()) for d in diffs):
                    return True
        return False

    def _put(self, i: int) -> None:
        old = self._where[i]
        if old is not None:
            lst = self._buckets[old[0]][old[1]]
            if i in lst:
                lst.remove(i)
        _, _, extra = self._kept[i]
        if not extra:
            self._where[i] = None
            return
        key = frozenset(extra)
        proj = tuple(self._kept[i][0][d] for d in range(len(self.names)) if d not in key)
        self._buckets[key][proj].append(i)
        self._where[i] = (key, proj)
        if key not in self._active:
            self._active.append(key)
