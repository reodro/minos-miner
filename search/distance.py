"""How different two configs are: by what they call, not by their settings (search/DESIGN.md, A4).

callset distance = sum over the score classes of weight x (1 - Jaccard of the two callsets in that class), divided by
the sum of the weights, with the score's own difficulty weights. Identical callsets are 0 whatever their settings;
two hills that call different variants are far apart even when their settings are close. The Jaccard per class comes
from the bottom-k sketches the evaluator attaches to every CallSet (evaluator/sketch.py, decoded here because search
does not import evaluator).

When a config has no sketch (a failed sketch, a simulator), the distance falls back to the mean normalised setting
difference (search.portfolio.distance's measure). The two measures are on different scales, so `far()` takes each
against its own threshold.
"""
from __future__ import annotations

import base64
import struct
from typing import Dict, List, Optional, Tuple

from common.scoring import DIFFICULTY_WEIGHTS
from common.types import FULL, Config

WEIGHTS: Dict[str, float] = dict(DIFFICULTY_WEIGHTS)


def _unpack(text: str) -> List[int]:
    raw = base64.b64decode(text)
    return list(struct.unpack(f">{len(raw) // 4}I", raw))


def jaccard(a: List[int], b: List[int], k: int) -> float:
    """Jaccard similarity of two sets from their bottom-k sketches (exact when the union holds k or fewer)."""
    if not a and not b:
        return 1.0
    sa, sb = set(a), set(b)
    union = sorted(sa | sb)[:k]
    both = sum(1 for x in union if x in sa and x in sb)
    return both / len(union)


def callset_distance(a: dict, b: dict) -> float:
    k = min(int(a.get("k", 128)), int(b.get("k", 128)))
    num = den = 0.0
    for cls, w in WEIGHTS.items():
        ha, hb = a["h"].get(cls, []), b["h"].get(cls, [])
        if not ha and not hb:
            continue
        num += w * (1.0 - jaccard(ha, hb, k))
        den += w
    return num / den if den else 0.0


def setting_distance(space, a: Config, b: Config) -> float:
    """Mean normalised difference over the searched settings (enums and bools: 0 or 1)."""
    names = space.names
    if not names:
        return 0.0
    total = 0.0
    for n in names:
        s = space.specs[n]
        x, y = a.params.get(n, s.default), b.params.get(n, s.default)
        if s.allowed or s.type == "bool":
            total += float(x != y)
        else:
            total += abs(space._to_unit(s, x) - space._to_unit(s, y))
    return total / len(names)


class Distance:
    """Distances between configs the client has evaluated; sketches decoded once per config."""

    def __init__(self, client, space):
        self.client, self.space = client, space
        self._decoded: Dict[str, Optional[dict]] = {}

    def sketch(self, c: Config) -> Optional[dict]:
        if c.key not in self._decoded:
            cs = self.client.callsets.get((c.key, FULL))
            sk = getattr(cs, "sketch", None) if cs is not None else None
            self._decoded[c.key] = None if not sk else {"k": sk.get("k", 128),
                                                         "h": {cls: _unpack(v) for cls, v in sk["h"].items()}}
        return self._decoded[c.key]

    def between(self, a: Config, b: Config) -> Tuple[float, str]:
        """(distance, "callset" | "settings")."""
        sa, sb = self.sketch(a), self.sketch(b)
        if sa is not None and sb is not None:
            return callset_distance(sa, sb), "callset"
        return setting_distance(self.space, a, b), "settings"

    def far(self, a: Config, b: Config, min_callset: float, min_settings: float) -> bool:
        d, kind = self.between(a, b)
        return d >= (min_callset if kind == "callset" else min_settings)
