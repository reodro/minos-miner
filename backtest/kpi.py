"""KPIs every part is judged by. Pure functions: easy to test, no I/O.

Legacy lesson: identical 63-minute search runs differed by 0.083 points. Any
comparison below ~0.1 needs replicated pairs before anyone believes it.
"""
from __future__ import annotations

from typing import Dict, List


def oracle_gap(ours: Dict[str, float], oracle: Dict[str, float]) -> float:
    """answer_set KPI. Scores of the same callsets under our answer set and the released truth.

    Returns oracle score of the oracle's pick minus oracle score of our pick (>= 0; 0 is perfect).
    """
    common_keys = ours.keys() & oracle.keys()
    if not common_keys:
        raise ValueError("no callsets scored under both answer sets")
    our_pick = max(common_keys, key=lambda k: ours[k])
    best = max(common_keys, key=lambda k: oracle[k])
    return oracle[best] - oracle[our_pick]


def rank_in_field(score: float, field: List[float]) -> int:
    """main KPI. 1-based rank our true score would have taken among the real submissions."""
    return 1 + sum(1 for s in field if s > score)


def would_win(score: float, field: List[float]) -> bool:
    return bool(field) and score > max(field)
