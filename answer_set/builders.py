"""Which answer-set builder runs. Owner: answer_set - main only calls builder(); adding or switching builders is a
change here, never in main."""
from __future__ import annotations

from typing import Optional

from answer_set.base import AnswerSetProvider

BUILDER_NAMES = ("provider",)          # the first one is the default


def builder(name: Optional[str] = None) -> AnswerSetProvider:
    """The answer-set builder `name` (None: the default)."""
    name = name or BUILDER_NAMES[0]
    if name == "provider":
        from answer_set.provider import Provider
        return Provider()
    raise ValueError(f"unknown answer-set builder {name!r}; known: {', '.join(BUILDER_NAMES)}")
