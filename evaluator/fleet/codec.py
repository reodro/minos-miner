"""JSON-safe dicts for what crosses the network between main and the nodes: jobs out, callsets (+ scores) back."""
from __future__ import annotations

from dataclasses import asdict
from typing import Optional

from common.types import CallSet, ClassCounts, Config, EvalJob, Round, ScoreResult


def round_to(r: Round) -> dict:
    return asdict(r)


def round_from(d: dict) -> Round:
    return Round(**d)


def job_to(j: EvalJob) -> dict:
    return {"job_id": j.job_id, "round": round_to(j.round), "config": j.config.to_dict(), "fidelity": j.fidelity,
            "priority": j.priority, "tag": j.tag, "lane": j.lane}


def job_from(d: dict) -> EvalJob:
    return EvalJob(d["job_id"], round_from(d["round"]), Config.from_dict(d["config"]), d["fidelity"], d["priority"],
                   d["tag"], d.get("lane", "throughput"))


def score_to(s: Optional[ScoreResult]) -> Optional[dict]:
    return None if s is None else asdict(s)


def score_from(d: Optional[dict]) -> Optional[ScoreResult]:
    if d is None:
        return None
    d = dict(d)
    d["counts"] = {k: ClassCounts(**v) for k, v in d["counts"].items()}
    return ScoreResult(**d)


def callset_to(c: CallSet) -> dict:
    d = asdict(c)
    d["score"] = score_to(c.score)
    return d


def callset_from(d: dict) -> CallSet:
    d = dict(d)
    d["score"] = score_from(d.get("score"))
    return CallSet(**d)
