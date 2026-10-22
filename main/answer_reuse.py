"""Reuse a round's finished answer set after a restart (owner: main; answer_set/ is not changed).

The k2 builder always deletes the round's build folder and builds again (answer_set/k2/provider.py), so every
restart of the live miner inside a round cost a whole build (71-179 s measured, 2026-10-21; twice in round 02:26).
`reusing(builder)` wraps the builder: when the round's folder holds a COMPLETE build from an earlier process (report.json
is written last by answer_set/k2/build.py, so it is the completeness mark, with the three output files beside it), the
answer set is read from it with the builder's own reader (read_output), which also redoes the N estimate as a build
would. Anything else (no folder, a partial build, another builder) builds as before.
"""
from __future__ import annotations

import logging
import os
import time
from pathlib import Path

log = logging.getLogger(__name__)
DONE_FILES = ("report.json", "truth.vcf.gz", "mutations.vcf.gz", "scores.tsv.gz")


class ReusingBuilder:
    def __init__(self, inner):
        self.inner = inner

    def build(self, round):
        out = self.finished(round)
        if out is not None:
            from answer_set.k2.provider import read_output
            try:
                answer = read_output(round, out)
                log.warning("%s answer set: reused the finished build in %s (no rebuild after a restart)",
                            round.tag, out)
                delay = float(os.environ.get("MINOS_ANSWER_DELAY_S") or 0)
                if delay > 0:            # a replay under live conditions: ready after the round's real build time
                    log.warning("%s answer set: held for %.0f s (the live build time)", round.tag, delay)
                    time.sleep(delay)
                return answer
            except Exception:  # noqa: BLE001 - an unreadable folder: build as before
                log.exception("%s answer set: the finished build could not be read; building again", round.tag)
        return self.inner.build(round)

    def finished(self, round):
        work_dir = getattr(self.inner, "work_dir", None)
        if work_dir is None or type(self.inner).__name__ != "K2Provider":
            return None
        out = Path(work_dir) / round.tag / "v1"
        return out if all((out / f).exists() for f in DONE_FILES) else None

    def update(self, answer_set, new_callsets):
        return self.inner.update(answer_set, new_callsets)

    def __getattr__(self, name):
        return getattr(self.inner, name)


BUILD_THREADS = 12       # live 02:26 (2026-10-21): 345 s with the builder's default 8 threads on a busy main; main has 16
                         # cores and its own loop needs few, so the build gets 12 (a setting of the builder, not its code)


def reusing(builder):
    if hasattr(builder, "threads") and int(getattr(builder, "threads") or 0) < BUILD_THREADS:
        builder.threads = BUILD_THREADS
    return ReusingBuilder(builder)
