"""Process limits a big node or main needs.

A node with hundreds of GATK workers holds a few files per worker at once (each persistent GATK server's log, the
scoring and cold-run pipes, the job files). The usual soft limit of 1024 open files is far below that: on
20261008T175300Z the 324-worker node hit "Too many open files" from minute 21 on, and 518 configs failed.
"""
from __future__ import annotations

import logging

log = logging.getLogger(__name__)


def raise_open_files() -> tuple:
    """Raise this process's soft open-file limit to its hard limit. Returns (old soft, new soft); a no-op where the
    platform has no such limit."""
    try:
        import resource
    except ImportError:                                 # not on Unix
        return (None, None)
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    want = hard if hard != resource.RLIM_INFINITY else max(soft, 1 << 20)
    if soft == resource.RLIM_INFINITY or soft >= want:
        return (soft, soft)
    try:
        resource.setrlimit(resource.RLIMIT_NOFILE, (want, hard))
    except (ValueError, OSError) as e:
        log.warning("open-file limit stays at %d: %s", soft, e)
        return (soft, soft)
    log.info("open-file limit raised from %d to %d", soft, want)
    return (soft, want)


# The validator kills a GATK run after this many seconds (minos_subnet PR #40, 2026-10-08:
# VALIDATOR_VARIANT_CALLING_TIMEOUT_SECONDS, default 900; it was 3600). Our runs may take 1800 s, so a config can score
# here and fail on the platform. Validators use 2-8 threads and we use 1, so a run over the limit here is a warning,
# not a certain failure there.
VALIDATOR_RUN_S = 900


def over_validator_limit(runtime_s) -> bool:
    """True when a run of ours took longer than the validator allows one."""
    try:
        return runtime_s is not None and float(runtime_s) > VALIDATOR_RUN_S
    except (TypeError, ValueError):
        return False
