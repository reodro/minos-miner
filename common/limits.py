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
