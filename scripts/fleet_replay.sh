#!/usr/bin/env bash
# A timed fleet replay of one stored round - how every engine test on the fleet is run:
#   1. sync this checkout's code to every node in fleet.json (copies the repo, pushes the self-test bundle,
#      rebuilds the GATK patch jar there)
#   2. wipe the round's caches on every node (callsets, scores, decision cache), so nothing is reused from an
#      earlier run of the same round
#   3. python -m main replay <tag> --fleet --minutes <min>: the search with dry-run submissions, then the truth
#      review and the rank, printed at the end of the log
# usage: scripts/fleet_replay.sh <round_tag> [minutes=72]          log: test/runs/replay_<tag>_<HHMM>.log
set -u
cd "$(dirname "$0")/.."
TAG=$1; MIN=${2:-72}
mkdir -p test/runs; LOG=test/runs/replay_${TAG}_$(date +%H%M).log
echo "[$(date +%H:%M:%S)] $(git branch --show-current) $(git log --oneline -1 | cut -c1-7) tag=$TAG minutes=$MIN log=$LOG"
echo "[$(date +%H:%M:%S)] sync this code to every node"
.venv/bin/python -m evaluator.fleet.control sync | cut -c1-160
echo "[$(date +%H:%M:%S)] wipe cache/{dcache,fleet,callsets}/$TAG on every node"
.venv/bin/python - "$TAG" <<'PY'
import sys
sys.path.insert(0, ".")
from concurrent.futures import ThreadPoolExecutor
from evaluator.fleet.control import sh
from evaluator.fleet.inventory import ensure_agent, load
tag = sys.argv[1]
hosts = load(); ensure_agent(hosts)
def one(h):
    r = sh(h, f"pgrep -f 'evaluator.fleet.node --stdi[o]' >/dev/null && echo AGENT-RUNNING; "
              f"rm -rf {h.repo}/cache/dcache/{tag} {h.repo}/cache/fleet/{tag} {h.repo}/cache/callsets/{tag} && echo wiped",
           timeout=120)
    return f"  {h.name}: {r.stdout.strip() or r.stderr.strip()[-200:]}"
with ThreadPoolExecutor(len(hosts)) as pool:
    for line in pool.map(one, hosts):
        print(line, flush=True)
PY
echo "[$(date +%H:%M:%S)] replay $TAG for $MIN min on the fleet (watch it on the dashboard)"
.venv/bin/python -m main replay "$TAG" --fleet --minutes "$MIN" > "$LOG" 2>&1
echo "[$(date +%H:%M:%S)] replay exit $?"
grep -E "answer set vs|  (final|baseline|top) |final TRUE|^rank|searched" "$LOG" | cut -c1-220
