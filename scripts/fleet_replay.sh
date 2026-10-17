#!/usr/bin/env bash
# A timed fleet replay of one stored round - how every engine test on the fleet is run:
#   1. sync this checkout's code to every node in fleet.json (copies the repo, pushes the self-test bundle,
#      rebuilds the GATK patch jar there)
#   2. wipe the round's caches on every node (callsets, scores, decision cache), so nothing is reused from an
#      earlier run of the same round. The replay itself also deletes this machine's saved results for the tag,
#      so the dashboard starts from this run.
#   3. python -m main replay <tag> --fleet --minutes <min> --no-review: the search with dry-run submissions
#   4. python -m main review <tag> in the background here (TRUE scores and rank); the fleet is already free
# usage: scripts/fleet_replay.sh <round_tag> [minutes=72]          log: test/runs/replay_<tag>_<HHMM>.log
set -u
cd "$(dirname "$0")/.."
TAG=$1; MIN=${2:-72}
mkdir -p test/runs; LOG=test/runs/replay_${TAG}_$(date +%H%M).log
echo "[$(date +%H:%M:%S)] $(git branch --show-current) $(git log --oneline -1 | cut -c1-7) tag=$TAG minutes=$MIN log=$LOG"
# the code is synced only when it changed since the last sync (a 72-minute test cadence has no 30 s to spare);
# a node with other code is rejected when it connects, so a stale stamp shows up at once
STAMP=cache/fleet_synced.stamp
CODE="$(git rev-parse HEAD) $( (git diff HEAD; git ls-files --others --exclude-standard -- '*.py' | xargs -r cat) | sha1sum | cut -c1-16) $(sha1sum fleet.json | cut -c1-16)"
if [ -f "$STAMP" ] && [ "$(cat "$STAMP")" = "$CODE" ]; then
    echo "[$(date +%H:%M:%S)] code unchanged since the last sync; not synced"
else
    echo "[$(date +%H:%M:%S)] sync this code to every node"
    .venv/bin/python -m evaluator.fleet.control sync | cut -c1-160
    mkdir -p cache; echo "$CODE" > "$STAMP"
fi
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
.venv/bin/python -m main replay "$TAG" --fleet --minutes "$MIN" --no-review > "$LOG" 2>&1
echo "[$(date +%H:%M:%S)] replay exit $?"
echo "[$(date +%H:%M:%S)] review $TAG here (fleet is free)"
nohup nice -n 19 .venv/bin/python -m main review "$TAG" >> "$LOG" 2>&1 &   # lowest CPU priority: the next round runs meanwhile
echo "[$(date +%H:%M:%S)] review pid $!"
grep -E "answer set vs|  (final|baseline|top) |final TRUE|^rank|searched|evaluated |local best" "$LOG" | cut -c1-220
