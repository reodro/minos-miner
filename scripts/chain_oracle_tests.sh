#!/usr/bin/env bash
# Back-to-back oracle fleet tests. The live search is never killed. When mining ends the fleet is free;
# review runs here. The next tag is the newest revealed round that already has bam and truth, and that
# has not already been mined.
# Engine succession: the next launch is v9.30 (compass on the highest score).
# usage: nohup bash scripts/chain_oracle_tests.sh >> test/runs/chain_oracle.log 2>&1 &
set -u
cd "$(dirname "$0")/.."
export PYTHONPATH=/minos/minos-miner
export ANSWER_SET_BUILDER=oracle
unset SEARCH_START
LOG=test/runs/chain_oracle.log
PHASE=test/runs/engine_chain_phase
mkdir -p test/runs
say() { echo "[$(date +%H:%M:%S)] $*"; }
if [ ! -f "$PHASE" ]; then
  echo v9.30 > "$PHASE"
fi

while true; do
  if pgrep -f "python -m main replay" >/dev/null; then
    sleep 30
    continue
  fi
  TAG=$(.venv/bin/python - <<'PY'
from __future__ import annotations
import glob, os
from common.platform import Platform
from common.settings import SETTINGS
from common.types import Round
from backtest.data import fetch, FILES

# standing: do not relaunch these tags
skip = {"20261009T093000Z", "20261009T104200Z"}

def mined(tag: str) -> bool:
    for path in glob.glob(f"test/runs/replay_{tag}_*.log"):
        try:
            text = open(path, errors="replace").read()
        except OSError:
            continue
        if "local best" in text or "final TRUE" in text or "\nevaluated " in text or " evaluated " in text:
            return True
    return False

p = Platform()
chosen = None
for r in p.rounds(limit=20):
    rid = r.get("round_id") or r.get("id")
    if not rid:
        continue
    tag = Round.tag_of(rid)
    if tag in skip:
        continue
    try:
        ver = p.verification(rid)
    except Exception:
        continue
    names = {f.get("name") for f in (ver.get("reveal") or {}).get("files") or []}
    if "input.bam" not in names or "truth.vcf.gz" not in names:
        continue
    if mined(tag):
        break  # newest with files already mined: wait for a newer tag
    region = r.get("region") or ""
    chrom = region.split(":")[0] if region else ""
    if chrom and not os.path.exists(SETTINGS.reference_fasta(chrom)):
        continue
    fetch(rid)
    d = os.path.join(SETTINGS.rounds_dir, tag)
    if not all(os.path.exists(os.path.join(d, f)) for f in ("input.bam", "truth.vcf.gz")):
        continue
    chosen = tag
    break
print(chosen or "")
PY
)
  if [ -z "$TAG" ]; then
    sleep 60
    continue
  fi
  VER=$(tr -d '[:space:]' < "$PHASE")
  case "$VER" in
    v9.16|v9.17|v9.18|v9.19|v9.20|v9.21|v9.22|v9.23|v9.24|v9.25|v9.26|v9.27|v9.28|v9.29) ;;
    *) VER=v9.29 ;;
  esac
  .venv/bin/python scripts/set_engine_tag.py "$VER"
  say "next $TAG engine $VER"
  unset SEARCH_START
  export ANSWER_SET_BUILDER=oracle
  bash scripts/fleet_replay.sh "$TAG" 72 >/tmp/chain_launch.log 2>&1
  say "fleet_replay $TAG $VER exit $?"
  if [ "$VER" = v9.16 ] || [ "$VER" = v9.17 ] || [ "$VER" = v9.18 ] || [ "$VER" = v9.19 ] || [ "$VER" = v9.20 ] || [ "$VER" = v9.21 ] || [ "$VER" = v9.22 ] || [ "$VER" = v9.23 ] || [ "$VER" = v9.24 ] || [ "$VER" = v9.25 ] || [ "$VER" = v9.26 ] || [ "$VER" = v9.27 ] || [ "$VER" = v9.28 ]; then
    echo v9.29 > "$PHASE"
    say "next engine v9.29"
  fi
done
