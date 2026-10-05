#!/usr/bin/env bash
# GRCh38 reference per chromosome, exactly as the validator uses it: FASTA + .fai + .dict + RTG SDF (for hap.py).
# Source: the platform's reference endpoint (302 -> storage; needs an explicit User-Agent).
#   bash scripts/fetch_genome.sh                    # chr14-chr22: every round so far (176/176) was on one of them
#   CHROMS="chr20 chr21" bash scripts/fetch_genome.sh
# A node fetches any other chromosome on demand when a round needs it (evaluator/fleet/node.py).
set -euo pipefail
DATA_DIR="${DATA_DIR:-/minos/data}"
BASE="https://api.theminos.ai/reference"
UA="minos-installer/0.1 (+https://github.com/minos-protocol/minos_subnet)"
CHROMS="${CHROMS:-$(echo chr{14..22})}"
SDF_FILES="done mainIndex nameIndex0 namedata0 namepointer0 progress seqdata0 seqpointer0 sequenceIndex0 summary.txt"

get() {   # url dest
  [ -s "$2" ] && return 0
  mkdir -p "$(dirname "$2")"
  curl -sSfL --retry 4 -A "$UA" -o "$2.part" "$1" && mv "$2.part" "$2"
}

for c in $CHROMS; do
  d="$DATA_DIR/reference/$c"
  for f in "$c.fa" "$c.fa.fai" "$c.dict"; do get "$BASE/$c/$f" "$d/$f"; done
  for f in $SDF_FILES; do get "$BASE/$c/$c.sdf/$f" "$d/$c.sdf/$f"; done
  echo "$(date -u +%FT%TZ) $c done"
done
