#!/usr/bin/env bash
# Download the reference data answer_set needs, into $DATA_DIR (default /minos/data).
#
#   1kg/phased/  2022 phased panel (1000G person lookup)
#   giab/        NIST v4.2.1 calls get a .tbi index for region queries (the files themselves are not fetched here)
#
#   bash scripts/fetch_reference_data.sh                 # chr14-22
#   CHROMS="1 2 3" bash scripts/fetch_reference_data.sh  # others
#
# Needs curl, bcftools. Safe to re-run: finished files are skipped.
set -uo pipefail

DATA_DIR="${DATA_DIR:-/minos/data}"
CHROMS="${CHROMS:-14 15 16 17 18 19 20 21 22}"
KG=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections
PHASED=$KG/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
OUT="$DATA_DIR/1kg"
mkdir -p "$OUT"/{phased,logs}

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$OUT/logs/fetch.log"; }

remote_size() { curl -sfI "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}'; }

fetch_plain() {  # url dst: resumable download of a small file with size check
  local url=$1 dst=$2 want
  want=$(remote_size "$url")
  [[ -s "$dst" && $(stat -c %s "$dst") == "$want" ]] && return 0
  curl -sfS --retry 5 -C - -o "$dst" "$url" && [[ $(stat -c %s "$dst") == "$want" ]] \
    && log "got $(basename "$dst")" || { log "FAILED $url"; return 1; }
}

for c in $CHROMS; do
  p="1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  fetch_plain "$PHASED/$p" "$OUT/phased/chr$c.vcf.gz"
  fetch_plain "$PHASED/$p.tbi" "$OUT/phased/chr$c.vcf.gz.tbi"
done
log "panels done"

# GIAB NIST v4.2.1 calls must be indexed for region queries (answer_set.k2)
GIAB_DIR="$(dirname "$OUT")/giab"
if [ -d "$GIAB_DIR" ]; then
  for v in "$GIAB_DIR"/HG00?.vcf.gz; do [ -f "$v.tbi" ] || bcftools index -t -f "$v"; done
fi
log "all done"
