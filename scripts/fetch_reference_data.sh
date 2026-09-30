#!/usr/bin/env bash
# Download the reference data answer_set needs, into $DATA_DIR (default /minos/data).
#
#   1kg/raw/     NYGC 30x raw callset (3,202 samples), every record and sample kept,
#                FORMAT reduced to GT:AD:DP:GQ (PL and AB dropped, ~53% of original).
#                Streamed through bcftools: the ~395 GB original is never stored.
#   1kg/phased/  2022 phased panel (donor fingerprinting)
#   1kg/mask/    1000G StrictMask (synthetic sites fall 100% inside "P")
#
#   bash scripts/fetch_reference_data.sh                 # chr14-22
#   CHROMS="1 2 3" bash scripts/fetch_reference_data.sh  # others
#
# Needs curl, bcftools, tabix. Safe to re-run: finished files are skipped.
set -uo pipefail

DATA_DIR="${DATA_DIR:-/minos/data}"
CHROMS="${CHROMS:-14 15 16 17 18 19 20 21 22}"
PARALLEL="${PARALLEL:-4}"
KG=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections
RAW=$KG/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot
PHASED=$KG/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
MASK=$KG/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask
OUT="$DATA_DIR/1kg"
mkdir -p "$OUT"/{raw,phased,mask,logs}

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$OUT/logs/fetch.log"; }

remote_size() { curl -sfI "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}'; }

fetch_raw() {
  local c=$1 url="$RAW/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr$1.recalibrated_variants.vcf.gz"
  local dst="$OUT/raw/chr$c.vcf.gz" tmp="$OUT/raw/chr$c.vcf.gz.part"
  [[ -s "$dst" && -s "$dst.tbi" ]] && { log "raw chr$c: already done"; return 0; }
  local want; want=$(remote_size "$url")
  for attempt in 1 2 3; do
    log "raw chr$c: attempt $attempt, $((want / 1000000)) MB remote"
    # No curl --retry here: a retry mid-stream would resend bytes into the pipe.
    # A broken stream is caught by the byte count and the whole chromosome is redone.
    local sizef="$OUT/logs/raw_chr$c.size" got rc
    curl -sfS -w '%{stderr}%{size_download}\n' "$url" 2>"$sizef" \
      | bcftools annotate -x FORMAT/PL,FORMAT/AB --threads 2 -Oz -o "$tmp" - 2>>"$OUT/logs/raw_chr$c.err"
    rc=("${PIPESTATUS[@]}")
    got=$(tail -n1 "$sizef")
    if [[ ${rc[0]} == 0 && ${rc[1]} == 0 && "$got" == "$want" ]] && tabix -f -p vcf "$tmp"; then
      mv "$tmp.tbi" "$dst.tbi" && mv "$tmp" "$dst"
      log "raw chr$c: done $(du -h "$dst" | cut -f1)"
      return 0
    fi
    log "raw chr$c: FAILED (curl rc=${rc[0]} bcftools rc=${rc[1]}, got $got of $want bytes); retrying"
    rm -f "$tmp" "$tmp.tbi"
  done
  log "raw chr$c: GAVE UP"; return 1
}

fetch_plain() {  # url dst: resumable download of a small file with size check
  local url=$1 dst=$2 want
  want=$(remote_size "$url")
  [[ -s "$dst" && $(stat -c %s "$dst") == "$want" ]] && return 0
  curl -sfS --retry 5 -C - -o "$dst" "$url" && [[ $(stat -c %s "$dst") == "$want" ]] \
    && log "got $(basename "$dst")" || { log "FAILED $url"; return 1; }
}

export -f fetch_raw log remote_size
export OUT RAW

# EBI throughput is per connection (3-25 MB/s), so panels and raw streams run concurrently.
(
  for c in $CHROMS; do
    p="1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    fetch_plain "$PHASED/$p" "$OUT/phased/chr$c.vcf.gz"
    fetch_plain "$PHASED/$p.tbi" "$OUT/phased/chr$c.vcf.gz.tbi"
    fetch_plain "$MASK/20160622.chr$c.mask.fasta.gz" "$OUT/mask/chr$c.strict.fa.gz"
  done
  log "panels and masks done"
) &

printf '%s\n' $CHROMS | xargs -P "$PARALLEL" -I{} bash -c 'fetch_raw {}'
wait

# 1000G 30x CRAM locations (answer_set.parent reads the donor's CRAM from the S3 mirror by region)
KGHC=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
mkdir -p "$OUT/cram_index"
for f in 1000G_2504_high_coverage.sequence.index 1000G_698_related_high_coverage.sequence.index; do
  fetch_plain "$KGHC/$f" "$OUT/cram_index/$f"
done

# GIAB NIST v4.2.1 calls must be indexed for region queries (answer_set.giab)
GIAB_DIR="$(dirname "$OUT")/giab"
if [ -d "$GIAB_DIR" ]; then
  for v in "$GIAB_DIR"/HG00?.vcf.gz; do [ -f "$v.tbi" ] || bcftools index -t -f "$v"; done
fi
log "all done"
