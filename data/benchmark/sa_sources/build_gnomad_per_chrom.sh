#!/bin/bash
# Per-chromosome tabix-extract + sa-build for gnomAD exomes v4.1.
# Deletes intermediate VCF after building .osa to keep disk bounded.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
SRC_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites"
EXTRACTS=$ROOT/data/benchmark/sa_sources/gnomad_extracts
SA_DB=$ROOT/data/benchmark/sa_db
REGIONS=$ROOT/data/benchmark/sa_sources/clinvar_2star_regions.bed
mkdir -p "$EXTRACTS" "$SA_DB"

build_one() {
  local chr=$1
  local url="${SRC_BASE}.chr${chr}.vcf.bgz"
  local vcf="${EXTRACTS}/gnomad_chr${chr}.vcf"
  local vcf_gz="${vcf}.gz"
  local osa="${SA_DB}/gnomad_chr${chr}.osa"

  if [ -s "$osa" ]; then
    echo "[chr${chr}] already built, skipping"
    return 0
  fi

  if [ ! -s "${vcf_gz}" ]; then
    local chr_regions
    chr_regions=$(awk -v c="chr${chr}" '$1==c {print $1":"$2+1"-"$3}' "$REGIONS")
    local n=$(echo "$chr_regions" | wc -l)
    echo "[chr${chr}] extracting $n regions..."
    # Header
    tabix -h "$url" "chr${chr}:1-1000" 2>/dev/null | grep "^#" > "$vcf"
    # Body
    echo "$chr_regions" | xargs -n 200 tabix "$url" >> "$vcf" 2>/dev/null || true
    local count=$(grep -vc "^#" "$vcf" || true)
    echo "[chr${chr}] extracted $count records"
    gzip -f "$vcf"
  fi

  echo "[chr${chr}] building .osa..."
  $ROOT/target/release/fastvep sa-build \
    --source gnomad \
    -i "${vcf_gz}" \
    -o "${SA_DB}/gnomad_chr${chr}" \
    --assembly GRCh38 2>&1 | tail -1
  # Delete the intermediate VCF.gz after build to save disk
  rm -f "${vcf_gz}"
  echo "[chr${chr}] done"
}

export -f build_one
export ROOT SRC_BASE EXTRACTS SA_DB REGIONS

CHROMS=${@:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 X Y}

# 4-way parallel
echo "$CHROMS" | tr ' ' '\n' | xargs -n 1 -P 4 bash -c 'build_one "$@"' _
echo "==> All done"
ls -la "$SA_DB"/gnomad_chr*.osa | head -30
