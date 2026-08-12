#!/bin/bash
# Per-chromosome tabix-extract + sa-build for gnomAD exomes v4.1.
# Deletes intermediate VCF after building to keep disk bounded.
#
# `sa-build --source gnomad` emits .osa2. The annotation loader picks up both
# .osa and .osa2 from --sa-dir, so a leftover v1 file beside a freshly built v2
# one would load the same source twice; the v1 copy is removed once the v2
# build has succeeded.
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
  local osa2="${SA_DB}/gnomad_chr${chr}.osa2"

  if [ -s "$osa2" ]; then
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

  echo "[chr${chr}] building .osa2..."
  $ROOT/target/release/fastvep sa-build \
    --source gnomad \
    -i "${vcf_gz}" \
    -o "${SA_DB}/gnomad_chr${chr}" \
    --assembly GRCh38 2>&1 | tail -1

  if [ ! -s "$osa2" ]; then
    echo "[chr${chr}] FAILED: $osa2 was not produced" >&2
    return 1
  fi

  # Retire the v1 database this run replaces, and drop the intermediate VCF.
  rm -f "${SA_DB}/gnomad_chr${chr}.osa" "${SA_DB}/gnomad_chr${chr}.osa.idx"
  rm -f "${vcf_gz}"
  echo "[chr${chr}] done"
}

export -f build_one
export ROOT SRC_BASE EXTRACTS SA_DB REGIONS

# gnomAD v4 releases no chrMT file for exomes (mitochondrial calls ship as a
# separate release with a different schema), so MT is not in this list even
# though clinvar_2star_regions.bed contains chrMT rows.
CHROMS=${@:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y}

# 4-way parallel
echo "$CHROMS" | tr ' ' '\n' | xargs -n 1 -P 4 bash -c 'build_one "$@"' _

# An earlier v1 run did produce an empty gnomad_chrMT.osa. Nothing above can
# ever replace it, and the loader accepts .osa and .osa2 side by side, so it
# would linger as a v1 gnomAD source in a directory this script exists to make
# v2. Retire it here rather than leaving it to be rediscovered.
if [ -e "${SA_DB}/gnomad_chrMT.osa" ]; then
  echo "Removing gnomad_chrMT.osa: gnomAD v4 has no chrMT exomes release to rebuild it from"
  rm -f "${SA_DB}/gnomad_chrMT.osa" "${SA_DB}/gnomad_chrMT.osa.idx"
fi

echo "==> All done"
ls -la "$SA_DB"/gnomad_chr*.osa2 | head -30
