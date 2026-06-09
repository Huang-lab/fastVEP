#!/usr/bin/env bash
# Per-chromosome: tabix-extract from gnomAD VCF → spliceai/phylop side-files →
# sa-build per chrom .osa. Skips chroms whose both .osa already exist.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
EXTRACTS=$ROOT/data/benchmark/sa_sources/spliceai_phylop_extracts
SA_DB=$ROOT/data/benchmark/sa_db
LOG_DIR=$ROOT/data/benchmark/sa_sources/spliceai_phylop_logs
mkdir -p "$EXTRACTS" "$SA_DB" "$LOG_DIR"

build_one() {
  local chr=$1
  local sa="$SA_DB/spliceai_chr${chr}.osa"
  local pp="$SA_DB/phylop_chr${chr}.osa"
  local log="$LOG_DIR/chr${chr}.log"

  if [ -s "$sa" ] && [ -s "$pp" ] && [ "${REBUILD:-0}" != "1" ]; then
    echo "[chr${chr}] both .osa exist, skipping"
    return 0
  fi

  echo "[chr${chr}] $(date +%H:%M:%S) extracting..." >>"$log"
  python3 "$ROOT/data/benchmark/sa_sources/extract_spliceai_phylop_from_gnomad.py" "$chr" >>"$log" 2>&1

  local sa_vcf="$EXTRACTS/spliceai_chr${chr}.vcf"
  local pp_tsv="$EXTRACTS/phylop_chr${chr}.tsv"

  if [ ! -s "$sa_vcf" ] || [ ! -s "$pp_tsv" ]; then
    echo "[chr${chr}] ERROR: extract empty"
    return 1
  fi

  "$ROOT/target/release/fastvep" sa-build --source spliceai \
    -i "$sa_vcf" -o "$SA_DB/spliceai_chr${chr}" --assembly GRCh38 >>"$log" 2>&1
  "$ROOT/target/release/fastvep" sa-build --source phylop \
    -i "$pp_tsv" -o "$SA_DB/phylop_chr${chr}" --assembly GRCh38 >>"$log" 2>&1

  rm -f "$sa_vcf" "$pp_tsv"
  echo "[chr${chr}] $(date +%H:%M:%S) done"
}

export -f build_one
export ROOT EXTRACTS SA_DB LOG_DIR

CHROMS=${@:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y}
PARALLEL=${PARALLEL:-4}

echo "$CHROMS" | tr ' ' '\n' | xargs -n 1 -P "$PARALLEL" bash -c 'build_one "$@"' _
echo "==> All done"
ls -la "$SA_DB"/spliceai_chr*.osa "$SA_DB"/phylop_chr*.osa | sort -V
