#!/bin/bash
# Per-chromosome tabix-extract + sa-build for one gnomAD v4.1 release: the
# exomes, the genomes, or the joint exomes + genomes callset.
#
#     build_gnomad_release_per_chrom.sh joint 21 22
#
# `build_gnomad_per_chrom.sh` does the same for the exomes only and keeps the
# whole INFO column. That is affordable for the exomes and not for the others:
# the joint release writes about 11 KB of INFO per record (every ancestry group
# by sex, for three callsets), and 200 of the benchmark's 24,350 regions on
# chr22 alone come to 470,664 records and 5.1 GB uncompressed. The builder reads
# a few dozen of those keys, so the filter below keeps only them and the
# extract that reaches the disk is about 3 % of that size.
#
# The keys kept are the ones `crates/fastvep-sa/src/sources/gnomad.rs` reads,
# under the names the chosen release spells them with:
#
#   exomes, genomes   AF AC AN nhomalt AF_<pop> AC_XY AN_XY faf95
#                     fafmax_faf95_max lcr segdup non_par
#   joint             the same with `_joint` appended (and `_joint_XY`,
#                     `AF_joint_<pop>`), plus exomes_filters / genomes_filters,
#                     which carry the per-callset FILTER names the joint release
#                     cannot put in its FILTER column
#
# Writes data/benchmark/sa_db_<release>/gnomad_chr<N>.osa2 and leaves the
# exomes database in sa_db/ alone, so a run can be pointed at either.
#
# A batch of regions that fails is retried, and a chromosome with a batch that
# still fails is not built: the exomes script ignored the exit status of tabix,
# which turns a dropped connection into a database with a hole in it, and a
# hole in a frequency database reads as "absent from gnomAD" - PM2.
set -euo pipefail

RELEASE="${1:-}"
case "$RELEASE" in
  exomes | genomes | joint) shift ;;
  *)
    echo "usage: $(basename "$0") <exomes|genomes|joint> [chrom ...]" >&2
    exit 2
    ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
FASTVEP="${FASTVEP:-$ROOT/target/release/fastvep}"
BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/${RELEASE}/gnomad.${RELEASE}.v4.1.sites"
REGIONS="${REGIONS:-$ROOT/data/benchmark/sa_sources/clinvar_2star_regions.bed}"
EXTRACTS="${EXTRACTS:-$ROOT/data/benchmark/sa_sources/gnomad_${RELEASE}_extracts}"
SA_DB="${SA_DB:-$ROOT/data/benchmark/sa_db_${RELEASE}}"
BATCH="${BATCH:-200}"
mkdir -p "$EXTRACTS" "$SA_DB"

# The INFO keys to keep, as an anchored alternation for awk.
if [ "$RELEASE" = joint ]; then
  KEEP='^(AC_joint|AN_joint|AF_joint|nhomalt_joint|AF_joint_(afr|amr|asj|eas|fin|mid|nfe|oth|remaining|sas)|AC_joint_XY|AN_joint_XY|faf95_joint|fafmax_faf95_max_joint|lcr|segdup|non_par|exomes_filters|genomes_filters)$'
else
  KEEP='^(AC|AN|AF|nhomalt|AF_(afr|amr|asj|eas|fin|mid|nfe|oth|remaining|sas)|AC_XY|AN_XY|faf95|fafmax_faf95_max|lcr|segdup|non_par)$'
fi

# Trim the INFO column of a VCF stream to the kept keys. Header lines pass
# untouched: the builder picks the key naming scheme from the `##INFO` lines.
trim_info() {
  awk -F'\t' -v OFS='\t' -v keep="$KEEP" '
    /^#/ { print; next }
    {
      n = split($8, kv, ";"); out = ""
      for (i = 1; i <= n; i++) {
        key = kv[i]; sub(/=.*/, "", key)
        if (key ~ keep) out = (out == "" ? kv[i] : out ";" kv[i])
      }
      $8 = (out == "" ? "." : out); print
    }'
}

build_one() {
  local chr=$1
  local url="${BASE}.chr${chr}.vcf.bgz"
  local vcf_gz="${EXTRACTS}/gnomad_chr${chr}.vcf.gz"
  local osa2="${SA_DB}/gnomad_chr${chr}.osa2"

  if [ -s "$osa2" ]; then
    echo "[${RELEASE} chr${chr}] already built, skipping"
    return 0
  fi

  if [ ! -s "$vcf_gz" ]; then
    local regions="${EXTRACTS}/regions_chr${chr}.txt"
    awk -v c="chr${chr}" '$1==c {print $1":"$2+1"-"$3}' "$REGIONS" > "$regions"
    local n
    n=$(wc -l < "$regions" | tr -d ' ')
    echo "[${RELEASE} chr${chr}] extracting $n regions in batches of $BATCH..."

    local part="${vcf_gz}.part"
    local start=1 batch_no=0
    # An interrupted run resumes after the last batch that finished, which matters for a chromosome of hours.
    # `.progress` holds the regions done, the bytes of the body they came to and the number of regions they were
    # counted against; the body is cut back to those bytes, so that a batch that was half written when the run
    # stopped is not kept twice, and a regions file that has changed starts again.
    if [ -s "$part.progress" ] && [ -f "$part.body" ]; then
      local regions_done bytes_done regions_of size
      read -r regions_done bytes_done regions_of < "$part.progress"
      size=$(wc -c < "$part.body" | tr -d ' ')
      if [ "$regions_of" = "$n" ] && [ "$size" -ge "$bytes_done" ]; then
        if [ "$size" -gt "$bytes_done" ]; then
          head -c "$bytes_done" "$part.body" > "$part.body.cut" && mv "$part.body.cut" "$part.body"
        fi
        start=$((regions_done + 1))
        batch_no=$(((regions_done + BATCH - 1) / BATCH))
        echo "[${RELEASE} chr${chr}] resuming after region ${regions_done} of $n"
      else
        : > "$part.body"
        rm -f "$part.progress"
      fi
    else
      : > "$part.body"
      rm -f "$part.progress"
    fi
    while [ "$start" -le "$n" ]; do
      batch_no=$((batch_no + 1))
      local ok=0 attempt
      for attempt in 1 2 3; do
        # shellcheck disable=SC2046 # one region per word is the point
        if sed -n "${start},$((start + BATCH - 1))p" "$regions" | xargs tabix "$url" 2> "$part.err" | trim_info > "$part.batch"; then
          # A pipeline's status is its last command's, so read tabix's own.
          if [ "${PIPESTATUS[1]:-0}" -eq 0 ] && [ "${PIPESTATUS[2]:-0}" -eq 0 ]; then
            ok=1
            break
          fi
        fi
        echo "[${RELEASE} chr${chr}] batch ${batch_no} attempt ${attempt} failed: $(head -c 300 "$part.err")" >&2
        sleep $((attempt * 20))
      done
      if [ "$ok" -ne 1 ]; then
        echo "[${RELEASE} chr${chr}] FAILED: batch ${batch_no} did not complete; not building a database with a hole in it" >&2
        return 1
      fi
      cat "$part.batch" >> "$part.body"
      rm -f "$part.batch"
      start=$((start + BATCH))
      local covered=$((start - 1))
      [ "$covered" -le "$n" ] || covered=$n
      printf '%s %s %s\n' "$covered" "$(wc -c < "$part.body" | tr -d ' ')" "$n" > "$part.progress.tmp" \
        && mv "$part.progress.tmp" "$part.progress"
    done

    # The header comes from the file itself, from a region that is certain to
    # exist, so the builder sees the same `##INFO` declarations the release has.
    { tabix -h "$url" "chr${chr}:1-1000" 2> /dev/null | grep '^#' || true; cat "$part.body"; } | bgzip -c > "$vcf_gz"
    rm -f "$part.body" "$part.err" "$part.progress"
    local count
    # `gzip -dc` and not `zcat`: on macOS zcat looks for a `.Z` file.
    count=$(gzip -dc "$vcf_gz" | grep -vc '^#' || true)
    echo "[${RELEASE} chr${chr}] extracted $count records"
  fi

  echo "[${RELEASE} chr${chr}] building .osa2..."
  # The builder says which region flags the release does not declare (the joint release has neither `lcr` nor
  # `segdup`) before it writes anything, so the log is kept, its warnings shown and its last line after them:
  # `| tail -1` alone cut the warning off.
  local build_log="${EXTRACTS}/sa_build_chr${chr}.log"
  "$FASTVEP" sa-build --source gnomad -i "$vcf_gz" -o "${SA_DB}/gnomad_chr${chr}" --assembly GRCh38 > "$build_log" 2>&1 \
    || { cat "$build_log" >&2; return 1; }
  { grep -E '^warning:' "$build_log" || true; } | sort -u
  tail -1 "$build_log"

  if [ ! -s "$osa2" ]; then
    echo "[${RELEASE} chr${chr}] FAILED: $osa2 was not produced" >&2
    return 1
  fi
  echo "[${RELEASE} chr${chr}] done"
}

export -f build_one trim_info
export RELEASE BASE REGIONS EXTRACTS SA_DB BATCH KEEP FASTVEP

CHROMS=${*:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y}
# Four at a time: the release's bucket is the bottleneck, not the CPU.
echo "$CHROMS" | tr ' ' '\n' | xargs -n 1 -P "${JOBS:-4}" bash -c 'build_one "$@"' _

echo "==> done"
ls -la "$SA_DB"/gnomad_chr*.osa2
