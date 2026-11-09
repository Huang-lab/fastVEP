#!/usr/bin/env bash
# The benchmark's databases with SpliceAI, PhyloP and gnomAD built from the whole of the gnomAD v4.1 exomes file.
#
#     bash analysis/acmg_benchmark/scripts/sa_sources/make_complete_stack.sh
#
# data/benchmark/sa_db/ holds the databases that the benchmark was first made with: SpliceAI and PhyloP as `.osa` files
# distilled in April from regions of the gnomAD file, and gnomAD as `.osa2` files of the same regions. They lack records
# that the gnomAD file carries (see results/complete_databases/README.md), so the benchmark now reads
# data/benchmark/sa_db_complete/ instead, which this script makes: a directory of symbolic links to every file of sa_db/
# except those three sources, and to the 72 complete files of data/benchmark/sa_sources/complete/ in their place. REVEL,
# ClinVar, the ClinVar protein index, ClinGen, gnomAD's gene constraints and RepeatMasker are the same files as before:
# their dates are part of the design of the benchmark.
#
# The 72 files are not built here: they were built on the production host of fastvep.org from the whole gnomAD file and
# copied (results/complete_databases/README.md says how), and this script refuses to link them unless every one of them has
# the checksum that results/complete_databases/SHA256SUMS records.
#
# Environment: BENCH_SA_DB_OLD (default sa_db) and BENCH_SA_DB_NEW (default sa_db_complete) name the two directories under
# data/benchmark/, and COMPLETE_DIR the directory of the 72 files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
BENCH="$ROOT/data/benchmark"
OLD="$BENCH/${BENCH_SA_DB_OLD:-sa_db}"
NEW="$BENCH/${BENCH_SA_DB_NEW:-sa_db_complete}"
COMPLETE_DIR="${COMPLETE_DIR:-$BENCH/sa_sources/complete}"
MANIFEST="$ROOT/analysis/acmg_benchmark/results/complete_databases/SHA256SUMS"

[ -d "$OLD" ] || { echo "missing $OLD" >&2; exit 1; }
[ -d "$COMPLETE_DIR" ] || { echo "missing $COMPLETE_DIR" >&2; exit 1; }
# The links name the files by their physical paths: in a worktree `data/` is a link to the main checkout's, and a link through
# it would break with the worktree.
OLD="$(cd "$OLD" && pwd -P)"
COMPLETE_DIR="$(cd "$COMPLETE_DIR" && pwd -P)"
[ -s "$MANIFEST" ] || { echo "missing $MANIFEST" >&2; exit 1; }

# Every file of the manifest is there and has its checksum, and nothing else of the three sources is.
if command -v sha256sum >/dev/null 2>&1; then check=(sha256sum -c --quiet); else check=(shasum -a 256 -c --quiet); fi
(cd "$COMPLETE_DIR" && "${check[@]}" "$MANIFEST") || { echo "a file of $COMPLETE_DIR does not have the checksum of $MANIFEST" >&2; exit 1; }
expected=$(wc -l < "$MANIFEST" | tr -d ' ')
present=$(find "$COMPLETE_DIR" -maxdepth 1 -name '*.osa2' | wc -l | tr -d ' ')
[ "$expected" = "$present" ] || { echo "$COMPLETE_DIR holds $present .osa2 files and the manifest lists $expected" >&2; exit 1; }

rm -rf "$NEW"
mkdir -p "$NEW"
kept=0
for f in "$OLD"/*; do
  name="$(basename "$f")"
  case "$name" in
    spliceai_chr*|phylop_chr*|gnomad_chr*) continue ;;   # replaced by the complete files
    .*) continue ;;
  esac
  ln -s "$f" "$NEW/$name"
  kept=$((kept + 1))
done
for f in "$COMPLETE_DIR"/*.osa2; do
  ln -s "$f" "$NEW/$(basename "$f")"
done

echo "$NEW: $kept links to files of $OLD and $present complete files"
ls "$NEW" | sed -e 's/chr[0-9XY]*//' | sort | uniq -c
