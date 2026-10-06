#!/usr/bin/env bash
# End-to-end ClinVar 2-star+ ACMG concordance benchmark.
#
# Replaces the run_full_benchmark_v9.sh ... _v14.sh family, which were seven
# byte-identical copies differing only in the version string. Pass the version
# instead:
#
#     bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15
#
# Code and documentation live under analysis/acmg_benchmark/ and are tracked;
# everything this script reads or writes lives under data/benchmark/ and is not.
# That split is the whole reason for the layout - see ../README.md.
#
# Inputs (all under data/benchmark/):
#   clinvar_2star.vcf         ClinVar 2-star+ VCF, ~673k SNV/small-indel
#   clinvar_2star_truth.tsv   truth table for concordance
#   sa_db/                    built .osa/.osa2/.oga databases
#   plus test_data/organisms/human/ for the FASTA + GFF3 + cache
#
# Outputs (data/benchmark/output_<version>/):
#   clinvar_2star.fastvep.vcf.gz   annotation + ACMG, bgzipped
#   concordance_summary.txt        headline metrics
#   concordance_matrix.csv         truth x predicted
#   concordance_by_chrom.csv, concordance_by_consequence.csv
#   criterion_firing_rates.csv, rule_distribution.csv
#   discrepancies.tsv              every non-exact call
#
# Output is VCF-bgzipped (~10x smaller than pretty-printed JSON and ~100x
# smaller than tab format with all fields). VCF is the only format that carries
# ACMG/ACMG_CRITERIA in the per-transcript CSQ entry.
#
# An existing annotated VCF is reused rather than regenerated, so re-running to
# recompute concordance after a scoring change is cheap. Delete it to force a
# fresh annotation.

set -euo pipefail

VERSION="${1:-}"
if [ -z "$VERSION" ]; then
  echo "usage: $(basename "$0") <version>    e.g. $(basename "$0") v15" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

INPUT="${INPUT:-$ROOT/data/benchmark/clinvar_2star.vcf}"
TRUTH="${TRUTH:-$ROOT/data/benchmark/clinvar_2star_truth.tsv}"
GFF3="${GFF3:-$ROOT/test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3}"
FASTA="${FASTA:-$ROOT/test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"
SA_DIR="${SA_DIR:-$ROOT/data/benchmark/sa_db}"
FASTVEP="${FASTVEP:-$ROOT/target/release/fastvep}"
OUT_DIR="${OUT_DIR:-$ROOT/data/benchmark/output_$VERSION}"

# An extra --acmg-config, to run a variant of the same version.
#
# Expanded below as ${CONFIG_ARGS[@]+...} rather than plain "${CONFIG_ARGS[@]}":
# under `set -u`, bash 3.2 - which is what macOS ships - treats an empty array
# as an unbound variable and aborts. An `if` rather than `[ ... ] && ...` for
# the same class of reason: a failing test as the last command in an && list
# is a `set -e` trap.
ACMG_CONFIG="${ACMG_CONFIG:-}"
CONFIG_ARGS=()
if [ -n "$ACMG_CONFIG" ]; then
  CONFIG_ARGS=(--acmg-config "$ACMG_CONFIG")
fi

# Any other `fastvep annotate` flags, for a variant of the same version - for
# instance `EXTRA_ARGS="--pick-order rank,mane_select"`. Split on whitespace,
# newlines included, and never globbed, so a flag and its value stay two words.
EXTRA_ARGS="${EXTRA_ARGS:-}"
EXTRA_ARR=()
if [ -n "$EXTRA_ARGS" ]; then
  set -f
  # shellcheck disable=SC2206 # splitting is the point
  EXTRA_ARR=($EXTRA_ARGS)
  set +f
fi

for f in "$INPUT" "$TRUTH" "$GFF3" "$FASTA" "$FASTVEP"; do
  [ -e "$f" ] || { echo "missing required input: $f" >&2; exit 1; }
done
[ -d "$SA_DIR" ] || { echo "missing SA database directory: $SA_DIR" >&2; exit 1; }

mkdir -p "$OUT_DIR"
VCFGZ="$OUT_DIR/clinvar_2star.fastvep.vcf.gz"
PROVENANCE="$OUT_DIR/run_provenance.txt"

# What produced the annotation, recorded beside it. The version name is the only
# part of the output's key, so without this a second run under the same name
# with a different config or flags would find the first run's annotation and
# score it as if it were its own.
ARGS_LINE="annotate --acmg --pick ${CONFIG_ARGS[*]+${CONFIG_ARGS[*]}} ${EXTRA_ARR[*]+${EXTRA_ARR[*]}}"
CONFIG_HASH="none"
if [ -n "$ACMG_CONFIG" ]; then
  CONFIG_HASH="$(shasum -a 256 "$ACMG_CONFIG" | cut -d' ' -f1)"
fi
CURRENT="args: $ARGS_LINE
acmg_config_sha256: $CONFIG_HASH"

if [ -s "$VCFGZ" ] && [ -f "$PROVENANCE" ]; then
  RECORDED="$(grep -E '^(args|acmg_config_sha256):' "$PROVENANCE")"
  if [ "$RECORDED" != "$CURRENT" ]; then
    echo "error: $VCFGZ was produced with different arguments than this run asks for." >&2
    echo "--- recorded ($PROVENANCE)" >&2
    echo "$RECORDED" >&2
    echo "--- requested" >&2
    echo "$CURRENT" >&2
    echo "Delete $OUT_DIR, or give this run its own version name." >&2
    exit 3
  fi
fi

if [ ! -s "$VCFGZ" ]; then
  echo "==> Annotating $(grep -vc '^#' "$INPUT") variants with --acmg ($VERSION)..."
  T0=$(date +%s)
  # Written under a temporary name and moved into place only once the whole
  # pipeline has succeeded. A failed annotate - a mistyped flag, say - used to
  # leave a valid, empty bgzip under the final name, which the next run found,
  # reused and scored as zero classified variants.
  rm -f "$VCFGZ.partial"
  "$FASTVEP" annotate \
    -i "$INPUT" \
    -o - \
    --gff3 "$GFF3" \
    --fasta "$FASTA" \
    --sa-dir "$SA_DIR" \
    --acmg --pick \
    ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} \
    ${EXTRA_ARR[@]+"${EXTRA_ARR[@]}"} \
    --output-format vcf \
    | bgzip -c > "$VCFGZ.partial"
  mv "$VCFGZ.partial" "$VCFGZ"
  T1=$(date +%s)
  echo "==> Annotation took $((T1-T0)) seconds"
  {
    echo "$CURRENT"
    echo "fastvep: $FASTVEP"
    echo "fastvep_version: $("$FASTVEP" --version 2>&1 | head -1)"
    echo "repo_commit: $(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "repo_dirty_files: $(git -C "$ROOT" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
    echo "annotation_seconds: $((T1-T0))"
    echo "recorded_at: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "$PROVENANCE"
elif [ ! -f "$PROVENANCE" ]; then
  echo "==> Reusing existing annotation at $VCFGZ, skipping. It predates run_provenance.txt, so what produced it is not recorded."
else
  echo "==> Reusing existing annotation at $VCFGZ, skipping. Delete it to re-annotate."
fi
ls -la "$VCFGZ"

echo "==> Computing concordance..."
python3 "$SCRIPT_DIR/03_evaluate_concordance.py" \
  --truth "$TRUTH" \
  --predictions "$VCFGZ" \
  --out "$OUT_DIR"

echo "==> Done. See $OUT_DIR/"
ls -la "$OUT_DIR/"
echo
echo "Record what changed in analysis/acmg_benchmark/RUN_VERSIONS.md."
