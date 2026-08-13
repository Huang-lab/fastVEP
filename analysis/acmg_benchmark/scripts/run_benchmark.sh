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

for f in "$INPUT" "$TRUTH" "$GFF3" "$FASTA" "$FASTVEP"; do
  [ -e "$f" ] || { echo "missing required input: $f" >&2; exit 1; }
done
[ -d "$SA_DIR" ] || { echo "missing SA database directory: $SA_DIR" >&2; exit 1; }

mkdir -p "$OUT_DIR"
VCFGZ="$OUT_DIR/clinvar_2star.fastvep.vcf.gz"

if [ ! -s "$VCFGZ" ]; then
  echo "==> Annotating $(grep -vc '^#' "$INPUT") variants with --acmg ($VERSION)..."
  T0=$(date +%s)
  "$FASTVEP" annotate \
    -i "$INPUT" \
    -o - \
    --gff3 "$GFF3" \
    --fasta "$FASTA" \
    --sa-dir "$SA_DIR" \
    --acmg --pick \
    ${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"} \
    --output-format vcf \
    | bgzip -c > "$VCFGZ"
  T1=$(date +%s)
  echo "==> Annotation took $((T1-T0)) seconds"
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
