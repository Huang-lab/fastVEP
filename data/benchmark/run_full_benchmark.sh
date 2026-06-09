#!/usr/bin/env bash
# End-to-end ClinVar 2-star+ benchmark.
#
# Inputs:
#   - data/benchmark/clinvar_2star.vcf            (ClinVar 2-star+ VCF, ~673k SNV/small-indel)
#   - data/benchmark/clinvar_2star_truth.tsv      (truth table for concordance)
#   - data/benchmark/sa_db/                       (built .osa/.oga databases)
#   - test_data/organisms/human/                  (FASTA + GFF3 + cache)
#
# Outputs:
#   - data/benchmark/output_v7/clinvar_2star.fastvep.vcf.gz   (annotation + ACMG, bgzipped)
#   - data/benchmark/output_v7/concordance_summary.txt
#   - data/benchmark/output_v7/concordance_matrix.csv
#   - data/benchmark/output_v7/concordance_by_chrom.csv
#   - data/benchmark/output_v7/concordance_by_consequence.csv
#   - data/benchmark/output_v7/criterion_firing_rates.csv
#   - data/benchmark/output_v7/discrepancies.tsv
#
# Output is VCF-bgzipped (~10× smaller than the prior pretty-printed JSON
# and ~100× smaller than tab format with all fields). VCF is the only
# format that includes ACMG/ACMG_CRITERIA in the per-transcript CSQ entry.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
INPUT=$ROOT/data/benchmark/clinvar_2star.vcf
GFF3=$ROOT/test_data/organisms/human/Homo_sapiens.GRCh38.115.gff3
FASTA=$ROOT/test_data/organisms/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
SA_DIR=$ROOT/data/benchmark/sa_db
OUT_DIR=$ROOT/data/benchmark/output_v7
mkdir -p "$OUT_DIR"
VCFGZ="$OUT_DIR/clinvar_2star.fastvep.vcf.gz"

if [ ! -s "$VCFGZ" ]; then
  echo "==> Annotating $(grep -vc '^#' $INPUT) variants with --acmg..."
  T0=$(date +%s)
  $ROOT/target/release/fastvep annotate \
    -i "$INPUT" \
    -o - \
    --gff3 "$GFF3" \
    --fasta "$FASTA" \
    --sa-dir "$SA_DIR" \
    --acmg --pick \
    --output-format vcf \
    | bgzip -c > "$VCFGZ"
  T1=$(date +%s)
  echo "==> Annotation took $((T1-T0)) seconds"
  ls -la "$VCFGZ"
else
  echo "==> Annotation already done at $VCFGZ ($(stat -c %s $VCFGZ) bytes), skipping."
fi

echo "==> Computing concordance..."
python3 $ROOT/analysis/acmg_benchmark/real_data/03_evaluate_concordance.py \
  --truth $ROOT/data/benchmark/clinvar_2star_truth.tsv \
  --predictions "$VCFGZ" \
  --out "$OUT_DIR"

echo "==> Done. See $OUT_DIR/"
ls -la "$OUT_DIR/"
