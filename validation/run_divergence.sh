#!/usr/bin/env bash
set -euo pipefail

# Re-measure every number in docs/VEP_DIVERGENCE.md against real Ensembl VEP.
#
# Three inputs, because they exercise different code paths: a stratified ClinVar
# 2-star+ sample built to be hard, a 1-in-200 systematic sample of the GIAB
# HG002 callset, which is what an ordinary WGS run looks like, and the 400
# ClinVar in-frame deletions, which are the only ones here that reach
# protein-level 3'-normalisation at a terminus. The two samples are cut by
# validation/sample_variants.py, which is deterministic - the numbers are
# reproducible from the inputs alone, with no seed to carry around.
#
# Prerequisites:
#   - cargo build --release
#   - Docker, and ensemblorg/ensembl-vep:release_115.1 (pulled if missing)
#   - bgzip and tabix (htslib): VEP rejects an unindexed GFF3 in --gff mode
#   - test_data/organisms/human/ (benchmarks/download_data.sh --human)
#   - data/benchmark/clinvar_2star.vcf
#     (analysis/acmg_benchmark/scripts/01_extract_clinvar_2star.py)
#
# Usage:
#   ./validation/run_divergence.sh [outdir]      # default: validation/results

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/results}"
VEP_IMAGE="ensemblorg/ensembl-vep:release_115.1"
HUMAN_DIR="$PROJECT_DIR/test_data/organisms/human"
GFF3="$HUMAN_DIR/Homo_sapiens.GRCh38.115.gff3"
FASTA="$HUMAN_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
SORTED_GFF3="$HUMAN_DIR/Homo_sapiens.GRCh38.115.sorted.gff3.gz"
CLINVAR="$PROJECT_DIR/data/benchmark/clinvar_2star.vcf"
HG002="$HUMAN_DIR/human_giab_hg002_full.vcf"
FASTVEP="$PROJECT_DIR/target/release/fastvep"

cd "$PROJECT_DIR"
mkdir -p "$OUT_DIR"

for f in "$FASTVEP" "$GFF3" "$FASTA"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

# VEP reads the GFF3 through tabix and rejects an unindexed one outright, so it
# gets a bgzipped copy of the same file fastVEP reads. Sorting is required by
# the index, not by the content.
if [[ ! -f "$SORTED_GFF3.tbi" ]]; then
    echo "--- preparing $SORTED_GFF3 (one time, a few minutes) ---"
    grep -v '^#' "$GFF3" | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -@ 4 -c > "$SORTED_GFF3"
    tabix -p gff "$SORTED_GFF3"
fi

docker image inspect "$VEP_IMAGE" >/dev/null 2>&1 || docker pull "$VEP_IMAGE"

run_vep() {
    local input="$1" output="$2"
    # VEP_DATABASE=0 rather than --offline: --offline sets `cache => 1` and VEP
    # then demands a species cache directory that --gff mode has no use for.
    # --allele_number because VEP does not emit CSQ entries in ALT order on a
    # multi-allelic record; the comparison needs to know which ALT each is.
    docker run --rm -e VEP_DATABASE=0 \
        -v "$PROJECT_DIR:/work" -v "$OUT_DIR:/out" \
        "$VEP_IMAGE" vep \
        --input_file "/out/$(basename "$input")" \
        --gff "/work/${SORTED_GFF3#"$PROJECT_DIR"/}" \
        --fasta "/work/${FASTA#"$PROJECT_DIR"/}" \
        --output_file "/out/$(basename "$output")" \
        --vcf --force_overwrite --no_stats \
        --hgvs --symbol --canonical --allele_number 2>&1 |
        grep -v '^WARNING: Parent entries' | tail -3
}

measure() {
    local name="$1" sample="$2"
    echo ""
    echo "================================================================"
    echo "  $name"
    echo "================================================================"
    "$FASTVEP" annotate -i "$sample" --gff3 "$GFF3" --fasta "$FASTA" \
        --hgvs --symbol --canonical -o "$OUT_DIR/fastvep_$name.vcf" 2>&1 | tail -1
    run_vep "$sample" "$OUT_DIR/vep_$name.vcf"
    python3 "$SCRIPT_DIR/compare_rows.py" \
        "$OUT_DIR/fastvep_$name.vcf" "$OUT_DIR/vep_$name.vcf" \
        --dump "$OUT_DIR/$name"
}

if [[ -f "$CLINVAR" ]]; then
    python3 "$SCRIPT_DIR/sample_variants.py" "$CLINVAR" "$OUT_DIR/clinvar_sample.vcf" --stratified
    measure clinvar "$OUT_DIR/clinvar_sample.vcf"
else
    echo "SKIP: $CLINVAR not found"
fi

if [[ -f "$HG002" ]]; then
    python3 "$SCRIPT_DIR/sample_variants.py" "$HG002" "$OUT_DIR/hg002_sample.vcf" --every 200
    measure hg002 "$OUT_DIR/hg002_sample.vcf"
else
    echo "SKIP: $HG002 not found"
fi

# The only input that exercises protein-level 3'-normalisation at a terminus,
# which neither sample above reaches often enough to measure.
INFRAME="$SCRIPT_DIR/human/clinvar_inframe_deletions.vcf"
if [[ -f "$INFRAME" ]]; then
    cp "$INFRAME" "$OUT_DIR/inframe_deletions.vcf"
    measure inframe_deletions "$OUT_DIR/inframe_deletions.vcf"
else
    echo "SKIP: $INFRAME not found"
fi

echo ""
echo "Done. Per-field disagreements are in $OUT_DIR/<sample>.<field>.tsv"
