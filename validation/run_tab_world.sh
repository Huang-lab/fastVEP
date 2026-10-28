#!/usr/bin/env bash
set -euo pipefail

# Run Ensembl VEP and fastVEP over one generated world and compare their tab files cell by
# cell. The measurements in docs/VEP_DIVERGENCE.md ("Part 3") are made this way.
#
# Three runs, each with VEP and with fastVEP:
#   default  the 17 columns a run without any flag has
#   flags    every optional column fastVEP has (--symbol ... --hgvsg), the world as it is
#   alone    the same, with every ALT of a multi-allelic record made a record of its own, so
#            that what is left is about the alleles and not about how the record was written
# then, for each tool, whether an allele is described the same in its record as alone.
#
# The world is crates/fastvep-cli/fixtures/vep_tab_world/gen_world.py's, with --all so that
# nothing the committed fixture leaves out is left out here.
#
# Prerequisites: cargo build --release, Docker, bgzip and tabix (htslib).
#
# Usage:
#   ./validation/run_tab_world.sh [N [SEED [outdir]]]     # 3000 records, seed 139

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
N="${1:-3000}"
SEED="${2:-139}"
OUT="${3:-$(mktemp -d)}"
# release_116.3, by digest, as the fixtures are made
VEP_IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FASTVEP="${FASTVEP:-$PROJECT_DIR/target/release/fastvep}"
GEN="$PROJECT_DIR/crates/fastvep-cli/fixtures/vep_tab_world/gen_world.py"

[[ -x "$FASTVEP" ]] || { echo "missing $FASTVEP: cargo build --release" >&2; exit 1; }
mkdir -p "$OUT"
python3 "$GEN" "$OUT" "$N" "$SEED" --all
python3 "$SCRIPT_DIR/split_alts.py" "$OUT/in.vcf" "$OUT/alone.vcf"
(grep '^#' "$OUT/genes.gff3"; grep -v '^#' "$OUT/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$OUT/genes.sorted.gff3.gz"
tabix -f -p gff "$OUT/genes.sorted.gff3.gz"
docker image inspect "$VEP_IMAGE" >/dev/null 2>&1 || docker pull "$VEP_IMAGE"

vep_flags=(--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg)
our_flags=(--symbol --biotype --canonical --numbers --variant-class --show-ref-allele --hgvs --hgvsg)

# run NAME INPUT.vcf [flags...]: both tools, writing vep_NAME.tab and fastvep_NAME.tab
run() {
    local name="$1" input="$2" with_flags="$3"
    local extra_vep=() extra_ours=()
    if [[ $with_flags == yes ]]; then extra_vep=("${vep_flags[@]}"); extra_ours=("${our_flags[@]}"); fi
    docker run --rm -e VEP_DATABASE=0 -v "$OUT:/data" "$VEP_IMAGE" vep \
        --input_file "/data/$input" --gff /data/genes.sorted.gff3.gz --fasta /data/genome.fa \
        --output_file "/data/vep_$name.tab" --tab --force_overwrite --no_stats \
        ${extra_vep[@]+"${extra_vep[@]}"} >/dev/null 2>&1
    "$FASTVEP" annotate -i "$OUT/$input" --gff3 "$OUT/genes.gff3" --fasta "$OUT/genome.fa" \
        --output-format tab -o "$OUT/fastvep_$name.tab" ${extra_ours[@]+"${extra_ours[@]}"} >/dev/null 2>&1
}

banner() { echo; echo "================================================================"; echo "  $1"; echo "================================================================"; }

# SOURCE and the column named after the GFF3 are two columns VEP adds for --gff
IGNORE=SOURCE,genes.sorted.gff3.gz

run default in.vcf no
banner "default: the 17 columns"
python3 "$SCRIPT_DIR/compare_tab.py" "$OUT/vep_default.tab" "$OUT/fastvep_default.tab" \
    --ignore "$IGNORE" --out "$OUT/diff_default"

run flags in.vcf yes
banner "flags: every optional column"
python3 "$SCRIPT_DIR/compare_tab.py" "$OUT/vep_flags.tab" "$OUT/fastvep_flags.tab" \
    --ignore "$IGNORE" --out "$OUT/diff_flags"

if [[ -f "$OUT/diff_flags/diff_HGVSg.tsv" ]]; then
    banner "HGVSg where the two differ: is the description of the variant?"
    python3 "$SCRIPT_DIR/check_hgvsg_roundtrip.py" "$OUT" "$OUT/diff_flags/diff_HGVSg.tsv"
fi

run alone alone.vcf yes
banner "alone: every ALT its own record, paired on record and transcript"
python3 "$SCRIPT_DIR/compare_tab.py" "$OUT/vep_alone.tab" "$OUT/fastvep_alone.tab" \
    --ignore "$IGNORE" --key Uploaded_variation,Feature --out "$OUT/diff_alone"

banner "an allele in its record against the same allele alone, Ensembl VEP"
python3 "$SCRIPT_DIR/compare_alone.py" "$OUT/in.vcf" "$OUT/vep_flags.tab" "$OUT/vep_alone.tab"
banner "an allele in its record against the same allele alone, fastVEP"
python3 "$SCRIPT_DIR/compare_alone.py" "$OUT/in.vcf" "$OUT/fastvep_flags.tab" "$OUT/fastvep_alone.tab"
echo
echo "Per-column disagreements are in $OUT/diff_<run>/diff_<column>.tsv"
