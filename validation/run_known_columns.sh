#!/usr/bin/env bash
# Run Ensembl VEP 116.3 and fastVEP over one VCF with the same known-variant flags, against the same
# VEP cache, and compare the columns they write for the known variants cell by cell.
#
#   validation/run_known_columns.sh CACHE_DIR GFF3 FASTA VCF OUTDIR [FLAG ...]
#
# CACHE_DIR  what VEP calls --dir_cache: the directory that holds homo_sapiens/116_GRCh38 (or
#            homo_sapiens_refseq/116_GRCh38 with --refseq). Its {chr}/all_vars.gz files are the known
#            variants for both tools; VEP reads its transcripts from the same cache and fastVEP
#            reads GFF3 ones, which is why the comparison pairs rows on (variant, allele) and not on
#            the transcript (compare_known_columns.py).
# GFF3, FASTA  fastVEP's gene models, and VEP's reference sequence.
# FLAG ...   the flags of the run, spelled as VEP spells them (`--af_1kg --check_existing`): fastVEP
#            accepts that spelling, so the same list serves both. Flags only VEP has (`--refseq
#            --exclude_predicted`) go to VEP alone through VEP_ONLY, and the other way round through
#            FASTVEP_ONLY, each a string of flags.
#
# Needs docker and python3, and a `fastvep` on the PATH (or FASTVEP=path). The image is VEP's
# release_116.3, by digest, so that a moved tag cannot change the answers unnoticed.
#
# Example, the case of issue #139:
#   validation/run_known_columns.sh ~/.vep Homo_sapiens.GRCh38.115.gff3 GRCh38.fa calls.vcf out \
#       --af_1kg --check_existing
set -euo pipefail

IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
here="$(cd "$(dirname "$0")" && pwd)"

if [[ $# -lt 5 ]]; then
  sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi
cache="$(cd "$1" && pwd)"
gff3="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
fasta="$(cd "$(dirname "$3")" && pwd)/$(basename "$3")"
vcf="$(cd "$(dirname "$4")" && pwd)/$(basename "$4")"
mkdir -p "$5"
out="$(cd "$5" && pwd)"
shift 5
flags=("$@")

species_dir="homo_sapiens"
if [[ " ${VEP_ONLY:-} " == *" --refseq "* ]]; then species_dir="homo_sapiens_refseq"; fi

# shellcheck disable=SC2086
docker run --rm \
  -v "$(dirname "$fasta")":/ref:ro -v "$cache":/cache:ro -v "$(dirname "$vcf")":/in:ro -v "$out":/out \
  "$IMAGE" vep --cache --offline --dir_cache /cache --assembly GRCh38 \
  --fasta "/ref/$(basename "$fasta")" --input_file "/in/$(basename "$vcf")" \
  --output_file /out/vep.tab --tab --force_overwrite --no_stats --fork 3 \
  "${flags[@]}" ${VEP_ONLY:-}

# shellcheck disable=SC2086
"${FASTVEP:-fastvep}" annotate -i "$vcf" --gff3 "$gff3" \
  --cache-dir "$cache/$species_dir/116_GRCh38" --output-format tab --no-progress \
  -o "$out/fv.tab" "${flags[@]}" ${FASTVEP_ONLY:-}

python3 "$here/compare_known_columns.py" "$out/vep.tab" "$out/fv.tab" --out "$out/diff"
