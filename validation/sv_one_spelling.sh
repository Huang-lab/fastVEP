#!/usr/bin/env bash
# Ensembl VEP 116.3 against fastVEP on a structural variant VCF whose contig is spelled one way throughout.
#
#   validation/sv_one_spelling.sh WORKDIR VCF CONTIG
#
# Ensembl VEP compares contig names letter for letter (`_compare_seq_region_names`), so a VCF that spells a
# contig two ways (gnomAD-SV: `21` in the contig column, `chr21` in `CHR2`) gets the rows of one spelling and
# not of the other, whichever the GFF3 uses (docs/STRUCTURAL_VARIANTS.md, section 10), and its --gff reader pairs
# a variant with transcripts of other contigs. This makes both go away: it cuts the GFF3 and the FASTA down to
# CONTIG, names it `chr` + CONTIG in them and in the VCF's contig column, and runs both tools over that.
# Every record whose mate is on another contig then has a mate that the reference does not have.
#
# Read from the environment:
#   FASTVEP       the fastVEP binary
#   GFF3          the gene models, bgzipped and tabix-indexed, with the contigs named without `chr`
#   FASTA         the reference, indexed, with the same names
#   VEP_IMAGE     the image of Ensembl VEP 116.3, by digest
# Needs docker, bgzip, tabix and samtools. VEP is run as it is by default, with its statistics.
#
# Writes WORKDIR/ref/in.tab (Ensembl VEP) and WORKDIR/out/in.tab (fastVEP), and prints the strict F1 of the
# vep-rs paper for them (validation/score_strict_key.py).
set -euo pipefail
export LC_ALL=C

work="${1:?usage: sv_one_spelling.sh WORKDIR VCF CONTIG}"
vcf="${2:?}"
contig="${3:?}"
: "${FASTVEP:?set FASTVEP}" "${GFF3:?set GFF3}" "${FASTA:?set FASTA}" "${VEP_IMAGE:?set VEP_IMAGE}"
here="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$work/ref" "$work/out"
work="$(cd "$work" && pwd)"
name="chr$contig"

tabix -h "$GFF3" "$contig" | awk -v n="$name" 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$1=n; print}' | bgzip -c > "$work/ref.gff3.gz"
tabix -p gff "$work/ref.gff3.gz"
samtools faidx "$FASTA" "$contig" | sed "1s/.*/>$name/" > "$work/ref.fa"
samtools faidx "$work/ref.fa"
awk -v n="$name" 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$1=n; print}' "$vcf" > "$work/in.vcf"

docker run --rm -e VEP_DATABASE=0 -v "$work:/data" "$VEP_IMAGE" vep --input_file /data/in.vcf \
  --gff /data/ref.gff3.gz --fasta /data/ref.fa --output_file /data/ref/in.tab --tab --force_overwrite \
  > "$work/vep.log" 2>&1
"$FASTVEP" annotate -i "$work/in.vcf" --gff3 "$work/ref.gff3.gz" --fasta "$work/ref.fa" \
  --output-format tab --no-progress -o "$work/out/in.tab" > "$work/fastvep.log" 2>&1

python3 -I "$here/score_strict_key.py" "$work/ref" "$work/out"
