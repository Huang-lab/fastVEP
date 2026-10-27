#!/usr/bin/env bash
# Ensembl VEP 116.3 on its own cache, and fastVEP over the matching GFF3, compared cell by cell.
#
#   run_cache_oracle.sh ensembl|refseq INPUT.vcf OUTDIR
#
# Environment (all required):
#   VEP_CACHE   the directory that holds homo_sapiens[_refseq]/116_GRCh38 (--dir_cache); chr21 and chr22 of it are
#               ~700 MB each of the two caches, which is what the numbers in docs/VEP_DIVERGENCE.md were measured on
#   GFF3        Ensembl's GFF3 for the same release, or RefSeq's with its contigs named as the FASTA names them
#               (RefSeq: keep NM_ and NR_ only, which is what `--exclude_predicted` does to the cache)
#   FASTA       the GRCh38 primary assembly the cache was built against
#   FASTVEP     the fastvep binary (default: fastvep)
#
# A RefSeq cache holds transcripts that VEP edits with an alignment (`_rna_edit`, BAM_EDIT=OK), which no GFF3 can
# reproduce; they are left out of the comparison. An Ensembl cache is built without readthrough transcripts and
# transcripts of biotype `artifact`, so fastVEP is told to leave them out too.
set -euo pipefail
# `${array[@]+"${array[@]}"}` below, for the arrays that are empty for one of the kinds: macOS's
# /bin/bash is 3.2, where "${empty[@]}" is an unbound variable under `set -u`.
kind="$1"; vcf="$2"; out="$3"
: "${VEP_CACHE:?}" "${GFF3:?}" "${FASTA:?}"
FASTVEP="${FASTVEP:-fastvep}"
IMAGE=ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3
here="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$out"
columns=(--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg)
case "$kind" in
  ensembl) vep_extra=(); fv_extra=(--exclude-readthrough --exclude-artifact); drop=() ;;
  refseq)  vep_extra=(--refseq --exclude_predicted); fv_extra=(--exclude-predicted); drop=(--drop-tx-where BAM_EDIT=OK) ;;
  *) echo "kind is ensembl or refseq" >&2; exit 2 ;;
esac
docker run --rm -e VEP_DATABASE=0 \
  -v "$(dirname "$FASTA"):/ref:ro" -v "$VEP_CACHE:/cache:ro" -v "$(cd "$(dirname "$vcf")" && pwd):/in:ro" -v "$(cd "$out" && pwd):/out" \
  "$IMAGE" vep --fork 2 --cache --offline --assembly GRCh38 --dir_cache /cache ${vep_extra[@]+"${vep_extra[@]}"} \
  --fasta "/ref/$(basename "$FASTA")" --input_file "/in/$(basename "$vcf")" --output_file /out/vep.tab \
  --tab --force_overwrite --no_stats "${columns[@]}"
"$FASTVEP" annotate -i "$vcf" --gff3 "$GFF3" --fasta "$FASTA" --output-format tab -o "$out/fastvep.tab" \
  --symbol --biotype --canonical --numbers --variant-class --show-ref-allele --hgvs --hgvsg "${fv_extra[@]}"
python3 "$here/compare_cache_run.py" "$out/vep.tab" "$out/fastvep.tab" "$out/diff" ${drop[@]+"${drop[@]}"}
