#!/usr/bin/env bash
# Regenerates this directory: a genome of repeats and edges (genome.fa), its variants
# (in.vcf), a one-gene-per-contig GFF3, and golden.tsv, which is the HGVS g. description
# Ensembl's own VEP writes for each allele of each record.
#
#   ./regen.sh [SEED]
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3): `vep --tab --hgvsg --allele_number` in --gff mode. Needs docker,
# python3, bgzip and tabix (htslib).
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
SEED="${1:-139}"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$here/gen_world.py" "$tmp" "$SEED"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
  --input_file /data/in.vcf --output_file /data/vep.tab --gff /data/genes.sorted.gff3.gz \
  --fasta /data/genome.fa --force_overwrite --no_stats --tab --hgvsg --allele_number >/dev/null 2>&1

cp "$tmp/genome.fa" "$tmp/genes.gff3" "$tmp/in.vcf" "$here/"

python3 - "$tmp/vep.tab" "$here/golden.tsv" "$IMAGE" "$SEED" <<'PY'
import sys
tab, out, image, seed = sys.argv[1:5]
header, version, golden = None, "?", {}
for line in open(tab):
    line = line.rstrip("\n")
    if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
        version = line.split()[-1]
    elif line.startswith("##"):
        continue
    elif line.startswith("#"):
        header = line[1:].split("\t")
    else:
        row = dict(zip(header, line.split("\t")))
        golden[(row["Uploaded_variation"], int(row["ALLELE_NUM"]))] = row["HGVSg"]
with open(out, "w") as f:
    f.write(f"# HGVSg of each allele (ALLELE_NUM, 1 is the first alt) of each record, as Ensembl VEP {version} writes it\n")
    f.write(f"# image {image}; gen_world.py SEED={seed}; `-` is VEP writing nothing\n")
    for (rid, n), hgvsg in sorted(golden.items(), key=lambda kv: (int(kv[0][0][1:]), kv[0][1])):
        f.write(f"{rid}\t{n}\t{hgvsg}\n")
print(f"{len(golden)} alleles, VEP {version}")
PY
