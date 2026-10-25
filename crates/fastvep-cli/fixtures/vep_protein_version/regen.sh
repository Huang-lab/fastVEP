#!/usr/bin/env bash
# Regenerates golden.tsv: the HGVSp Ensembl's own VEP writes for in.vcf over genes.gff3.
#
#   ./regen.sh            # keeps the committed world
#   ./regen.sh remake     # also rebuilds genome.fa, genes.gff3 and in.vcf from gen_world.py
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3, the release of the cache in issue #139): `vep --tab --hgvs --gff --fasta`.
# VEP 116 reads a protein's version from the `version` of its first CDS record, and
# 115 did not (it gave every translation version 1). Needs docker, bgzip and tabix.
#
# golden.tsv keeps the rows that have an HGVSp: (Uploaded_variation, Feature, HGVSp).
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"

if [[ "${1:-}" == remake ]]; then
    python3 "$here/gen_world.py" "$here"
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cp "$here/in.vcf" "$here/genes.gff3" "$here/genome.fa" "$tmp/"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
  --input_file /data/in.vcf --output_file /data/vep.tab --tab --force_overwrite --no_stats \
  --gff /data/genes.sorted.gff3.gz --fasta /data/genome.fa --hgvs >/dev/null

python3 - "$tmp/vep.tab" "$here/golden.tsv" "$IMAGE" <<'PY'
import sys
tab, out, image = sys.argv[1:4]
version = "?"
header, rows = None, []
for line in open(tab):
    line = line.rstrip("\n")
    if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
        version = line.split()[-1]
    elif line.startswith("##"):
        continue
    elif line.startswith("#"):
        header = line[1:].split("\t")
    else:
        r = dict(zip(header, line.split("\t")))
        if r["HGVSp"] != "-":
            rows.append((r["Uploaded_variation"], r["Feature"], r["HGVSp"]))
with open(out, "w") as f:
    f.write(f"# HGVSp of Ensembl VEP {version} for in.vcf over genes.gff3: record, transcript, HGVSp\n")
    f.write(f"# image {image}\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
print(f"{len(rows)} rows, VEP {version}")
PY
