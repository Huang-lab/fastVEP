#!/usr/bin/env bash
# Regenerates golden.tsv: what Ensembl's own VEP writes for in.vcf over a RefSeq-style genes.gff3.
#
#   ./regen.sh            # keeps the committed world
#   ./regen.sh remake     # also rebuilds genome.fa, genes.gff3 and in.vcf from gen_world.py
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3, the release of the cache in issue #139): `vep --vcf --gff --fasta`.
# VEP's GFF reader names a transcript, a gene and a biotype by rules of its own, and
# a RefSeq GFF3 is where they show (see gen_world.py). Needs docker, bgzip and tabix.
#
# golden.tsv keeps, for every CSQ entry VEP writes, the fields that name things and the
# ones that depend on the name: record, Allele, Feature, Gene, SYMBOL, BIOTYPE,
# Consequence, HGVSc, HGVSp.
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
  --input_file /data/in.vcf --output_file /data/vep.vcf --vcf --force_overwrite --no_stats \
  --gff /data/genes.sorted.gff3.gz --fasta /data/genome.fa \
  --symbol --biotype --hgvs >/dev/null

python3 - "$tmp/vep.vcf" "$here/golden.tsv" "$IMAGE" <<'PY'
import re, sys
vcf, out, image = sys.argv[1:4]
version = "?"
fields = None
rows = []
for line in open(vcf):
    if line.startswith("##VEP="):
        version = re.search(r'VEP="?(v?[0-9.]+)', line).group(1)
    if line.startswith("##INFO=<ID=CSQ"):
        fields = re.search(r"Format: ([^\"]+)\"", line).group(1).split("|")
    if line.startswith("#"):
        continue
    cols = line.rstrip("\n").split("\t")
    csq = [kv[4:] for kv in cols[7].split(";") if kv.startswith("CSQ=")][0]
    for entry in csq.split(","):
        f = dict(zip(fields, entry.split("|")))
        rows.append((cols[2], f["Allele"], f["Feature"], f["Gene"], f["SYMBOL"], f["BIOTYPE"],
                     f["Consequence"], f["HGVSc"], f["HGVSp"]))
with open(out, "w") as o:
    o.write(f"# Ensembl VEP {version} for in.vcf over genes.gff3: record, Allele, Feature, Gene, SYMBOL, BIOTYPE, Consequence, HGVSc, HGVSp\n")
    o.write(f"# image {image}\n")
    for r in rows:
        o.write("\t".join(r) + "\n")
print(f"{len(rows)} rows, VEP {version}")
PY
