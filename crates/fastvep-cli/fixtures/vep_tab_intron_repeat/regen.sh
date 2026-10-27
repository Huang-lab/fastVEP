#!/usr/bin/env bash
# Regenerates this directory: the world (genome.fa, genes.gff3) and its changes (in.vcf), from
# gen_world.py, and what Ensembl's own VEP writes for them with every optional column fastVEP has
# (golden_flags.tab.gz).
#
#   ./regen.sh [SAMPLE]            keep every SAMPLE-th change (default 3); 1 keeps all of them
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3, the release of the cache in issue #139), in --gff mode over this world's gene
# models and sequence, with --distance 0 so that the transcripts of the world do not annotate each
# other's changes as up and downstream. Needs docker, python3, bgzip and tabix (htslib). The golden
# keeps the columns the test reads and no others.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
SAMPLE="${1:-3}"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FLAGS="--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg --distance 0"

python3 "$here/gen_world.py" "$here" "$SAMPLE"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cp "$here/genome.fa" "$here/genes.gff3" "$here/in.vcf" "$tmp/"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

# shellcheck disable=SC2086
docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
  --input_file /data/in.vcf --output_file /data/flags.tab --gff /data/genes.sorted.gff3.gz \
  --fasta /data/genome.fa --force_overwrite --no_stats --tab $FLAGS >/dev/null 2>&1

python3 - "$tmp" "$here" "$IMAGE" "$SAMPLE" <<'PY'
import gzip, sys
tmp, here, image, sample = sys.argv[1:5]
out, version = [], "?"
for line in open(f"{tmp}/flags.tab"):
    line = line.rstrip("\n")
    if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
        version = line.split()[-1]
    elif line.startswith("##"):
        continue
    else:
        out.append(line)
header = out[0][1:].split("\t")
# only the columns the test reads: the rest of the table is the same for every record
keep = {"Uploaded_variation", "Location", "Allele", "Feature", "Consequence", "IMPACT",
        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
        "HGVSc", "HGVSp", "HGVS_OFFSET", "STRAND", "VARIANT_CLASS"}
drop = {i for i, name in enumerate(header) if name not in keep}
with gzip.GzipFile(f"{here}/golden_flags.tab.gz", "wb", mtime=0) as f:
    f.write(f"## the tab file of Ensembl VEP {version}; image {image}; gen_world.py SAMPLE={sample}\n".encode())
    for row in out:
        cells = [c for i, c in enumerate(row.split("\t")) if i not in drop]
        f.write(("\t".join(cells) + "\n").encode())
print(f"VEP {version}: {len(out) - 1} rows")
PY
