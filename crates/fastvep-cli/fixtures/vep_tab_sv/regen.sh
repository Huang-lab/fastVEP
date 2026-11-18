#!/usr/bin/env bash
# Regenerates this directory: the world (genome.fa, genes.gff3, from gen_world.py), the structural
# variants and other changes in it (in.vcf), and what Ensembl's own VEP writes for them with every
# optional column fastVEP has (golden_flags.tab.gz).
#
#   ./regen.sh
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below (release_116.3),
# in --gff mode: its gene models are genes.gff3 and its sequence genome.fa. Needs docker, python3,
# bgzip and tabix (htslib). The two columns VEP adds for --gff (SOURCE and the file name) are dropped
# from the golden.
#
# VEP is run as it is by default, with its statistics (the HTML file lands in the temporary
# directory and is dropped), and not with --no_stats as the other fixtures do: the two runs write
# different rows for a breakend that has no mate, and the default run's are the ones the code is
# written to give (docs/STRUCTURAL_VARIANTS.md).
#
# --distance is 1000, not VEP's 5000: the world is 14 kb with ten transcripts, and at 5000 every
# variant is a neighbour of nearly every transcript, so most of what it asks would be answered by
# the flanks.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FLAGS="--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg --total_length"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$here/gen_world.py" "$tmp"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

# shellcheck disable=SC2086
docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
  --input_file /data/in.vcf --output_file /data/flags.tab --gff /data/genes.sorted.gff3.gz \
  --fasta /data/genome.fa --distance 1000 --force_overwrite --tab $FLAGS >/dev/null 2>&1

cp "$tmp/genome.fa" "$tmp/genes.gff3" "$tmp/in.vcf" "$here/"

python3 - "$tmp" "$here" "$IMAGE" <<'PY'
import gzip, sys
tmp, here, image = sys.argv[1:4]
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
drop = {i for i, name in enumerate(header) if name == "SOURCE" or name.endswith(".gz")}
with gzip.GzipFile(f"{here}/golden_flags.tab.gz", "wb", mtime=0) as f:
    f.write(f"## the tab file of Ensembl VEP {version}; image {image}; gen_world.py\n".encode())
    for row in out:
        cells = [c for i, c in enumerate(row.split("\t")) if i not in drop]
        f.write(("\t".join(cells) + "\n").encode())
print(f"VEP {version}: {len(out) - 1} rows")
PY
