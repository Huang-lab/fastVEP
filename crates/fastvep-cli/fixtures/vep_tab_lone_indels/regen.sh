#!/usr/bin/env bash
# Regenerates this directory: lone indels written with bases the REF and ALT share beyond the
# anchor (in.vcf, from gen_lone.py), and what Ensembl's own VEP writes for them with every
# optional column fastVEP has (golden_flags.tab.gz).
#
#   ./regen.sh [N [SEED]]          N records (default 1500), seed 139
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below (release_116.3,
# the release of the cache in issue #139), in --gff mode over ../vep_tab_world's genes and genome.
# Needs docker, python3, bgzip and tabix (htslib). The two columns VEP adds for --gff (SOURCE and
# the file name) are dropped from the golden.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
world="$here/../vep_tab_world"
N="${1:-1500}"
SEED="${2:-139}"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FLAGS="--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$here/gen_lone.py" "$tmp/in.vcf" "$N" "$SEED"
cp "$world/genome.fa" "$world/genes.gff3" "$tmp/"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

# shellcheck disable=SC2086
docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
  --input_file /data/in.vcf --output_file /data/flags.tab --gff /data/genes.sorted.gff3.gz \
  --fasta /data/genome.fa --force_overwrite --no_stats --tab $FLAGS >/dev/null 2>&1

cp "$tmp/in.vcf" "$here/in.vcf"

python3 - "$tmp" "$here" "$IMAGE" "$N" "$SEED" <<'PY'
import gzip, sys
tmp, here, image, n, seed = sys.argv[1:6]
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
    f.write(f"## the tab file of Ensembl VEP {version}; image {image}; gen_lone.py N={n} SEED={seed}\n".encode())
    for row in out:
        cells = [c for i, c in enumerate(row.split("\t")) if i not in drop]
        f.write(("\t".join(cells) + "\n").encode())
print(f"VEP {version}: {len(out) - 1} rows")
PY
