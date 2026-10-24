#!/usr/bin/env bash
# Regenerates this directory: the VCF, the VEP-format variation cache, a one-gene
# GFF3 and golden.tsv, which is what Ensembl's own VEP answers for them.
#
#   ./regen.sh [N [SEED]]          N pairs (default 1200), default seed 139
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3, the release the cache of issue #139 is for): `vep --tab
# --check_existing --cache --offline` over the cache this script builds, with the
# one-gene GFF3 standing in for transcripts. Needs docker, python3, bgzip and tabix
# (htslib). The genome VEP insists on is generated, used and thrown away; only the
# inputs fastVEP needs and the answers are kept.
#
# The cache is built the way VEP's own converter builds it
# (ensembl-vep/convert_cache.pl): bgzip the file, then
#
#     tabix -C -s 1 -b 5 -e 5 all_vars.gz
#
# so it has a .csi and no .tbi.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
N="${1:-1200}"
SEED="${2:-139}"
# release_116.3, by digest, so a moved tag cannot change the answers unnoticed
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$here/gen_world.py" "$tmp" "$N" "$SEED"
cache="$tmp/cache/homo_sapiens/116_GRCh38"
bgzip -c "$tmp/cache_10.txt" > "$cache/10/all_vars.gz"
tabix -C -s 1 -b 5 -e 5 "$cache/10/all_vars.gz"
# VEP reads a GFF3 through tabix and rejects an unindexed one
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

# VEP_DATABASE=0 rather than --offline alone: see validation/run_divergence.sh
docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
  --input_file /data/in.vcf --output_file /data/vep.tab --tab --force_overwrite --no_stats \
  --offline --cache --dir_cache /data/cache --species homo_sapiens --assembly GRCh38 --cache_version 116 \
  --gff /data/genes.sorted.gff3.gz --fasta /data/genome.fa --check_existing >/dev/null

rm -rf "$here/cache"
cp -R "$tmp/cache" "$here/cache"
cp "$tmp/in.vcf" "$tmp/genes.gff3" "$here/"

python3 - "$tmp/vep.tab" "$here/golden.tsv" "$IMAGE" "$N" "$SEED" <<'PY'
import collections, sys
tab, out, image, n, seed = sys.argv[1:6]
header, rows, version = None, collections.defaultdict(list), "?"
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
        rows[row["Uploaded_variation"]].append(row["Existing_variation"])
with open(out, "w") as f:
    f.write(f"# Existing_variation per input record, as Ensembl VEP {version} prints it\n")
    f.write(f"# image {image}\n# gen_world.py N={n} SEED={seed}; vep --tab --check_existing --cache --offline\n")
    f.write("# a record with several rows (alts) lists their cells sorted, joined with '|'\n")
    for rid in sorted(rows, key=lambda r: int(r[2:])):
        f.write(f"{rid}\t{'|'.join(sorted(rows[rid]))}\n")
print(f"{len(rows)} records, VEP {version}")
PY
