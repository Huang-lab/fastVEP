#!/usr/bin/env bash
# Regenerates this directory: a RefSeq-style world (genome.fa, genes.gff3, in.vcf, from gen_world.py)
# and what Ensembl's own VEP writes for it, and for the edge fixture's world, under the flags
# that change how the text is written:
#
#   golden_no_escape.tsv.gz     the Feature, SYMBOL and HGVSp of every CSQ entry of vep --vcf over
#                               this world, as it is and with --no_escape
#   golden_total_length.tsv.gz  the three position cells of every row of vep --tab --total_length
#                               over ../vep_tab_edges/in.vcf in the world of ../vep_tab_world
#
#   ./regen.sh
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below (release_116.3,
# the release of the cache in issue #139), in --gff mode. Needs docker, python3, bgzip and tabix
# (htslib).
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FLAGS="--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg"

python3 "$here/gen_world.py" "$here"

prepare() { # <directory>: sort, compress and index the GFF3 for --gff
  (grep '^#' "$1/genes.gff3"; grep -v '^#' "$1/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
    | bgzip -c > "$1/genes.sorted.gff3.gz"
  tabix -p gff "$1/genes.sorted.gff3.gz"
}
vep() { # <directory> <output> <vep arguments...>
  local dir="$1" out="$2"; shift 2
  docker run --rm -e VEP_DATABASE=0 -v "$dir:/data" "$IMAGE" vep \
    --input_file /data/in.vcf --output_file "/data/$out" --gff /data/genes.sorted.gff3.gz \
    --fasta /data/genome.fa --force_overwrite --no_stats "$@" >/dev/null 2>&1
}

compat="$(mktemp -d)"
edges="$(mktemp -d)"
trap 'rm -rf "$compat" "$edges"' EXIT
cp "$here/genome.fa" "$here/genes.gff3" "$here/in.vcf" "$compat/"
cp "$here/../vep_tab_world/genome.fa" "$here/../vep_tab_world/genes.gff3" "$here/../vep_tab_edges/in.vcf" "$edges/"
prepare "$compat"
prepare "$edges"

# shellcheck disable=SC2086
vep "$compat" escaped.vcf --vcf $FLAGS
# shellcheck disable=SC2086
vep "$compat" unescaped.vcf --vcf $FLAGS --no_escape
vep "$edges" total.tab --tab --total_length

python3 - "$compat" "$edges" "$here" "$IMAGE" <<'PY'
import gzip, re, sys
compat, edges, here, image = sys.argv[1:5]

def tab(path):
    header, rows, version = None, [], "?"
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
            version = line.split()[-1]
        elif line.startswith("##"):
            continue
        elif line.startswith("#"):
            header = line[1:].split("\t")
        else:
            rows.append(line.split("\t"))
    return header, rows, version

def write(name, lines, note):
    with gzip.GzipFile(f"{here}/{name}", "wb", mtime=0) as f:
        f.write(f"## {note}; image {image}\n".encode())
        for line in lines:
            f.write((line + "\n").encode())

header, rows, version = tab(f"{edges}/total.tab")

# the CSQ of each record, as it is and with --no_escape
def csq(path):
    out, fields = [], None
    for line in open(path):
        if line.startswith("##INFO=<ID=CSQ"):
            fields = re.search(r'Format: ([^"]+)"', line).group(1).split("|")
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        info = [kv[4:] for kv in cols[7].split(";") if kv.startswith("CSQ=")]
        for entry in (info[0].split(",") if info else []):
            f = dict(zip(fields, entry.split("|")))
            out.append((cols[2], f["Allele"], f["Feature"], f["SYMBOL"], f["HGVSp"]))
    return out
lines = ["#mode\tUploaded_variation\tAllele\tFeature\tSYMBOL\tHGVSp"]
for mode, path in (("escaped", "escaped.vcf"), ("no_escape", "unescaped.vcf")):
    for r in csq(f"{compat}/{path}"):
        lines.append("\t".join((mode,) + r))
write("golden_no_escape.tsv.gz", lines, f"the CSQ of Ensembl VEP {version} --vcf over this world, as it is and with --no_escape")

# the position cells of --total_length
cols = ["Uploaded_variation", "Allele", "Feature", "cDNA_position", "CDS_position", "Protein_position"]
at = [header.index(c) for c in cols]
write("golden_total_length.tsv.gz",
      ["#" + "\t".join(cols)] + ["\t".join(r[i] for i in at) for r in rows],
      f"the position cells of Ensembl VEP {version} --tab --total_length over ../vep_tab_edges/in.vcf in the world of ../vep_tab_world")
print(f"VEP {version}")
PY
