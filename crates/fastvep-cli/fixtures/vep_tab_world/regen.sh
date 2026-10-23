#!/usr/bin/env bash
# Regenerates this directory: a generated world (genome.fa, genes.gff3, in.vcf) and what
# Ensembl's own VEP writes for it, in three forms:
#
#   golden_default.tab.gz   vep --tab, no other flag
#   golden_flags.tab.gz     vep --tab with every optional column fastVEP writes
#   golden_flags.vcf.gz     vep --vcf with the same flags (the CSQ of each record)
#
#   ./regen.sh [N [SEED]]          N records (default 500), seed 139
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3, the release of the cache in issue #139), in --gff mode: its gene models
# are genes.gff3, its sequence is genome.fa. Needs docker, python3, bgzip and tabix
# (htslib). Two columns VEP adds for --gff (SOURCE and the file name) are dropped from
# the tab goldens; the same two fields are dropped from the CSQ.
#
# What the world leaves out, and why, is in gen_world.py.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
N="${1:-500}"
SEED="${2:-139}"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FLAGS="--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$here/gen_world.py" "$tmp" "$N" "$SEED"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

vep() { # <output> <vep arguments...>
  local out="$1"; shift
  docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
    --input_file /data/in.vcf --output_file "/data/$out" --gff /data/genes.sorted.gff3.gz \
    --fasta /data/genome.fa --force_overwrite --no_stats "$@" >/dev/null 2>&1
}
vep default.tab --tab
# shellcheck disable=SC2086
vep flags.tab --tab $FLAGS
# shellcheck disable=SC2086
vep flags.vcf --vcf $FLAGS

cp "$tmp/genome.fa" "$tmp/genes.gff3" "$tmp/in.vcf" "$here/"

python3 - "$tmp" "$here" "$IMAGE" "$N" "$SEED" <<'PY'
import gzip, sys
tmp, here, image, n, seed = sys.argv[1:6]
note = f"image {image}; gen_world.py N={n} SEED={seed}"

def tab(src, dst):
    out, version = [], "?"
    for line in open(f"{tmp}/{src}"):
        line = line.rstrip("\n")
        if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
            version = line.split()[-1]
        elif line.startswith("##"):
            continue
        else:
            out.append(line)
    header = out[0][1:].split("\t")
    drop = {i for i, name in enumerate(header) if name == "SOURCE" or name.endswith(".gz")}
    with gzip.GzipFile(f"{here}/{dst}", "wb", mtime=0) as f:
        f.write(f"## the tab file of Ensembl VEP {version}; {note}\n".encode())
        for row in out:
            cells = [c for i, c in enumerate(row.split("\t")) if i not in drop]
            f.write(("\t".join(cells) + "\n").encode())
    return len(out) - 1, version

rows, version = tab("default.tab", "golden_default.tab.gz")
rows_f, _ = tab("flags.tab", "golden_flags.tab.gz")

fields, records = None, 0
with gzip.GzipFile(f"{here}/golden_flags.vcf.gz", "wb", mtime=0) as f:
    f.write(f"##the CSQ of Ensembl VEP {version}; {note}\n".encode())
    for line in open(f"{tmp}/flags.vcf"):
        line = line.rstrip("\n")
        if line.startswith("##INFO=<ID=CSQ"):
            fields = line.split("Format: ")[1].rstrip('">').split("|")
            keep = [i for i, name in enumerate(fields) if name != "SOURCE" and not name.endswith(".gz")]
            f.write(("##CSQ-fields: " + "|".join(fields[i] for i in keep) + "\n").encode())
        elif line.startswith("#"):
            continue
        else:
            c = line.split("\t")
            info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
            entries = []
            for entry in info.get("CSQ", "").split(","):
                v = entry.split("|")
                v += [""] * (len(fields) - len(v))
                entries.append("|".join(v[i] for i in keep))
            f.write(("\t".join(c[:5]) + "\t" + ",".join(entries) + "\n").encode())
            records += 1
print(f"VEP {version}: {rows} default rows, {rows_f} flag rows, {records} records")
PY
