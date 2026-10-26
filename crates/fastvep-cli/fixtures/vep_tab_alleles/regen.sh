#!/usr/bin/env bash
# Regenerates this directory: records whose ALT or REF is not a run of A, C, G and T (in.vcf, from
# gen_alleles.py), and what Ensembl's own VEP writes for them with every optional column fastVEP
# has: golden_flags.tab.gz (the tab file) and golden_flags.vcf.gz (the CSQ of each record).
#
#   ./regen.sh
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below (release_116.3,
# the release of the cache in issue #139), in --gff mode over ../vep_tab_world's genes and genome.
# Needs docker, python3, bgzip and tabix (htslib). The two columns VEP adds for --gff (SOURCE and
# the file name) are dropped from the goldens; the same two fields are dropped from the CSQ.
#
# non_ref.vcf is written by hand and is not run through VEP: VEP skips a record with a
# <NON_REF> ALT whole, and writes a `coding_sequence_variant` row for <*> that says nothing.
# What fastVEP does with them instead is in docs/VEP_DIVERGENCE.md and in the test.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
world="$here/../vep_tab_world"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FLAGS="--symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

python3 "$here/gen_alleles.py" "$tmp/in.vcf"
cp "$world/genome.fa" "$world/genes.gff3" "$tmp/"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

vep() { # <output> <vep arguments...>
  local out="$1"; shift
  docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
    --input_file /data/in.vcf --output_file "/data/$out" --gff /data/genes.sorted.gff3.gz \
    --fasta /data/genome.fa --force_overwrite --no_stats "$@" >/dev/null 2>&1
}
# shellcheck disable=SC2086
vep flags.tab --tab $FLAGS
# shellcheck disable=SC2086
vep flags.vcf --vcf $FLAGS

cp "$tmp/in.vcf" "$here/in.vcf"

python3 - "$tmp" "$here" "$IMAGE" <<'PY'
import gzip, sys
tmp, here, image = sys.argv[1:4]
note = f"image {image}; gen_alleles.py"

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
    f.write(f"## the tab file of Ensembl VEP {version}; {note}\n".encode())
    for row in out:
        cells = [c for i, c in enumerate(row.split("\t")) if i not in drop]
        f.write(("\t".join(cells) + "\n").encode())

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
                if not entry:
                    continue
                v = entry.split("|")
                v += [""] * (len(fields) - len(v))
                entries.append("|".join(v[i] for i in keep))
            f.write(("\t".join(c[:5]) + "\t" + ",".join(entries) + "\n").encode())
            records += 1
print(f"VEP {version}: {len(out) - 1} tab rows, {records} records in the VCF")
PY
