#!/usr/bin/env bash
# Regenerates golden.tab: the tab file Ensembl's own VEP writes for in.vcf and genes.gff3.
#
#   ./regen.sh
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3, the release of the cache in issue #139): `vep --tab` with no other
# flag, so what is compared is the output format and not any option. Needs docker,
# bgzip and tabix (htslib).
#
# VEP's --gff mode skips a variant on a contig no feature of the GFF3 sits on, and
# in.vcf has a chrM record that has no gene near it, the case a cache run gives an
# intergenic row for. So VEP is given a cache as well: one that holds a single
# dummy record per chromosome and so never reports a known variant, and exists only
# to make `10` and `MT` the contigs VEP knows. A random genome stands in for the FASTA
# --gff mode wants (VEP indexes it itself); none of the rows reads a base of it.
#
# Only the columns of the default tab file are kept (the first 17). The two VEP adds
# for --gff and --check_existing (SOURCE, the name of the GFF3, CLIN_SIG, SOMATIC,
# PHENO) are not the output format of the issue.
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cp "$here/in.vcf" "$here/genes.gff3" "$tmp/"

python3 - "$tmp" <<'PY'
import random, sys
tmp = sys.argv[1]
rng = random.Random(1)
with open(f"{tmp}/genome.fa", "w") as f:
    for name, length in (("chr10", 20000), ("chrM", 1000)):
        f.write(f">{name}\n")
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        for i in range(0, length, 60):
            f.write(seq[i:i + 60] + "\n")
PY

cache="$tmp/cache/homo_sapiens/116_GRCh38"
mkdir -p "$cache/10" "$cache/MT"
printf '%s\n' '# CACHE UPDATED 2026-04-08 00:00:00' 'species	homo_sapiens' 'assembly	GRCh38' \
  'cache_region_size	1000000' \
  'variation_cols	chr,variation_name,failed,somatic,start,end,allele_string,strand,minor_allele,minor_allele_freq,clin_sig,phenotype_or_disease,clin_sig_allele,pubmed,var_synonyms,AFR,AMR,EAS,EUR,SAS' \
  'var_type	tabix' > "$cache/info.txt"
printf '10\tchr10\nMT\tchrM\n' > "$cache/chr_synonyms.txt"
printf '10\tdummy1\t0\t0\t1\t1\tA/C\t1\t.\t.\t.\t0\t.\t.\t.\t.\t.\t.\t.\t.\n' | bgzip -c > "$cache/10/all_vars.gz"
printf 'MT\tdummy2\t0\t0\t1\t1\tA/C\t1\t.\t.\t.\t0\t.\t.\t.\t.\t.\t.\t.\t.\n' | bgzip -c > "$cache/MT/all_vars.gz"
tabix -C -s 1 -b 5 -e 5 "$cache/10/all_vars.gz"
tabix -C -s 1 -b 5 -e 5 "$cache/MT/all_vars.gz"
(grep '^#' "$tmp/genes.gff3"; grep -v '^#' "$tmp/genes.gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
  | bgzip -c > "$tmp/genes.sorted.gff3.gz"
tabix -p gff "$tmp/genes.sorted.gff3.gz"

docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
  --input_file /data/in.vcf --output_file /data/vep.tab --tab --force_overwrite --no_stats \
  --offline --cache --dir_cache /data/cache --species homo_sapiens --assembly GRCh38 --cache_version 116 \
  --gff /data/genes.sorted.gff3.gz --fasta /data/genome.fa --check_existing >/dev/null

python3 - "$tmp/vep.tab" "$here/golden.tab" "$IMAGE" <<'PY'
import sys
tab, out, image = sys.argv[1:4]
version = "?"
rows = []
for line in open(tab):
    line = line.rstrip("\n")
    if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
        version = line.split()[-1]
    elif not line.startswith("##"):
        rows.append("\t".join(line.split("\t")[:17]))
with open(out, "w") as f:
    f.write(f"## the default tab file of Ensembl VEP {version}, the first 17 columns\n")
    f.write(f"## image {image}\n")
    f.write("\n".join(rows) + "\n")
print(f"{len(rows) - 1} rows, VEP {version}")
PY
