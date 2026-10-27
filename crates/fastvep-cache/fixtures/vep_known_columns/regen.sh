#!/usr/bin/env bash
# Regenerates this directory: the VCF, the VEP-format variation cache with the real
# release-116 column set, a two-transcript GFF3, and one golden file per set of
# flags: what Ensembl's own VEP answers for them in the known-variant columns.
#
#   ./regen.sh [N [SEED]]          N generated pairs (default 400), default seed 139
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below
# (release_116.3, the release the cache of issue #139 is for): `vep --tab --cache
# --offline` over the cache this script builds, with the GFF3 standing in for
# transcripts. Needs docker, python3, bgzip and tabix (htslib). The genome VEP
# insists on is generated, used and thrown away; only the inputs fastVEP needs and
# the answers are kept.
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
N="${1:-400}"
SEED="${2:-139}"
# release_116.3, by digest, so a moved tag cannot change the answers unnoticed
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"

# name  VEP flags: the sets differ in what they read (`--af_1kg` does not read a known
# variant that matched another allele of the site; `--af` does; `--max_af` reads every
# population whether or not it is written)
SETS=(
  "reporter|--af_1kg --check_existing"
  "af|--af"
  "gnomad|--af_gnomade --af_gnomadg"
  "max_af|--max_af"
  "pubmed|--pubmed"
  "existing|--check_existing"
  "csa0|--check_existing --pubmed --clin_sig_allele 0"
  "all|--af --af_1kg --af_gnomade --af_gnomadg --max_af --pubmed --check_existing"
)

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

for set in "${SETS[@]}"; do
  name="${set%%|*}"
  flags="${set#*|}"
  # VEP_DATABASE=0 rather than --offline alone: see validation/run_divergence.sh
  # (flags are split on purpose: they are a fixed list of plain words)
  # shellcheck disable=SC2086
  docker run --rm -e VEP_DATABASE=0 -v "$tmp:/data" "$IMAGE" vep \
    --input_file /data/in.vcf --output_file "/data/vep_$name.tab" --tab --force_overwrite --no_stats \
    --offline --cache --dir_cache /data/cache --species homo_sapiens --assembly GRCh38 --cache_version 116 \
    --gff /data/genes.sorted.gff3.gz --fasta /data/genome.fa $flags >/dev/null
  echo "ran VEP: $name ($flags)" >&2
done

rm -rf "$here/cache"
cp -R "$tmp/cache" "$here/cache"
cp "$tmp/in.vcf" "$tmp/genes.gff3" "$here/"
rm -f "$here"/golden_*.tsv

python3 - "$tmp" "$here" "$IMAGE" "$N" "$SEED" "${SETS[@]}" <<'PY'
import sys
tmp, here, image, n, seed, *sets = sys.argv[1:]
# the columns this directory is about: everything else VEP writes is the business of other goldens
KNOWN = ("Existing_variation AF AFR_AF AMR_AF EAS_AF EUR_AF SAS_AF "
         "gnomADe_AF gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF "
         "gnomADe_MID_AF gnomADe_NFE_AF gnomADe_REMAINING_AF gnomADe_SAS_AF "
         "gnomADg_AF gnomADg_AFR_AF gnomADg_AMI_AF gnomADg_AMR_AF gnomADg_ASJ_AF gnomADg_EAS_AF "
         "gnomADg_FIN_AF gnomADg_MID_AF gnomADg_NFE_AF gnomADg_REMAINING_AF gnomADg_SAS_AF "
         "MAX_AF MAX_AF_POPS CLIN_SIG SOMATIC PHENO PUBMED").split()
for entry in sets:
    name, flags = entry.split("|", 1)
    header, rows, version = None, [], "?"
    for line in open(f"{tmp}/vep_{name}.tab"):
        line = line.rstrip("\n")
        if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
            version = line.split()[-1]
        elif line.startswith("##"):
            continue
        elif line.startswith("#"):
            header = line[1:].split("\t")
        else:
            rows.append(dict(zip(header, line.split("\t"))))
    cols = [c for c in header if c in KNOWN]
    with open(f"{here}/golden_{name}.tsv", "w") as f:
        f.write(f"# the known-variant columns of every row, as Ensembl VEP {version} prints them\n")
        f.write(f"# image {image}\n# gen_world.py N={n} SEED={seed}; vep --tab --cache --offline {flags}\n")
        f.write("# columns, in VEP's order; each row is Uploaded_variation, Allele, Feature and then the\n")
        f.write("# cells that are not `-`, as COLUMN=value, tab separated\n")
        f.write("#columns\t" + "\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["Uploaded_variation"], r["Allele"], r["Feature"])):
            cells = [f"{c}={r[c]}" for c in cols if r[c] != "-"]
            f.write("\t".join([r["Uploaded_variation"], r["Allele"], r["Feature"], *cells]) + "\n")
    print(f"golden_{name}.tsv: {len(rows)} rows, VEP {version}")
PY
