#!/usr/bin/env bash
# Regenerates the goldens of this directory: what Ensembl's own VEP writes for each damaged GFF3
# next to its VCF (golden_<case>.tab.gz).
#
# The inputs are not generated: each is a small RefSeq-style gene model (no `rank` on its exons)
# with one record damaged, kept as it was found by fuzzing the GFF3 reader of fastVEP:
#
#   phase_dot      the first CDS record of NM_090001.3 has the phase `.`
#   phase_minus1   the same, `-1`
#   overlap        the two exons of NM_090002.2 (reverse strand) overlap: 3-6899 and 6000-6299
#   base1          a reverse-strand CDS of one base, at base 1 of its contig
#
# A CDS phase of `x` or `3` has no golden: VEP stops on it (`Bad value (x) for exon phase. Should
# only be -1,0,1,2`, exit 2) and writes no row, so the test only holds fastVEP to saying so.
#
# The oracle is VEP itself, the program, run in Docker from the image pinned below (release_116.3)
# in --gff mode. Needs docker, bgzip and tabix (htslib).
set -euo pipefail
export LC_ALL=C

here="$(cd "$(dirname "$0")" && pwd)"
IMAGE="ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
FLAGS="--symbol --biotype --canonical --numbers --hgvs"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# case, gff3, vcf, fasta
for spec in "phase_dot phase_dot.gff3 phase.vcf genome.fa" \
            "phase_minus1 phase_minus1.gff3 phase.vcf genome.fa" \
            "overlap overlap.gff3 overlap.vcf genome.fa" \
            "base1 base1.gff3 base1.vcf base1.fa"; do
  # shellcheck disable=SC2086
  set -- $spec
  name="$1"; gff3="$2"; vcf="$3"; fasta="$4"
  work="$tmp/$name"
  mkdir -p "$work"
  cp "$here/$vcf" "$work/in.vcf"
  cp "$here/$fasta" "$work/genome.fa"
  (grep '^#' "$here/$gff3" || true; grep -v '^#' "$here/$gff3" | sort -k1,1 -k4,4n -k5,5n -t$'\t') \
    | bgzip -c > "$work/genes.sorted.gff3.gz"
  tabix -p gff "$work/genes.sorted.gff3.gz"
  # shellcheck disable=SC2086
  docker run --rm -e VEP_DATABASE=0 -v "$work:/data" "$IMAGE" vep \
    --input_file /data/in.vcf --output_file /data/flags.tab --gff /data/genes.sorted.gff3.gz \
    --fasta /data/genome.fa --force_overwrite --no_stats --tab $FLAGS >/dev/null 2>&1
  python3 - "$work" "$here" "$IMAGE" "$name" <<'PY'
import gzip, sys
work, here, image, name = sys.argv[1:5]
out, version = [], "?"
for line in open(f"{work}/flags.tab"):
    line = line.rstrip("\n")
    if line.startswith("## ENSEMBL VARIANT EFFECT PREDICTOR"):
        version = line.split()[-1]
    elif line.startswith("##"):
        continue
    else:
        out.append(line)
header = out[0][1:].split("\t")
drop = {i for i, n in enumerate(header) if n == "SOURCE" or n.endswith(".gz")}
with gzip.GzipFile(f"{here}/golden_{name}.tab.gz", "wb", mtime=0) as f:
    f.write(f"## the tab file of Ensembl VEP {version}; image {image}; case {name}\n".encode())
    for row in out:
        cells = [c for i, c in enumerate(row.split("\t")) if i not in drop]
        f.write(("\t".join(cells) + "\n").encode())
print(f"{name}: VEP {version}: {len(out) - 1} rows")
PY
done
