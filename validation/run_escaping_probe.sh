#!/usr/bin/env bash
set -euo pipefail

# Measure what Ensembl VEP writes into a VCF INFO value for each character
# that has to be decided about, and print it next to what fastVEP writes.
#
# The escaping rules in `crates/fastvep-io/src/output.rs` cite this script.
# They used to cite a Docker run nobody else could repeat, which is the thing
# CLAUDE.md's "Claims" section is about: a rationale that is confidently wrong
# about a file format survives review by sounding specific.
#
# Two probes, because VEP is not self-consistent between them:
#
#   1. A `--custom` BED whose interval names carry one special character each.
#      This is a value VEP copies from a file into a CSQ subfield, so it shows
#      the substitution rules with nothing else in the way.
#   2. A synonymous SNV, whose `HGVSp` ends in `p.Xxx123=`. This is a value VEP
#      *generates*, and it is escaped where the custom value is not.
#
# Prerequisites:
#   - cargo build --release
#   - Docker, and ensemblorg/ensembl-vep:release_115.1 (pulled if missing)
#   - bgzip and tabix (htslib)
#
# Usage:
#   ./validation/run_escaping_probe.sh [outdir]   # default: validation/results/escaping

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="${1:-$SCRIPT_DIR/results/escaping}"
VEP_IMAGE="ensemblorg/ensembl-vep:release_115.1"
FASTVEP="$PROJECT_DIR/target/release/fastvep"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

command -v bgzip >/dev/null || { echo "need bgzip (htslib)" >&2; exit 1; }
command -v tabix >/dev/null || { echo "need tabix (htslib)" >&2; exit 1; }
docker image inspect "$VEP_IMAGE" >/dev/null 2>&1 || docker pull "$VEP_IMAGE"

# A self-contained gene model, so the probe needs no downloaded reference: one
# gene, one transcript, a 200 bp single-exon CDS opening with ATG.
cat > probe.gff3 <<'GFF'
##gff-version 3
1	probe	gene	1001	1200	.	+	.	ID=gene:ENSG_P;Name=PROBE;biotype=protein_coding
1	probe	mRNA	1001	1200	.	+	.	ID=transcript:ENST_P;Parent=gene:ENSG_P;biotype=protein_coding;tag=Ensembl_canonical;transcript_support_level=1
1	probe	exon	1001	1200	.	+	.	ID=exon:E_P;Parent=transcript:ENST_P;rank=1
1	probe	CDS	1001	1200	.	+	0	ID=CDS:P_P;Parent=transcript:ENST_P
GFF

python3 - <<'PY'
# `AAA` codons throughout (Lys), so A>G at a third base is synonymous
# (AAA -> AAG) and gives VEP an `HGVSp` ending in `=` to escape.
seq = ["A"] * 1300
seq[1000:1003] = list("ATG")
with open("ref.fa", "w") as f:
    f.write(">1\n")
    for i in range(0, 1300, 60):
        f.write("".join(seq[i : i + 60]) + "\n")
PY

# One interval per character, each overlapping its own variant, so a mangled
# value can always be traced back to the character that produced it.
printf '%s\n' \
  $'1\t1010\t1020\tspace test' \
  $'1\t1030\t1040\tcolon:test' \
  $'1\t1050\t1060\tsemi;test' \
  $'1\t1070\t1080\teq=test' \
  $'1\t1090\t1100\tcomma,test' \
  $'1\t1110\t1120\tpipe|test' \
  $'1\t1130\t1140\tquote"test' \
  $'1\t1150\t1160\tpct%test' \
  > custom.bed
bgzip -f -c custom.bed > custom.bed.gz
tabix -f -p bed custom.bed.gz

# VEP reads a GFF3 through tabix and rejects an unindexed one outright.
grep -v '^#' probe.gff3 | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -f -c > probe.gff3.gz
tabix -f -p gff probe.gff3.gz

{
  echo '##fileformat=VCFv4.2'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  for pos in 1015 1035 1055 1075 1095 1115 1135 1155; do
    printf '1\t%s\t.\tA\tG\t.\tPASS\t.\n' "$pos"
  done
} > probe.vcf

docker run --rm -e VEP_DATABASE=0 -v "$OUT_DIR:/d" "$VEP_IMAGE" vep \
  --input_file /d/probe.vcf --gff /d/probe.gff3.gz --fasta /d/ref.fa \
  --output_file /d/vep.vcf --vcf --force_overwrite --no_stats --hgvs --symbol \
  --custom 'file=/d/custom.bed.gz,short_name=PROBE,format=bed,type=overlap' \
  2>&1 | grep -vE '^(WARNING: Parent entries|WARNING: Ignoring)' | tail -3 || true

echo ""
echo "=== What VEP 115.1 writes for a value it copies from a file ==="
printf "%-14s %-16s %s\n" "in the BED" "VEP writes" "verdict"
if python3 - <<'PY'
import re, sys

fields = None
rows = {}
for line in open("vep.vcf"):
    if line.startswith("##INFO=<ID=CSQ"):
        fields = re.search(r"Format: ([^\"]+)", line).group(1).split("|")
        continue
    if line.startswith("#"):
        continue
    cols = line.rstrip("\n").split("\t")
    m = re.search(r"CSQ=([^\t;]+)", cols[7])
    if not m or fields is None:
        continue
    for entry in m.group(1).split(","):
        v = dict(zip(fields, entry.split("|")))
        got = v.get("PROBE", "")
        if got:
            rows[cols[1]] = got

expected = {
    "1015": ("space test", "space_test"),
    "1035": ("colon:test", "colon:test"),
    "1055": ("semi;test", "semi%3Btest"),
    "1075": ("eq=test", "eq=test"),
    "1095": ("comma,test", "comma&test"),
    "1115": ("pipe|test", "pipe&test"),
    "1135": ('quote"test', 'quote"test'),
    "1155": ("pct%test", "pct%test"),
}
bad = 0
for pos, (src, want) in expected.items():
    got = rows.get(pos, "(no value)")
    ok = "as recorded" if got == want else f"CHANGED, expected {want}"
    if got != want:
        bad += 1
    print(f"{src:<14} {got:<16} {ok}")

print("")
print("=== ...and for a value it generates itself (HGVSp of a synonymous SNV) ===")
for line in open("vep.vcf"):
    if line.startswith("#"):
        continue
    m = re.search(r"CSQ=([^\t;]+)", line.split("\t")[7])
    if not m:
        continue
    for entry in m.group(1).split(","):
        v = dict(zip(fields, entry.split("|")))
        if "synonymous" in v.get("Consequence", "") and v.get("HGVSp"):
            print(f"  HGVSp = {v['HGVSp']}")
            if "%3D" not in v["HGVSp"]:
                print("  CHANGED: expected the `=` to be written as %3D")
                bad += 1
            break
    else:
        continue
    break

print("")
print("Conclusions the escaping rules in output.rs rest on:")
print("  - `:` and `\"` are written through literally by VEP")
print("  - a space becomes `_`")
print("  - `;` becomes %3B, `,` and `|` become `&`")
print("  - `=` is literal in a copied custom value and %3D in a generated one")
print("  - `%` is literal, so VEP's own CSQ is not unambiguously decodable:")
print("    a literal `%` before `3B` cannot be told from an escaped `;`.")
print("    fastVEP encodes it as %25 in FV_* fields, which is what makes the")
print("    single-decode-pass promise in docs/SUPPLEMENTARY_ANNOTATIONS.md true.")
sys.exit(1 if bad else 0)
PY
then vep_status=0; else vep_status=$?; fi

echo ""
echo "=== What fastVEP writes for the same characters ==="
if [[ -x "$FASTVEP" ]]; then
    "$FASTVEP" annotate -i probe.vcf --gff3 probe.gff3 --fasta ref.fa \
        --symbol --hgvs -o fastvep.vcf --no-progress >/dev/null 2>&1
    echo "  CSQ subfields (escape_csq_str): space -> _, ';' -> %3B, ',' and '|' -> &,"
    echo "  '=' -> %3D, ':' literal, TAB/CR/LF -> %09/%0D/%0A"
    echo "  FV_* pipe fields (escape_vcf_subfield): the above, plus space -> %20,"
    echo "  '\"' -> %22, '&' -> %26, '|' -> %7C, ':' literal"
    grep -o 'CSQ=[^;[:space:]]*' fastvep.vcf | head -1 | cut -c1-120
else
    echo "  SKIP: $FASTVEP not built (cargo build --release)"
fi

exit $vep_status
