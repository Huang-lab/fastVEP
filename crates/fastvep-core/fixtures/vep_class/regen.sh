#!/usr/bin/env bash
# Regenerates golden.tsv: the class Ensembl VEP's own `SO_variation_class` gives
# each allele string from gen_alleles.py. That function is what `--variant_class`
# reports (`VariationFeature::class_SO_term`).
#
#   ./regen.sh                       rewrite golden.tsv (the committed sample)
#   ./regen.sh big N OUT.tsv         N strings with VEP's answers, for the
#                                    env-gated test in tests/vep_variant_class_golden.rs
#
# Needs perl, python3, curl and shasum. The oracle is VEP's code run as it is:
# SO_variation_class is extracted by name from ensembl-variation and the
# SO_TERM_* constants it returns come from the same commit's Constants.pm, both
# pinned to the release/115 commit below, fetched here and checked against the
# SHA-256 recorded below. Nothing from them is kept in this repository, only the
# answers.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"

VARIATION_SHA=23c76f60b1592e4df86159cf5530bdc326120c3d
SEQUENCE_SHA256=15619ca7453899439b92a519bc9ed9cff2c5d5226e35ac157b5649591c1e8328
CONSTANTS_SHA256=98021460cfada22118c6d6b7865bbed3b25c501ca484ba2399380662f1012051

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

fetch() { # <path in the repo> <dest> <sha256>
  curl -fsS -m 60 -o "$2" "https://raw.githubusercontent.com/Ensembl/ensembl-variation/$VARIATION_SHA/$1"
  local got
  got="$(shasum -a 256 "$2" | cut -d' ' -f1)"
  if [ "$got" != "$3" ]; then
    echo "SHA-256 of $1 is $got, expected $3: the pinned file changed" >&2
    exit 1
  fi
}

fetch modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm "$tmp/sequence.pm" "$SEQUENCE_SHA256"
fetch modules/Bio/EnsEMBL/Variation/Utils/Constants.pm "$tmp/constants.pm" "$CONSTANTS_SHA256"

extract() { awk -v name="$2" '$0 ~ "^sub " name " *\\{" {p=1} p {print} p && /^}/ {p=0}' "$1"; }

{
  echo 'use strict; use warnings;'
  grep '^use constant SO_TERM_' "$tmp/constants.pm"
  extract "$tmp/sequence.pm" SO_variation_class
  cat <<'PERL'
while (my $line = <STDIN>) {
  chomp $line;
  print SO_variation_class($line), "\n";
}
PERL
} >"$tmp/oracle.pl"
perl -c "$tmp/oracle.pl" 2>/dev/null

run() { perl "$tmp/oracle.pl" <"$1" >"$2"; }

if [ "${1:-}" = "big" ]; then
  n="${2:?string count}"; out="${3:?output file}"
  python3 "$here/gen_alleles.py" --count "$n" >"$tmp/alleles.txt"
  run "$tmp/alleles.txt" "$tmp/classes.txt"
  paste "$tmp/alleles.txt" "$tmp/classes.txt" >"$out"
  echo "wrote $(wc -l <"$out") strings to $out"
  exit 0
fi

# The committed sample: every handwritten string, then up to 150 of each class
# VEP answers with, so no class can go untested.
python3 "$here/gen_alleles.py" >"$tmp/alleles.txt"
run "$tmp/alleles.txt" "$tmp/classes.txt"
paste "$tmp/alleles.txt" "$tmp/classes.txt" >"$tmp/all.tsv"
handwritten="$(python3 -c "
import sys; sys.dont_write_bytecode = True; sys.path.insert(0, '$here')
import gen_alleles; print(len(gen_alleles.handwritten()))")"

{
  echo "# Ensembl VEP's own answers, from running SO_variation_class as it is."
  echo "# ensembl-variation@$VARIATION_SHA (release/115) Utils/Sequence.pm sha256 $SEQUENCE_SHA256"
  echo "# and Utils/Constants.pm sha256 $CONSTANTS_SHA256 (the SO_TERM_* values)."
  echo "# Regenerate with regen.sh; do not edit."
  echo "# allele_string<TAB>class"
  awk -v n="$handwritten" 'NR <= n' "$tmp/all.tsv"
  awk -F'\t' -v n="$handwritten" 'NR > n && ++seen[$2] <= 150' "$tmp/all.tsv"
} >"$here/golden.tsv"
echo "wrote $(grep -vc '^#' "$here/golden.tsv") strings to $here/golden.tsv"
