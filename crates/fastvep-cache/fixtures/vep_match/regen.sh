#!/usr/bin/env bash
# Regenerates golden.tsv: what Ensembl VEP's own allele matching answers for
# every pair in gen_cases.py's output.
#
#   ./regen.sh                       rewrite golden.tsv (the committed sample)
#   ./regen.sh big N OUT.tsv         N-pair set with VEP's answers, for the
#                                    env-gated test in tests/vep_match_golden.rs
#
# Needs perl, python3, curl and shasum. The oracle is VEP's code, run as-is:
# trim_sequences, get_matched_variant_alleles and _get_trim_directions are
# extracted by name from ensembl-variation, reverse_comp from ensembl core, both
# pinned to the release/115 commits below. They are fetched here, checked
# against the SHA-256 recorded below, and not kept in this repository; only the
# answers are. Four helpers (throw, warning, assert_ref and the sub-name
# plumbing) are stubbed so the file runs without the Ensembl API installed.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"

VARIATION_SHA=23c76f60b1592e4df86159cf5530bdc326120c3d
CORE_SHA=266b84da60126688c6182b7739b9c43700432cfe
VARIATION_SHA256=15619ca7453899439b92a519bc9ed9cff2c5d5226e35ac157b5649591c1e8328
CORE_SHA256=91ae9a63bf72a649f927cd84877929bb6db901caeacac95b5bf083a7f3ce0bfd

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

fetch() { # <url> <dest> <sha256>
  curl -fsS -m 60 -o "$2" "$1"
  local got
  got="$(shasum -a 256 "$2" | cut -d' ' -f1)"
  if [ "$got" != "$3" ]; then
    echo "SHA-256 of $1 is $got, expected $3: the pinned file changed" >&2
    exit 1
  fi
}

fetch "https://raw.githubusercontent.com/Ensembl/ensembl-variation/$VARIATION_SHA/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm" \
  "$tmp/variation.pm" "$VARIATION_SHA256"
fetch "https://raw.githubusercontent.com/Ensembl/ensembl/$CORE_SHA/modules/Bio/EnsEMBL/Utils/Sequence.pm" \
  "$tmp/core.pm" "$CORE_SHA256"

# a sub, from `sub NAME {` to the first line that is just `}`
extract() { awk -v name="$2" '$0 ~ "^sub " name " *\\{" {p=1} p {print} p && /^}/ {p=0}' "$1"; }

{
  cat <<'PERL'
use strict;
use warnings;
sub throw   { die join(' ', @_) . "\n" }
sub warning { warn join(' ', @_) . "\n" }
sub assert_ref { my ($r, $t) = @_; die "expected $t ref\n" unless ref($r) eq $t; }
PERL
  extract "$tmp/core.pm" reverse_comp
  extract "$tmp/variation.pm" trim_sequences
  extract "$tmp/variation.pm" _get_trim_directions
  extract "$tmp/variation.pm" get_matched_variant_alleles
  cat <<'PERL'
while (my $line = <STDIN>) {
  chomp $line;
  my ($as, $ap, $astr, $bs, $bp, $bstr) = split /\t/, $line;
  my $m = get_matched_variant_alleles(
    { allele_string => $as, pos => $ap, strand => $astr },
    { allele_string => $bs, pos => $bp, strand => $bstr },
  );
  print join(';', map { "$_->{a_allele},$_->{b_allele}" } @$m) || '-', "\n";
}
PERL
} >"$tmp/oracle.pl"
perl -c "$tmp/oracle.pl" 2>/dev/null

# run <cases file> <expected file>
run() { perl "$tmp/oracle.pl" <"$1" >"$2"; }

if [ "${1:-}" = "big" ]; then
  n="${2:?pair count}"; out="${3:?output file}"
  share=$((n / 3))
  python3 "$here/gen_cases.py" --structured $((n - 2 * share)) --multiallelic "$share" --random "$share" >"$tmp/cases.tsv"
  run "$tmp/cases.tsv" "$tmp/expected.txt"
  paste "$tmp/cases.tsv" "$tmp/expected.txt" >"$out"
  echo "wrote $(wc -l <"$out") pairs to $out"
  exit 0
fi

# The committed sample: every handwritten pair, then the pairs VEP matches and
# the pairs it does not, in a fixed ratio so a matcher that always says "no"
# (or "yes") cannot pass.
python3 "$here/gen_cases.py" >"$tmp/cases.tsv"
run "$tmp/cases.tsv" "$tmp/expected.txt"
paste "$tmp/cases.tsv" "$tmp/expected.txt" >"$tmp/all.tsv"

handwritten="$(python3 -c "
import sys; sys.path.insert(0, '$here')
import gen_cases; print(len(gen_cases.handwritten()))")"

{
  echo "# Ensembl VEP's own answers, from running get_matched_variant_alleles as-is."
  echo "# ensembl-variation@$VARIATION_SHA (release/115) Utils/Sequence.pm sha256 $VARIATION_SHA256"
  echo "# ensembl@$CORE_SHA (release/115) Utils/Sequence.pm sha256 $CORE_SHA256"
  echo "# Regenerate with regen.sh; do not edit."
  echo "# a_allele_string<TAB>a_pos<TAB>a_strand<TAB>b_allele_string<TAB>b_pos<TAB>b_strand<TAB>matches"
  echo "# matches: '-' for none, else a_allele,b_allele pairs joined by ';'"
  # awk with a counter rather than `head`: head closing the pipe is a SIGPIPE,
  # which `pipefail` turns into a failed script and a truncated golden.tsv
  awk -v n="$handwritten" 'NR <= n' "$tmp/all.tsv"
  awk -F'\t' -v n="$handwritten" 'NR > n && $7 != "-" && ++kept <= 700' "$tmp/all.tsv"
  awk -F'\t' -v n="$handwritten" 'NR > n && $7 == "-" && ++kept <= 1300' "$tmp/all.tsv"
} >"$here/golden.tsv"
echo "wrote $(grep -vc '^#' "$here/golden.tsv") pairs to $here/golden.tsv"
