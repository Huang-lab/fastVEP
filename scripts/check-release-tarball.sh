#!/usr/bin/env bash
#
# Run a release archive the way a user would, before a user does.
#
#   scripts/check-release-tarball.sh fastVEP-v0.5.0-aarch64-unknown-linux-gnu.tar.gz [v0.5.0]
#
# Unpacks the archive into a scratch directory, starts both binaries, and
# annotates the README's example. The release workflow runs it on every archive
# it builds. A Linux archive is run inside Rocky Linux 9 (glibc 2.34), the
# oldest system the release claims to run on, because a binary that links a
# newer glibc builds and packs without a word and then fails at the first
# `fastvep --version` on the user's machine.
#
# Given the version (the tag, on a release), it also requires `fastvep` to
# report it. A tag pushed before Cargo.toml is bumped would otherwise publish an
# archive named v0.5.0 whose binary says 0.4.0.
#
# Plain bash 3.2 and coreutils only: it runs on a stock macOS and in a minimal
# container.
set -euo pipefail

archive=${1:?usage: check-release-tarball.sh ARCHIVE [VERSION]}
want_version=${2:-}
repo=$(cd "$(dirname "$0")/.." && pwd)
archive=$(cd "$(dirname "$archive")" && pwd)/$(basename "$archive")

fail() {
  echo "check-release-tarball: $*" >&2
  exit 1
}

# The checksum is published beside the archive, so a file that does not verify
# is a defect of the release as much as a bad binary is.
if [ -f "$archive.sha256" ]; then
  (
    cd "$(dirname "$archive")"
    if command -v sha256sum >/dev/null; then
      sha256sum -c "$(basename "$archive").sha256"
    else
      shasum -a 256 -c "$(basename "$archive").sha256"
    fi
  ) || fail "$(basename "$archive") does not match its .sha256 file"
fi

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

tar -xzf "$archive" -C "$scratch"
unpacked=("$scratch"/fastVEP-*/)
if [ "${#unpacked[@]}" -ne 1 ] || [ ! -d "${unpacked[0]}" ]; then
  fail "$(basename "$archive") does not unpack to a single fastVEP-* directory"
fi
dir=${unpacked[0]%/}

for file in fastvep fastvep-web README.md LICENSE.md CHANGELOG.md; do
  [ -e "$dir/$file" ] || fail "$(basename "$archive") has no $file"
done

version=$("$dir/fastvep" --version)
echo "$version"
if [ -n "$want_version" ] && [ "$version" != "fastvep ${want_version#v}" ]; then
  fail "the tag is $want_version but the binary reports '$version': Cargo.toml was not bumped before the tag"
fi

# fastvep-web has no --version; --help is enough to prove it loads and starts.
"$dir/fastvep-web" --help >/dev/null || fail "fastvep-web does not start"

# The README's try-it example. The inputs are copied first because annotate
# writes its transcript cache next to the GFF3, and the checkout may be mounted
# read-only.
cp "$repo/tests/test.vcf" "$repo/tests/test.gff3" "$scratch/"
(
  cd "$scratch"
  "$dir/fastvep" annotate -i test.vcf --gff3 test.gff3 --hgvs --output-format tab \
    >annotated.tab 2>annotate.log
) || {
  cat "$scratch/annotate.log" >&2
  fail "fastvep annotate failed on the README example"
}

# Two rows that show the annotation went through the GFF3 parser, the exon
# geometry and the strand. Not the whole output, and nothing a release changes
# on purpose (v0.4.0 and master write different HGVSc for these variants):
# this checks that the binary works on its platform, and the test suite checks
# the rest. The rows go through a variable because `cut | grep -q` can die of
# SIGPIPE under pipefail.
rows=$(cut -f1,7 "$scratch/annotated.tab")
for want in $'rs_5utr_brca1\t5_prime_UTR_variant' $'rs_intron_brca1\tintron_variant'; do
  grep -qxF "$want" <<<"$rows" || fail "no '${want/$'\t'/ }' row in the annotation of the README example"
done

echo "ok: $(basename "$archive") starts, reports '$version', and annotates the README example"
