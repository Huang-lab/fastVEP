#!/usr/bin/env bash
#
# Cut a fastVEP release.
#
#   scripts/release.sh            # auto-compute next version from conventional commits
#   scripts/release.sh 0.4.0      # force an explicit version
#
# What it does (nothing is pushed; you review, then tag):
#   1. Determine the target version (git-cliff --bumped-version, or the arg).
#   2. Draft the changelog entry for the new release from commits since the
#      last tag and print it — you paste/edit the curated prose into CHANGELOG.md.
#   3. Bump the workspace version in Cargo.toml and refresh Cargo.lock.
#
# After running: edit CHANGELOG.md, then:
#   git commit -am "release: v<VERSION>" && git tag -a v<VERSION> -m "v<VERSION>"
#   git push && git push --tags
#
set -euo pipefail
cd "$(dirname "$0")/.."

command -v git-cliff >/dev/null || { echo "git-cliff not found (brew install git-cliff)"; exit 1; }

if [[ $# -ge 1 ]]; then
  VERSION="${1#v}"
else
  # git-cliff picks major/minor/patch from the conventional commits since the last tag.
  VERSION="$(git-cliff --bumped-version 2>/dev/null | sed 's/^v//')"
fi
[[ -n "${VERSION}" ]] || { echo "could not determine version"; exit 1; }

CURRENT="$(grep -m1 '^version' Cargo.toml | sed -E 's/.*"([^"]+)".*/\1/')"
echo "Current: v${CURRENT}   ->   Target: v${VERSION}"
echo

echo "=== Draft changelog for v${VERSION} (edit the prose in CHANGELOG.md) ==="
git-cliff --unreleased --tag "v${VERSION}" 2>/dev/null | tail -n +5
echo "=========================================================================="
echo

# Bump the single workspace version; all crates inherit version.workspace = true.
sed -i.bak -E "0,/^version = \"[^\"]+\"/s//version = \"${VERSION}\"/" Cargo.toml && rm -f Cargo.toml.bak
cargo update --workspace >/dev/null 2>&1 || true

echo "Bumped Cargo.toml -> ${VERSION} and refreshed Cargo.lock."
echo "Next: finalize CHANGELOG.md, then:"
echo "  git commit -am \"release: v${VERSION}\" && git tag -a \"v${VERSION}\" -m \"v${VERSION}\""
echo "  git push && git push --tags"
