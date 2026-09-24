#!/usr/bin/env bash
set -euo pipefail

# This file is the source of truth for recipes/fastvep/build.sh in
# bioconda-recipes. The two copies have to stay identical: the bioconda one is
# what CI proves, and a fix applied to only one of them is a fix that reaches
# nobody. Sync both when you touch either.

# Keep cargo caches inside the build sandbox.
export CARGO_HOME="${SRC_DIR}/.cargo"

# The workspace release profile sets `strip = "symbols"`. conda-build does its
# own Mach-O post-processing (llvm-otool / install_name_tool) before stripping
# binaries itself, and llvm-otool aborts with SIGABRT on the stripped
# fastvep-web binary on osx-64. Hand conda-build an unstripped binary instead;
# it strips the result during packaging, so the shipped package is unaffected.
export CARGO_PROFILE_RELEASE_STRIP=none

# Generate a bundled third-party license manifest (bioconda convention for Rust).
cargo-bundle-licenses --format yaml --output THIRDPARTY.yml

# Build and install both workspace binaries into $PREFIX/bin/.
cargo install --no-track --locked --verbose \
    --path crates/fastvep-cli \
    --root "${PREFIX}"

cargo install --no-track --locked --verbose \
    --path crates/fastvep-web \
    --root "${PREFIX}"
