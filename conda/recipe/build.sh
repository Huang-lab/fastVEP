#!/usr/bin/env bash
set -xeuo pipefail

# Keep cargo's registry/cache inside the build sandbox.
export CARGO_HOME="${BUILD_PREFIX}/.cargo"

# Bioconda policy for Rust packages: vendor the dependency licences so they can
# be shipped next to ours. `about.license_file` in meta.yaml picks this up.
cargo-bundle-licenses --format yaml --output THIRDPARTY.yml

# Build and install both workspace binaries into $PREFIX/bin/.
cargo install --locked --no-track \
    --path crates/fastvep-cli \
    --root "${PREFIX}"

cargo install --locked --no-track \
    --path crates/fastvep-web \
    --root "${PREFIX}"
