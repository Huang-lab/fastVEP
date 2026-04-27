#!/usr/bin/env bash
set -euo pipefail

# Keep cargo's registry/cache inside the build sandbox.
export CARGO_HOME="${SRC_DIR}/.cargo"

# Build and install both workspace binaries into $PREFIX/bin/.
cargo install --locked --no-track \
    --path crates/fastvep-cli \
    --root "${PREFIX}"

cargo install --locked --no-track \
    --path crates/fastvep-web \
    --root "${PREFIX}"
