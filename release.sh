#!/bin/sh
set -e
cargo build
cargo test
rm -rf target/doc
cargo doc --no-deps
cp -r target/doc docs

# Site is published at BASE/accumulator.
# Docs are at BASE/accumulator/docs.
mv docs/accumulator docs/docs

# Add Jekyll redirect.
echo '0a
---
redirect_from: "/index.html"
---
.
w' | ed docs/docs/index.html
