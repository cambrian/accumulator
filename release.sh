#!/bin/sh
set -e
cargo build
cargo test
rm -rf target/doc
cargo doc --no-deps
cp -r target/doc docs
