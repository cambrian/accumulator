#!/bin/sh
set -e

# Sanity check the code.
cargo build
cargo test

# Clone docs branch in submodule.
git submodule update --init docs
cd docs
git checkout gh-pages
cd ..

# Backup manually-placed files.
mkdir -p docs_tmp
cp docs/404.html docs_tmp
cp docs/_config.yml docs_tmp
cp docs/README.md docs_tmp

# Generate docs.
rm -rf target/doc docs/*
cargo doc --no-deps
cp -r target/doc/* docs

# Restore manually-placed files.
cp docs_tmp/* docs
rm -rf docs_tmp

# Site is published at BASE/accumulator.
# Docs are at BASE/accumulator/docs.
mv docs/accumulator docs/docs

# Add Jekyll redirect for BASE/accumulator.
echo '0a
---
redirect_from: "/index.html"
---
.
w' | ed docs/docs/index.html

# Commit regenerated docs to submodule.
cd docs
git add --all
git commit -m "Regenerating docs via script."
cd ..
