# crypto
Crypto primitives in Rust.

## Setup
1. Clone this repo.
2. Run `./setup.sh`.
3. Open VS Code and install the recommended extensions.
4. Restart VS Code.

## Development
1. Leave `cargo watch` running to type-check on save.
2. Run `cargo build [--package PACKAGE] [--release]` to build package(s).
3. Run `cargo run [--package PACKAGE]` to run package executables.

## Troubleshooting
- If your RLS hangs at `building`, run `cargo clean && rustup update`.
- If you get unexpected build errors, delete `Cargo.lock`, run `cargo update`, and re-build.
