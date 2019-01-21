# NOTE: if your rls hangs at "building", quit VS Code and run:
#   cargo clean && rustup update
# After this, reopen VS Code and you should be good.

# Set up the Rust toolchain and modify PATH.
curl https://sh.rustup.rs -sSf | sh
source $HOME/.cargo/env
rustup update

# Install useful extensions.
cargo install cargo-watch
cargo install cargo-edit
rustup component add rls rust-analysis rust-src clippy
