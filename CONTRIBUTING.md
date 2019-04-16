# Contributing
Contributions are encouraged from anyone!

## Guidelines
- We adhere to the Rust [Code of Conduct](https://www.rust-lang.org/conduct.html).
- For bugs and feature requests, please
  [submit an issue](https://github.com/cambrian/accumulator/issues).
- For code contributions, take a look at the information below.

## Recommended Setup
1. Clone this repository.
2. Run `./setup.sh` (if you already have the Rust toolchain installed, take a look at the setup
   script and add only the tools you are missing).
3. Open VS Code and install the recommended extensions. **Note**: You are welcome to use your
   favorite IDE, but the original authors of this library have found VS Code extremely handy.
4. Restart VS Code.

## Development
- Leave `cargo watch -x "clippy -- -W clippy::pedantic"` running to type-check on save. Please
  adhere to the `pedantic` option in Clippy, as any contribution must pass those checks (without
  errors or warnings) to be included.
- Ensure that your code is formatted with `rustfmt`. If you use the recommended VS Code setup, this
  should happen whenever you save a file.
- Write tests! The repository has many examples of tests; run `cargo test` early and often.
- The command `cargo bench` uses [Criterion](https://crates.io/crates/criterion) benchmarks.
- When you are ready to submit your branch, create a pull request to `master`. A code owner will
  shepherd your PR through a review process prior to merge.

## Starting Points
- Comments in the code labeled `REVIEW` or `TODO` need attention of some sort.
- You are also welcome to triage and resolve
  [issues](https://github.com/cambrian/accumulator/issues), particularly those labeled `good first
  issue`.
- So little time, so much [code to review](https://github.com/cambrian/accumulator/pulls)...

## Notes
- This repository falls under the MIT License. Ensure that your contributions are MIT-compatible (or
  ask for help if you are not sure).
- To become a code owner or collaborator on this repository, please message one of the code owners.

## Troubleshooting
- If your Rust Language Server (RLS) hangs at `building` in VS Code, run
  `cargo clean && rustup update`.
- If you get unexpected build errors, delete `Cargo.lock`, run `cargo update`, and re-build.
