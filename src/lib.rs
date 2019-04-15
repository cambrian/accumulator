//! Cryptographic accumulator and vector commitment library, written and maintained by Cambrian
//! Labs.
//!
//! Most users will simply need to write
//!
//! `use accumulator::{Accumulator, VectorCommitment};`
//!
//! or similar to use the library. However, internal modules are exported as well in case you would
//! like to use them.
//!
//! We have spent significant time optimizing performance. Most accumulator or vector-commitment
//! functions will bottleneck in hashing to large primes. To alleviate this, we created a
//! zero-allocation U256 type that uses the low-level `mpn_` functions in GMP, which is used in
//! `hash_to_prime`. TODO: Benchmark our U256 vs. 256-bit `rug::Integer` vs. Parity U256.
#![allow(clippy::unknown_clippy_lints)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::empty_enum)]

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate arrayref;

mod accumulator;
pub use crate::accumulator::*;
mod vector_commitment;
pub use vector_commitment::*;

pub mod group;
pub mod hash;
pub mod proof;
pub mod uint;
pub mod util;
