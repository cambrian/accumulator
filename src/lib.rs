//! Fast cryptographic accumulator and vector commitment library, written and maintained by Cambrian
//! Labs.
//!
//! This library uses Rug Integers.
//!
//! # Important Usage Points for the Accumulator
//!
//! In order for accumulators to work and be secure in your application, you need to take care
//! of two invariants:
//!
//! 1. The cryptographic accumulator implemented here is meant to accumulate **odd primes only**.
//! To help you enforce this invariant we provide the function `hash_to_prime` in the `hash` module.
//!
//! 2. Cryptographic accumulators as implemented here assume that **no element is accumulated
//! twice**.
//!
//! # What is an accumulator?
//!
//! An accumulator is a cryptographic primitive that “accumulates”
//! values and lets you provide efficiently-verifiable proofs to other parties that
//! certain values have or have not been accumulated without needing to provide
//! the entire accumulated set.
//!
//! An accumulator is similar to a set but it is _constant-sized_,
//! meaning it can accumulate O(n) elements while remaining an O(1) data structure. For the
//! mathematics behind accumulators, see _Batching Techniques for Accumulators with
//! Applications to IOPs and Stateless Blockchains_ by Boneh, Bünz, and Fisch.
//!
//! # Examples
//!
//! Correct Usage:
//! ```
//! let acc = Accumulator::<ClassGroup>::new();
//!
//! let prime1 = Integer::from(7);
//! let prime2 = Integer::from(11);
//!
//! // Accumulate prime1 and prime2. The add method returns the new accumulator state and a
//! // proof you accumulated prime1 and prime2.
//! let (acc, proof) = acc.add(&[prime1.clone(), prime2.clone()]);
//!
//! // A network participant who sees (acc, proof, and [prime1, prime2]) can verify
//! // the update was made correctly...
//! assert!(acc.verify_membership(&[prime1, prime2], &proof));
//!
//! // ... and trying to verify something that has not been accumulated will not work.
//! assert!(!acc.verify_membership(&[Integer::from(15)], &proof));
//! ```
//!
//! Incorrect Usage:
//!
//! ```
//! let acc = Accumulator::<ClassGroup>::new();
//!
//! let prime1 = Integer::from(7);
//! let prime2 = Integer::from(7);
//!
//! // The same element should not be accumulated twice.
//! let (acc, proof) = acc.add(&[prime1, prime2]);
//!
//! let nonprime = Integer::from(10);
//!
//! // You should not accumulate composite numbers or 2.
//! let (acc, proof) = acc.add(&[nonprime]);
//!
//! ```
//!
//! Typical users of this library will only need to access `hash_to_prime` and public-facing
//! routines on `accumulator` and `vector_commitment`, but we also export several internal modules
//! for special applications. **Use at your own risk**.
//!
//! # Groups
//!
//! Cryptographic accumulators and vector commitments are powered by group operations. We provide
//! an RSA group implementation with the [RSA-2048 modulus](https://en.wikipedia.org/wiki/RSA_numbers#RSA-2048) and a form class group implementation with a fixed discriminant.
//!
//! The RSA group is fast but relies on the security of the RSA-2048 modulus and needs trusted
//! setup if using a different modulus. The class group implementation is slower but eliminates
//! the need for a trusted setup.
//!
//! # Disclaimer
//!
//! This implementation is research-grade and has not been audited. Use with caution.
//!
//! # References
//!
//! This library is based on the description of Accumulators and Vector Commitments in
//! _Batching Techniques for Accumulators with Applications to IOPs and Stateless Blockchains_
//! by Boneh, Bünz, and Fisch
//!
//! # Remarks on Performance
//!
//! Most accumulator or vector-commitment functions will bottleneck in hashing to large primes.
//! To alleviate this, we created a zero-allocation U256 type that uses the low-level
//! `mpn_` functions in GMP, which is used in `hash_to_prime`.

// TODO: Benchmark our U256 vs. 256-bit `rug::Integer` vs. Parity U256.
#![allow(clippy::unknown_clippy_lints)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::empty_enum)]

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate arrayref;

mod accumulator;
pub use accumulator::*;
mod vector_commitment;
pub use vector_commitment::*;

#[warn(missing_docs)]
pub mod group;

pub mod hash;
pub mod proof;
pub mod uint;
pub mod util;
