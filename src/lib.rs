//! Fast cryptographic accumulator and vector commitment library, originally written and maintained
//! by Cambrian Labs.
//!
//! **Disclaimer**: This library is intended to be production-quality code, but it has not been
//! independently-audited for correctness or tested to our satisfaction. As such, please **treat
//! this library as research-grade**.
//!
//! # Important Note
//!
//! In order for accumulators to work and be secure in your application, you need to take care of an
//! important invariant: Cryptographic accumulators as implemented here assume that **no element is
//! accumulated twice**.
//!
//! # What is an accumulator?
//!
//! An accumulator is a cryptographic primitive that "accumulates" values and lets you provide
//! efficiently-verifiable proofs to other parties that certain values have or have not been
//! accumulated without needing to provide the entire accumulated set.
//!
//! An accumulator is similar to a set or Merkle Tree but it is stored in _constant-size_ and
//! generates inclusion and exclusion proofs in constant-size as well. For the mathematics behind
//! accumulators, see _Batching Techniques for Accumulators with Applications to IOPs and Stateless
//! Blockchains_ (Boneh, Bünz, and Fisch 2018) [[Link]](https://eprint.iacr.org/2018/1188.pdf).
//!
//! Throughout our code, we refer to this paper as `BBF`. We also refer to _Universal Accumulators
//! with Efficient Nonmembership Proofs_ (Li, Li, Xue 2007)
//! [[Link]](https://link.springer.com/content/pdf/10.1007/978-3-540-72738-5_17.pdf), which we
//! abbreviate as `LLX`.
//!
//! # What is a vector commitment?
//!
//! A vector commitment is like an accumulator, except that it is position binding. You can not
//! only claim that an element has been accumulated, but also at what _position_ it has been
//! accumulated.
//!
//! Our vector commitment implementation is very much a work-in-progress (WIP), and should be
//! treated with even more skepticism than our accumulators.
//!
//! # Usage
//! ```
//! // A very basic example.
//! use accumulator::Accumulator;
//! use accumulator::group::Rsa2048;
//!
//! let acc = Accumulator::<Rsa2048, &'static str>::empty();
//!
//! // Accumulate "dog" and "cat". The `add_with_proof` method returns the new accumulator state
//! // and a proof that you accumulated "dog" and "cat".
//! let (acc, proof) = acc.add_with_proof(&["dog", "cat"]);
//!
//! // A network participant who sees (acc, proof, and ["dog", "cat"]) can verify that the update
//! // was formed correctly ...
//! assert!(acc.verify_membership_batch(&["dog", "cat"], &proof));
//!
//! // ... and trying to verify something that has not been accumulated will not work.
//! assert!(!acc.verify_membership(&"cow", &proof));
//! ```
//!
//! Typical users of this library will mostly access public-facing routines on `accumulator` and
//! `vector_commitment`, but we also export internal modules for useful traits, types (such as the
//! `Rsa2048` group), and specialized procedures. **Use internal components at your own risk**.
//!
//! You can find a more interesting application of our library
//! [here](https://github.com/cambrian/accumulator-demo), where we create a proof-of-concept for
//! stateless Bitcoin nodes!
//!
//! # Groups
//!
//! Cryptographic accumulators and vector commitments are powered by group operations. We provide
//! an RSA group implementation with the
//! [RSA-2048 modulus](https://en.wikipedia.org/wiki/RSA_numbers#RSA-2048) and a form class group
//! implementation with a fixed discriminant, generated by OpenSSL.
//!
//! The RSA group is fast but relies on the security of the RSA-2048 modulus and needs trusted
//! setup if using a different modulus. The class group implementation is slower but eliminates
//! the need for a trusted setup.
//!
//! # Performance
//!
//! Most accumulator or vector-commitment functions will bottleneck in hashing to large primes. To
//! alleviate this, we created a zero-allocation `U256` type that uses the low-level `mpn_`
//! functions in [GMP](https://gmplib.org). Our `hash_to_prime` uses this type internally.
//!
//! Class groups are currently not performant for any meaningful use case. A PR is in the
//! works to drastically speed them up, using techniques learned from the
//! [Chia VDF competition](https://github.com/Chia-Network/vdf-competition).
#![allow(clippy::unknown_clippy_lints)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::empty_enum)]
#![warn(missing_docs)]

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
#[allow(missing_docs)]
pub mod uint;
pub mod util;
