//! This module contains implementations for different succinct proofs for hidden-order groups.
//! They are used as building blocks for many of the cryptographic primitives in this library.
//! Use with caution.
//!
//! Implementations are based on section 3 of _Batching Techniques for Accumulators
//! with Applications to IOPs and Stateless Blockchains_.
mod poe;
pub use poe::Poe;
mod pokcr;
pub use pokcr::Pokcr;
mod poke2;
pub use poke2::Poke2;
