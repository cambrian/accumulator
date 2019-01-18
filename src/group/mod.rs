use num::BigInt;
use num::BigUint;
use std::marker::Sized;

pub mod class;
pub mod rsa;

/// Abstract group interface.
/// TODO: See which references can/should be made consuming to work with ring.
pub trait Group: Sized + Eq + Clone {
  fn identity() -> Self;
  fn op(&self, rhs: &Self) -> Self;

  /// May be overridden for performant specializations.
  /// (E.g. To do Montgomery multiplication for the RSA group)
  /// TODO: Default implementation with repeated squaring.
  fn exp(&self, n: &BigUint) -> Self;
}

/// Groups that support efficient inversions.
/// E.g. RSA groups using Fermat's Little Theorem
pub trait InvertibleGroup: Group {
  fn inv(&self) -> Self;
  fn exp_signed(&self, n: &BigInt) -> Self {
    // None case implies a negative value.
    // TODO: Make inv() take an exponent and change this accordingly.
    match n.to_biguint() {
      Some(value) => self.exp(&value),
      None => self
        .inv()
        .exp(&(-n).to_biguint().expect("negative BigInt expected")),
    }
  }
}

/// Groups that have generators.
/// TODO: Remove and move into Group.
pub trait CyclicGroup: Group {
  fn generator() -> Self;
}
