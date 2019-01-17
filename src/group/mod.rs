use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigInt;
use num::BigUint;
use std::marker::Sized;

pub mod class;
pub mod rsa;

/// Abstract group interface
/// Not using alga because we need to include some performance optimizations
/// TODO: see which references can/should be made consuming to work with ring.
pub trait Group: Sized + Eq {
  fn identity() -> Self;
  fn op(&self, rhs: &Self) -> Self;

  /// May be overridden for performant specializations
  /// (e.g. to do montgomery multiplication for rsa group)
  fn exp(&self, n: &BigUint) -> Self;
}

/// Groups that support efficient inversions.
/// E.g. RSA groups using Fermat's Little Theorem
pub trait InvertibleGroup: Group {
  fn inv(&self) -> Self;
  fn exp_signed(&self, n: &BigInt) -> Self {
    // REVIEW: check sign to avoid converting bigint->biguint twice in the negative case.
    match n.to_biguint() {
      Some(value) => self.exp(&value),
      None => self.inv().exp(&(-n).to_biguint().expect("negative BigInt expected"))
    }
  }
}

/// Groups that have generators
pub trait CyclicGroup: Group {
  fn generator() -> Self;
}
