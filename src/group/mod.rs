use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigUint;

pub mod class;
pub mod rsa;

/// Selection of an element in the generating set.
pub trait Generator<O: Operator>: AbstractGroup<O> {
  fn generator() -> Self;
}

/// Efficient computation of group inverses.
pub trait Inverse<O: Operator>: AbstractGroup<O> {
  fn inverse(&self) -> Self;
}

/// Efficient exponentiation in a group.
pub trait Pow<O: Operator>: AbstractGroup<O> {
  // TODO: Write default impl using repeated squaring.
  fn pow(&self, exp: BigUint) -> Self;
}
