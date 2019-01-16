use alga::general::AbstractGroup;
use alga::general::Operator;
use bigint::uint::U256;

pub mod class;
pub mod rsa;

/// Picks an element from the set generating the relevant group.
pub trait Generator<O: Operator>: AbstractGroup<O> {
  fn generator() -> Self;
}

/// Efficiently computes the inverse of a group element.
pub trait Inverse<O: Operator>: AbstractGroup<O> {
  fn inverse(&self) -> Self;
}

/// Exponentiation (repeated operation) in a group.
pub trait Pow<O: Operator>: AbstractGroup<O> {
  // TODO: Write default impl using repeated squaring.
  fn pow(&self, exp: U256) -> Self;
}
