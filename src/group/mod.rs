use num::BigInt;
use num::BigUint;
use std::marker::Sized;
use serde::ser::Serialize;
use num::integer::Integer;

pub mod class;
pub mod rsa;
pub mod dummy;

/// We need a runtime representation for the group itself because reading in group parameters
/// (i.e. RSA modulus) is infeasible to do at the type-level in Rust.
/// We treat Groups as singletons to get as close to mimicking type-level programming as possible.
/// TODO: Would be nice to have this implement a trait to get dot-notation for group operations.
/// Not sure how to do that tho
pub trait Group: Sized {
  type Elem: Eq + Serialize + Clone + Sized;
  /// This function should return either a const or lazy_static group representation.
  /// This should be replaced by a const fn when they are added to Rust.
  fn get() -> Self;
  fn op_(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
  fn id_(&self) -> Self::Elem;
  /// E.g. 2, for RSA groups.
  fn base_elem_(&self) -> Self::Elem;

  // Convenience functions below
  fn id() -> Self::Elem {
    Self::id_(&Self::get())
  }
  fn base_elem() -> Self::Elem {
    Self::base_elem_(&Self::get())
  }
  fn op(a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
    Self::op_(&Self::get(), a, b)
  }
  fn exp(a: &Self::Elem, n: &BigUint) -> Self::Elem {
    if n == &From::<u32>::from(0) { Self::id() }
    else if n == &From::<u32>::from(1) { a.clone() }
    else if n.is_odd() { Self::op(a, &Self::exp(&Self::op(a, a), &(n >> 1))) }
    else { Self::exp(&Self::op(a, a), &(n >> 1)) }
  }
}

pub trait InvertibleGroup: Group where {
  fn inv_(&self, a: &Self::Elem) -> Self::Elem;
  fn inv(a: &Self::Elem) -> Self::Elem {
    Self::inv_(&Self::get(), a)
  }
  fn exp_signed(a: &Self::Elem, n: &BigInt) -> Self::Elem {
    // REVIEW: check sign to avoid converting bigint->biguint twice in the negative case.
    match n.to_biguint() {
      Some(value) => Self::exp(a, &value),
      None => Self::exp(&Self::inv(&a), &(-n).to_biguint().expect("negative BigInt expected"))
    }
  }
}
