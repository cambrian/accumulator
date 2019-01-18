use num::BigInt;
use num::BigUint;
use num_traits::identities::Zero;
use num_traits::identities::One;
use serde::ser::Serialize;
use num::integer::Integer;
use std::marker::Sized;

pub mod class;
pub mod dummy;
pub mod rsa;

/// We need a runtime representation for the group itself because reading in group parameters
/// (i.e. RSA modulus) is infeasible to do at the type-level in Rust.
///
/// We treat Groups as singletons to get as close to mimicking type-level programming as possible.
///
/// TODO: Would be nice to have this implement a trait to get dot-notation for group operations.
/// Not sure how to do that tho
pub trait Group: Sized {
  // Consider: Add an operator trait (e.g. Mul for *) to Elem that defines the op itself?
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
  /// Repeated squaring algorithm. Implementations may override this (e.g. Montgomery multiplication
  /// for RSA groups) for performance reasons.
  fn exp(a: &Self::Elem, n: &BigUint) -> Self::Elem {
    if *n == BigUint::zero() { Self::id() }
    else if *n == BigUint::one() { a.clone() }
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
    // After further discussion: Writing a specialized inv() that takes an exponent is only a
    // marginal speedup over inv() then exp() in the negative exponent case. (That is, the
    // complexity does not change.)
    if *n >= BigInt::zero() {
      Self::exp(a, &n.to_biguint().expect("positive BigInt expected"))
    } else {
      Self::exp(
        &Self::inv(a),
        &(-n).to_biguint().expect("negative BigInt expected"),
      )
    }
  }
}
