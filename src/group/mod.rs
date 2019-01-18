use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigInt;
use num::BigUint;
use std::marker::Sized;
use serde::ser::Serialize;

pub mod class;
pub mod rsa;
pub mod dummy;

/// We need a runtime representation for the group itself because reading in group parameters
/// (i.e. RSA modulus) is infeasible to do at the type-level in Rust.
/// We treat Groups as singletons to get as close to mimicking type-level programming as possible.
pub trait Group {
  type Elem; // Would be nice to have this implement a trait so we get dot-notation for group operations. Not sure how to do that tho
  fn init(/* TODO */) -> Self;
  fn op_(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
  fn id_(&self) -> Self::Elem;
  /// E.g. 2, for RSA groups.
  fn base_elem_(&self) -> Self::Elem;

  // Convenience functions below
  fn get() -> Self; // read singleton, or call init if not yet initialized
  fn id() -> Self::Elem {
    Self::id_(&Self::get())
  }
  fn base_elem() -> Self::Elem {
    Self::base_elem_(&Self::get())
  }
  fn op(a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
    Self::op_(&Self::get(), a, b)
  }
  fn exp(a: &Self::Elem, n: &BigUint) -> Self::Elem; // TODO
}

pub trait InvertibleGroup: Group {
  fn inv_(&self, a: &Self::Elem) -> Self::Elem;
  fn inv(a: &Self::Elem) -> Self::Elem {
    Self::inv_(&Self::get(), a)
  }
  fn exp_signed(a: &Self::Elem, n: &BigInt) -> Self::Elem {
    // REVIEW: check sign to avoid converting bigint->biguint twice in the negative case.
    // For Pablo: Sanjay mentioned an optimization that suggests changing one of these fn
    // signatures?
    match n.to_biguint() {
      Some(value) => Self::exp(a, &value),
      None => Self::exp(&Self::inv(&a), &(-n).to_biguint().expect("negative BigInt expected"))
    }
  }
}
