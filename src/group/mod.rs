use num::BigInt;
use num::BigUint;
use num_traits::identities::Zero;
use serde::ser::Serialize;
use std::marker::Sized;

pub mod class;
pub mod dummy;
pub mod rsa;

/// We need a runtime representation for the group itself because reading in group parameters
/// (i.e. RSA modulus) is infeasible to do at the type-level in Rust.
/// We treat Groups as singletons to get as close to mimicking type-level programming as possible.
pub trait Group: Sized {
  // Consider: Add an operator trait (e.g. Mul for *) to Elem that defines the op itself?
  type Elem: Eq + Serialize + Clone + Sized;

  fn init() -> Self;
  fn op_(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
  fn id_(&self) -> Self::Elem;
  /// E.g. 2, for RSA groups.
  fn base_elem_(&self) -> Self::Elem;

  // Convenience functions below:
  fn get() -> Self; // Read singleton, or call `init` if not yet initialized.
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
