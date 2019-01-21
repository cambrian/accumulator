use num::integer::Integer;
use num::BigInt;
use num::BigUint;
use num_traits::identities::One;
use num_traits::identities::Zero;
use serde::ser::Serialize;
use std::marker::Sized;

pub mod class;
pub mod dummy;
pub mod rsa;

/// We need a runtime representation for the group itself because reading in group parameters
/// (i.e. RSA modulus) is infeasible to do at the type-level in Rust.
///
/// We mimic type-level programming by using the singleton pattern here. For each group, there
/// should be a single, constant, static instance of the group representation accessible at all
/// times. This way, we can "reflect" information about the group type by accessing the singleton.
/// Refer to dummy.rs for an example.
pub trait Group: Sized + 'static {
  /// In theory the association Group::Elem is bijective, such that it makes sense to write
  /// something like Elem::Group::get(). This would let us define op, exp, inv, etc on the Elem
  /// type and avoid using prefix notation for all of our group operations.
  /// But afaik bijective associated types are not supported by Rust.
  type Elem: Eq + Serialize + Clone + Sized;

  // TODO: possible to make private??
  fn get() -> &'static Self;

  fn id_(&'static self) -> Self::Elem;

  /// Returns an element with unknown order.
  /// REVIEW: Consider renaming Group such that it's clear that the group has elements of unknown
  /// order.
  /// REVIEW: Consider renaming to unknown_order_elem (maybe small_unknown_order_elem).
  /// E.g. 2, for RSA groups.
  fn base_elem_(&'static self) -> Self::Elem;

  fn op_(&'static self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

  /// Default implementation of exponentiation via repeated squaring.
  /// Group implementations may provide more performant specializations
  /// (e.g. Montgomery multiplication for RSA groups).
  fn exp_(&'static self, a: &Self::Elem, n: &BigUint) -> Self::Elem {
    if *n == BigUint::zero() {
      Self::id()
    } else if *n == BigUint::one() {
      a.clone()
    } else if n.is_odd() {
      Self::op(a, &Self::exp_(&self, &Self::op(a, a), &(n >> 1)))
    } else {
      Self::exp_(&self, &Self::op(a, a), &(n >> 1))
    }
  }

  // -------------------
  // END OF REQUIRED FNS
  // -------------------

  fn id() -> Self::Elem {
    Self::id_(Self::get())
  }

  fn base_elem() -> Self::Elem {
    Self::base_elem_(Self::get())
  }

  fn op(a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
    Self::op_(Self::get(), a, b)
  }

  fn exp(a: &Self::Elem, n: &BigUint) -> Self::Elem {
    Self::exp_(Self::get(), a, n)
  }
}

/// Trait for groups that support efficient inverse calculations.
/// NOT used to mean a cyclic group (where every element has an inverse).
pub trait InvertibleGroup: Group {
  fn inv_(&self, a: &Self::Elem) -> Self::Elem;

  // -------------------
  // END OF REQUIRED FNS
  // -------------------

  fn inv(a: &Self::Elem) -> Self::Elem {
    Self::inv_(Self::get(), a)
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
