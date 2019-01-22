use super::util::bi;
use super::util::Singleton;
use num::integer::Integer;
use num::{BigInt, BigUint};
use num_traits::identities::{One, Zero};
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
pub trait Group: Singleton {
  /// In theory the association Group::Elem is bijective, such that it makes sense to write
  /// something like Elem::Group::get(). This would let us define op, exp, inv, etc on the Elem
  /// type and avoid using prefix notation for all of our group operations.
  /// But afaik bijective associated types are not supported by Rust.
  type Elem: Eq + Serialize + Clone + Sized;

  fn id_(rep: &Self::Rep) -> Self::Elem;

  /// Returns an element with unknown order.
  /// REVIEW: Consider renaming Group such that it's clear that the group has elements of unknown
  /// order.
  /// REVIEW: Consider renaming to unknown_order_elem (maybe small_unknown_order_elem).
  /// E.g. 2, for RSA groups.
  fn base_elem_(rep: &Self::Rep) -> Self::Elem;

  fn op_(rep: &Self::Rep, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

  /// Default implementation of exponentiation via repeated squaring.
  /// Group implementations may provide more performant specializations
  /// (e.g. Montgomery multiplication for RSA groups).
  fn exp_(rep: &Self::Rep, a: &Self::Elem, n: &BigUint) -> Self::Elem {
    if *n == BigUint::zero() {
      Self::id()
    } else if *n == BigUint::one() {
      a.clone()
    } else if n.is_odd() {
      Self::op(a, &Self::exp_(rep, &Self::op(a, a), &(n >> 1)))
    } else {
      Self::exp_(rep, &Self::op(a, a), &(n >> 1))
    }
  }

  // -------------------
  // END OF REQUIRED FNS
  // -------------------

  fn id() -> Self::Elem {
    Self::id_(Self::rep())
  }

  fn base_elem() -> Self::Elem {
    Self::base_elem_(Self::rep())
  }

  fn op(a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
    Self::op_(Self::rep(), a, b)
  }

  fn exp(a: &Self::Elem, n: &BigUint) -> Self::Elem {
    Self::exp_(Self::rep(), a, n)
  }
}

/// Trait for groups that support efficient inverse calculations.
/// NOT used to mean a cyclic group (where every element has an inverse).
pub trait InvertibleGroup: Group {
  fn inv_(rep: &Self::Rep, a: &Self::Elem) -> Self::Elem;

  // -------------------
  // END OF REQUIRED FNS
  // -------------------

  fn inv(a: &Self::Elem) -> Self::Elem {
    Self::inv_(Self::rep(), a)
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

/// Not tested thoroughly, auditing/review welcome.
pub fn multi_exp<G: Group>(n: usize, alphas: &[G::Elem], x: &[BigInt]) -> G::Elem {
  if n == 1 {
    return alphas[0].clone();
  }
  let n_half: usize = n / 2;
  let alpha_l = &alphas[..n_half];
  let alpha_r = &alphas[n_half..];
  let x_l = &x[..n_half];
  let x_r = &x[n_half..];
  // G::op expects a BigUint.
  let x_star_l = (x_l.iter().fold(bi(1), |a, b| a * b)).to_biguint().unwrap();
  let x_star_r = (x_r.iter().fold(bi(1), |a, b| a * b)).to_biguint().unwrap();
  let l = multi_exp::<G>(n_half, alpha_l, x_l);
  let r = multi_exp::<G>(n - n_half, alpha_r, x_r);
  G::op(&G::exp(&l, &x_star_r), &G::exp(&r, &x_star_l))
}

#[cfg(test)]
mod tests {
  use super::dummy::DummyRSA;
  use super::*;

  #[test]
  fn test_multi_exp() {
    // TODO: Build more general testing framework.
    let alpha_1 = DummyRSA::elem_of(2);
    let alpha_2 = DummyRSA::elem_of(3);
    let x_1 = bi(3);
    let x_2 = bi(2);
    let res = multi_exp::<DummyRSA>(
      2,
      &[alpha_1.clone(), alpha_2.clone()],
      &[x_1.clone(), x_2.clone()],
    );
    assert!(res == DummyRSA::elem_of(108));
    let alpha_3 = DummyRSA::elem_of(5);
    let x_3 = bi(1);
    let res_2 = multi_exp::<DummyRSA>(3, &[alpha_1, alpha_2, alpha_3], &[x_1, x_2, x_3]);
    assert!(res_2 == DummyRSA::elem_of(1_687_500));
  }
}
