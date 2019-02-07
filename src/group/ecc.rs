use super::Group;
use crate::util::{int, TypeRep};
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::Identity;
use rug::integer::Order;
use rug::ops::Pow;
use rug::Integer;
use std::hash::{Hash, Hasher};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Ed25519 {}

lazy_static! {
  pub static ref MAX_SAFE_EXPONENT: Integer = int(2).pow(255) - 1;
  pub static ref MAX_SAFE_SCALAR: Scalar = {
    let mut digits: [u8; 32] = [0; 32];
    MAX_SAFE_EXPONENT.write_digits(&mut digits, Order::LsfLe);
    Scalar::from_bytes_mod_order(digits)
  };
}

impl Ed25519 {
  fn max_safe_exponent() -> &'static Integer {
    &MAX_SAFE_EXPONENT
  }
}

/// REVIEW: Ideally we'd just use RistrettoPoint here, but only traits defined in this crate can
/// be implemented for arbitrary types. How to fix without wrapping?
///
/// It may make sense to fork curve25519-dalek to add the Hash impl. Then we won't need to wrap it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ed25519Elem(RistrettoPoint);

#[allow(clippy::derive_hash_xor_eq)]
impl Hash for Ed25519Elem {
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.0.compress().as_bytes().hash(state);
  }
}

impl TypeRep for Ed25519 {
  type Rep = ();
  fn rep() -> &'static Self::Rep {
    &()
  }
}

impl Group for Ed25519 {
  type Elem = Ed25519Elem;

  fn op_(_: &(), a: &Ed25519Elem, b: &Ed25519Elem) -> Ed25519Elem {
    Ed25519Elem(a.0 + b.0)
  }

  fn id_(_: &()) -> Ed25519Elem {
    Ed25519Elem(RistrettoPoint::identity())
  }

  fn inv_(_: &(), x: &Ed25519Elem) -> Ed25519Elem {
    Ed25519Elem(-x.0)
  }

  fn exp_(_: &(), x: &Ed25519Elem, n: &Integer) -> Ed25519Elem {
    let mut remaining = n.clone();
    let mut result = Ed25519::id();

    while remaining > *MAX_SAFE_EXPONENT {
      result = Ed25519Elem(result.0 + x.0 * (*MAX_SAFE_SCALAR));
      remaining -= Ed25519::max_safe_exponent();
    }

    let mut digits: [u8; 32] = [0; 32];
    remaining.write_digits(&mut digits, Order::LsfLe);
    let factor = Scalar::from_bytes_mod_order(digits);
    Ed25519Elem(result.0 + x.0 * factor)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::util::int;
  use curve25519_dalek::constants;

  #[test]
  fn test_inv() {
    let bp = Ed25519Elem(constants::RISTRETTO_BASEPOINT_POINT);
    let bp_inv = Ed25519::inv(&bp);
    assert!(Ed25519::op(&bp, &bp_inv) == Ed25519::id());
    assert_ne!(bp, bp_inv);
  }

  #[test]
  fn test_exp() {
    let bp = Ed25519Elem(constants::RISTRETTO_BASEPOINT_POINT);
    let exp_a = Ed25519::exp(&bp, &int(2).pow(258));
    let exp_b = Ed25519::exp(&bp, &int(2).pow(257));
    let exp_b_2 = Ed25519::exp(&exp_b, &int(2));
    assert_eq!(exp_a, exp_b_2);
  }
}
