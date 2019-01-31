/// TODO
use super::Group;
use crate::util::{int, TypeRep};
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::Identity;
use rug::integer::Order;
use rug::Integer;
use std::hash::{Hash, Hasher};

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum Ed25519 {}

lazy_static! {
  pub static ref MAX_SAFE_INTEGER: Integer = int(2 ^ 256);
}
/// Derive copy?
#[derive(Clone, Debug, Eq)]
pub struct Ed25519Elem(RistrettoPoint);

// TODO
impl Hash for Ed25519Elem {
  fn hash<H: Hasher>(&self, _state: &mut H) {}
}

// TODO
impl PartialEq for Ed25519Elem {
  fn eq(&self, other: &Ed25519Elem) -> bool {
    self.0 == other.0
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
    let mut digits: [u8; 32] = [0; 32];
    let mut result = x.clone();
    while remaining > *MAX_SAFE_INTEGER {
      MAX_SAFE_INTEGER.write_digits(&mut digits, Order::LsfLe);
      let factor = Scalar::from_bytes_mod_order(digits);
      result = Ed25519Elem(result.0 * factor);
      // TODO: find way to avoid clone
      remaining -= MAX_SAFE_INTEGER.clone();
    }
    remaining.write_digits(&mut digits, Order::LsfLe);
    let factor = Scalar::from_bytes_mod_order(digits);
    Ed25519Elem(result.0 * factor)
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
    let exp_512 = Ed25519::exp(&bp, &int(256));
    let exp_256 = Ed25519::exp(&bp, &int(128));
    let exp_256 = Ed25519::exp(&exp_256, &int(2));
    assert_eq!(exp_512, exp_256);
  }
}
