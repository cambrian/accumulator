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
  pub static ref MAX_SAFE_EXPONENT: Integer = int(2 ^ 256);
}

impl Ed25519 {
  fn max_safe_exponent() -> &'static Integer {
    &MAX_SAFE_EXPONENT
  }
}

/// Derive copy?
#[derive(Clone, Debug, Eq)]
pub struct Ed25519Elem(RistrettoPoint);

impl Hash for Ed25519Elem {
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.0.compress().as_bytes().hash(state);
  }
}

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
    while remaining > *Ed25519::max_safe_exponent() {
      Ed25519::max_safe_exponent().write_digits(&mut digits, Order::LsfLe);
      let factor = Scalar::from_bytes_mod_order(digits);
      result = Ed25519Elem(result.0 * factor);
      remaining -= Ed25519::max_safe_exponent();
    }
    remaining.write_digits(&mut digits, Order::LsfLe);
    let factor = Scalar::from_bytes_mod_order(digits);
    Ed25519Elem(result.0 * factor)
  }
}
