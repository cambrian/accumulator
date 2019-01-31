/// TODO
use super::Group;
use crate::util::TypeRep;
use curve25519_dalek::ristretto::RistrettoPoint;
// use curve25519_dalek::subtle::ConstantTimeEq;
use curve25519_dalek::traits::Identity;
use rug::Integer;
use std::hash::{Hash, Hasher};

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum Ed25519 {}

lazy_static! {
  pub static ref rp: RistrettoPoint = RistrettoPoint::identity();
}

/// Derive copy?
#[derive(Clone, Debug, Eq)]
pub struct Ed25519Elem(RistrettoPoint);

// TODO
impl Hash for Ed25519Elem {
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.0.compress().as_bytes().hash(state);
  }
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
    // assert(not point at infinity)
    Ed25519Elem(-x.0)
  }

  fn exp_(_: &(), _x: &Ed25519Elem, _n: &Integer) -> Ed25519Elem {
    // Need to implement Integer -> Scalar, or x * Integer
    // Ed25519(x.0 * n)
    unimplemented!();
  }
}
