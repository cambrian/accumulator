/// TODO
use super::Group;
use crate::util::TypeRep;
use curve25519_dalek::constants::BASEPOINT_ORDER;
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::Identity;
use rug::Integer;
use std::hash::{Hash, Hasher};

// REVIEW: rename to Ed25519 since the "Ed" is not an acronym (it's short for "Edwards" iirc)
#[derive(Debug, PartialEq, Eq, Clone)]
pub enum ED25519 {}

lazy_static! {
  pub static ref rp: RistrettoPoint = RistrettoPoint::identity();
}

/// Derive copy?
// REVIEW: rename to Ed25519 since the "Ed" is not an acronym (it's short for "Edwards" iirc)
#[derive(Clone, Debug, Eq)]
pub struct ED25519Elem(RistrettoPoint);

// TODO
impl Hash for ED25519Elem {
  fn hash<H: Hasher>(&self, _state: &mut H) {}
}

// TODO
impl PartialEq for ED25519Elem {
  fn eq(&self, _other: &ED25519Elem) -> bool {
    true
  }
}

/// Unsure of correct rep for this group. BASEPOINT_ORDER is described as the order of the Ristretto
/// Group. I don't believe this will actually be used in the ops since it is likely used internally
/// but it seemed like the best fit for TypeRep. Review welcome.
/// REVIEW: If you never use the type rep, just make Rep = ()
impl TypeRep for ED25519 {
  type Rep = Scalar;
  fn rep() -> &'static Self::Rep {
    &BASEPOINT_ORDER
  }
}

impl Group for ED25519 {
  type Elem = ED25519Elem;

  fn op_(_: &Scalar, a: &ED25519Elem, b: &ED25519Elem) -> ED25519Elem {
    ED25519Elem(a.0 + b.0)
  }

  fn id_(_: &Scalar) -> ED25519Elem {
    ED25519Elem(RistrettoPoint::identity())
  }

  fn inv_(_: &Scalar, _x: &ED25519Elem) -> ED25519Elem {
    // assert(not point at infinity)
    // Convert x.Y to -x.Y, but I believe this requires control of FieldElems, which are in the
    // private crate field.
    unimplemented!();
  }

  fn exp_(_: &Scalar, _x: &ED25519Elem, _n: &Integer) -> ED25519Elem {
    // Need to implement Integer -> Scalar, or x * Integer
    // ED25519(x.0 * n)
    unimplemented!();
  }
}
