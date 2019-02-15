use crate::group::Group;
use rug::{Assign, Integer};

/// Poor man's type-level programming.
/// This trait allows us to reflect "type-level" (i.e. static) information at runtime.
pub trait TypeRep: 'static {
  type Rep: 'static;
  fn rep() -> &'static Self::Rep;
}

pub fn int<T>(val: T) -> Integer
where
  Integer: From<T>,
{
  Integer::from(val)
}

/// Computes the `(xy)`th root of `g` given the `x`th and `y`th roots of `g` and `(x, y)` coprime.
/// Consider moving this to accumulator?
#[allow(clippy::similar_names)]
pub fn shamir_trick<G: Group>(
  xth_root: &G::Elem,
  yth_root: &G::Elem,
  x: &Integer,
  y: &Integer,
  gcd: &mut Integer,
  a: &mut Integer,
  b: &mut Integer,
) -> Option<G::Elem> {
  if G::exp(xth_root, x) != G::exp(yth_root, y) {
    return None;
  }

  (&mut *gcd, &mut *a, &mut *b).assign(x.gcd_cofactors_ref(&y));

  if *gcd != int(1) {
    return None;
  }

  Some(G::op(&G::exp(xth_root, &b), &G::exp(yth_root, &a)))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, Rsa2048, UnknownOrderGroup};
  use crate::util::int;

  #[test]
  fn test_shamir_trick() {
    let (x, y, z) = (int(13), int(17), &int(19));
    let xth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(y.clone() * z));
    let yth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(x.clone() * z));
    let xyth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), z);
    let (mut gcd, mut a, mut b) = (Integer::new(), Integer::new(), Integer::new());
    assert!(
      shamir_trick::<Rsa2048>(&xth_root, &yth_root, &x, &y, &mut gcd, &mut a, &mut b)
        == Some(xyth_root)
    );
  }

  // #[test]
  // fn test_shamir_trick_failure() {
  //   let (x, y, z) = (&int(7), &int(14), &int(19)); // Inputs not coprime.
  //   let xth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(y * z));
  //   let yth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(x * z));
  //   assert!(shamir_trick::<Rsa2048>(&xth_root, &yth_root, x, y) == None);
  // }
}
