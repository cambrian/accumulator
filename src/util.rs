use crate::group::Group;
use rug::Integer;

/// We use the singleton pattern to fake type-level programming.
/// Self::Rep stores info that we would like to "reflect" from the type-level at runtime.
/// We use a separate type Self::Rep from Self so that Self can be an uninhabitable type and exist
/// purely at the type-level.
pub trait Singleton {
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
pub fn shamir_trick<G: Group>(
  xth_root: &G::Elem,
  yth_root: &G::Elem,
  x: &Integer,
  y: &Integer,
) -> Option<G::Elem> {
  if G::exp(xth_root, x) != G::exp(yth_root, y) {
    return None;
  }

  let (gcd, a, b) = x.clone().gcd_cofactors(y.clone(), Integer::new());

  if gcd != int(1) {
    return None;
  }

  Some(G::op(&G::exp(xth_root, &b), &G::exp(yth_root, &a)))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, UnknownOrderGroup, RSA2048};
  use crate::util::int;

  #[test]
  fn test_shamir_trick() {
    let (x, y, z) = (&int(13), &int(17), &int(19));
    let xth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(y * z));
    let yth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(x * z));
    let xyth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), z);
    assert!(shamir_trick::<RSA2048>(&xth_root, &yth_root, x, y) == Some(xyth_root));
  }
}
