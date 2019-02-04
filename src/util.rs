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

  // Could use `gcd_cofactors_mut` here, but it looks messy and doesn't do much.
  let (gcd, a, b) = x.clone().gcd_cofactors(y.clone(), Integer::new());

  if gcd != int(1) {
    return None;
  }

  Some(G::op(&G::exp(xth_root, &b), &G::exp(yth_root, &a)))
}

// Solve a linear congruence of form `ax = b mod m` for the set of solutions x,
// characterized by integers mu and v such that x = mu + vn where n is any integer.
// REVIEW: instead of panicking, return Result<(Integer, Integer), ErrType>
pub fn solve_linear_congruence(a: &Integer, b: &Integer, m: &Integer) -> (Integer, Integer) {
  // g = gcd(a, m) => da + em = g
  let (g, d, _) = a.clone().gcd_cofactors(m.clone(), Integer::new());

  // q = floor_div(b, g)
  // r = b % g
  let (q, r) = b.clone().div_rem_floor(g.clone());
  if r != Integer::from(0) {
    panic!("No solution to linear congruence.");
  }

  // mu = (q * d) % m
  // v = m / g
  let mu = (q * d) % m;
  let (v, _) = m.clone().div_rem_floor(g);
  (mu, v)
}

// Compute GCD of three integers.
// REVIEW: you only use this once so remove this fn and just use the body inline
#[inline]
pub fn three_gcd(a: &Integer, b: &Integer, c: &Integer) -> Integer {
  a.clone().gcd(&b).gcd(&c)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, UnknownOrderGroup, RSA2048};
  use crate::util::int;

  #[test]
  fn test_linear_congruence_solver_basic() {
    assert_eq!(
      (Integer::from(-2), Integer::from(4)),
      solve_linear_congruence(&Integer::from(3), &Integer::from(2), &Integer::from(4))
    );

    assert_eq!(
      (Integer::from(-2), Integer::from(4)),
      solve_linear_congruence(&Integer::from(3), &Integer::from(2), &Integer::from(4))
    );

    assert_eq!(
      (Integer::from(1), Integer::from(2)),
      solve_linear_congruence(&Integer::from(5), &Integer::from(1), &Integer::from(2))
    );

    assert_eq!(
      (Integer::from(-3), Integer::from(5)),
      solve_linear_congruence(&Integer::from(2), &Integer::from(4), &Integer::from(5))
    );
  }

  #[test]
  fn test_shamir_trick() {
    let (x, y, z) = (&int(13), &int(17), &int(19));
    let xth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(y * z));
    let yth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(x * z));
    let xyth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), z);
    assert!(shamir_trick::<RSA2048>(&xth_root, &yth_root, x, y) == Some(xyth_root));
  }

  #[test]
  fn test_shamir_trick_failure() {
    let (x, y, z) = (&int(7), &int(14), &int(19)); // Inputs not co-prime.
    let xth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(y * z));
    let yth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(x * z));
    assert!(shamir_trick::<RSA2048>(&xth_root, &yth_root, x, y) == None);
  }
}
