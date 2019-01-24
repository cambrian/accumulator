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

/// Returns `(a, b, GCD(x, y))` s.t. `ax + by = GCD(x, y)`.
pub fn bezout(x: &Integer, y: &Integer) -> (Integer, Integer, Integer) {
  let (mut s, mut old_s) = (int(0), int(1));
  let (mut t, mut old_t) = (int(1), int(0));
  let (mut r, mut old_r) = (y.clone(), x.clone());

  while r != int(0) {
    let quotient = old_r.clone() / &r;
    let (temp_r, temp_s, temp_t) = (r, s, t);

    r = old_r.clone() - &quotient * &temp_r;
    s = old_s.clone() - &quotient * &temp_s;
    t = old_t.clone() - &quotient * &temp_t;

    old_r = temp_r;
    old_s = temp_s;
    old_t = temp_t;
  }

  (old_s, old_t, old_r)
}

// /// A homebrew version of the future Rust function mod_euc.
// pub fn mod_euc_big<T: Clone>(x: &Integer, m: &T) -> Integer
// where
//   Integer: From<T>,
// {
//   let m_big = Integer::from(m.clone());
//   (x % &m_big + &m_big) % &m_big
// }

/// Computes the `(xy)`th root of `g` given the `x`th and `y`th roots of `g` and `(x, y)` coprime.
pub fn shamir_trick<G: Group>(
  xth_root: &G::Elem,
  yth_root: &G::Elem,
  x: &Integer,
  y: &Integer,
) -> Option<G::Elem> {
  if G::exp(xth_root, x) != G::exp(yth_root, y) {
    return None;
  }

  let (a, b, gcd) = bezout(x, y);

  if gcd != int(1) {
    return None;
  }

  Some(G::op(
    &G::exp_signed(xth_root, &b),
    &G::exp_signed(yth_root, &a),
  ))
}

#[cfg(test)]
mod tests {
  use super::*;
  use num_traits::identities::One;

  #[test]
  fn test_bezout() {
    let x = int(7);
    let y = int(165);
    let (a, b, gcd) = bezout(&x, &y);
    assert!(gcd.is_one());
    assert!(a == int(-47));
    assert!(b == int(2));
  }

  #[test]
  fn test_mod_euc_big() {
    let r = mod_euc_big(&int(-8), &(3u8));
    assert!(r == BigUint::one());
  }

  #[test]
  fn test_product() {
    let elems = [int(2), int(3), int(4), int(5), int(6), int(7)];
    assert!(product(&elems) == int(5040));
  }
}
