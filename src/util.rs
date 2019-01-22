use super::group::{Group, InvertibleGroup};
use num::{BigInt, BigUint, Unsigned};
use num_bigint::Sign::Plus;
use num_traits::identities::{One, Zero};

pub trait Singleton: 'static {
  // TODO: possible to make private??
  fn get() -> &'static Self;
}

/// Returns `(a, b, GCD(x, y))` s.t. `ax + by = GCD(x, y)`.
pub fn bezout(x: &BigUint, y: &BigUint) -> (BigInt, BigInt, BigInt) {
  let (mut s, mut old_s): (BigInt, BigInt) = (num::zero(), num::one());
  let (mut t, mut old_t): (BigInt, BigInt) = (num::one(), num::zero());
  let (mut r, mut old_r) = (
    BigInt::from_biguint(Plus, y.clone()),
    BigInt::from_biguint(Plus, x.clone()),
  );

  while !r.is_zero() {
    let quotient = &old_r / &r;
    let (temp_r, temp_s, temp_t) = (r, s, t);

    r = &old_r - &quotient * &temp_r;
    s = &old_s - &quotient * &temp_s;
    t = &old_t - &quotient * &temp_t;

    old_r = temp_r;
    old_s = temp_s;
    old_t = temp_t;
  }

  (old_s, old_t, old_r)
}

/// A homebrew version of the future Rust function mod_euc.
pub fn mod_euc_big<U: Unsigned + Clone>(x: &BigInt, m: &U) -> BigUint
where
  BigInt: From<U>,
{
  let m_big = BigInt::from(m.clone());
  ((x % &m_big + &m_big) % &m_big)
    .to_biguint()
    .expect("positive BigInt expected")
}

/// Trait definition to enable equality checks between group elements and other custom
/// data structures by means of byte array comparison.
pub trait ConvertBytes: Group {
  fn to_le_bytes(x: &Self::Elem) -> Vec<u8>;
  fn to_be_bytes(x: &Self::Elem) -> Vec<u8>;
}

/// REVIEW: generalize
pub fn product(elems: &[&BigUint]) -> BigUint {
  elems.iter().fold(num::one(), |a, b| a * *b)
}

/// Computes the `(xy)`th root of `g` given the `x`th and `y`th roots of `g` and `(x, y)` coprime.
pub fn shamir_trick<G: InvertibleGroup>(
  xth_root: &G::Elem,
  yth_root: &G::Elem,
  x: &BigUint,
  y: &BigUint,
) -> Option<G::Elem> {
  if G::exp(xth_root, x) != G::exp(yth_root, y) {
    return None;
  }

  let (a, b, gcd) = bezout(x, y);

  if !gcd.is_one() {
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
    let x = BigUint::from(7 as u16);
    let y = BigUint::from(165 as u16);
    let (a, b, gcd) = bezout(&x, &y);
    assert!(gcd.is_one());
    assert!(a == BigInt::from(-47 as i16));
    assert!(b == BigInt::from(2 as i16));
  }

  #[test]
  fn test_mod_euc_big() {
    let r = mod_euc_big(&BigInt::from(-8), &(3 as u8));
    assert!(r == BigUint::one());
  }

  #[test]
  fn test_product() {
    let elems = [
      &BigUint::from(2u32),
      &BigUint::from(3u32),
      &BigUint::from(4u32),
      &BigUint::from(5u32),
      &BigUint::from(6u32),
      &BigUint::from(7u32),
    ];
    assert!(product(&elems) == BigUint::from(5040u32));
  }
}
