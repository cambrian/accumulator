use crate::group::{Group, InvertibleGroup};
use num::{BigInt, BigUint, Unsigned};
use num_traits::identities::{One, Zero};
use std::ops::Mul;

/// Trait definition to enable equality checks between group elements and other custom data
/// structures by means of byte array comparison.
pub trait ConvertBytes: Group {
  fn to_le_bytes(x: &Self::Elem) -> Vec<u8>;
  fn to_be_bytes(x: &Self::Elem) -> Vec<u8>;
}

/// We use the singleton pattern to fake type-level programming.
/// Self::Rep stores info that we would like to "reflect" from the type-level at runtime.
/// We use a separate type Self::Rep from Self so that Self can be an uninhabitable type and exist
/// purely at the type-level.
/// TODO: can we enforce Self to be uninhabitable?
pub trait Singleton {
  type Rep: 'static;
  fn rep() -> &'static Self::Rep;
}

pub fn bi<T>(val: T) -> BigInt
where
  BigInt: From<T>,
{
  BigInt::from(val)
}

pub fn bu<U: Unsigned>(val: U) -> BigUint
where
  BigUint: From<U>,
{
  BigUint::from(val)
}

/// Returns `(a, b, GCD(x, y))` s.t. `ax + by = GCD(x, y)`.
pub fn bezout(x: &BigUint, y: &BigUint) -> (BigInt, BigInt, BigInt) {
  let (mut s, mut old_s) = (bi(0), bi(1));
  let (mut t, mut old_t) = (bi(1), bi(0));
  let (mut r, mut old_r) = (bi(y.clone()), bi(x.clone()));

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
  let m_big = bi(m.clone());
  ((x % &m_big + &m_big) % &m_big)
    .to_biguint()
    .expect("positive BigInt expected")
}

pub fn product<T: Mul + One + Clone>(elems: &[&T]) -> T {
  elems.iter().fold(num::one(), |a, &b| a * b.clone())
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
    let x = bu(7u16);
    let y = bu(165u16);
    let (a, b, gcd) = bezout(&x, &y);
    assert!(gcd.is_one());
    assert!(a == bi(-47));
    assert!(b == bi(2));
  }

  #[test]
  fn test_mod_euc_big() {
    let r = mod_euc_big(&bi(-8), &(3u8));
    assert!(r == BigUint::one());
  }

  #[test]
  fn test_product() {
    let elems = [
      &bu(2u32),
      &bu(3u32),
      &bu(4u32),
      &bu(5u32),
      &bu(6u32),
      &bu(7u32),
    ];
    assert!(product(&elems) == bu(5040u32));
  }
}
