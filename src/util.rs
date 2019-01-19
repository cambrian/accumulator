use num::{BigInt, BigUint, Unsigned};
use num_bigint::Sign::Plus;
use num_traits::identities::Zero;

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

// pub fn mod_euc<S: Signed, U: Unsigned + Clone>(x: &S, m: &U) -> U
// where
//   S: From<U>,
//   U: TryFrom<S>,
// {
//   let m_big = S::from(m.clone());
//   U::from((*x % m_big + m_big) % m_big)
// }

pub fn mod_euc_big<U: Unsigned + Clone>(x: &BigInt, m: &U) -> BigUint
where
  BigInt: From<U>,
{
  let m_big = BigInt::from(m.clone());
  ((x % &m_big + &m_big) % &m_big)
    .to_biguint()
    .expect("positive BigInt expected")
}

#[cfg(test)]
mod tests {
  use super::*;
  use num_traits::identities::One;

  #[test]
  fn test_bezout_simple() {
    let x = BigUint::from(7 as u16);
    let y = BigUint::from(165 as u16);
    let (a, b, gcd) = bezout(&x, &y);
    assert!(gcd.is_one());
    assert!(a == BigInt::from(-47 as i16));
    assert!(b == BigInt::from(2 as i16));
  }
}
