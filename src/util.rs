use super::group::Group;
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

/// Not tested thoroughly, auditing/review welcome.
pub fn multi_exp<G: Group>(n: u16, alphas: &[G::Elem], x: &[BigInt]) -> G::Elem {
  if n == 1 {
    return alphas[0].to_owned();
  }
  let n_half = n / 2;
  let alpha_l = &alphas[..n_half as usize];
  let alpha_r = &alphas[n_half as usize..];
  let x_l = &x[..n_half as usize];
  let x_r = &x[n_half as usize..];
  let x_star_l = (x_l.iter().fold(BigInt::from(1 as u8), |a, b| a * b))
    .to_biguint()
    .unwrap();
  let x_star_r = (x_r.iter().fold(BigInt::from(1 as u8), |a, b| a * b))
    .to_biguint()
    .unwrap();
  let l = multi_exp::<G>(n_half, alpha_l, x_l);
  let r = multi_exp::<G>(n - n_half, alpha_r, x_r);
  G::op(&G::exp(&l, &x_star_r), &G::exp(&r, &x_star_l))
}

#[cfg(test)]
mod tests {
  use super::super::group::dummy::{DummyRSA, DummyRSAElem};
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
  fn test_multi_exp() {
    // TODO: Build more general testing framework
    let alpha_1 = DummyRSAElem::of(2);
    let alpha_2 = DummyRSAElem::of(3);
    let x_1 = BigInt::from(3 as u8);
    let x_2 = BigInt::from(2 as u8);
    let res = multi_exp::<DummyRSA>(
      2,
      &[alpha_1.clone(), alpha_2.clone()],
      &[x_1.clone(), x_2.clone()],
    );
    assert!(res == DummyRSAElem::of(108));
    let alpha_3 = DummyRSAElem::of(5);
    let x_3 = BigInt::from(1 as u8);
    let res_2 = multi_exp::<DummyRSA>(3, &[alpha_1, alpha_2, alpha_3], &[x_1, x_2, x_3]);
    assert!(res_2 == DummyRSAElem::of(1_687_500));
  }

  #[test]
  fn test_mod_euc_big() {
    let r = mod_euc_big(&BigInt::from(-8), &(3 as u8));
    assert!(r == BigUint::one());
  }
}
