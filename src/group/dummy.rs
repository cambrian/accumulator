//! Dummy RSA group for 64-bit numbers.
//! Use this group for testing while we figure out ring integration.

use super::super::util;
use super::{Group, InvertibleGroup};
use num::BigUint;
use num_traits::cast::ToPrimitive;
use num_traits::identities::One;
use std::u64;

pub struct DummyRSA {
  modulus: u64,
}

const P: u64 = 226_022_213;
const Q: u64 = 12_364_769;

const DUMMY_RSA: DummyRSA = DummyRSA { modulus: P*Q };

impl Group for DummyRSA {
  type Elem = u64;
  fn get() -> Self {
    DUMMY_RSA
  }
  fn op_(&self, a: &u64, b: &u64) -> u64 {
    ((*a as u128 * *b as u128) % self.modulus as u128) as u64
  }
  fn id_(&self) -> u64 {
    1
  }
  fn base_elem_(&self) -> u64 {
    2
  }
}

impl InvertibleGroup for DummyRSA {
  // TODO: a potentially faster algorithm exists via Euler's theorem
  fn inv_(&self, x: &u64) -> u64 {
    let x_big = BigUint::from(*x);
    let mod_big = BigUint::from(DummyRSA::get().modulus);
    let (a, _, gcd) = util::bezout(&x_big, &mod_big);
    dbg!(&a);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    a.to_u64().expect("should be in u64 range")
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_op() {
    let a = DummyRSA::op(&2, &3);
    assert!(a == 6);
    let b = DummyRSA::op(&(P+1), &(Q+1));
    assert!(b == P+Q+1);
  }

  #[test]
  fn test_exp() {
    let a = DummyRSA::exp(&2, &From::<u64>::from(3));
    assert!(a == 8);
    let b = DummyRSA::exp(&2, &From::<u64>::from(128));
    assert!(b == 782_144_413_693_680);
  }

  #[test]
  fn test_inv() {
    let r = DummyRSA::inv(&2);
    assert!(r == 1397356226306899);
  }
}
