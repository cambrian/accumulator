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

// const P: u64 = 2_413_575_613;
// const Q: u64 = 1_090_574_917;

const DUMMY_RSA: DummyRSA = DummyRSA {
  modulus: 2_632_185_023_820_699_121,
};

#[derive(PartialEq, Eq, Clone, Serialize)]
struct DummyRSAElem(pub u64);

impl Group for DummyRSA {
  type Elem = u64;
  fn get() -> Self {
    DUMMY_RSA
  }
  fn op_(&self, a: &u64, b: &u64) -> u64 {
    (a * b) % self.modulus
  }
  fn id_(&self) -> u64 {
    1
  }
  fn base_elem_(&self) -> u64 {
    2
  }
}

impl InvertibleGroup for DummyRSA {
  fn inv_(&self, x: &u64) -> u64 {
    let x_big = BigUint::from(*x);
    let mod_big = BigUint::from(DummyRSA::get().modulus);
    let (a, _, gcd) = util::bezout(&x_big, &mod_big);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    a.to_u64().expect("should be in u64 range")
  }
}

#[cfg(test)]
mod tests {
  use super::DummyRSA;
  use super::InvertibleGroup;
  // TODO: Test op and solve overflow.

  #[test]
  fn test_inverse() {
    let inv_result = DummyRSA::inv(&32_416_188_490);
    assert!(inv_result == 1_312_590_163_415_190_100);
  }
}
