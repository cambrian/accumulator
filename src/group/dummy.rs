//! Dummy RSA group for 64-bit numbers
//! Use this group for testing while we figure out ring integration

use std::u64;
use super::{Group, InvertibleGroup};

// you might need to add something to support inversion
// (e.g. p or q for use with euler's theorem)
pub struct DummyRSA {
  modulus: u64
}

const P: u64 = 3_728_871_397;
const Q: u64 = 2_235_920_447;

const DUMMY_RSA: DummyRSA = DummyRSA { modulus: P*Q };

impl Group for DummyRSA {
  type Elem = u64;
  fn get() -> Self {
    DUMMY_RSA
  }
  fn op_(&self, a: &u64, b: &u64) -> u64{
    ((*a as u128 * *b as u128) % self.modulus as u128) as u64
  }
  fn id_(&self) -> u64 {
    1
  }
  fn base_elem_(&self) -> u64 {
    2
  }
}

// TODO
// impl InvertibleGroup for DummyRSA {
//   fn inv_(&self, a: &u64) {

//   }
// }

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_op() {
    let a = 2;
    let b = 3;
    let c = DummyRSA::op(&a, &b);
    assert!(c == 6);
    let x = P+1;
    let y = Q+1;
    let z = DummyRSA::op(&x, &y);
    assert!(z == P+Q+1);
  }

  #[test]
  fn test_exp() {
    let a = 2;
    let na = From::<u64>::from(3);
    let ra = DummyRSA::exp(&a, &na);
    assert!(ra == 8);
    let b = 2;
    let nb = From::<u64>::from(128);
    let rb = DummyRSA::exp(&b, &nb);
    assert!(rb == 3_902_709_244_137_443_127);
  }
}