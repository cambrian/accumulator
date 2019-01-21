//! Dummy RSA group for 64-bit numbers.
//! Use this group for testing while we figure out ring integration.

use super::super::util;
use super::super::util::ConvertBytes;
use super::{Group, InvertibleGroup};
use num::{BigInt, BigUint};
use num_traits::cast::ToPrimitive;
use num_traits::identities::One;
use std::u64;

pub struct DummyRSA {
  modulus: u64,
}

const P: u64 = 226_022_213;
const Q: u64 = 12_364_769;

const DUMMY_RSA: DummyRSA = DummyRSA { modulus: P * Q };

#[derive(Clone, Eq, PartialEq, Serialize)]
pub struct DummyRSAElem {
  val: u64,
}

impl DummyRSA {
  pub fn elem_of(val_unbounded: u64) -> DummyRSAElem {
    let val = val_unbounded % Self::get().modulus;
    if val > DUMMY_RSA.modulus / 2 {
      DummyRSAElem {
        val: util::mod_euc_big(&-BigInt::from(val), &Self::get().modulus)
          .to_u64()
          .expect("positive BigInt expected"),
      }
    } else {
      DummyRSAElem { val }
    }
  }
}

impl Group for DummyRSA {
  type Elem = DummyRSAElem;
  fn get() -> &'static Self {
    &DUMMY_RSA
  }
  fn op_(&self, a_elem: &DummyRSAElem, b_elem: &DummyRSAElem) -> DummyRSAElem {
    // Note: This is a pretty naive implementation of op.
    let (a, b) = (a_elem.val, b_elem.val);
    let op_result = ((u128::from(a) * u128::from(b)) % u128::from(self.modulus)) as u64;
    DummyRSA::elem_of(op_result)
  }
  fn id_(&self) -> DummyRSAElem {
    DummyRSA::elem_of(1)
  }
  fn base_elem_(&self) -> DummyRSAElem {
    DummyRSA::elem_of(2)
  }
}

impl InvertibleGroup for DummyRSA {
  fn inv_(&self, x: &DummyRSAElem) -> DummyRSAElem {
    let x_big = BigUint::from(x.val);
    let mod_big = BigUint::from(self.modulus);
    let (a, _, gcd) = util::bezout(&x_big, &mod_big);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    DummyRSA::elem_of(
      util::mod_euc_big(&a, &self.modulus)
        .to_u64()
        .expect("u64-sized BigInt expected"),
    )
  }
}

impl ConvertBytes for DummyRSA {
  fn to_le_bytes(x: &DummyRSAElem) -> Vec<u8> {
    x.val.to_le_bytes().to_vec()
  }

  fn to_be_bytes(x: &DummyRSAElem) -> Vec<u8> {
    x.val.to_be_bytes().to_vec()
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_op() {
    let a = DummyRSA::op(&DummyRSA::elem_of(2), &DummyRSA::elem_of(3));
    assert!(a == DummyRSA::elem_of(6));
    let b = DummyRSA::op(&DummyRSA::elem_of(P + 1), &DummyRSA::elem_of(Q + 1));
    assert!(b == DummyRSA::elem_of(P + Q + 1));
  }

  /// Tests that -x and x are treated as the same element.
  #[test]
  fn test_cosets() {
    assert!(DummyRSA::elem_of(2) == DummyRSA::elem_of(2_794_712_452_613_795));
    let r = DummyRSA::op(
      &DummyRSA::elem_of(931_570_817_537_932),
      &DummyRSA::elem_of(2),
    );
    assert!(r == DummyRSA::elem_of(931_570_817_537_933));
  }

  #[test]
  fn test_exp() {
    let a = DummyRSA::exp(&DummyRSA::elem_of(2), &From::<u64>::from(3));
    assert!(a == DummyRSA::elem_of(8));
    let b = DummyRSA::exp(&DummyRSA::elem_of(2), &From::<u64>::from(128));
    assert!(b == DummyRSA::elem_of(782_144_413_693_680));
  }

  #[test]
  fn test_inv() {
    let r = DummyRSA::inv(&DummyRSA::elem_of(2));
    assert!(r == DummyRSA::elem_of(1_397_356_226_306_899));
    let r = DummyRSA::inv(&DummyRSA::elem_of(32_416_188_490));
    assert!(r == DummyRSA::elem_of(173_039_603_491_119));
  }
}
