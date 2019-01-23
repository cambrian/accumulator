//! Dummy RSA group for 64-bit numbers.
//! Use this group for testing while we figure out ring integration.

use super::{Group, InvertibleGroup};
use crate::util;
use crate::util::{bi, bu, ConvertBytes, Singleton};
use num_traits::cast::ToPrimitive;
use num_traits::identities::One;
use std::u64;

#[derive(PartialEq, Eq)]
pub enum DummyRSA {}

pub struct DummyRSAModulus {
  modulus: u64,
}

const P: u64 = 226_022_213;
const Q: u64 = 12_364_769;

const DUMMY_RSA_MODULUS: DummyRSAModulus = DummyRSAModulus { modulus: P * Q };

#[derive(Clone, Eq, PartialEq, Serialize)]
pub struct DummyRSAElem {
  val: u64,
}

impl Singleton for DummyRSA {
  type Rep = DummyRSAModulus;
  fn rep() -> &'static DummyRSAModulus {
    &DUMMY_RSA_MODULUS
  }
}

impl DummyRSA {
  pub fn elem_of(val_unbounded: u64) -> DummyRSAElem {
    let n = Self::rep().modulus;
    let val = val_unbounded % n;
    if val > n / 2 {
      DummyRSAElem {
        val: util::mod_euc_big(&-bi(val), &n)
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
  fn op_(
    DummyRSAModulus { modulus }: &DummyRSAModulus,
    a_elem: &DummyRSAElem,
    b_elem: &DummyRSAElem,
  ) -> DummyRSAElem {
    // Note: This is a pretty naive implementation of op.
    let (a, b) = (a_elem.val, b_elem.val);
    let op_result = ((u128::from(a) * u128::from(b)) % u128::from(*modulus)) as u64;
    DummyRSA::elem_of(op_result)
  }
  fn id_(_: &DummyRSAModulus) -> DummyRSAElem {
    DummyRSA::elem_of(1)
  }
  fn base_elem_(_: &DummyRSAModulus) -> DummyRSAElem {
    DummyRSA::elem_of(2)
  }
}

impl InvertibleGroup for DummyRSA {
  fn inv_(DummyRSAModulus { modulus }: &DummyRSAModulus, x: &DummyRSAElem) -> DummyRSAElem {
    let x_big = bu(x.val);
    let mod_big = bu(*modulus);
    let (a, _, gcd) = util::bezout(&x_big, &mod_big);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    DummyRSA::elem_of(
      util::mod_euc_big(&a, modulus)
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
