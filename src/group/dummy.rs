//! Dummy RSA group for 64-bit numbers.
//! Use this group for testing while we figure out ring integration.

use super::super::util;
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

#[derive(Clone, Serialize)]
pub struct DummyRSAElem {
  val: u64,
}

impl DummyRSAElem {
  pub fn of(val: u64) -> DummyRSAElem {
    if val > DUMMY_RSA.modulus / 2 {
      DummyRSAElem {
        val: util::mod_euc_big(&-BigInt::from(val), &DUMMY_RSA.modulus)
          .to_u64()
          .expect("positive BigInt expected"),
      }
    } else {
      DummyRSAElem { val }
    }
  }
}

impl PartialEq for DummyRSAElem {
  fn eq(&self, other: &DummyRSAElem) -> bool {
    self.val == other.val
  }
}

impl Eq for DummyRSAElem {}

impl Group for DummyRSA {
  type Elem = DummyRSAElem;
  fn get() -> &'static Self {
    &DUMMY_RSA
  }
  fn op_(&self, a_elem: &DummyRSAElem, b_elem: &DummyRSAElem) -> DummyRSAElem {
    let (a, b) = (a_elem.val, b_elem.val);
    let op_result = ((u128::from(a) * u128::from(b)) % u128::from(self.modulus)) as u64;
    DummyRSAElem::of(op_result)
  }
  fn id_(&self) -> DummyRSAElem {
    DummyRSAElem::of(1)
  }
  fn base_elem_(&self) -> DummyRSAElem {
    DummyRSAElem::of(2)
  }
}

/// Trait for groups that support efficient inverse calculations.
/// NOT used to mean a cyclic group (where every element has an inverse).
impl InvertibleGroup for DummyRSA {
  fn inv_(&self, x: &DummyRSAElem) -> DummyRSAElem {
    let x_big = BigUint::from(x.val);
    let mod_big = BigUint::from(self.modulus);
    let (a, _, gcd) = util::bezout(&x_big, &mod_big);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    DummyRSAElem::of(
      util::mod_euc_big(&a, &self.modulus)
        .to_u64()
        .expect("u64-sized BigInt expected"),
    )
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_op() {
    let a = DummyRSA::op(&DummyRSAElem::of(2), &DummyRSAElem::of(3));
    assert!(a == DummyRSAElem::of(6));
    let b = DummyRSA::op(&DummyRSAElem::of(P + 1), &DummyRSAElem::of(Q + 1));
    assert!(b == DummyRSAElem::of(P + Q + 1));
  }

  /// Tests that -x and x are treated as the same element.
  #[test]
  fn test_cosets() {
    assert!(DummyRSAElem::of(2) == DummyRSAElem::of(2_794_712_452_613_795));
    let r = DummyRSA::op(&DummyRSAElem::of(931_570_817_537_932), &DummyRSAElem::of(2));
    assert!(r == DummyRSAElem::of(931_570_817_537_933));
  }

  #[test]
  fn test_exp() {
    let a = DummyRSA::exp(&DummyRSAElem::of(2), &From::<u64>::from(3));
    assert!(a == DummyRSAElem::of(8));
    let b = DummyRSA::exp(&DummyRSAElem::of(2), &From::<u64>::from(128));
    assert!(b == DummyRSAElem::of(782_144_413_693_680));
  }

  #[test]
  fn test_inv() {
    let r = DummyRSA::inv(&DummyRSAElem::of(2));
    assert!(r == DummyRSAElem::of(1_397_356_226_306_899));
    let r = DummyRSA::inv(&DummyRSAElem::of(32_416_188_490));
    assert!(r == DummyRSAElem::of(173_039_603_491_119));
  }
}
