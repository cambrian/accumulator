//! Dummy RSA group for 64-bit numbers
//! Use this group for testing while we figure out ring integration

use std::u64;
use super::{Group, InvertibleGroup};

// you might need to add something to support inversion
// (e.g. p or q for use with euler's theorem)
pub struct DummyRSA {
  modulus: u64
}

const P: u64 = 2_413_575_613;
const Q: u64 = 1_090_574_917;

const DUMMY_RSA: DummyRSA = DummyRSA { modulus: 2_632_185_023_820_699_121};

#[derive(PartialEq, Eq, Clone, Serialize)]
struct DummyRSAElem (pub u64);

impl Group for DummyRSA {
  type Elem = u64;
  fn get() -> Self {
    DUMMY_RSA
  }
  fn op_(&self, a: &u64, b: &u64) -> u64{
    (a * b) % self.modulus
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
