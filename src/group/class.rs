//! Class Group implementation (draft)
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::util::{int, Singleton};
use rug::Integer;

#[derive(Debug, PartialEq, Eq)]
pub enum ClassGroup {}

const CLASS_GROUP_DISCRIMINANT_: i32 = 7;

lazy_static! {
  pub static ref CLASS_GROUP_DISCRIMINANT: Integer = int(CLASS_GROUP_DISCRIMINANT_);
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ClassElem {
  a: Integer,
  b: Integer,
  c: Integer
}

impl Singleton for ClassGroup {
  type Rep = Integer;
  fn rep() -> &'static Self::Rep {
    &CLASS_GROUP_DISCRIMINANT
  }
}

impl Group for ClassGroup {
  type Elem = ClassElem;
  fn op_(discriminant: &Integer, a: &ClassElem, b: &ClassElem) -> ClassElem {
    unimplemented!();
  }

  fn id_(_: &Integer) -> ClassElem {
    unimplemented!();
  }

  fn inv_(discriminant: &Integer, x: &ClassElem) -> ClassElem {
    unimplemented!();
  }

  fn exp_(discriminant: &Integer, x: &ClassElem, n: &Integer) -> ClassElem {
    unimplemented!();
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(_: &Integer) -> ClassElem {
    unimplemented!();
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_init() {
    assert!(false);
  }

  #[test]
  fn test_op() {
    assert!(false);
  }

  #[test]
  fn test_id() {
    assert!(false);
  }

  #[test]
  fn test_inv() {
    assert!(false);
  }

  #[test]
  fn test_exp() {
    assert!(false);
  }
}
