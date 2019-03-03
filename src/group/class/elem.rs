//! Defines the ClassElem struct and associated traits.
use crate::num::mpz::Mpz;
use std::hash::{Hash, Hasher};

#[allow(clippy::stutter)]
#[derive(Debug)]
pub struct ClassElem {
  pub a: Mpz, // TODO: CAN WE MAKE THIS PRIVATE?
  pub b: Mpz,
  pub c: Mpz,
}

impl Default for ClassElem {
  fn default() -> Self {
    ClassElem {
      a: Mpz::default(),
      b: Mpz::default(),
      c: Mpz::default(),
    }
  }
}

impl Clone for ClassElem {
  fn clone(&self) -> Self {
    let mut ret = ClassElem::default();
    ret.a = self.a.clone();
    ret.b = self.b.clone();
    ret.c = self.c.clone();
    ret
  }
}

impl PartialEq for ClassElem {
  fn eq(&self, other: &ClassElem) -> bool {
    // ClassElems only ever exist in reduced form, unless manually
    // instantiated with ClassElem {...}
    self.a == other.a && self.b == other.b && self.c == other.c
  }
}

impl Hash for ClassElem {
  // Assumes `ClassElem` is reduced and normalized, which will be the case unless a struct is
  // instantiated manually in this module.
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.a.hash(state);
    self.b.hash(state);
    self.c.hash(state);
  }
}

impl Eq for ClassElem {}
unsafe impl Send for ClassElem {}
unsafe impl Sync for ClassElem {}
