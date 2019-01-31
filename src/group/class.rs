//! Class Group implementation (draft)
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::util::{int, Singleton};
use rug::{Assign, Integer};
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub enum ClassGroup {}

// 2048-bit prime, negated.  Generated using OpenSSL.
//
// According to "A Survey of IQ Cryptography" (Buchmann & Hamdy) Table 1,
// IQ-MPQS for computing discrete logarithms in class groups with a 2048-bit discriminant is
// comparable in complexity to GNFS for factoring a 4096-bit integer.
const DISCRIMINANT2048_DECIMAL: &str = "-3061606903480752394709365751632081521549287\
                                        63761650679027169886578024000373319144482182\
                                        51590830110218951921584943041318477665819248\
                                        19762767207780092618088326303048417113668721\
                                        61223643645001916696949342349722487050631171\
                                        04912335573294798164577233813687887340799331\
                                        65653042145718668727765268057567320767851636\
                                        96501234808269893879755485983099594863614250\
                                        21860161020248607833276306314923730985457097\
                                        27023505674117797343725737548405701383103177\
                                        54359137013512655926325773048926718050691092\
                                        94533717273440872863614264045883351603859982\
                                        80988603297435639020911295652025967761702701\
                                        7014711623966286152805654229445219531956098223";

lazy_static! {
  pub static ref CLASS_GROUP_DISCRIMINANT: Integer =
    Integer::from_str(DISCRIMINANT2048_DECIMAL).unwrap();
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ClassElem {
  a: Integer,
  b: Integer,
  c: Integer,
}

impl ClassElem {
  fn normalize(&mut self) {
    if self.b > Integer::from(-&self.a) && self.b <= self.a {
      return;
    }
    // r = (a - b) / 2a, floor division
    // a = a
    // b = b + 2ra
    // c = ar^2 + br + c
    let aa = Integer::from(2 * &self.a);
    let r = Integer::from(&self.a - &self.b) / &aa;

    self.c = Integer::from(&self.a * &r) * &r + &self.b * &r + &self.c;
    self.b = &self.b + 2 * r * &self.a;
  }

  fn reduce(&mut self) {
    self.normalize();
    while self.a > self.c || (self.a == self.c && self.b < 0) {
      // s = (c + b) / 2c, floor division
      // a = c
      // b = -b + 2sc
      // c = cs^2 - bs + a

      let cc = Integer::from(2 * &self.c);
      let s = Integer::from(&self.c + &self.b) / &cc;
      let old_c = Integer::from(&self.c);

      self.c = Integer::from(&self.c * &s) * &s - &self.b * &s + &self.a;
      self.a.assign(&old_c);
      self.b = Integer::from(-&self.b) + 2 * Integer::from(&s * &old_c);
    }
    self.normalize();
  }
}

pub fn solve_linear_congruence(a: &Integer, b: &Integer, m: &Integer) -> (Integer, Integer) {
  // g = gcd(a, m), da + em = g
  let (g, d, _) = a.clone().gcd_cofactors(m.clone(), Integer::new());
  // q = b/g, r = b % g
  // mu = (q * d) % m
  // v = m / g

  let q = b / g.clone();
  let mu = (q * d) % m;
  let v = m / g;
  (mu, v)
}

#[inline]
pub fn three_gcd(a: &Integer, b: &Integer, c: &Integer) -> Integer {
  a.clone().gcd(&b).gcd(&c)
}

impl Singleton for ClassGroup {
  type Rep = Integer;
  fn rep() -> &'static Self::Rep {
    &CLASS_GROUP_DISCRIMINANT
  }
}

impl Group for ClassGroup {
  type Elem = ClassElem;
  fn op_(_: &Integer, x: &ClassElem, y: &ClassElem) -> ClassElem {
    // g = (b1 + b2) / 2
    let g = (Integer::from(&x.b) + &y.b) / 2;
    // h = (b2 - b1) / 2
    let h = (&y.b - Integer::from(&x.b)) / 2;
    // w = gcd(a1, a2, g)
    let w = three_gcd(&x.a, &y.a, &g);
    // j = w
    let j = Integer::from(&w);
    // r = 0
    let r = Integer::from(0);
    // s = a1 / w
    let s = Integer::from(&x.a / &w);
    // t = a2 / w
    let t = Integer::from(&y.a / &w);
    // u = g / w
    let u = Integer::from(&g / &w);
    // a = tu
    let a = Integer::from(&t * &u);
    // b = hu + sc1
    let b = Integer::from(&h * &u) + (&s * &x.c);
    // m = st
    let mut m = Integer::from(&s * &t);

    let (mu, v) = solve_linear_congruence(&a, &b, &m);
    // a = tv
    let a = Integer::from(&t * &v);
    // b = h - t * mu
    let b = &h - Integer::from(&t * &mu);
    // m = s
    m.assign(&s);

    let (lambda, _) = solve_linear_congruence(&a, &b, &m);
    // k = mu + v*lambda
    let k = &mu + Integer::from(&v * &lambda);
    // l = (k*t - h) / s
    let l = (Integer::from(&k * &t) - &h) / &s;
    // m = (tuk - hu - c1s) / st
    let m = (Integer::from(&t * &u) * &k - &h * &u - &x.c * &s) / Integer::from(&s * &t);
    // A = st - ru
    let a = Integer::from(&s * &t) - Integer::from(&r * &u);
    // B = ju + mr - (kt + ls)
    let b = Integer::from(&j * &u) + Integer::from(&m * &r)
      - (Integer::from(&k * &t) + Integer::from(&l * &s));
    // C = kl - jm
    let c = Integer::from(&k * &l) - Integer::from(&j * &m);
    let mut ret = ClassElem { a, b, c };
    ret.reduce();
    ret
  }

  fn id_(d: &Integer) -> ClassElem {
    let a = Integer::from(1);
    let b = Integer::from(1);
    // c = (b*b - d) / 4a
    let c = Integer::from(-d) / 4;
    ClassElem { a, b, c }
  }

  fn inv_(_: &Integer, x: &ClassElem) -> ClassElem {
    ClassElem {
      a: Integer::from(&x.a),
      b: Integer::from(-(&x.b)),
      c: Integer::from(&x.c),
    }
  }

  fn exp_(_: &Integer, x: &ClassElem, n: &Integer) -> ClassElem {
    unimplemented!();
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Integer) -> ClassElem {
    let a = Integer::from(2);
    let b = Integer::from(1);
    // c = (b*b - d) / 4a
    let c = Integer::from(1 - d) / Integer::from(8);
    let mut ret = ClassElem { a, b, c };
    ret.reduce();
    ret
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
    let g1 = ClassGroup::unknown_order_elem();
    let g2 = ClassGroup::unknown_order_elem();
    let x = ClassGroup::op(&g1, &g2);
    dbg!(x);
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
