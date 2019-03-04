//! Class Group implementation
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::mpz::Mpz;
use crate::util::{int, log_2, TypeRep};
use rug::Integer;
use std::cell::RefCell;
use std::cmp::{max, min};
use std::fmt;
use std::fmt::Display;
use std::fmt::Formatter;
use std::hash::{Hash, Hasher};
use std::str::FromStr;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ClassGroup {}

const EXP_THRESH: i64 = 31;
const THRESH: i64 = ((1 as u64) << 31) as i64;

// REVIEW: It appears that this group implementation always expects its elements to be reduced.
// However, the type signature fn reduce(&mut self) implies that there is a valid representation
// for un-reduced elements. This is suboptimal. Can you restructure the code in such a way that
// ClassElem instances are only created for valid elements? Likewise for normalize.
// Hint: change reduce to look something like
// fn reduce(a: Integer, b: Integer, c: Integer) -> ClassElem { ... }

// 2048-bit prime, negated, congruent to 3 mod 4.  Generated using OpenSSL.

// According to "A Survey of IQ Cryptography" (Buchmann & Hamdy) Table 1,
// IQ-MPQS for computing discrete logarithms in class groups with a 2048-bit discriminant is
// comparable in complexity to GNFS for factoring a 4096-bit integer.
const DISCRIMINANT2048_DECIMAL: &str =
  "-30616069034807523947093657516320815215492876376165067902716988657802400037331914448218251590830\
  1102189519215849430413184776658192481976276720778009261808832630304841711366872161223643645001916\
  6969493423497224870506311710491233557329479816457723381368788734079933165653042145718668727765268\
  0575673207678516369650123480826989387975548598309959486361425021860161020248607833276306314923730\
  9854570972702350567411779734372573754840570138310317754359137013512655926325773048926718050691092\
  9453371727344087286361426404588335160385998280988603297435639020911295652025967761702701701471162\
  3966286152805654229445219531956098223";

lazy_static! {
  pub static ref CLASS_GROUP_DISCRIMINANT: Mpz = Mpz::from_str(DISCRIMINANT2048_DECIMAL).unwrap();
}

thread_local! {
  static CTX: RefCell<ClassCtx> = Default::default();
}

#[derive(Debug)]
pub struct ClassElem {
  a: Mpz,
  b: Mpz,
  c: Mpz,
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

impl Display for ClassElem {
  fn fmt(&self, f: &mut Formatter) -> fmt::Result {
    write!(f, "\na: {}\nb: {}\nc: {}", self.a, self.b, self.c)
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
    let mut r_self = self.clone();
    let mut r_other = other.clone();
    r_self.reduce();
    r_other.reduce();
    r_self.a == r_other.a && r_self.b == r_other.b && r_self.c == r_other.c
  }
}

impl Eq for ClassElem {}

unsafe impl Send for ClassElem {}
unsafe impl Sync for ClassElem {}

pub struct LinCongruenceCtx {
  g: Mpz,
  d: Mpz,
  e: Mpz,
  q: Mpz,
  r: Mpz,
}

impl Default for LinCongruenceCtx {
  fn default() -> Self {
    Self {
      g: Mpz::default(),
      d: Mpz::default(),
      e: Mpz::default(),
      q: Mpz::default(),
      r: Mpz::default(),
    }
  }
}

impl LinCongruenceCtx {
  pub fn solve_linear_congruence(&mut self, mu: &mut Mpz, v: &mut Mpz, a: &Mpz, b: &Mpz, m: &Mpz) {
    // g = gcd(a, m) => da + em = g
    self.g.gcd_cofactors(&mut self.d, &mut self.e, a, m);
    // q = floor_div(b, g)
    // r = b % g
    self.q.floor_div_rem(&mut self.r, b, &self.g);
    // mu = (q * d) % m
    mu.mul(&self.q, &self.d);
    mu.modulo_mut(m);
    // v = m / g
    v.floor_div(m, &self.g);
  }
}

pub struct ClassCtx {
  negative_a: Mpz,
  r: Mpz,
  denom: Mpz,
  old_b: Mpz,
  ra: Mpz,
  s: Mpz,
  x: Mpz,
  old_a: Mpz,
  g: Mpz,
  w: Mpz,
  u: Mpz,
  a: Mpz,
  b: Mpz,
  m: Mpz,
  k: Mpz,
  mu: Mpz,
  v: Mpz,
  sigma: Mpz,
  lambda: Mpz,
  h: Mpz,
  t: Mpz,
  l: Mpz,
  j: Mpz,
  sctx: LinCongruenceCtx,
}

impl Default for ClassCtx {
  fn default() -> Self {
    Self {
      negative_a: Mpz::default(),
      r: Mpz::default(),
      denom: Mpz::default(),
      old_b: Mpz::default(),
      ra: Mpz::default(),
      s: Mpz::default(),
      x: Mpz::default(),
      old_a: Mpz::default(),
      g: Mpz::default(),
      w: Mpz::default(),
      u: Mpz::default(),
      a: Mpz::default(),
      b: Mpz::default(),
      m: Mpz::default(),
      k: Mpz::default(),
      mu: Mpz::default(),
      v: Mpz::default(),
      sigma: Mpz::default(),
      lambda: Mpz::default(),
      h: Mpz::default(),
      t: Mpz::default(),
      l: Mpz::default(),
      j: Mpz::default(),
      sctx: LinCongruenceCtx::default(),
    }
  }
}

pub fn signed_shift(op: u64, shift: i64) -> u64 {
  match shift {
    x if x > 0 => op << shift,
    x if x <= -64 => 0,
    _ => op >> (-shift),
  }
}

pub fn mpz_get_si_2exp(op: &Mpz) -> (i64, i64) {
  let size = op.size();
  let last = op.limbn((size - 1) as i64);
  let lg2 = log_2(last) + 1;
  let mut exp = lg2 as i64;
  let mut ret = signed_shift(last, 63 - exp);
  if size > 1 {
    exp += ((size - 1) * 64) as i64;
    let prev = op.limbn((size - 2) as i64);
    ret += signed_shift(prev, (-1 - lg2 as i32) as i64);
  }
  if op.is_neg() {
    return (-(ret as i64), exp);
  }
  (ret as i64, exp)
}

pub fn test_reduction(x: &mut ClassElem) -> bool {
  let a_b = x.a.cmp_abs(&x.b);
  let c_b = x.c.cmp_abs(&x.b);

  if a_b < 0 || c_b < 0 {
    return false;
  }

  let a_c = x.a.cmp(&x.c);
  if a_c > 0 {
    x.a.swap(&mut x.c);
    x.b.neg_mut();
  }
  if a_c == 0 && x.b.is_neg() {
    x.b.neg_mut();
  }
  true
}

impl ClassCtx {
  fn normalize(&mut self, x: &mut ClassElem) {
    self.mu.add(&x.b, &x.c);
    self.s.mul_ui(&x.c, 2); // a2 = s

    self.denom.floor_div(&self.mu, &self.s);

    self.a.set(&x.c);

    self.s.mul_ui(&self.denom, 2);
    self.b.neg(&x.b);
    self.b.addmul(&x.c, &self.s);

    self.r.set(&x.a);
    self.r.submul(&x.b, &self.denom);
    self.denom.square_mut();
    self.r.addmul(&x.c, &self.denom);

    x.a.set(&self.a);
    x.b.set(&self.b);
    x.c.set(&self.r);
  }

  fn reduce(&mut self, x: &mut ClassElem) {
    while !test_reduction(x) {
      let (mut a, a_exp) = mpz_get_si_2exp(&x.a);
      let (mut b, b_exp) = mpz_get_si_2exp(&x.b);
      let (mut c, c_exp) = mpz_get_si_2exp(&x.c);

      let mut max_exp = a_exp;
      let mut min_exp = a_exp;

      max_exp = max(max_exp, b_exp);
      max_exp = max(max_exp, c_exp);
      min_exp = min(min_exp, b_exp);
      min_exp = min(min_exp, c_exp);

      println!("about to check normalize");

      if max_exp - min_exp > EXP_THRESH {
        self.normalize(x);
        continue;
      }
      //  println!("a: {}", x.a);
      //  println!("b: {}", x.b);
      //  println!("c: {}", x.c);
      max_exp += 1; // for overflow safety
      a >>= max_exp - a_exp;
      b >>= max_exp - b_exp;
      c >>= max_exp - c_exp;

      let mut u_ = 1;
      let mut v_ = 0;
      let mut w_ = 0;
      let mut y_ = 1;

      let mut u;
      let mut v;
      let mut w;
      let mut y;

      //    println!("starting do-while loop");
      loop {
        println!("start of loop");
        u = u_;
        v = v_;
        w = w_;
        y = y_;
        let delta = if b >= 0 {
          //      println!("top");
          (b + c) / (c << 1)
        } else {
          //      println!("bottom");
          -(-b + c) / (c << 1)
        };
        let a_ = c;
        let mut c_ = c * delta;
        let b_ = -b + (c_ << 1);
        let gamma = b - c_;
        //    println!("a: {}", a);
        //    println!("delta: {}", delta);
        //    println!("gamma: {}", gamma);
        c_ = a - delta * gamma;

        a = a_;
        b = b_;
        c = c_;

        u_ = v;
        v_ = -u + delta * v;
        w_ = y;
        y_ = -w + delta * y;
        if !((v_.abs() | y_.abs()) <= THRESH && a > c && c > 0) {
          break;
        }
      }
      println!("finished loop");
      if (v_.abs() | y_.abs()) <= THRESH {
        u = u_;
        v = v_;
        w = w_;
        y = y_;
      }
      let aa = u * u;
      println!("aa: {}", aa);
      let ab = u * w;
      println!("ab: {}", ab);
      let ac = w * w;
      println!("ac: {}", ac);
      let ba = (u * v) << 1;
      println!("ba: {}", ba);
      let bb = u * y + v * w;
      println!("bb: {}", bb);
      let bc = (w * y) << 1;
      println!("bc: {}", bc);
      let ca = v * v;
      println!("ca: {}", ca);
      let cb = v * y;
      println!("cb: {}", cb);
      let cc = y * y;
      println!("cc: {}", cc);

      self.a.mul_si(&x.a, aa); // a = faa
      self.b.mul_si(&x.b, ab); // b = fab
      self.h.mul_si(&x.c, ac); // h = fac

      self.g.mul_si(&x.a, ba); // g = fba
      self.j.mul_si(&x.b, bb); // j = fbb
      self.k.mul_si(&x.c, bc); // k = fbc

      self.s.mul_si(&x.a, ca); // s = fca
      self.w.mul_si(&x.b, cb); // w = fcb
      self.l.mul_si(&x.c, cc); // l = fcc

      x.a.add(&self.a, &self.b);
      x.a.add_mut(&self.h);

      x.b.add(&self.g, &self.j);
      x.b.add_mut(&self.k);

      x.c.add(&self.s, &self.w);
      x.c.add_mut(&self.l);
    }
  }

  fn square(&mut self, x: &mut ClassElem) {
    // Solve `bk = c mod a` for k, represented by mu, v and any integer n s.t. k = mu + v * n
    self
      .sctx
      .solve_linear_congruence(&mut self.mu, &mut self.v, &x.b, &x.c, &x.a);

    // tmp = (b * mu) / a
    self.m.mul(&x.b, &self.mu);
    self.m.sub_mut(&x.c);
    self.m.floor_div_mut(&x.a);

    // A = a^2
    self.old_a.set(&x.a);
    x.a.square_mut();

    // B = b - 2a * mu
    self.a.mul(&self.mu, &self.old_a);
    self.a.mul_ui_mut(2);
    x.b.sub_mut(&self.a);

    // C = mu^2 - tmp
    x.c.mul(&self.mu, &self.mu);
    x.c.sub_mut(&self.m);

    self.reduce(x);
  }

  // Expects discriminants of x and y to be equal
  fn op(&mut self, x: &ClassElem, y: &ClassElem) -> ClassElem {
    // g = (b1 + b2) / 2
    self.g.add(&x.b, &y.b);
    self.g.floor_div_ui_mut(2);
    // h = (b2 - b1) / 2
    self.h.sub(&y.b, &x.b);
    self.h.floor_div_ui_mut(2);
    // w = gcd(a1, a2, g)
    self.w.gcd(&x.a, &y.a);
    self.w.gcd_mut(&self.g);
    // j = w
    self.j.set(&self.w);
    // r = 0
    self.r.set_ui(0);
    // s = a1 / w
    self.s.floor_div(&x.a, &self.w);
    // t = a2 / w
    self.t.floor_div(&y.a, &self.w);
    // u = g / w
    self.u.floor_div(&self.g, &self.w);
    // a = tu
    self.a.mul(&self.t, &self.u);
    // b = hu + sc
    self.b.mul(&self.h, &self.u);
    self.m.mul(&self.s, &x.c);
    self.b.add_mut(&self.m);
    // m = st
    self.m.mul(&self.s, &self.t);
    // Solve linear congruence `(tu)k = hu + sc mod st` or `ak = b mod m` for solutions k.
    self
      .sctx
      .solve_linear_congruence(&mut self.mu, &mut self.v, &self.a, &self.b, &self.m);

    // a = tv
    self.a.mul(&self.t, &self.v);
    // b = h - t * mu
    self.m.mul(&self.t, &self.mu);
    self.b.sub(&self.h, &self.m);
    // m = s
    self.m.set(&self.s);
    // Solve linear congruence `(tv)k = h - t * mu mod s` or `ak = b mod m` for solutions k
    self
      .sctx
      .solve_linear_congruence(&mut self.lambda, &mut self.sigma, &self.a, &self.b, &self.m);

    // k = mu + v * lambda
    self.a.mul(&self.v, &self.lambda);
    self.k.add(&self.mu, &self.a);
    // l = (k * t - h) / s
    self.l.mul(&self.k, &self.t);
    self.l.sub_mut(&self.h);
    self.l.floor_div_mut(&self.s);
    // m = (tuk - hu - c1s) / st
    self.m.mul(&self.t, &self.u);
    self.m.mul_mut(&self.k);
    self.a.mul(&self.h, &self.u);
    self.m.sub_mut(&self.a);
    self.a.mul(&x.c, &self.s);
    self.m.sub_mut(&self.a);
    self.a.mul(&self.s, &self.t);
    self.m.floor_div_mut(&self.a);

    let mut ret = ClassElem::default();
    // A = st - ru
    ret.a.mul(&self.s, &self.t);
    self.a.mul(&self.r, &self.u);
    ret.a.sub_mut(&self.a);

    // B = ju - kt + ls
    ret.b.mul(&self.j, &self.u);
    self.a.mul(&self.m, &self.r);
    ret.b.add_mut(&self.a);
    self.a.mul(&self.k, &self.t);
    ret.b.sub_mut(&self.a);
    self.a.mul(&self.l, &self.s);
    ret.b.sub_mut(&self.a);

    // C = kl - jm
    ret.c.mul(&self.k, &self.l);
    self.a.mul(&self.j, &self.m);
    ret.c.sub_mut(&self.a);
    self.reduce(&mut ret);
    ret
  }

  fn id(&mut self, d: &Mpz) -> ClassElem {
    let mut ret = ClassElem::default();
    ret.a.set_ui(1);
    ret.b.set_ui(1);
    // c = (b * b - d) / 4a
    self.a.sub(&ret.b, &d); // b == b*b
    ret.c.floor_div_ui(&self.a, 4);
    ret
  }
}

// ClassElem and ClassGroup ops based on Chia's fantastic doc explaining applied class groups:
//  https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf.

impl ClassElem {
  fn normalize(&mut self) {
    CTX.with(|ctx| ctx.borrow_mut().normalize(self))
  }

  fn reduce(&mut self) {
    CTX.with(|ctx| ctx.borrow_mut().reduce(self))
  }

  fn square(&mut self) {
    CTX.with(|ctx| ctx.borrow_mut().square(self))
  }

  fn discriminant(&self, d: &mut Mpz) -> Mpz {
    let mut tmp = Mpz::default();
    d.mul(&self.b, &self.b);
    tmp.mul(&self.a, &self.c);
    tmp.mul_ui_mut(4);
    d.sub_mut(&tmp);
    d.clone()
  }

  fn validate(&self) -> bool {
    let mut d = Mpz::default();
    self.discriminant(&mut d);
    d == *ClassGroup::rep()
  }

  // expects normalized element
  fn is_reduced(&self) -> bool {
    !(self.a > self.c || (self.a == self.c && self.b.cmp_si(0) < 0))
  }

  fn is_normalized(&self, ctx: &mut ClassCtx) -> bool {
    ctx.negative_a.neg(&self.a);
    self.b > ctx.negative_a && self.b <= self.a
  }

  #[inline]
  #[cfg(feature = "benchmark")]
  pub fn bench_square(&mut self) {
    self.square()
  }

  #[inline]
  #[cfg(feature = "benchmark")]
  pub fn bench_normalize(&mut self) {
    self.normalize()
  }

  #[inline]
  #[cfg(feature = "benchmark")]
  pub fn bench_reduce(&mut self) {
    self.reduce()
  }

  #[cfg(feature = "benchmark")]
  pub fn assign(&mut self, a: &Mpz, b: &Mpz, c: &Mpz) {
    self.a.set(&a);
    self.b.set(&b);
    self.c.set(&c);
  }
}

impl TypeRep for ClassGroup {
  type Rep = Mpz;
  fn rep() -> &'static Self::Rep {
    &CLASS_GROUP_DISCRIMINANT
  }
}

impl Group for ClassGroup {
  type Elem = ClassElem;

  fn op_(_: &Mpz, x: &ClassElem, y: &ClassElem) -> ClassElem {
    CTX.with(|ctx| ctx.borrow_mut().op(x, y))
  }

  fn id_(d: &Mpz) -> ClassElem {
    CTX.with(|ctx| ctx.borrow_mut().id(d))
  }

  fn inv_(_: &Mpz, x: &ClassElem) -> ClassElem {
    let mut ret = ClassElem::default();
    ret.a.set(&x.a);
    ret.b.neg(&x.b);
    ret.c.set(&x.c);
    ret
  }

  fn exp_(d: &Mpz, a: &ClassElem, n: &Integer) -> ClassElem {
    CTX.with(|ctx| {
      let mut ctx = ctx.borrow_mut();
      let (mut val, mut a, mut n) = {
        if *n < int(0) {
          (ctx.id(d), Self::inv(a), int(-n))
        } else {
          (ctx.id(d), a.clone(), n.clone())
        }
      };
      loop {
        if n == int(0) {
          return val;
        }
        if n.is_odd() {
          val = ctx.op(&val, &a);
        }
        ctx.square(&mut a);
        n >>= 1;
      }
    })
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Mpz) -> ClassElem {
    // a = 2
    // b = 1
    // c = (b * b - d) / 4a
    let mut ret = ClassElem::default();
    ret.a.set_ui(2);
    ret.b.set_ui(1);

    ret.c.set_ui(1);
    ret.c.sub_mut(&d);
    ret.c.floor_div_ui_mut(8);

    ret.reduce();
    ret
  }
}

impl Hash for ClassElem {
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.a.hash(state);
    self.b.hash(state);
    self.c.hash(state);
  }
}

impl<T> ElemFrom<(T, T, T)> for ClassGroup
where
  Mpz: From<T>,
{
  fn elem(t: (T, T, T)) -> ClassElem {
    let mut class_elem = ClassElem {
      a: Mpz::from(t.0),
      b: Mpz::from(t.1),
      c: Mpz::from(t.2),
    };

    // Ideally, this should return an error and the
    // return type of ElemFrom should be Result<Self::Elem, Self:err>,
    // but this would require a lot of ugly "unwraps" in the accumulator
    // library. Besides, users should not need to create new class group
    // elements, so an invalid ElemFrom here should signal a severe internal error.
    assert!(class_elem.validate());

    class_elem.reduce();
    class_elem
  }
}

// Caveat: tests that use "ground truth" use outputs from
//  Chia's sample implementation in python:
//    https://github.com/Chia-Network/vdf-competition/blob/master/inkfish/classgroup.py.
#[cfg(test)]
mod tests {
  use super::*;
  use std::collections::hash_map::DefaultHasher;

  // Makes a class elem tuple but does not reduce.
  fn construct_raw_elem_from_strings(a: &str, b: &str, c: &str) -> ClassElem {
    ClassElem {
      a: Mpz::from_str(a).unwrap(),
      b: Mpz::from_str(b).unwrap(),
      c: Mpz::from_str(c).unwrap(),
    }
  }

  #[should_panic]
  #[test]
  fn test_bad_elem() {
    let _ = ClassGroup::elem((1, 2, 3));
  }

  #[test]
  fn test_elem_from() {
    let a1 = Mpz::from_str("16").unwrap();
    let b1 = Mpz::from_str("105").unwrap();
    let c1 = Mpz::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",
    )
    .unwrap();

    let a2 = Mpz::from_str("16").unwrap();
    let b2 = Mpz::from_str("9").unwrap();
    let c2 = Mpz::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036",
    )
    .unwrap();

    let reduced_elem = ClassGroup::elem((a1, b1, c1));
    let also_reduced_elem = ClassGroup::elem((a2, b2, c2));
    assert_eq!(reduced_elem, also_reduced_elem);
  }
  #[test]
  fn test_equality() {
    let not_reduced = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    let diff_elem = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    assert!(not_reduced == reduced_ground_truth);
    assert!(not_reduced == not_reduced.clone());
    assert!(reduced_ground_truth == reduced_ground_truth.clone());
    assert!(not_reduced != diff_elem);
    assert!(reduced_ground_truth != diff_elem);
  }

  #[test]
  fn test_hash() {
    let mut not_reduced = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    let mut reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    let diff_elem = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    let mut hasher_lh = DefaultHasher::new();
    let mut hasher_rh = DefaultHasher::new();
    not_reduced.reduce();
    reduced_ground_truth.reduce();
    not_reduced.hash(&mut hasher_lh);
    reduced_ground_truth.hash(&mut hasher_rh);
    assert!(hasher_lh.finish() == hasher_rh.finish());
    assert!(hasher_lh.finish() == hasher_lh.finish());
    assert!(hasher_rh.finish() == hasher_rh.finish());

    hasher_lh = DefaultHasher::new();
    hasher_rh = DefaultHasher::new();
    not_reduced.hash(&mut hasher_lh);
    diff_elem.hash(&mut hasher_rh);
    assert!(hasher_lh.finish() != hasher_rh.finish());
  }

  #[test]
  fn test_reduce_basic() {
    let to_reduce = construct_raw_elem_from_strings(
      "59162244921619725812008939143220718157267937427074598447911241410131470159247784852210767449\
      675610037288729551814191198624164179866076352187405442496568188988272422133088755036699145362\
      385840772236403043664778415471196678638241785773530531198720497580622741709880533724904220122\
      358854068046553219863419609777498761804625479650772123754523807001976654588225908928022367436\
      8",
      "18760351095004839755193532164856605650590306627169248964100884295652838905828158941233738613\
      175821849253748329102319504958410190952820220503570113920576542676928659211807590199941027958\
      195895385446372444261885022800653454209101497963588809819572703579484085278913354621371362285\
      341138299691587953249270188429393417132110841259813122945626515477865766896056280729710478647\
      13",
      "14872270891432803054791175727694631095755964943358394411314110783404577714102170379700365256\
      599679049493824862742803590079461712691146098397470840896560034332315858221821103076776907123\
      277315116632337385101204055232891361405428635972040596205450316747012080794838691280547894128\
      246741601088755087359234554141346980837292342320288111397175220296098629890108459305643419353\
      36"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "26888935961824081232597112540509824504614070059776273347136888921115497522070287009841688662\
      983066376019079593372296556420848446780369918809384119124783870290778875424468497961559643807\
      918398860928578027038014112641529893817109240852544158309292025321122680747989987560029531021\
      808743313150630063377037854944",
      "14529985196481999393995154363327100184407232892559561136140792409262328867440167480822808496\
      853924547751298342980606034124112579835255733824790020119078588372593288210628255956605240171\
      744703418426092073347584357826862813733154338737148962212641444735717023402201569115323580814\
      54099903972209626147819759991",
      "28467266502267127591420289007165819749231433586093061478772560429058231137856046130384492811\
      816456933286039468940950129263300933723839212086399375780796041634531383342902918719073416087\
      614456845205980227091403964285870107268917183244016635907926846271829374679124848388403486656\
      1564478239095738726823372184204"
    );
    // Equality checking should be defined to reduce
    assert_eq!(to_reduce, reduced_ground_truth);
    let already_reduced = reduced_ground_truth.clone();
    assert_eq!(already_reduced, reduced_ground_truth);
  }

  #[test]
  fn test_normalize_basic() {
    let mut unnormalized = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",

    );

    let normalized_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "4783760786688675616733383986925127377420761933775791859799529477781625005833111632534101811\
       06720472171123774764735020601213528425753087932376215639471576300984851315174010737751911943\
       19531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526\
       85419887133231339948938699768182757831793879217091871179468485931169743972659665650159413844\
       97394942286170683296647767144847422761580905834957146491938390841109871491186151583613524884\
       88402038894799695420483272708933239751363849397287571692736881031223140446926522431859701738\
       9945629057462766047140854869124473221137588347335081555186814036",
    );

    unnormalized.normalize();
    assert_eq!(normalized_ground_truth, unnormalized);
  }

  #[test]
  fn test_discriminant_basic() {
    let g = ClassGroup::unknown_order_elem();
    let mut d = Mpz::default();
    g.discriminant(&mut d);
    assert!(d == *ClassGroup::rep())
  }

  #[test]
  fn test_discriminant_across_ops() {
    let id = ClassGroup::id();
    let g1 = ClassGroup::unknown_order_elem();
    let g2 = ClassGroup::op(&g1, &g1);
    let g3 = ClassGroup::op(&id, &g2);
    let g3_inv = ClassGroup::inv(&g3);

    assert!(id.validate());
    assert!(g1.validate());
    assert!(g2.validate());
    assert!(g3.validate());
    assert!(g3_inv.validate());
  }

  #[test]
  fn test_op_single() {
    let a = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    let b = construct_raw_elem_from_strings(
      "16",
      "41",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814061"
    );

    let ground_truth = construct_raw_elem_from_strings(
      "64",
      "9",
      "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      50972369992385512081817723330993784096234932189292318422025780578511173163060796492543474864\
      07264365691511785213717281118305284397086833770388796703509"
    );

    let x = ClassGroup::op(&a, &b);

    assert_eq!(x, ground_truth);
  }

  #[test]
  fn test_op_alternating() {
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();
    let mut g_star = ClassGroup::id();

    // g
    g = ClassGroup::op(&g_anchor, &g);

    // g^2, g^* = g^2
    g = ClassGroup::op(&g_anchor, &g);
    g_star = ClassGroup::op(&g, &g_star);

    // g^3
    g = ClassGroup::op(&g_anchor, &g);

    // g^4, g^* = g^2 * g^4 = g^6
    g = ClassGroup::op(&g_anchor, &g);
    g_star = ClassGroup::op(&g, &g_star);

    let ground_truth = construct_raw_elem_from_strings(
      "64",
      "9",
      "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      509723699923855120818177233309937840962349321892923184220257805785111731630607964925434748640\
      7264365691511785213717281118305284397086833770388796703509"
    );

    assert_eq!(ground_truth, g_star);
  }

  #[test]
  fn test_op_complex() {
    // 1. Take g^100, g^200, ..., g^1000.
    // 2. Compute g^* = g^100 * ... * g^1000
    // 3. For each of g^100, g^200, ..., g^1000
    //    compute the inverse of that element and
    //    assert that g^* * current_inverse = product of g^100, g^200, ..., g^1000
    //    without the inversed-out element.
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    let mut gs = vec![];
    let mut gs_invs = vec![];

    let mut g_star = ClassGroup::id();
    for i in 1..=1000 {
      g = ClassGroup::op(&g_anchor, &g);
      assert!(g.validate());
      if i % 100 == 0 {
        gs.push(g.clone());
        gs_invs.push(ClassGroup::inv(&g));
        g_star = ClassGroup::op(&g, &g_star);
        println!("{}", i);
        assert!(g_star.validate());
      }
    }

    let elems_n_invs = gs.iter().zip(gs_invs.iter());
    for (g_elem, g_inv) in elems_n_invs {
      assert!(g_elem.validate());
      assert!(g_inv.validate());
      let mut curr_prod = ClassGroup::id();
      for elem in &gs {
        if elem != g_elem {
          curr_prod = ClassGroup::op(&curr_prod, &elem);
          assert!(curr_prod.validate());
        }
      }
      assert_eq!(ClassGroup::id(), ClassGroup::op(&g_inv, &g_elem));
      assert_eq!(curr_prod, ClassGroup::op(&g_inv, &g_star));
    }
  }

  #[test]
  fn test_id_basic() {
    let g = ClassGroup::unknown_order_elem();
    let id = ClassGroup::id();
    assert_eq!(g, ClassGroup::op(&g, &id));
    assert_eq!(g, ClassGroup::op(&id, &g));
    assert_eq!(id, ClassGroup::op(&id, &id));
  }

  #[test]
  fn test_id_repeated() {
    let mut id = ClassGroup::id();
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();
    for _ in 0..1000 {
      id = ClassGroup::op(&id, &id);
      assert_eq!(id, ClassGroup::id());
      g = ClassGroup::op(&g, &ClassGroup::id());
      assert_eq!(g, g_anchor);
    }
  }

  #[test]
  fn test_inv() {
    let id = ClassGroup::id();
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();

    for i in 0..1000 {
      g = ClassGroup::op(&g, &g_anchor);
      let g_inv = ClassGroup::inv(&g);
      //  println!("i: {}", i);
      println!("g_inv: {}", g_inv);
      let x1 = ClassGroup::op(&g_inv, &g);
      let x2 = ClassGroup::op(&g, &g_inv);
      println!("g_inv * g: {}", x1);
      println!("g * g_inv: {}", x2);
      println!("id: {}", id);
      println!("g: {}", g);
      assert_eq!(id, x1);
      assert_eq!(id, x2);
      assert_eq!(g, ClassGroup::inv(&g_inv));
    }
  }

  #[test]
  fn test_exp_basic() {
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    for i in 1..=1000 {
      g = ClassGroup::op(&g, &g_anchor);
      assert_eq!(&g, &ClassGroup::exp(&g_anchor, &Integer::from(i)));
    }
  }

  #[test]
  fn test_square_basic() {
    let g = ClassGroup::unknown_order_elem();
    let mut g4 = ClassGroup::id();

    // g^4
    for _ in 0..4 {
      g4 = ClassGroup::op(&g, &g4);
    }

    // g^2
    let mut g2 = g.clone();
    // g^4
    g2.square();
    g2.square();

    assert_eq!(&g2, &g4);
  }
}
