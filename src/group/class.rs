//! Class Group implementation, with optimizations.
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::num::flint;
use crate::num::flint::fmpz;
use crate::num::mpz::{flint_mpz_struct, Mpz};
use crate::util::{int, LinCongruenceCtx, TypeRep};
use rug::Integer;
use std::cell::RefCell;
use std::hash::{Hash, Hasher};
use std::str::FromStr;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ClassGroup {}

// 2048-bit prime, negated, congruent to 3 mod 4.  Generated using OpenSSL.
// According to "A Survey of IQ Cryptography" (Buchmann & Hamdy) Table 1, IQ-MPQS for computing
// discrete logarithms in class groups with a 2048-bit discriminant is comparable in complexity to
// GNFS for factoring a 4096-bit integer.
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

// ClassElem and ClassGroup ops based on Chia's fantastic doc explaining applied class groups:
//  https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf,
//  hereafter refered to as "Binary Quadratic Forms"

#[allow(clippy::stutter)]
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

impl Eq for ClassElem {}
unsafe impl Send for ClassElem {}
unsafe impl Sync for ClassElem {}

#[allow(non_snake_case)]
pub struct ClassCtx {
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
  G_sq_op: Mpz,
  x_sq_op: Mpz,
  y_sq_op: Mpz,
  z_sq_op: Mpz,
  By_sq_op: Mpz,
  Dy_sq_op: Mpz,
  bx_sq_op: Mpz,
  by_sq_op: Mpz,
  L_sq_op: Mpz,
  dx_sq_op: Mpz,
  dy_sq_op: Mpz,
  t_sq_op: Mpz,
  q_sq_op: Mpz,
  ax_sq_op: Mpz,
  ay_sq_op: Mpz,
  Q1_sq_op: Mpz,
  scratch: Mpz,
  sctx: LinCongruenceCtx,
}

impl Default for ClassCtx {
  fn default() -> Self {
    let mut s = Self {
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

      // Used in squaring ops
      G_sq_op: Mpz::default(),
      x_sq_op: Mpz::default(),
      y_sq_op: Mpz::default(),
      z_sq_op: Mpz::default(),
      By_sq_op: Mpz::default(),
      Dy_sq_op: Mpz::default(),
      bx_sq_op: Mpz::default(),
      by_sq_op: Mpz::default(),
      L_sq_op: Mpz::default(),
      dx_sq_op: Mpz::default(),
      dy_sq_op: Mpz::default(),
      t_sq_op: Mpz::default(),
      q_sq_op: Mpz::default(),
      ax_sq_op: Mpz::default(),
      ay_sq_op: Mpz::default(),
      Q1_sq_op: Mpz::default(),
      scratch: Mpz::default(),

      sctx: LinCongruenceCtx::default(),
    };

    // Precomputation needed for NUDULP squaring.
    s.L_sq_op.abs(&CLASS_GROUP_DISCRIMINANT.clone());
    s.L_sq_op.root_mut(4);
    s
  }
}

impl ClassCtx {
  fn normalize_(&mut self, a: &mut Mpz, b: &mut Mpz, c: &mut Mpz) {
    // Binary Quadratic Forms, 5.1.1
    self.r.sub(&a, &b);
    self.denom.mul_ui(&a, 2);
    self.r.fdiv_q_mut(&self.denom);

    self.old_b.set(&b);

    self.ra.mul(&self.r, &a);
    b.add_mut(&self.ra);
    b.add_mut(&self.ra);

    self.ra.mul_mut(&self.r);
    c.add_mut(&self.ra);

    self.ra.set(&self.r);
    self.ra.mul_mut(&self.old_b);
    c.add_mut(&self.ra);
  }

  fn normalize_mut(&mut self, x: &mut ClassElem) {
    if ClassGroup::is_normal(&x.a, &x.b, &x.c) {
      return;
    }
    self.normalize_(&mut x.a, &mut x.b, &mut x.c);
  }

  fn normalize(&mut self, mut a: Mpz, mut b: Mpz, mut c: Mpz) -> (Mpz, Mpz, Mpz) {
    if ClassGroup::is_normal(&a, &b, &c) {
      return (a, b, c);
    }
    self.normalize_(&mut a, &mut b, &mut c);
    (a, b, c)
  }

  fn reduce_(&mut self, elem: &mut ClassElem) {
    // Binary Quadratic Forms, 5.2.1
    while !ClassGroup::is_reduced(&elem.a, &elem.b, &elem.c) {
      self.s.add(&elem.c, &elem.b);
      self.x.mul_ui(&elem.c, 2);
      self.s.fdiv_q_mut(&self.x);
      self.old_a.set(&elem.a);
      self.old_b.set(&elem.b);
      elem.a.set(&elem.c);
      elem.b.neg_mut();
      self.x.mul(&self.s, &elem.c);
      self.x.mul_ui_mut(2);
      elem.b.add_mut(&self.x);

      elem.c.mul_mut(&self.s);
      elem.c.mul_mut(&self.s);
      self.x.mul(&self.old_b, &self.s);
      elem.c.sub_mut(&self.x);
      elem.c.add_mut(&self.old_a);
    }
  }

  fn reduce(&mut self, a: Mpz, b: Mpz, c: Mpz) -> (Mpz, Mpz, Mpz) {
    let (a, b, c) = self.normalize(a, b, c);
    let mut elem = ClassElem { a, b, c };
    self.reduce_(&mut elem);
    self.normalize(elem.a, elem.b, elem.c)
  }

  fn reduce_mut(&mut self, x: &mut ClassElem) {
    self.normalize_mut(x);
    self.reduce_(x);
    self.normalize_mut(x);
  }

  fn square(&mut self, x: &mut ClassElem) {
    // Binary Quadratic Forms, 6.3.1
    self
      .sctx
      .solve_linear_congruence(&mut self.mu, &mut self.v, &x.b, &x.c, &x.a)
      .unwrap();

    self.m.mul(&x.b, &self.mu);
    self.m.sub_mut(&x.c);
    self.m.fdiv_q_mut(&x.a);

    self.old_a.set(&x.a);
    x.a.square_mut();

    self.a.mul(&self.mu, &self.old_a);
    self.a.mul_ui_mut(2);
    x.b.sub_mut(&self.a);

    x.c.mul(&self.mu, &self.mu);
    x.c.sub_mut(&self.m);

    self.reduce_mut(x);
  }

  #[allow(non_snake_case)]
  fn square_nudulp_help_flint(ctx: &mut ClassCtx) {
    unsafe {
      let mut fy: fmpz = 0;
      let mut fx: fmpz = 0;
      let mut fby: fmpz = 0;
      let mut fbx: fmpz = 0;
      let mut fL: fmpz = 0;

      let mut y_square_clone = flint_mpz_struct::from(&ctx.y_sq_op);
      let mut x_square_clone = flint_mpz_struct::from(&ctx.x_sq_op);
      let mut by_square_clone = flint_mpz_struct::from(&ctx.by_sq_op);
      let mut bx_square_clone = flint_mpz_struct::from(&ctx.bx_sq_op);
      let mut L_square_clone = flint_mpz_struct::from(&ctx.L_sq_op);

      flint::fmpz_set_mpz(&mut fy, &mut y_square_clone);
      flint::fmpz_set_mpz(&mut fx, &mut x_square_clone);
      flint::fmpz_set_mpz(&mut fby, &mut by_square_clone);
      flint::fmpz_set_mpz(&mut fbx, &mut bx_square_clone);
      flint::fmpz_set_mpz(&mut fL, &mut L_square_clone);

      // Flint Lehmer partial extended GCD.
      flint::fmpz_xgcd_partial(&mut fy, &mut fx, &mut fby, &mut fbx, &mut fL);

      flint::fmpz_get_mpz(&mut y_square_clone, &mut fy);
      flint::fmpz_get_mpz(&mut x_square_clone, &mut fx);
      flint::fmpz_get_mpz(&mut by_square_clone, &mut fby);
      flint::fmpz_get_mpz(&mut bx_square_clone, &mut fbx);

      ctx.y_sq_op = Mpz::from(y_square_clone);
      ctx.x_sq_op = Mpz::from(x_square_clone);
      ctx.by_sq_op = Mpz::from(by_square_clone);
      ctx.bx_sq_op = Mpz::from(bx_square_clone);
    }

    ctx.x_sq_op.neg_mut();
    if ctx.x_sq_op.sgn() > 0 {
      ctx.y_sq_op.neg_mut();
    } else {
      ctx.by_sq_op.neg_mut();
    }
  }

  fn square_nudulp_help(ctx: &mut ClassCtx) {
    ctx.x_sq_op.set_ui(1);
    ctx.y_sq_op.set_ui(0);
    while ctx.by_sq_op.cmpabs(&ctx.L_sq_op) > 0 && ctx.bx_sq_op.sgn() == 0 {
      ctx.q_sq_op.fdiv_q(&ctx.by_sq_op, &ctx.bx_sq_op);
      ctx.t_sq_op.fdiv_r(&ctx.by_sq_op, &ctx.bx_sq_op);

      ctx.by_sq_op.set(&ctx.bx_sq_op);
      ctx.bx_sq_op.set(&ctx.t_sq_op);
      ctx.y_sq_op.submul(&ctx.q_sq_op, &ctx.x_sq_op);
      ctx.t_sq_op.set(&ctx.y_sq_op);
      ctx.y_sq_op.set(&ctx.x_sq_op);
      ctx.x_sq_op.set(&ctx.t_sq_op);
      ctx.z_sq_op.add_ui_mut(1);
    }

    if ctx.z_sq_op.odd() != 0 {
      ctx.by_sq_op.neg_mut();
      ctx.y_sq_op.neg_mut();
    }
  }

  fn square_nudulp(&mut self, x: &mut ClassElem) {
    // Jacobson, Michael J., and Alfred J. Van Der Poorten. "Computational aspects of NUCOMP."
    // Algorithm 2 (Alg 2).

    // Step 1 in Alg 2.
    self
      .G_sq_op
      .gcdext(&mut self.y_sq_op, &mut self.scratch, &x.b, &x.a);
    self.By_sq_op.divexact(&x.a, &self.G_sq_op);
    self.Dy_sq_op.divexact(&x.b, &self.G_sq_op);

    // Step 2 in Alg 2.
    self.bx_sq_op.mul(&self.y_sq_op, &x.c);
    self.bx_sq_op.modulo_mut(&self.By_sq_op);
    self.by_sq_op.set(&self.By_sq_op);

    if self.by_sq_op.cmpabs(&self.L_sq_op) <= 0 {
      // Step 4 in Alg 2.
      self.dx_sq_op.mul(&self.bx_sq_op, &self.Dy_sq_op);
      self.dx_sq_op.sub_mut(&x.c);
      self.dx_sq_op.divexact_mut(&self.By_sq_op);
      x.a.mul(&self.by_sq_op, &self.by_sq_op);
      x.c.mul(&self.bx_sq_op, &self.bx_sq_op);
      self.t_sq_op.add(&self.bx_sq_op, &self.by_sq_op);
      self.t_sq_op.square_mut();
      x.b.sub_mut(&self.t_sq_op);
      x.b.add_mut(&x.a);
      x.b.add_mut(&x.c);
      self.t_sq_op.mul(&self.G_sq_op, &self.dx_sq_op);
      x.c.sub_mut(&self.t_sq_op);
      return;
    }

    // Most of Step 3 in Alg 2.
    if cfg!(feature = "flint") {
      // Subroutine as handled by top entry to the Chia VDF competition "bulaiden."
      ClassCtx::square_nudulp_help_flint(self);
    } else {
      // Subroutine as presented in "Computational aspects of NUCOMP.", Algorithm 2, most of step 3.
      ClassCtx::square_nudulp_help(self);
    }

    self.ax_sq_op.mul(&self.G_sq_op, &self.x_sq_op);
    self.ay_sq_op.mul(&self.G_sq_op, &self.y_sq_op);

    // Step 5 in Alg 2.
    self.t_sq_op.mul(&self.Dy_sq_op, &self.bx_sq_op);
    self.t_sq_op.submul(&x.c, &self.x_sq_op);
    self.dx_sq_op.divexact(&self.t_sq_op, &self.By_sq_op);
    self.Q1_sq_op.mul(&self.y_sq_op, &self.dx_sq_op);
    self.dy_sq_op.add(&self.Q1_sq_op, &self.Dy_sq_op);
    x.b.add(&self.dy_sq_op, &self.Q1_sq_op);
    x.b.mul_mut(&self.G_sq_op);
    self.dy_sq_op.divexact_mut(&self.x_sq_op);
    x.a.mul(&self.by_sq_op, &self.by_sq_op);
    x.c.mul(&self.bx_sq_op, &self.bx_sq_op);
    self.t_sq_op.add(&self.bx_sq_op, &self.by_sq_op);
    x.b.submul(&self.t_sq_op, &self.t_sq_op);
    x.b.add_mut(&x.a);
    x.b.add_mut(&x.c);
    x.a.submul(&self.ay_sq_op, &self.dy_sq_op);
    x.c.submul(&self.ax_sq_op, &self.dx_sq_op);
    self.reduce_mut(x);
  }

  fn op(&mut self, x: &ClassElem, y: &ClassElem) -> ClassElem {
    // Binary Quadratic Forms, 6.1.1
    self.g.add(&x.b, &y.b);
    self.g.fdiv_q_ui_mut(2);
    self.h.sub(&y.b, &x.b);
    self.h.fdiv_q_ui_mut(2);
    self.w.gcd(&x.a, &y.a);
    self.w.gcd_mut(&self.g);
    self.j.set(&self.w);
    self.r.set_ui(0);
    self.s.fdiv_q(&x.a, &self.w);
    self.t.fdiv_q(&y.a, &self.w);
    self.u.fdiv_q(&self.g, &self.w);
    self.a.mul(&self.t, &self.u);
    self.b.mul(&self.h, &self.u);
    self.m.mul(&self.s, &x.c);
    self.b.add_mut(&self.m);
    self.m.mul(&self.s, &self.t);
    self
      .sctx
      .solve_linear_congruence(&mut self.mu, &mut self.v, &self.a, &self.b, &self.m)
      .unwrap();

    self.a.mul(&self.t, &self.v);
    self.m.mul(&self.t, &self.mu);
    self.b.sub(&self.h, &self.m);
    self.m.set(&self.s);
    self
      .sctx
      .solve_linear_congruence(&mut self.lambda, &mut self.sigma, &self.a, &self.b, &self.m)
      .unwrap();

    self.a.mul(&self.v, &self.lambda);
    self.k.add(&self.mu, &self.a);
    self.l.mul(&self.k, &self.t);
    self.l.sub_mut(&self.h);
    self.l.fdiv_q_mut(&self.s);
    self.m.mul(&self.t, &self.u);
    self.m.mul_mut(&self.k);
    self.a.mul(&self.h, &self.u);
    self.m.sub_mut(&self.a);
    self.a.mul(&x.c, &self.s);
    self.m.sub_mut(&self.a);
    self.a.mul(&self.s, &self.t);
    self.m.fdiv_q_mut(&self.a);

    let mut ret = ClassElem::default();

    ret.a.mul(&self.s, &self.t);
    self.a.mul(&self.r, &self.u);
    ret.a.sub_mut(&self.a);

    ret.b.mul(&self.j, &self.u);
    self.a.mul(&self.m, &self.r);
    ret.b.add_mut(&self.a);
    self.a.mul(&self.k, &self.t);
    ret.b.sub_mut(&self.a);
    self.a.mul(&self.l, &self.s);
    ret.b.sub_mut(&self.a);

    ret.c.mul(&self.k, &self.l);
    self.a.mul(&self.j, &self.m);
    ret.c.sub_mut(&self.a);
    let (a, b, c) = self.reduce(ret.a, ret.b, ret.c);
    ClassElem { a, b, c }
  }

  fn id(&mut self, d: &Mpz) -> ClassElem {
    // Binary Quadratic Forms, Definition 5.4
    let mut ret = ClassElem::default();
    ret.a.set_ui(1);
    ret.b.set_ui(1);
    self.a.sub(&ret.b, &d);
    ret.c.fdiv_q_ui(&self.a, 4);
    ret
  }
}

impl ClassElem {
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
}

impl TypeRep for ClassGroup {
  type Rep = Mpz;
  fn rep() -> &'static Self::Rep {
    &CLASS_GROUP_DISCRIMINANT
  }
}

impl ClassGroup {
  pub fn normalize(a: Mpz, b: Mpz, c: Mpz) -> (Mpz, Mpz, Mpz) {
    CTX.with(|ctx| ctx.borrow_mut().normalize(a, b, c))
  }

  pub fn reduce(a: Mpz, b: Mpz, c: Mpz) -> (Mpz, Mpz, Mpz) {
    CTX.with(|ctx| ctx.borrow_mut().reduce(a, b, c))
  }

  pub fn square(x: &mut ClassElem) {
    CTX.with(|ctx| ctx.borrow_mut().square(x))
  }

  // expects normalized element
  fn is_reduced(a: &Mpz, b: &Mpz, c: &Mpz) -> bool {
    ClassGroup::is_normal(a, b, c) && (a <= c && !(a == c && b.cmp_si(0) < 0))
  }

  fn is_normal(a: &Mpz, b: &Mpz, _c: &Mpz) -> bool {
    let mut neg_a = Mpz::default();
    neg_a.neg(a);
    neg_a < *b && b <= a
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

        if cfg!(feature = "nudulp") {
          // NUDULP optimization, maybe using FLINT bindings.
          ctx.square_nudulp(&mut a);
        } else {
          ctx.square(&mut a);
        }
        n >>= 1;
      }
    })
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Mpz) -> ClassElem {
    // Binary Quadratic Forms, Definition 5.4
    let mut ret = ClassElem::default();
    ret.a.set_ui(2);
    ret.b.set_ui(1);
    ret.c.set_ui(1);
    ret.c.sub_mut(&d);
    ret.c.fdiv_q_ui_mut(8);

    let (a, b, c) = ClassGroup::reduce(ret.a, ret.b, ret.c);
    ClassElem { a, b, c }
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

impl<A, B, C> ElemFrom<(A, B, C)> for ClassGroup
where
  Mpz: From<A>,
  Mpz: From<B>,
  Mpz: From<C>,
{
  fn elem(abc: (A, B, C)) -> ClassElem {
    let (a, b, c) = ClassGroup::reduce(Mpz::from(abc.0), Mpz::from(abc.1), Mpz::from(abc.2));
    let class_elem = ClassElem { a, b, c };

    // Ideally, this should return an error and the
    // return type of ElemFrom should be Result<Self::Elem, Self:err>,
    // but this would require a lot of ugly "unwraps" in the accumulator
    // library. Besides, users should not need to create new class group
    // elements, so an invalid ElemFrom here should signal a severe internal error.
    assert!(class_elem.validate());

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

    assert!(not_reduced != reduced_ground_truth);
    assert!(not_reduced == not_reduced.clone());
    assert!(reduced_ground_truth == reduced_ground_truth.clone());
    assert!(not_reduced != diff_elem);
    assert!(reduced_ground_truth != diff_elem);

    let not_reduced = ClassGroup::reduce(not_reduced.a, not_reduced.b, not_reduced.c);
    assert!(
      not_reduced
        == (
          reduced_ground_truth.a,
          reduced_ground_truth.b,
          reduced_ground_truth.c
        )
    );
  }

  #[test]
  fn test_hash() {
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

    let mut hasher_lh = DefaultHasher::new();
    let mut hasher_rh = DefaultHasher::new();
    not_reduced.hash(&mut hasher_lh);
    reduced_ground_truth.hash(&mut hasher_rh);

    assert!(hasher_lh.finish() != hasher_rh.finish());
    assert!(hasher_lh.finish() == hasher_lh.finish());
    assert!(hasher_rh.finish() == hasher_rh.finish());

    hasher_lh = DefaultHasher::new();
    hasher_rh = DefaultHasher::new();
    let reduced = ClassGroup::reduce(not_reduced.a, not_reduced.b, not_reduced.c);
    reduced.hash(&mut hasher_lh);
    reduced_ground_truth.hash(&mut hasher_rh);
    assert!(hasher_lh.finish() == hasher_rh.finish());

    hasher_lh = DefaultHasher::new();
    hasher_rh = DefaultHasher::new();
    reduced.hash(&mut hasher_lh);
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

    let already_reduced = reduced_ground_truth.clone();
    assert_eq!(already_reduced, reduced_ground_truth);

    assert_ne!(to_reduce, reduced_ground_truth);
    let reduced = ClassGroup::reduce(to_reduce.a, to_reduce.b, to_reduce.c);
    assert_eq!(
      reduced,
      (
        reduced_ground_truth.a,
        reduced_ground_truth.b,
        reduced_ground_truth.c
      )
    );
  }

  #[test]
  // REVIEW: this test should be restructured to not construct ClassElems but it will do for now
  fn test_normalize_basic() {
    let unnormalized = construct_raw_elem_from_strings(
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

    let elem_tuple = ClassGroup::normalize(unnormalized.a, unnormalized.b, unnormalized.c);
    assert_eq!(
      (
        normalized_ground_truth.a,
        normalized_ground_truth.b,
        normalized_ground_truth.c
      ),
      elem_tuple
    );
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

    assert_eq!(ClassGroup::op(&a, &b), ground_truth);
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
    // 2. Compute g^* = g^100 * ... * g^1000.
    // 3. For each of g^100, g^200, ..., g^1000 compute the inverse of that element and assert that
    //    g^* * current_inverse = product of g^100, g^200, ..., g^1000 without the inversed-out
    //    element.
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

    for _ in 0..1000 {
      g = ClassGroup::op(&g, &g_anchor);
      let g_inv = ClassGroup::inv(&g);
      assert_eq!(id, ClassGroup::op(&g_inv, &g));
      assert_eq!(id, ClassGroup::op(&g, &g_inv));
      assert_eq!(g, ClassGroup::inv(&g_inv));
    }
  }

  #[test]
  fn test_exp_basic() {
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    for i in 1..=1000 {
      g = ClassGroup::op(&g, &g_anchor);
      assert_eq!(&g, &ClassGroup::exp(&g_anchor, &int(i)));
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
    ClassGroup::square(&mut g2);
    ClassGroup::square(&mut g2);

    assert_eq!(&g2, &g4);
  }

  #[test]
  fn test_square_repeated() {
    let g = ClassGroup::unknown_order_elem();
    let mut g1024 = ClassGroup::id();

    // g^1024
    for _ in 0..4096 {
      g1024 = ClassGroup::op(&g, &g1024);
    }

    // g^2
    let mut g2 = g.clone();

    // g^1024
    for _ in 0..12 {
      ClassGroup::square(&mut g2);
    }
    assert_eq!(g2, g1024);
  }
}
