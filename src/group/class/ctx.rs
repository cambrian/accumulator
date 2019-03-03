//! Defines mpz context struct used for no-rellocation class group computations.
use super::ClassElem;
use super::ClassGroup;
use super::CLASS_GROUP_DISCRIMINANT;
use crate::num::flint;
use crate::num::flint::fmpz;
use crate::num::mpz::{flint_mpz_struct, Mpz};
use crate::util::int;
use crate::util::{LinCongruenceCtx, TypeRep};
use rug::Integer;

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

//  Class group operations based on Chia's fantastic doc explaining applied class groups:
//  https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf,
//  hereafter refered to as "Binary Quadratic Forms".  Includes optimizations from Chia's VDF
//  competition and from Jacobson, Michael J., and Alfred J. Van Der Poorten.
//  "Computational aspects of NUCOMP."
impl ClassCtx {
  pub fn normalize(&mut self, mut a: Mpz, mut b: Mpz, mut c: Mpz) -> (Mpz, Mpz, Mpz) {
    if self.elem_is_normal(&a, &b, &c) {
      return (a, b, c);
    }
    self.normalize_(&mut a, &mut b, &mut c);
    (a, b, c)
  }

  pub fn reduce(&mut self, a: Mpz, b: Mpz, c: Mpz) -> (Mpz, Mpz, Mpz) {
    let (a, b, c) = self.normalize(a, b, c);
    let mut elem = ClassElem { a, b, c };
    self.reduce_(&mut elem);
    self.normalize(elem.a, elem.b, elem.c)
  }

  pub fn validate(&mut self, a: &Mpz, b: &Mpz, c: &Mpz) -> bool {
    self.discriminant(a, b, c) == *ClassGroup::rep()
  }

  pub fn square(&mut self, x: &mut ClassElem) {
    if cfg!(feature = "nudulp") {
      // NUDULP optimization, maybe using FLINT bindings.
      self.square_nudulp(x);
    } else {
      self.square_regular(x);
    }
  }

  pub fn op(&mut self, x: &ClassElem, y: &ClassElem) -> ClassElem {
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

  pub fn id(&mut self, d: &Mpz) -> ClassElem {
    // Binary Quadratic Forms, Definition 5.4
    let mut ret = ClassElem::default();
    ret.a.set_ui(1);
    ret.b.set_ui(1);
    self.a.sub(&ret.b, &d);
    ret.c.fdiv_q_ui(&self.a, 4);
    ret
  }

  pub fn exp(&mut self, d: &Mpz, a: &ClassElem, n: &Integer) -> ClassElem {
    let (mut val, mut a, mut n) = {
      if *n < int(0) {
        (self.id(d), self.inv(&a), int(-n))
      } else {
        (self.id(d), a.clone(), n.clone())
      }
    };
    loop {
      if n == int(0) {
        return val;
      }
      if n.is_odd() {
        val = self.op(&val, &a);
      }

      self.square(&mut a);
      n >>= 1;
    }
  }

  pub fn inv(&mut self, x: &ClassElem) -> ClassElem {
    let mut ret = ClassElem::default();
    ret.a.set(&x.a);
    ret.b.neg(&x.b);
    ret.c.set(&x.c);
    ret
  }

  pub fn generator(&mut self, d: &Mpz) -> ClassElem {
    // Binary Quadratic Forms, Definition 5.4
    let mut ret = ClassElem::default();
    ret.a.set_ui(2);
    ret.b.set_ui(1);
    ret.c.set_ui(1);
    ret.c.sub_mut(&d);
    ret.c.fdiv_q_ui_mut(8);

    let (a, b, c) = self.reduce(ret.a, ret.b, ret.c);
    ClassElem { a, b, c }
  }

  // Private functions

  fn discriminant(&mut self, a: &Mpz, b: &Mpz, c: &Mpz) -> Mpz {
    let mut d = Mpz::default();
    d.mul(&b, &b);
    self.scratch.mul(&a, &c);
    self.scratch.mul_ui_mut(4);
    d.sub_mut(&self.scratch);
    d
  }

  fn square_regular(&mut self, x: &mut ClassElem) {
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

  fn square_nudulp(&mut self, x: &mut ClassElem) {
    // Jacobson, Michael J., and Alfred J. Van Der Poorten. "Computational aspects of NUCOMP."
    // Algorithm 2 (Alg 2).

    // Step 1 in Alg 2.
    self
      .G_sq_op
      .gcdext(&mut self.scratch, &mut self.y_sq_op, &x.a, &x.b);
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
      self.reduce_mut(x);
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

  fn elem_is_reduced(&mut self, a: &Mpz, b: &Mpz, c: &Mpz) -> bool {
    self.elem_is_normal(a, b, c) && (a <= c && !(a == c && b.cmp_si(0) < 0))
  }

  fn elem_is_normal(&mut self, a: &Mpz, b: &Mpz, _c: &Mpz) -> bool {
    self.scratch.neg(&a);
    self.scratch < *b && b <= a
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

  fn reduce_(&mut self, elem: &mut ClassElem) {
    // Binary Quadratic Forms, 5.2.1
    while !self.elem_is_reduced(&elem.a, &elem.b, &elem.c) {
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

  fn reduce_mut(&mut self, x: &mut ClassElem) {
    self.normalize_mut(x);
    self.reduce_(x);
    self.normalize_mut(x);
  }

  fn normalize_mut(&mut self, x: &mut ClassElem) {
    if self.elem_is_normal(&x.a, &x.b, &x.c) {
      return;
    }
    self.normalize_(&mut x.a, &mut x.b, &mut x.c);
  }

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
}
