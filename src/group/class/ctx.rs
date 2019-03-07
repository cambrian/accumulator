//! Defines Mpz context struct used for no-reallocation class group computations.
use super::ClassElem;
use super::ClassGroup;
use super::CLASS_GROUP_DISCRIMINANT;
use crate::num::flint;
use crate::num::flint::fmpz;
use crate::num::mpz::{flint_mpz_struct, Mpz};
use crate::util::int;
use crate::util::TypeRep;
use rug::Integer;

macro_rules! mut_tuple_elems {
  ($tuple:expr, $($tpl_idx:tt),+) => {
    (
      $(
        &mut $tuple.$tpl_idx,
      )*
    )
  };
}

pub struct LinCongruenceCtx {
  ctx: (Mpz, Mpz, Mpz, Mpz, Mpz),
}

impl Default for LinCongruenceCtx {
  fn default() -> Self {
    Self {
      ctx: (
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
      ),
    }
  }
}

impl LinCongruenceCtx {
  pub fn solve_linear_congruence(
    &mut self,
    mu: &mut Mpz,
    v: &mut Mpz,
    a: &Mpz,
    b: &Mpz,
    m: &Mpz,
  ) -> Option<()> {
    let (g, d, e, q, r) = mut_tuple_elems!(self.ctx, 0, 1, 2, 3, 4);

    // Binary Quadratic Forms, 7.4.1
    g.gcdext(d, e, a, m);
    q.fdiv_qr(r, b, &g);

    if r.sgn() != 0 {
      // No solution.
      return None;
    }

    mu.mul(&q, &d);
    mu.modulo_mut(m);
    v.fdiv_q(m, &g);
    Some(())
  }
}

#[allow(non_snake_case)]
#[allow(clippy::type_complexity)]
pub struct ClassCtx {
  L_sq_op: Mpz,
  ctx: (
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
    Mpz,
  ),
  lin_cong_ctx: LinCongruenceCtx,
}

impl Default for ClassCtx {
  fn default() -> Self {
    let mut s = Self {
      L_sq_op: Mpz::default(),
      ctx: (
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
      ),
      lin_cong_ctx: LinCongruenceCtx::default(),
    };

    // Precomputation needed for NUDULP.
    s.L_sq_op.abs(&CLASS_GROUP_DISCRIMINANT);
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
    let (scratch,) = mut_tuple_elems!(self.ctx, 0);
    if Self::elem_is_normal(scratch, &a, &b, &c) {
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
    let (g, h, j, w, r, s, t, u, a, b, l, m, mut mu, mut v, mut lambda, mut sigma, k) =
      mut_tuple_elems!(self.ctx, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

    // Binary Quadratic Forms, 6.1.1
    g.add(&x.b, &y.b);
    g.fdiv_q_ui_mut(2);
    h.sub(&y.b, &x.b);
    h.fdiv_q_ui_mut(2);
    w.gcd(&x.a, &y.a);
    w.gcd_mut(&g);
    j.set(&w);
    r.set_ui(0);
    s.fdiv_q(&x.a, &w);
    t.fdiv_q(&y.a, &w);
    u.fdiv_q(&g, &w);
    a.mul(&t, &u);
    b.mul(&h, &u);
    m.mul(&s, &x.c);
    b.add_mut(&m);
    m.mul(&s, &t);
    self
      .lin_cong_ctx
      .solve_linear_congruence(&mut mu, &mut v, &a, &b, &m)
      .unwrap();

    a.mul(&t, &v);
    m.mul(&t, &mu);
    b.sub(&h, &m);
    m.set(&s);
    self
      .lin_cong_ctx
      .solve_linear_congruence(&mut lambda, &mut sigma, &a, &b, &m)
      .unwrap();

    a.mul(&v, &lambda);
    k.add(&mu, &a);
    l.mul(&k, &t);
    l.sub_mut(&h);
    l.fdiv_q_mut(&s);
    m.mul(&t, &u);
    m.mul_mut(&k);
    a.mul(&h, &u);
    m.sub_mut(&a);
    a.mul(&x.c, &s);
    m.sub_mut(&a);
    a.mul(&s, &t);
    m.fdiv_q_mut(&a);

    let mut ret = ClassElem::default();

    ret.a.mul(&s, &t);
    a.mul(&r, &u);
    ret.a.sub_mut(&a);

    ret.b.mul(&j, &u);
    a.mul(&m, &r);
    ret.b.add_mut(&a);
    a.mul(&k, &t);
    ret.b.sub_mut(&a);
    a.mul(&l, &s);
    ret.b.sub_mut(&a);

    ret.c.mul(&k, &l);
    a.mul(&j, &m);
    ret.c.sub_mut(&a);

    self.reduce_mut(&mut ret);
    ret
  }

  pub fn id(&mut self, d: &Mpz) -> ClassElem {
    let (a,) = mut_tuple_elems!(self.ctx, 0);

    // Binary Quadratic Forms, Definition 5.4
    let mut ret = ClassElem::default();
    ret.a.set_ui(1);
    ret.b.set_ui(1);
    a.sub(&ret.b, &d);
    ret.c.fdiv_q_ui(&a, 4);
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
    let (scratch,) = mut_tuple_elems!(self.ctx, 0);

    let mut d = Mpz::default();
    d.mul(&b, &b);
    scratch.mul(&a, &c);
    scratch.mul_ui_mut(4);
    d.sub_mut(&scratch);
    d
  }

  fn square_regular(&mut self, x: &mut ClassElem) {
    let (mut mu, mut v, m, a, old_a) = mut_tuple_elems!(self.ctx, 0, 1, 2, 3, 4);

    // Binary Quadratic Forms, 6.3.1
    self
      .lin_cong_ctx
      .solve_linear_congruence(&mut mu, &mut v, &x.b, &x.c, &x.a)
      .unwrap();

    m.mul(&x.b, &mu);
    m.sub_mut(&x.c);
    m.fdiv_q_mut(&x.a);

    old_a.set(&x.a);
    x.a.square_mut();

    a.mul(&mu, &old_a);
    a.mul_ui_mut(2);
    x.b.sub_mut(&a);

    x.c.mul(&mu, &mu);
    x.c.sub_mut(&m);

    self.reduce_mut(x);
  }

  #[allow(non_snake_case)]
  fn square_nudulp(&mut self, x: &mut ClassElem) {
    // Jacobson, Michael J., and Alfred J. Van Der Poorten. "Computational aspects of NUCOMP."
    // Algorithm 2 (Alg 2).

    let (
      G_sq_op,
      scratch,
      y_sq_op,
      By_sq_op,
      Dy_sq_op,
      bx_sq_op,
      by_sq_op,
      dx_sq_op,
      q_sq_op,
      t_sq_op,
      ax_sq_op,
      ay_sq_op,
      Q1_sq_op,
      x_sq_op,
      z_sq_op,
      dy_sq_op,
    ) = mut_tuple_elems!(self.ctx, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);

    let L_sq_op = &mut self.L_sq_op;

    // Step 1 in Alg 2.
    G_sq_op.gcdext(scratch, y_sq_op, &x.a, &x.b);
    By_sq_op.divexact(&x.a, &G_sq_op);
    Dy_sq_op.divexact(&x.b, &G_sq_op);

    // Step 2 in Alg 2.
    bx_sq_op.mul(&y_sq_op, &x.c);
    bx_sq_op.modulo_mut(&By_sq_op);
    by_sq_op.set(&By_sq_op);

    if by_sq_op.cmpabs(&L_sq_op) <= 0 {
      // Step 4 in Alg 2.
      dx_sq_op.mul(&bx_sq_op, &Dy_sq_op);
      dx_sq_op.sub_mut(&x.c);
      dx_sq_op.divexact_mut(&By_sq_op);
      x.a.mul(&by_sq_op, &by_sq_op);
      x.c.mul(&bx_sq_op, &bx_sq_op);
      t_sq_op.add(&bx_sq_op, &by_sq_op);
      t_sq_op.square_mut();

      x.b.sub_mut(&t_sq_op);
      x.b.add_mut(&x.a);
      x.b.add_mut(&x.c);
      t_sq_op.mul(&G_sq_op, &dx_sq_op);
      x.c.sub_mut(&t_sq_op);
      self.reduce_mut(x);
      return;
    }

    // Most of Step 3 in Alg 2.
    if cfg!(feature = "flint") {
      // Subroutine as handled by top entry to the Chia VDF competition "bulaiden."
      unsafe {
        let mut fy: fmpz = 0;
        let mut fx: fmpz = 0;
        let mut fby: fmpz = 0;
        let mut fbx: fmpz = 0;
        let mut fL: fmpz = 0;

        let mut y_square_clone = flint_mpz_struct::from(*y_sq_op);
        let mut x_square_clone = flint_mpz_struct::from(*x_sq_op);
        let mut by_square_clone = flint_mpz_struct::from(*by_sq_op);
        let mut bx_square_clone = flint_mpz_struct::from(*bx_sq_op);
        let mut L_square_clone = flint_mpz_struct::from(*L_sq_op);

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

        *y_sq_op = Mpz::from(y_square_clone);
        *x_sq_op = Mpz::from(x_square_clone);
        *by_sq_op = Mpz::from(by_square_clone);
        *bx_sq_op = Mpz::from(bx_square_clone);
      }

      x_sq_op.neg_mut();
      if x_sq_op.sgn() > 0 {
        y_sq_op.neg_mut();
      } else {
        by_sq_op.neg_mut();
      }
    } else {
      // Subroutine as presented in "Computational aspects of NUCOMP.", Algorithm 2, most of step 3.
      x_sq_op.set_ui(1);
      y_sq_op.set_ui(0);
      z_sq_op.set_ui(0);
      while by_sq_op.cmpabs(&L_sq_op) > 0 && bx_sq_op.sgn() == 0 {
        q_sq_op.fdiv_q(&by_sq_op, &bx_sq_op);
        t_sq_op.fdiv_r(&by_sq_op, &bx_sq_op);

        by_sq_op.set(&bx_sq_op);
        bx_sq_op.set(&t_sq_op);
        y_sq_op.submul(&q_sq_op, &x_sq_op);
        t_sq_op.set(&y_sq_op);
        y_sq_op.set(&x_sq_op);
        x_sq_op.set(&t_sq_op);
        z_sq_op.add_ui_mut(1);
      }

      if z_sq_op.odd() != 0 {
        by_sq_op.neg_mut();
        y_sq_op.neg_mut();
      }
    }

    ax_sq_op.mul(&G_sq_op, &x_sq_op);
    ay_sq_op.mul(&G_sq_op, &y_sq_op);

    // Step 5 in Alg 2.
    t_sq_op.mul(&Dy_sq_op, &bx_sq_op);
    t_sq_op.submul(&x.c, &x_sq_op);
    dx_sq_op.divexact(&t_sq_op, &By_sq_op);
    Q1_sq_op.mul(&y_sq_op, &dx_sq_op);
    dy_sq_op.add(&Q1_sq_op, &Dy_sq_op);
    x.b.add(&dy_sq_op, &Q1_sq_op);
    x.b.mul_mut(&G_sq_op);
    dy_sq_op.divexact_mut(&x_sq_op);
    x.a.mul(&by_sq_op, &by_sq_op);
    x.c.mul(&bx_sq_op, &bx_sq_op);
    t_sq_op.add(&bx_sq_op, &by_sq_op);
    x.b.submul(&t_sq_op, &t_sq_op);
    x.b.add_mut(&x.a);
    x.b.add_mut(&x.c);
    x.a.submul(&ay_sq_op, &dy_sq_op);
    x.c.submul(&ax_sq_op, &dx_sq_op);

    self.reduce_mut(x);
  }

  fn elem_is_reduced(scratch: &mut Mpz, a: &Mpz, b: &Mpz, c: &Mpz) -> bool {
    Self::elem_is_normal(scratch, a, b, c) && (a <= c && !(a == c && b.cmp_si(0) < 0))
  }

  fn elem_is_normal(scratch: &mut Mpz, a: &Mpz, b: &Mpz, _c: &Mpz) -> bool {
    scratch.neg(&a);
    let ret_val = *scratch < *b && b <= a;
    ret_val
  }

  fn reduce_(&mut self, elem: &mut ClassElem) {
    let (x, s, old_a, old_b, scratch) = mut_tuple_elems!(self.ctx, 0, 1, 2, 3, 4);

    // Binary Quadratic Forms, 5.2.1
    while !Self::elem_is_reduced(scratch, &elem.a, &elem.b, &elem.c) {
      s.add(&elem.c, &elem.b);
      x.mul_ui(&elem.c, 2);
      s.fdiv_q_mut(&x);
      old_a.set(&elem.a);
      old_b.set(&elem.b);
      elem.a.set(&elem.c);
      elem.b.neg_mut();
      x.mul(&s, &elem.c);
      x.mul_ui_mut(2);
      elem.b.add_mut(&x);

      elem.c.mul_mut(&s);
      elem.c.mul_mut(&s);
      x.mul(&old_b, &s);
      elem.c.sub_mut(&x);
      elem.c.add_mut(&old_a);
    }
  }

  fn reduce_mut(&mut self, x: &mut ClassElem) {
    self.normalize_mut(x);
    self.reduce_(x);
    self.normalize_mut(x);
  }

  fn normalize_mut(&mut self, x: &mut ClassElem) {
    let (scratch,) = mut_tuple_elems!(self.ctx, 0);
    if Self::elem_is_normal(scratch, &x.a, &x.b, &x.c) {
      return;
    }
    self.normalize_(&mut x.a, &mut x.b, &mut x.c);
  }

  fn normalize_(&mut self, a: &mut Mpz, b: &mut Mpz, c: &mut Mpz) {
    let (r, denom, old_b, ra) = mut_tuple_elems!(self.ctx, 0, 1, 2, 3);

    // Binary Quadratic Forms, 5.1.1
    r.sub(&a, &b);
    denom.mul_ui(&a, 2);
    r.fdiv_q_mut(&denom);

    old_b.set(&b);

    ra.mul(&r, &a);
    b.add_mut(&ra);
    b.add_mut(&ra);

    ra.mul_mut(&r);
    c.add_mut(&ra);

    ra.set(&r);
    ra.mul_mut(&old_b);
    c.add_mut(&ra);
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_linear_congruence_solver() {
    fn test_congruence_problem_w_solution(x: i64, y: i64, z: i64, mu_i64: u64, v_i64: u64) {
      let mut ctx: LinCongruenceCtx = LinCongruenceCtx::default();
      let (mut a, mut b, mut c) = (Mpz::default(), Mpz::default(), Mpz::default());
      let (mut mu, mut v) = (Mpz::default(), Mpz::default());
      let (mut mu_truth, mut v_truth) = (Mpz::default(), Mpz::default());

      a.set_si(x);
      b.set_si(y);
      c.set_si(z);

      mu_truth.set_ui(mu_i64);
      v_truth.set_ui(v_i64);

      ctx
        .solve_linear_congruence(&mut mu, &mut v, &a, &b, &c)
        .unwrap();
      assert_eq!((mu_truth, v_truth), (mu, v));
    }

    test_congruence_problem_w_solution(3, 2, 4, 2, 4);
    test_congruence_problem_w_solution(5, 1, 2, 1, 2);
    test_congruence_problem_w_solution(2, 4, 5, 2, 5);
    test_congruence_problem_w_solution(230, 1081, 12167, 2491, 529);
  }

  #[test]
  fn test_linear_congruence_solver_no_solution() {
    fn test_congruence_problem_no_solution(x: i64, y: i64, z: i64) {
      let mut ctx: LinCongruenceCtx = LinCongruenceCtx::default();
      let (mut a, mut b, mut c) = (Mpz::default(), Mpz::default(), Mpz::default());
      let (mut mu, mut v) = (Mpz::default(), Mpz::default());

      a.set_si(x);
      b.set_si(y);
      c.set_si(z);

      // Let g = gcd(a, m). If b is not divisible by g, there are no solutions. If b is divisible by
      // g, there are g solutions.
      let result = ctx.solve_linear_congruence(&mut mu, &mut v, &a, &b, &c);
      assert!(result.is_none());
    }

    test_congruence_problem_no_solution(33, 7, 143);
    test_congruence_problem_no_solution(13, 14, 39);
  }
}
