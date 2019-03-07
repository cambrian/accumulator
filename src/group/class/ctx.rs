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
  ($ctx:expr, $($tpl_idx:tt),+) => {
    (
      $(
        &mut $ctx.inner.$tpl_idx,
      )*
    )
  };
}

pub struct LinCongruenceCtx {
  pub inner: (Mpz, Mpz, Mpz, Mpz, Mpz),
}

impl Default for LinCongruenceCtx {
  fn default() -> Self {
    Self {
      inner: (
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
        Mpz::default(),
      ),
    }
  }
}

pub struct OpCtx {
  pub inner: (
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
}

impl Default for OpCtx {
  fn default() -> Self {
    Self {
      inner: (
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
    let (g, d, e, q, r) = mut_tuple_elems!(self, 0, 1, 2, 3, 4);

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
  pub L_sq_op: Mpz,
  pub op_ctx: OpCtx,
  pub lin_cong_ctx: LinCongruenceCtx,
}

impl Default for ClassCtx {
  fn default() -> Self {
    let mut s = Self {
      L_sq_op: Mpz::default(),
      op_ctx: OpCtx::default(),
      lin_cong_ctx: LinCongruenceCtx::default(),
    };

    // Precomputation needed for NUDULP.
    s.L_sq_op.abs(&CLASS_GROUP_DISCRIMINANT);
    s.L_sq_op.root_mut(4);
    s
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
