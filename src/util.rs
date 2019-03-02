use crate::group::Group;
use crate::num::mpz::Mpz;
use rug::Integer;

/// Poor man's type-level programming.
/// This trait allows us to reflect "type-level" (i.e. static) information at runtime.
pub trait TypeRep: 'static {
  type Rep: 'static;
  fn rep() -> &'static Self::Rep;
}

pub fn int<T>(val: T) -> Integer
where
  Integer: From<T>,
{
  Integer::from(val)
}

/// Computes the `(xy)`th root of `g` given the `x`th and `y`th roots of `g` and `(x, y)` coprime.
/// Consider moving this to accumulator?
#[allow(clippy::similar_names)]
pub fn shamir_trick<G: Group>(
  xth_root: &G::Elem,
  yth_root: &G::Elem,
  x: &Integer,
  y: &Integer,
) -> Option<G::Elem> {
  if G::exp(xth_root, x) != G::exp(yth_root, y) {
    return None;
  }

  let (gcd, a, b) = <(Integer, Integer, Integer)>::from(x.gcd_cofactors_ref(&y));

  if gcd != int(1) {
    return None;
  }

  Some(G::op(&G::exp(xth_root, &b), &G::exp(yth_root, &a)))
}

/// Folds over `xs` but in a divide-and-conquer fashion: Instead of `F(F(F(F(acc, a), b), c), d))`
/// this computes `F(acc, F(F(a, b), F(c, d)))`.
pub fn divide_and_conquer<F, T: Clone, E>(f: F, acc: T, xs: &[T]) -> Result<T, E>
where
  F: Fn(&T, &T) -> Result<T, E>,
{
  if xs.is_empty() {
    return Ok(acc);
  }

  Ok(f(&acc, &divide_and_conquer_(&f, xs)?)?)
}

fn divide_and_conquer_<F, T: Clone, E>(f: &F, xs: &[T]) -> Result<T, E>
where
  F: Fn(&T, &T) -> Result<T, E>,
{
  if xs.len() == 1 {
    return Ok(xs[0].clone());
  }

  let mid = xs.len() / 2;
  let left = &xs[..mid];
  let right = &xs[mid..];
  Ok(f(
    &divide_and_conquer_(f, left)?,
    &divide_and_conquer_(f, right)?,
  )?)
}

// Linear Congruence Mpz context, default implementation, and solver.
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
  pub fn solve_linear_congruence(
    &mut self,
    mu: &mut Mpz,
    v: &mut Mpz,
    a: &Mpz,
    b: &Mpz,
    m: &Mpz,
  ) -> Option<()> {
    // Binary Quadratic Forms, 7.4.1
    self.g.gcdext(&mut self.d, &mut self.e, a, m);
    self.q.fdiv_qr(&mut self.r, b, &self.g);

    if self.r.sgn() != 0 {
      // No solution.
      return None;
    }

    mu.mul(&self.q, &self.d);
    mu.modulo_mut(m);
    v.fdiv_q(m, &self.g);
    Some(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, Rsa2048, UnknownOrderGroup};
  use crate::util::int;

  #[derive(Debug)]
  enum Never {}

  /// Merge-based computation of Integer array products. Faster than  the iterative `iter.product()`
  /// for really large Integers.
  fn merge_product(xs: &[Integer]) -> Integer {
    divide_and_conquer(
      |a, b| -> Result<Integer, Never> { Ok(int(a * b)) },
      int(1),
      &xs,
    )
    .unwrap()
  }

  #[test]
  fn test_shamir_trick() {
    let (x, y, z) = (&int(13), &int(17), &int(19));
    let xth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(y * z));
    let yth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(x * z));
    let xyth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), z);
    assert!(shamir_trick::<Rsa2048>(&xth_root, &yth_root, x, y) == Some(xyth_root));
  }

  #[test]
  fn test_shamir_trick_failure() {
    let (x, y, z) = (&int(7), &int(14), &int(19)); // Inputs not coprime.
    let xth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(y * z));
    let yth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(x * z));
    assert!(shamir_trick::<Rsa2048>(&xth_root, &yth_root, x, y) == None);
  }

  #[test]
  fn test_merge_product() {
    let ints = vec![int(3), int(5), int(7), int(9), int(11)];
    assert!(merge_product(&ints) == int(10395));
  }

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
