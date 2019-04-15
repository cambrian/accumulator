//! Miscellaneous functions commonly used throughout the library.
use crate::group::Group;
use crate::hash::hash_to_prime;
use rug::Integer;
use std::hash::Hash;

// Get a tuple of mutable reference from a tuple.
#[macro_export]
macro_rules! mut_tuple_elems {
  ($ctx:expr, $($tpl_idx:tt),+) => {
    (
      $(
        &mut $ctx.inner.$tpl_idx,
      )*
    )
  };
}

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

pub fn prime_hash_product<T: Hash>(ts: &[T]) -> Integer {
  ts.iter().map(hash_to_prime).product()
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
}
