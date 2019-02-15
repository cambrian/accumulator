use crate::group::Group;
use rug::Integer;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Never {}

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

/// Merge-based computation of Integer array products. Faster than  the iterative `iter.product()`
/// for really large Integers.
pub fn merge_product(xs: &mut [Integer]) -> Integer {
  try_merge_reduce(
    |a, b| -> Result<Integer, Never> { Ok(int(a * b)) },
    int(1),
    xs,
  )
  .unwrap()
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

/// Folds over `xs` but in a fashion that repeatedly merges adjacent elements: Instead of
/// `F(F(F(F(acc, a), b), c), d))` this computes `F(acc, F(F(a, b), F(c, d)))`.
pub fn try_merge_reduce<F, T, E>(merge: F, acc: T, xs: &mut [T]) -> Result<T, E>
where
  F: Fn(&T, &T) -> Result<T, E>,
{
  // First iter is merge_skip == 2.
  let mut merge_step = 2_usize.pow(0);
  let mut merge_skip = 2_usize.pow(1);
  let num_xs = xs.len();
  loop {
    for i in (0..num_xs).step_by(merge_skip) {
      if i + merge_step >= num_xs {
        break;
      }

      xs[i] = merge(&xs[i], &xs[i + merge_step])?;
    }

    merge_step *= 2;
    merge_skip *= 2;

    if merge_skip >= xs.len() {
      break;
    }
  }

  if num_xs == 0 {
    Ok(acc)
  } else {
    Ok(merge(&acc, &xs[0])?)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, Rsa2048, UnknownOrderGroup};
  use crate::util::int;

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
}
