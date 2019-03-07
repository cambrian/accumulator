//! This module contains implementations for different mathematical groups, each of which satisfies
//! our `UnknownOrderGroup` trait. They can be used with the accumulator and vector commitment
//! structures, or standalone if you have a custom application.
//!
//! The preferred elliptic group implementation is the Ristretto group, which is a cyclic subset of
//! the Ed25519 group.
use crate::util::{int, TypeRep};
use rug::Integer;
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::Sized;

mod class;
pub use class::{ClassElem, ClassGroup};
mod ristretto;
pub use ristretto::{Ristretto, RistrettoElem};
mod rsa;
pub use rsa::{Rsa2048, Rsa2048Elem};

/// We avoid having to pass group objects around by using the `TypeRep` trait.
///
/// The other traits are only required here because Rust can't figure out how to do stuff with an
/// `Accumulator<G>` even though it's just a wrapped `G::Elem`. If possible we'd remove them.
pub trait Group: Clone + Debug + Eq + Hash + TypeRep + Send + Sync {
  /// In theory the association `Group::Elem` is bijective, such that it makes sense to write
  /// something like `Elem::Group::get()`. This would let us define `op`, `exp`, `inv`, etc. on the
  /// `Elem` type and avoid using prefix notation for all of our group operations.
  /// But AFAIK bijective associated types are not supported by Rust.
  type Elem: Clone + Debug + Eq + Hash + Sized + Send + Sync;

  fn id_(rep: &Self::Rep) -> Self::Elem;

  fn op_(rep: &Self::Rep, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

  /// Default implementation of exponentiation via repeated squaring.
  /// Implementations may provide more performant specializations
  /// (e.g. Montgomery multiplication for RSA groups).
  fn exp_(_rep: &Self::Rep, a: &Self::Elem, n: &Integer) -> Self::Elem {
    let (mut val, mut a, mut n) = {
      if *n < int(0) {
        (Self::id(), Self::inv(a), int(-n))
      } else {
        (Self::id(), a.clone(), n.clone())
      }
    };
    while n > int(0) {
      if n.is_odd() {
        val = Self::op(&val, &a);
      }
      a = Self::op(&a, &a);
      n >>= 1;
    }
    val
  }

  fn inv_(rep: &Self::Rep, a: &Self::Elem) -> Self::Elem;

  // -------------------
  // END OF REQUIRED FNS
  // -------------------

  fn id() -> Self::Elem {
    Self::id_(Self::rep())
  }

  fn op(a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
    Self::op_(Self::rep(), a, b)
  }

  fn exp(a: &Self::Elem, n: &Integer) -> Self::Elem {
    Self::exp_(Self::rep(), a, n)
  }

  fn inv(a: &Self::Elem) -> Self::Elem {
    Self::inv_(Self::rep(), a)
  }
}

/// We use this to mean a group containing elements of unknown order, not necessarily that the group
/// itself has unknown order. E.g. RSA groups.
#[allow(clippy::stutter)]
pub trait UnknownOrderGroup: Group {
  /// E.g. 2, for RSA groups.
  fn unknown_order_elem_(rep: &Self::Rep) -> Self::Elem;

  fn unknown_order_elem() -> Self::Elem {
    Self::unknown_order_elem_(Self::rep())
  }
}

/// Like `From<T>`, but implemented on the Group instead of the element type.
pub trait ElemFrom<T>: Group {
  fn elem(val: T) -> Self::Elem;
}

pub fn multi_exp<G: Group>(alphas: &[G::Elem], x: &[Integer]) -> G::Elem {
  if alphas.len() == 1 {
    return alphas[0].clone();
  }

  let n_half = alphas.len() / 2;
  let alpha_l = &alphas[..n_half];
  let alpha_r = &alphas[n_half..];
  let x_l = &x[..n_half];
  let x_r = &x[n_half..];
  let x_star_l = x_l.iter().product();
  let x_star_r = x_r.iter().product();
  let l = multi_exp::<G>(alpha_l, x_l);
  let r = multi_exp::<G>(alpha_r, x_r);
  G::op(&G::exp(&l, &x_star_r), &G::exp(&r, &x_star_l))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::util::int;

  #[test]
  fn test_multi_exp() {
    let alpha_1 = Rsa2048::elem(2);
    let alpha_2 = Rsa2048::elem(3);
    let x_1 = int(3);
    let x_2 = int(2);
    let res = multi_exp::<Rsa2048>(
      &[alpha_1.clone(), alpha_2.clone()],
      &[x_1.clone(), x_2.clone()],
    );
    assert!(res == Rsa2048::elem(108));
    let alpha_3 = Rsa2048::elem(5);
    let x_3 = int(1);
    let res_2 = multi_exp::<Rsa2048>(&[alpha_1, alpha_2, alpha_3], &[x_1, x_2, x_3]);
    assert!(res_2 == Rsa2048::elem(1_687_500));
  }
}
