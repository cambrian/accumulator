use crate::util::{int, Singleton};
use rug::Integer;
use std::hash::Hash;
use std::marker::Sized;

mod class;
mod rsa;
pub use rsa::{RSA2048Elem, RSA2048};

/// We need a runtime representation for the group itself because reading in group parameters
/// (i.e. RSA modulus) is infeasible to do at the type-level in Rust.
///
/// We mimic type-level programming by using the singleton pattern here. For each group, there
/// should be a single, constant, static instance of the group representation accessible at all
/// times. This way, we can "reflect" information about the group type by accessing the singleton.
/// Refer to dummy.rs for an example.
pub trait Group: Singleton {
  /// In theory the association Group::Elem is bijective, such that it makes sense to write
  /// something like Elem::Group::get(). This would let us define op, exp, inv, etc on the Elem
  /// type and avoid using prefix notation for all of our group operations.
  /// But afaik bijective associated types are not supported by Rust.
  type Elem: Eq + Clone + Sized + Hash;

  fn id_(rep: &Self::Rep) -> Self::Elem;

  fn op_(rep: &Self::Rep, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

  /// Default implementation of exponentiation via repeated squaring.
  /// Group implementations may provide more performant specializations
  /// (e.g. Montgomery multiplication for RSA groups).
  fn exp_(_rep: &Self::Rep, a: &Self::Elem, n: &Integer) -> Self::Elem {
    let (mut val, mut a, mut n) = {
      if *n < int(0) {
        (Self::id(), Self::inv(a), int(-n))
      } else {
        (Self::id(), a.clone(), n.clone())
      }
    };
    loop {
      if n == int(0) {
        return val;
      }
      if n.is_odd() {
        val = Self::op(&val, &a);
      }
      a = Self::op(&a, &a);
      n >>= 1;
    }
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
pub trait UnknownOrderGroup: Group {
  /// E.g. 2, for RSA groups.
  fn unknown_order_elem_(rep: &Self::Rep) -> Self::Elem;

  fn unknown_order_elem() -> Self::Elem {
    Self::unknown_order_elem_(Self::rep())
  }
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
  // G::op expects a Integer.
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
    let alpha_1 = RSA2048::elem(2);
    let alpha_2 = RSA2048::elem(3);
    let x_1 = int(3);
    let x_2 = int(2);
    let res = multi_exp::<RSA2048>(
      &[alpha_1.clone(), alpha_2.clone()],
      &[x_1.clone(), x_2.clone()],
    );
    assert!(res == RSA2048::elem(108));
    let alpha_3 = RSA2048::elem(5);
    let x_3 = int(1);
    let res_2 = multi_exp::<RSA2048>(&[alpha_1, alpha_2, alpha_3], &[x_1, x_2, x_3]);
    assert!(res_2 == RSA2048::elem(1_687_500));
  }
}

/// Like From<T>, but implemented on the Group instead of the element type.
pub trait ElemFrom<T>: Group {
  fn elem(val: T) -> Self::Elem;
}
