use crate::internal::group::Group;
use crate::internal::hash::hash_to_prime;
use crate::internal::util::int;
use rug::Integer;

#[allow(non_snake_case)]
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Poe<G: Group> {
  Q: G::Elem,
}

impl<G: Group> Poe<G> {
  /// See page 16 of B&B.
  pub fn prove(base: &G::Elem, exp: &Integer, result: &G::Elem) -> Self {
    let l = hash_to_prime(&(base, exp, result));
    let q = exp / l;
    Self {
      Q: G::exp(&base, &q),
    }
  }

  /// See page 16 of B&B.
  pub fn verify(base: &G::Elem, exp: &Integer, result: &G::Elem, proof: &Self) -> bool {
    let l = hash_to_prime(&(base, exp, result));
    let r = int(exp % &l);
    // w = Q^l * u^r
    let w = G::op(&G::exp(&proof.Q, &l), &G::exp(&base, &r));
    w == *result
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::internal::group::{ElemFrom, Rsa2048, UnknownOrderGroup};
  use crate::internal::util::int;

  #[test]
  fn test_poe_small_exp() {
    // 2^20 = 1048576
    let base = Rsa2048::unknown_order_elem();
    let exp = int(20);
    let result = Rsa2048::elem(1_048_576);
    let proof = Poe::<Rsa2048>::prove(&base, &exp, &result);
    assert!(Poe::verify(&base, &exp, &result, &proof));
    assert!(
      proof
        == Poe {
          Q: Rsa2048::elem(1)
        }
    );

    // 2^35 = 34359738368
    let exp_2 = int(35);
    let result_2 = Rsa2048::elem(34_359_738_368u64);
    let proof_2 = Poe::<Rsa2048>::prove(&base, &exp_2, &result_2);
    assert!(Poe::verify(&base, &exp_2, &result_2, &proof_2));
    assert!(
      proof_2
        == Poe {
          Q: Rsa2048::elem(1)
        }
    );
  }
}
