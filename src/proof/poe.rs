use crate::group::Group;
use crate::hash::{hash_to_prime, Blake2b};
use crate::util::int;
use rug::Integer;

#[allow(non_snake_case)]
#[derive(Debug, PartialEq, Eq)]
pub struct PoE<G: Group> {
  Q: G::Elem,
}

impl<G: Group> PoE<G> {
  /// See page 16 of B&B.
  pub fn prove(base: &G::Elem, exp: &Integer, result: &G::Elem) -> PoE<G> {
    let l = hash_to_prime(&Blake2b::default, &(base, exp, result));
    let q = exp / l;
    PoE {
      Q: G::exp(&base, &q),
    }
  }

  /// See page 16 of B&B.
  pub fn verify(base: &G::Elem, exp: &Integer, result: &G::Elem, proof: &PoE<G>) -> bool {
    let l = hash_to_prime(&Blake2b::default, &(base, exp, result));
    let r = int(exp % &l);
    // w = Q^l * u^r
    let w = G::op(&G::exp(&proof.Q, &l), &G::exp(&base, &r));
    w == *result
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{ElemFrom, UnknownOrderGroup, RSA2048};
  use crate::util::int;

  #[test]
  fn test_poe_small_exp() {
    // 2^20 = 1048576
    let base = RSA2048::unknown_order_elem();
    let exp = int(20);
    let result = RSA2048::elem(1_048_576);
    let proof = PoE::<RSA2048>::prove(&base, &exp, &result);
    assert!(PoE::verify(&base, &exp, &result, &proof));
    assert!(
      proof
        == PoE {
          Q: RSA2048::elem(1)
        }
    );

    // 2^35 = 34359738368
    let exp_2 = int(35);
    let result_2 = RSA2048::elem(34_359_738_368u64);
    let proof_2 = PoE::<RSA2048>::prove(&base, &exp_2, &result_2);
    assert!(PoE::verify(&base, &exp_2, &result_2, &proof_2));
    assert!(
      proof_2
        == PoE {
          Q: RSA2048::elem(1)
        }
    );
  }
}
