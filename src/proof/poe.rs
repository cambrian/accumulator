use crate::group::Group;
use crate::hash::{hash_to_prime, Blake2b};
use num::BigUint;
use num_integer::Integer;

#[allow(non_snake_case)]
#[derive(Debug, PartialEq, Eq)]
pub struct PoE<G: Group> {
  Q: G::Elem,
}

impl<G: Group> PoE<G> {
  /// See page 16 of B&B.
  pub fn prove(base: &G::Elem, exp: &BigUint, result: &G::Elem) -> PoE<G> {
    let l = hash_to_prime(&Blake2b::default, &(base, exp, result));
    let q = exp.div_floor(&l);
    PoE {
      Q: G::exp(&base, &q),
    }
  }

  /// See page 16 of B&B.
  pub fn verify(base: &G::Elem, exp: &BigUint, result: &G::Elem, proof: &PoE<G>) -> bool {
    let l = hash_to_prime(&Blake2b::default, &(base, exp, result));
    let r = exp % l.clone();
    // w = Q^l * u^r
    let w = G::op(&G::exp(&proof.Q, &l), &G::exp(&base, &r));
    w == *result
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{ElemFromUnsigned, UnknownOrderGroup, RSA2048};
  use crate::util::bu;

  #[test]
  fn test_poe_small_exp() {
    // 2^20 = 1048576
    let base = RSA2048::unknown_order_elem();
    let exp = bu(20u8);
    let result = RSA2048::elem_of(1_048_576u32);
    let proof = PoE::<RSA2048>::prove(&base, &exp, &result);
    assert!(PoE::verify(&base, &exp, &result, &proof));
    assert!(
      proof
        == PoE {
          Q: RSA2048::elem_of(1u8)
        }
    );

    // 2^35 = 34359738368
    let exp_2 = bu(35u8);
    let result_2 = RSA2048::elem_of(34_359_738_368u64);
    let proof_2 = PoE::<RSA2048>::prove(&base, &exp_2, &result_2);
    assert!(PoE::verify(&base, &exp_2, &result_2, &proof_2));
    assert!(
      proof_2
        == PoE {
          Q: RSA2048::elem_of(1u8)
        }
    );
  }
}
