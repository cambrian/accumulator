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
  use crate::group::{DummyRSA, UnknownOrderGroup};
  use crate::util::bu;

  // Current exponents are far smaller than generated primes, so all PoE proofs are producing
  // pretty useless Q = 1 proofs. This will be remedied when we implement RSA2048.
  #[test]
  fn test_poe() {
    // 2^20 = 1048576
    let base = DummyRSA::unknown_order_elem();
    let exp = bu(20u8);
    let result = DummyRSA::elem_of(1_048_576);
    let proof = PoE::<DummyRSA>::prove(&base, &exp, &result);
    assert!(PoE::verify(&base, &exp, &result, &proof));
    assert!(
      proof
        == PoE {
          Q: DummyRSA::elem_of(1)
        }
    );

    // 2^35 = 34359738368
    let exp_2 = bu(35u8);
    let result_2 = DummyRSA::elem_of(34_359_738_368);
    let proof_2 = PoE::<DummyRSA>::prove(&base, &exp_2, &result_2);
    assert!(PoE::verify(&base, &exp_2, &result_2, &proof_2));
    assert!(
      proof_2
        == PoE {
          Q: DummyRSA::elem_of(1)
        }
    );
    // Cannot verify wrong base/exp/result triple with wrong pair.
    // assert!(!verify_poe::<DummyRSA>(&base, &exp_2, &result_2, &proof));
  }
}
