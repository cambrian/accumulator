use super::super::group::Group;
use super::super::hash::hashes;
use num::BigUint;
use num_integer::Integer;
use serde::ser::Serialize;

#[allow(non_snake_case)]
#[derive(PartialEq, Eq)]
pub struct PoE<G: Group> {
  Q: G::Elem,
}

impl<G: Group> PoE<G> {
  /// See page 16 of B&B.
  pub fn prove(base: &G::Elem, exp: &BigUint, result: &G::Elem) -> PoE<G> {
    let l = hash_prime(exp, base, result);
    let q = exp.div_floor(&l);
    PoE {
      Q: G::exp(&base, &q),
    }
  }

  /// See page 16 of B&B.
  pub fn verify(base: &G::Elem, exp: &BigUint, result: &G::Elem, proof: &PoE<G>) -> bool {
    let l = hash_prime(exp, base, result);
    let r = exp % l.clone();
    // w = Q^l * u^r
    let w = G::op(&G::exp(&proof.Q, &l), &G::exp(&base, &r));
    w == *result
  }
}

<<<<<<< HEAD
fn hash_prime<G: Serialize>(_exp: &BigUint, _base: &G, _result: &G) -> BigUint {
=======
/// See page 16 of B&B.
pub fn verify_poe<G: Group>(
  base: &G::Elem,
  exp: &BigUint,
  result: &G::Elem,
  proof: &PoE<G>,
) -> bool {
  let l = hash_prime(exp, base, result);
  let r = exp % l.clone();
  // w = Q^l * u^r
  let w = G::op(&G::exp(&proof.Q, &l), &G::exp(&base, &r));
  w == *result
}

fn hash_prime<G: Serialize>(exp: &BigUint, base: &G, result: &G) -> BigUint {
>>>>>>> Using hash prime function
  // TODO: Replace with commented out when hash_prime is implemented.
  let mut hash_string = exp.to_str_radix(16);
  hash_string.push_str(&serde_json::to_string(&base).unwrap());
  hash_string.push_str(&serde_json::to_string(&result).unwrap());
  hashes::h_prime(&hashes::blake2, hash_string.as_bytes())
}

#[cfg(test)]
mod tests {
  use super::super::super::group::dummy::DummyRSA;
  use super::*;

  // Current exponents are far smaller than generated primes, so all PoE proofs are producing
  // Q = 1 proofs, which is pretty useless. This will be remedied when we implement RSA2048
  #[test]
  fn test_poe() {
    // 2^20 = 1048576
    let base = DummyRSA::base_elem();
    let exp = BigUint::from(20 as u8);
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
    let exp_2 = BigUint::from(35 as u8);
    let result_2 = DummyRSA::elem_of(34_359_738_368);
    let proof_2 = PoE::<DummyRSA>::prove(&base, &exp_2, &result_2);
    assert!(PoE::verify(&base, &exp_2, &result_2, &proof_2));
    assert!(
      proof_2
        == PoE {
          Q: DummyRSA::elem_of(1)
        }
    );
<<<<<<< HEAD
    // Cannot verify wrong base/exp/result triple with wrong pair.
    assert!(!PoE::verify(&base, &exp_2, &result_2, &proof));
=======
    // // Cannot verify wrong base/exp/result triple with wrong pair.
    // assert!(!verify_poe::<DummyRSA>(&base, &exp_2, &result_2, &proof));
>>>>>>> Using hash prime function
  }
}
