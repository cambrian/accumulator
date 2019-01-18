use super::super::group::Group;
// use super::super::hash::hashes;
use num::BigUint;
use serde::ser::Serialize;

pub struct PoE<T> {
  q: T,
}

/// See page 16 of B&B.
pub fn prove_poe<G: Group>(base: &G::Elem, exp: &BigUint, result: &G::Elem) -> PoE<G::Elem> {
  let l = hash_prime(exp, base, result);
  let q = exp / l;
  PoE {
    q: G::exp(&base, &q),
  }
}

/// See page 16 of B&B.
pub fn verify_poe<G: Group>(
  base: &G::Elem,
  exp: &BigUint,
  result: &G::Elem,
  proof: &PoE<G::Elem>,
) -> bool {
  let l = hash_prime(exp, base, result);
  let r = exp % l.clone();
  // w = Q^l * u^r
  let w = G::op(&G::exp(&proof.q, &l), &G::exp(&base, &r));
  w == *result
}

fn hash_prime<G: Serialize>(_exp: &BigUint, _base: &G, _result: &G) -> BigUint {
  // TODO: replace with commented out when hash_prime is implemented
  BigUint::from(13 as u8)
  // let mut hash_string = exp.to_str_radix(16);
  // hash_string.push_str(&serde_json::to_string(&base).unwrap());
  // hash_string.push_str(&serde_json::to_string(&result).unwrap());
  // hashes::h_prime(&hashes::blake2, hash_string.as_bytes())
}

#[cfg(test)]
mod tests {
  use super::super::super::group::dummy;
  use super::*;

  #[test]
  fn test_poe() {
    let dummy = dummy::DummyRSA::get();
    let base = dummy::DummyRSA::base_elem_(&dummy);
    let exp = BigUint::from(20 as u8);
    let result = dummy::DummyRSA::op(&1_048_576, &1);
    let proof = prove_poe::<dummy::DummyRSA>(&base, &exp, &result);
    // 2^20 = 1048576
    assert!(verify_poe::<dummy::DummyRSA>(&base, &exp, &result, &proof));
    let exp_2 = BigUint::from(35 as u8);
    let result_2 = dummy::DummyRSA::op(&34_359_738_368, &1);
    let proof_2 = prove_poe::<dummy::DummyRSA>(&base, &exp_2, &result_2);
    // 2^35 = 34359738368
    assert!(verify_poe::<dummy::DummyRSA>(
      &base, &exp_2, &result_2, &proof_2
    ));
    // Cannot verify wrong base/exp/result triple with wrong pair
    assert!(!verify_poe::<dummy::DummyRSA>(
      &base, &exp_2, &result_2, &proof
    ));
  }
}
