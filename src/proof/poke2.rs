use super::super::group::Group;
use super::super::group::InvertibleGroup;
use super::super::hash::hashes;
use num::BigInt;
use num::BigUint;
use num_bigint::Sign::Plus;
use serde::ser::Serialize;

pub struct PoKE2<T> {
  z: T,
  q: T,
  r: BigUint,
}

/// See page 16 of B&B.
pub fn prove_poke2<G: InvertibleGroup>(
  base: &G::Elem,
  exp: &BigInt,
  result: &G::Elem,
) -> PoKE2<G::Elem> {
  let g = G::base_elem();
  let z = G::exp_signed(&g, exp);
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let q = exp / BigInt::from_biguint(Plus, l.clone());
  // Should be a safe unwrap because Rem is implemented for BigInt as
  // `BigInt::from_biguint(self.sign, self.data % other)
  let r = (exp % BigInt::from_biguint(Plus, l)).to_biguint().unwrap();
  PoKE2 {
    z,
    q: G::exp_signed(&G::op(&base, &G::exp(&g, &alpha)), &q),
    r,
  }
}

/// See page 16 of B&B.
pub fn verify_poke2<G: Group>(base: &G::Elem, result: &G::Elem, proof: &PoKE2<G::Elem>) -> bool {
  let PoKE2 { z, q, r } = proof;
  let g = G::base_elem();
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let lhs = G::op(
    &G::exp(q, &l),
    &G::exp(&G::op(&base, &G::exp(&g, &alpha)), &r),
  );
  let rhs = G::op(result, &G::exp(&z, &alpha));
  lhs == rhs
}

fn hash_prime<G: Serialize>(_u: &G, _w: &G, _z: &G) -> BigUint {
  // TODO: replace with commented out when hash_prime is implemented
  BigUint::from(13 as u8)
  // let mut hash_string = serde_json::to_string(&u).unwrap();
  // hash_string.push_str(&serde_json::to_string(&w).unwrap());
  // hash_string.push_str(&serde_json::to_string(&z).unwrap());
  // hashes::h_prime(&hashes::blake2, hash_string.as_bytes())
}

fn hash_inputs<G: Serialize>(u: &G, w: &G, z: &G, l: &BigUint) -> BigUint {
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hash_string.push_str(&l.to_str_radix(16));
  hashes::blake2(hash_string.as_bytes(), None)
}

#[cfg(test)]
mod tests {
  use super::super::super::group::dummy;
  use super::*;
  #[test]
  fn test_poke2() {
    let dummy = dummy::DummyRSA::get();
    let base = dummy::DummyRSA::base_elem_(&dummy);
    let exp = BigInt::from(20 as u8);
    let result = dummy::DummyRSA::op(&1_048_576, &1);
    let proof = prove_poke2::<dummy::DummyRSA>(&base, &exp, &result);
    // 2^20 = 1048576
    assert!(verify_poke2::<dummy::DummyRSA>(&base, &result, &proof));
    let exp_2 = BigInt::from(35 as u8);
    let result_2 = dummy::DummyRSA::op(&34_359_738_368, &1);
    let proof_2 = prove_poke2::<dummy::DummyRSA>(&base, &exp_2, &result_2);
    // 2^35 = 34359738368
    assert!(verify_poke2::<dummy::DummyRSA>(&base, &result_2, &proof_2));
    // Cannot verify wrong base/exp/result triple with wrong pair
    assert!(!verify_poke2::<dummy::DummyRSA>(&base, &result_2, &proof));
  }
}
