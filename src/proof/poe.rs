use super::super::group::Group;
use super::super::hash::hashes;
use super::PoE;
use num::BigUint;
use serde::ser::Serialize;

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

fn hash_prime<G: Serialize>(exp: &BigUint, base: &G, result: &G) -> BigUint {
  let mut hash_string = exp.to_str_radix(16);
  hash_string.push_str(&serde_json::to_string(&base).unwrap());
  hash_string.push_str(&serde_json::to_string(&result).unwrap());
  hashes::h_prime(&hashes::blake2, hash_string.as_bytes())
}
