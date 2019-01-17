use super::super::group::Group;
use super::super::hash::hashes;
use super::PoE;
use num::BigUint;
use serde::ser::Serialize;

/// See page 16 of B&B.
pub fn prove_poe<G: Group + Serialize>(base: &G, exp: &BigUint, result: &G) -> PoE<G> {
  let l = hash_prime(exp, base, result);
  let q = exp / l;
  PoE { q: base.exp(&q) }
}

/// See page 16 of B&B.
pub fn verify_poe<G: Group + Serialize>(
  base: &G,
  exp: &BigUint,
  result: &G,
  proof: &PoE<G>,
) -> bool {
  let l = hash_prime(exp, base, result);
  let r = exp % l.clone();
  // w = Q^l * u^r
  let w = proof.q.exp(&l).op(&base.exp(&r));
  w == *result
}

fn hash_prime<G: Group + Serialize>(exp: &BigUint, base: &G, result: &G) -> BigUint {
  let mut hash_string = exp.to_str_radix(16);
  hash_string.push_str(&serde_json::to_string(&base).unwrap());
  hash_string.push_str(&serde_json::to_string(&result).unwrap());
  hashes::blake2_prime(hash_string.as_bytes())
}
