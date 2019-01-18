use super::super::group::{Group, GroupElem};
use super::super::hash::hashes;
use super::PoE;
use num::BigUint;

/// See page 16 of B&B.
/// REVIEW: rename to prove_poe
pub fn compute_poe<G: Group>(
  base: &GroupElem<G>,
  exp: &BigUint,
  result: &GroupElem<G>,
) -> PoE<GroupElem<G>> {
  let l = hash_prime(exp, base, result);
  let q = exp / l;
  PoE { q: base.exp(&q) }
}

/// See page 16 of B&B.
pub fn verify_poe<G: Group>(
  base: &GroupElem<G>,
  exp: &BigUint,
  result: &GroupElem<G>,
  proof: &PoE<GroupElem<G>>,
) -> bool {
  let l = hash_prime(exp, base, result);
  let r = exp % l.clone();
  // w = Q^l * u^r
  let w = proof.q.exp(&l).op(&base.exp(&r));
  w == *result
}

fn hash_prime<G: Group>(exp: &BigUint, base: &GroupElem<G>, result: &GroupElem<G>) -> BigUint
{
  let mut hash_string = exp.to_str_radix(16);
  hash_string.push_str(&serde_json::to_string(&base).unwrap());
  hash_string.push_str(&serde_json::to_string(&result).unwrap());
  hashes::blake2_prime(hash_string.as_bytes())
}
