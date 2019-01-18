use super::super::group::{Group, InvertibleGroup, GroupElem, InvertibleGroupElem};
use super::super::hash::hashes;
use super::PoKE2;
use num::BigInt;
use num::BigUint;
use num_bigint::Sign::Plus;

/// REVIEW: rename to prove_poke2
/// See page 16 of B&B.
pub fn compute_poke2<G: InvertibleGroup>(
  base: &InvertibleGroupElem<G>,
  exp: &BigInt,
  result: &InvertibleGroupElem<G>,
) -> PoKE2<InvertibleGroupElem<G>> {
  let g = G::base_elem();
  let z = g.exp_signed(exp);
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let q = exp / BigInt::from_biguint(Plus, l.clone());
  let r = exp % BigInt::from_biguint(Plus, l);
  PoKE2 {
    z,
    q: base.op(&g.exp(&alpha)).exp_signed(&q),
    r,
  }
}

/// See page 16 of B&B.
pub fn verify_poke2<G: InvertibleGroup>(
  base: &InvertibleGroupElem<G>,
  result: &InvertibleGroupElem<G>,
  proof: &PoKE2<InvertibleGroupElem<G>>,
) -> bool
{
  let PoKE2 { z, q, r } = proof;
  let g = G::base_elem();
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let lhs = q
    .exp(&l)
    .op(&(base.op(&g.exp(&alpha))).exp_signed(&r));
  let rhs = result.op(&z.exp(&alpha));
  lhs == rhs
}

fn hash_prime<G: Group>(u: &GroupElem<G>, w: &GroupElem<G>, z: &GroupElem<G>) -> BigUint
{
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hashes::blake2_prime(hash_string.as_bytes())
}

fn hash_inputs<G: Group>(u: &GroupElem<G>, w: &GroupElem<G>, z: &GroupElem<G>, l: &BigUint) -> BigUint
{
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hash_string.push_str(&l.to_str_radix(16));
  hashes::blake2(hash_string.as_bytes(), None)
}
