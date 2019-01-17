use super::super::group::{CyclicGroup, Group, InvertibleGroup};
use super::super::hash::hashes;
use super::PoKE2;
use num::BigInt;
use num::BigUint;
use num_bigint::Sign::Plus;
use serde::ser::Serialize;

/// See page 16 of B&B.
pub fn prove_poke2<G: InvertibleGroup + CyclicGroup + Serialize>(
  base: &G,
  exp: &BigInt,
  result: &G,
) -> PoKE2<G> {
  let g = G::generator();
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
pub fn verify_poke2<G: InvertibleGroup + CyclicGroup + Serialize>(
  base: &G,
  result: &G,
  proof: &PoKE2<G>,
) -> bool {
  let PoKE2 { z, q, r } = proof;
  let g = G::generator();
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let lhs = q.exp(&l).op(&(base.op(&g.exp(&alpha))).exp_signed(&r));
  let rhs = result.op(&z.exp(&alpha));
  lhs == rhs
}

fn hash_prime<G: Group + Serialize>(u: &G, w: &G, z: &G) -> BigUint {
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hashes::blake2_prime(hash_string.as_bytes())
}

fn hash_inputs<G: Group + Serialize>(u: &G, w: &G, z: &G, l: &BigUint) -> BigUint {
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hash_string.push_str(&l.to_str_radix(16));
  hashes::blake2(hash_string.as_bytes(), None)
}
