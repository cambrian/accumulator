use super::super::group::{Group, InvertibleGroup};
use super::super::hash::hashes;
use super::PoKE2;
use num::BigInt;
use num::BigUint;
use num_bigint::Sign::Plus;

/// See page 16 of B&B.
pub fn prove_poke2<G: InvertibleGroup>(
  base: &G::Elem,
  exp: &BigInt,
  result: &G::Elem,
) -> PoKE2<G::Elem> {
  let g = G::base_elem();
  let z = G::exp_signed(&g, exp);
  let l = hash_prime::<G>(base, result, &z);
  let alpha = hash_inputs::<G>(base, result, &z, &l);
  let q = exp / BigInt::from_biguint(Plus, l.clone());
  let r = exp % BigInt::from_biguint(Plus, l);
  PoKE2 {
    z,
    q: G::exp_signed(&G::op(&base, &G::exp(&g, &alpha)), &q),
    r,
  }
}

/// See page 16 of B&B.
pub fn verify_poke2<G: InvertibleGroup>(
  base: &G::Elem,
  result: &G::Elem,
  proof: &PoKE2<G::Elem>,
) -> bool {
  let PoKE2 { z, q, r } = proof;
  let g = G::base_elem();
  let l = hash_prime::<G>(base, result, &z);
  let alpha = hash_inputs::<G>(base, result, &z, &l);
  let lhs = G::op(&G::exp(q, &l), &G::exp_signed(&G::op(&base, &G::exp(&g, &alpha)), &r));
  let rhs = G::op(result, &G::exp(&z, &alpha));
  lhs == rhs
}

/// Review: simplify Group element requirement fromt u,w,z to Serialize, and remove unnecessary
/// group qualifications from callsites
fn hash_prime<G: Group>(u: &G::Elem, w: &G::Elem, z: &G::Elem) -> BigUint {
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hashes::blake2_prime(hash_string.as_bytes())
}

/// Review: simplify Group element requirement fromt u,w,z to Serialize, and remove unnecessary
/// group qualifications from callsites
fn hash_inputs<G: Group>(u: &G::Elem, w: &G::Elem, z: &G::Elem, l: &BigUint) -> BigUint {
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hash_string.push_str(&l.to_str_radix(16));
  hashes::blake2(hash_string.as_bytes(), None)
}
