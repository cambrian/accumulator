use super::super::group::Generator;
use super::super::group::Inverse;
use super::super::hash::hashes;
use super::PoKE2;
use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigInt;
use num::BigUint;
use num_bigint::Sign::Plus;
use serde::ser::Serialize;

pub fn compute_poke2<O, G: AbstractGroup<O> + Generator<O> + Inverse<O> + Serialize>(
  base: &G,
  exp: &BigInt,
  result: &G,
) -> PoKE2<G>
where
  O: Operator,
{
  let g = G::generator();
  let z = g.pow_signed(exp);
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let q = exp / BigInt::from_biguint(Plus, l.clone());
  let r = exp % BigInt::from_biguint(Plus, l);
  PoKE2 {
    z,
    q: base.operate(&g.pow(&alpha)).pow_signed(&q),
    r,
  }
}

pub fn verify_poke2<O, G: AbstractGroup<O> + Generator<O> + Inverse<O> + Serialize>(
  base: &G,
  result: &G,
  proof: &PoKE2<G>,
) -> bool
where
  O: Operator,
{
  let PoKE2 { z, q, r } = proof;
  let g = G::generator();
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let lhs = q
    .pow(&l)
    .operate(&(base.operate(&g.pow(&alpha))).pow_signed(&r));
  let rhs = result.operate(&z.pow(&alpha));
  lhs == rhs
}

fn hash_prime<O, G: AbstractGroup<O> + Serialize>(u: &G, w: &G, z: &G) -> BigUint
where
  O: Operator,
{
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hashes::blake2_prime(hash_string.as_bytes())
}

fn hash_inputs<O, G: AbstractGroup<O> + Serialize>(u: &G, w: &G, z: &G, l: &BigUint) -> BigUint
where
  O: Operator,
{
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hash_string.push_str(&l.to_str_radix(16));
  hashes::blake2(hash_string.as_bytes(), None)
}
