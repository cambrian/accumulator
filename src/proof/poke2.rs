// TODO
use super::super::group::Generator;
use super::super::group::Pow;
// use super::super::hash::hashes;
use super::PoKE2;
use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigUint;
use serde::ser::Serialize;

pub fn compute_poke2<O, G: AbstractGroup<O> + Generator<O> + Pow<O> + Serialize>(
  base: &G,
  exp: &BigUint,
  result: &G,
) -> PoKE2<G>
where
  O: Operator,
{
  let g = G::generator();
  let z = g.pow(exp);
  let l = hash_prime(base, result, &z);
  let alpha = hash_inputs(base, result, &z, &l);
  let q = exp / l.clone();
  let r = exp % l;
  PoKE2 {
    z,
    q: base.operate(&g.pow(&alpha)).pow(&q),
    r,
  }
}

pub fn verify_poke2<O, G: AbstractGroup<O> + Generator<O> + Pow<O> + Serialize>(
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
  let lhs = q.pow(&l).operate(&(base.operate(&g.pow(&alpha))).pow(&r));
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
  let mut target = vec![0u8; 32];
  // TODO: Use HPrime function when defined
  // let _ = hashes::blake2(hash_string.as_bytes(), None).to_big_endian(&mut target);
  BigUint::from_bytes_be(&target[..])
}

fn hash_inputs<O, G: AbstractGroup<O> + Serialize>(u: &G, w: &G, z: &G, l: &BigUint) -> BigUint
where
  O: Operator,
{
  let mut hash_string = serde_json::to_string(&u).unwrap();
  hash_string.push_str(&serde_json::to_string(&w).unwrap());
  hash_string.push_str(&serde_json::to_string(&z).unwrap());
  hash_string.push_str(&l.to_str_radix(16));
  let mut target = vec![0u8; 32];
  // let _ = hashes::blake2(hash_string.as_bytes(), None).to_big_endian(&mut target);
  BigUint::from_bytes_be(&target[..])
}
