// TODO
use super::super::group::Pow;
use super::super::hash::hashes;
use super::PoE;
use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigUint;
use serde::ser::Serialize;

pub fn compute_poe<O, G: AbstractGroup<O> + Pow<O> + Serialize>(
  base: &G,
  exp: &BigUint,
  result: &G,
) -> PoE<G>
where
  O: Operator,
{
  let l = hash_prime(exp, base, result);
  let q = exp / l;
  PoE { q: base.pow(&q) }
}

pub fn verify_poe<O, G: AbstractGroup<O> + Pow<O> + Serialize>(
  base: &G,
  exp: &BigUint,
  result: &G,
  proof: &PoE<G>,
) -> bool
where
  O: Operator,
{
  let l = hash_prime(exp, base, result);
  let r = exp % l.clone();
  // w = Q^l * u^r
  let w = proof.q.pow(&l).operate(&base.pow(&r));
  w == *result
}

fn hash_prime<O, G: AbstractGroup<O> + Serialize>(exp: &BigUint, base: &G, result: &G) -> BigUint
where
  O: Operator,
{
  let mut hash_string = exp.to_str_radix(16);
  hash_string.push_str(&serde_json::to_string(&base).unwrap());
  hash_string.push_str(&serde_json::to_string(&result).unwrap());
  let mut target = vec![0u8; 32];
  // TODO: Use HPrime function when defined
  let _ = hashes::blake2(hash_string.as_bytes(), None).to_big_endian(&mut target);
  BigUint::from_bytes_be(&target[..])
}
