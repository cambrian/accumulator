// TODO
use super::super::group::Pow;
use super::super::hash;
use super::PoE;
use alga::general::AbstractGroup;
use alga::general::Multiplicative;
use alga::general::Operator;
use bigint::uint::U256;
use num::BigUint;
use serde_json;

// x, u, w: u^x = w
pub fn compute_poe<O, G: AbstractGroup<O> + Pow<O>>(base: &G, exp: &BigUint, result: &G) -> PoE<G>
where
  O: Operator,
{
  // TODO: Use HPrime function when defined
  // let z = hash::blake2(base, None);
  let l = BigUint::from(1 as u8); // HPrime(exp, base, result);
  let q = exp / l;
  PoE { q: base.pow(q) }
}

pub fn verify_poe<O, G: AbstractGroup<O> + Pow<O>>(
  base: &G,
  exp: &BigUint,
  result: &G,
  proof: &PoE<G>,
) -> bool
where
  O: Operator,
{
  // TODO: Use HPrime function when defined
  let l = BigUint::from(1 as u8); // HPrime(exp, base, result);
  let r = exp % l;
  let l = BigUint::from(23 as u8); // BigInt::from(l as U256);
  let w = proof.q.pow(l);
  // let w = w.op(Multiplicative, &base.pow(l));
  false
}
