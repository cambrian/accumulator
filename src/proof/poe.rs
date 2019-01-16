// TODO
use super::super::group::Pow;
use super::PoE;
use alga::general::AbstractGroup;
use alga::general::Operator;
use bigint::uint::U256;
use num::BigUint;

// x, u, w: u^x = w
pub fn compute_poe<O, G: AbstractGroup<O> + Pow<O>>(base: G, exp: U256, result: G) -> PoE<G>
where
  O: Operator,
{
  // TODO: Use HPrime function when defined
  let l = U256::from(1); // HPrime(exp, base, result);
  let q = exp / l;
  // TODO: Convert q U256 result into BigUint
  let q = BigUint::from(23 as u8);
  PoE { q: base.pow(q) }
}

pub fn verify_poe<O, G: AbstractGroup<O> + Pow<O>>(
  base: G,
  exp: U256,
  result: G,
  proof: PoE<G>,
) -> bool
where
  O: Operator,
{
  // TODO: Use HPrime function when defined
  let l = U256::from(1); // HPrime(exp, base, result);
  let r = exp % l;
  let l = BigUint::from(23 as u8); // BigInt::from(l as U256);
  let w = proof.q.pow(l);
  // let w = w.mul(base.pow(r));
  w == result
}

// Prove (x, u, w)
// l <- Hprime (x, u, w)
// q <- Int divide x / l
// r <- x mod l
// return PoE <- u^q

// Verify (x, u, w, Q)
// l <- Hprime (x, u, w)
// r <- x mod l
// Check Q^lu^r = w
// 1 if yes, else 0
