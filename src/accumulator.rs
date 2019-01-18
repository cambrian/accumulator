use super::group::{CyclicGroup, Group, InvertibleGroup};
use super::proof::{poe, poke2, PoE, PoKE2};
use num;
use num::BigInt;
use num::BigUint;
use num_bigint::Sign::Plus;
use num_traits::identities::One;
use num_traits::identities::Zero;
use serde::ser::Serialize;

#[derive(Debug)]
pub enum AccError {
  BadWitness,
  InputsNotCoPrime,
}

/// Initializes the accumulator to a group element.
pub fn setup<G: CyclicGroup>() -> G {
  G::generator()
}

/// Adds `elems` to the accumulator `acc`.
pub fn add<G: Group + Serialize>(acc: &G, elems: &[BigUint]) -> (G, PoE<G>) {
  let x = product(elems);
  let new_acc = acc.exp(&x);
  let poe_proof = poe::prove_poe(acc, &x, &new_acc);
  (new_acc, poe_proof)
}

/// Removes the elements in `elem_witnesses` from the accumulator `acc`.
pub fn delete<G: InvertibleGroup + Serialize>(
  acc: &G,
  elem_witnesses: &[(BigUint, G)],
) -> Result<(G, PoE<G>), AccError> {
  if elem_witnesses.is_empty() {
    let poe_proof = poe::prove_poe(acc, &BigUint::zero(), acc);
    return Ok((acc.clone(), poe_proof));
  }

  let (mut elem_aggregate, mut acc_next) = elem_witnesses[0].clone();

  for (elem, witness) in elem_witnesses
    .split_first() // Chop off first entry.
    .expect("unexpected witnesses")
    .1
  {
    if &witness.exp(elem) != acc {
      return Err(AccError::BadWitness);
    }

    let acc_next_option = shamir_trick(&acc_next, witness, &elem_aggregate, elem);
    match acc_next_option {
      Some(acc_next_value) => acc_next = acc_next_value,
      None => return Err(AccError::InputsNotCoPrime),
    };

    elem_aggregate = elem_aggregate * elem;
  }

  let poe_proof = poe::prove_poe(&acc_next, &elem_aggregate, acc);
  Ok((acc_next, poe_proof))
}

/// See `delete`.
pub fn prove_membership<G: InvertibleGroup + Serialize>(
  acc: &G,
  elem_witnesses: &[(BigUint, G)],
) -> Result<(G, PoE<G>), AccError> {
  delete(acc, elem_witnesses)
}

/// Verifies the PoE returned by `prove_membership` s.t. `witness` ^ `elems` = `result`.
pub fn verify_membership<G: Group + Serialize>(
  witness: &G,
  elems: &[BigUint],
  result: &G,
  proof: &PoE<G>,
) -> bool {
  let exp = product(elems);
  poe::verify_poe(witness, &exp, result, proof)
}

/// Returns a proof (and associated variables) that `elems` are not in `acc_set`.
pub fn prove_nonmembership<G: InvertibleGroup + CyclicGroup + Serialize>(
  acc: &G,
  acc_set: &[BigUint],
  elems: &[BigUint],
) -> Result<(G, G, G, PoKE2<G>, PoE<G>), AccError> {
  let x = product(elems);
  let s = product(acc_set);
  let (a, b, gcd) = bezout(&x, &s);

  if !gcd.is_one() {
    return Err(AccError::InputsNotCoPrime);
  }

  let g = G::generator();
  let d = g.exp_signed(&a);
  let v = acc.exp_signed(&b);
  let gv_inverse = g.op(&v.inv());

  let poke2_proof = poke2::prove_poke2(acc, &b, &v);
  let poe_proof = poe::prove_poe(&d, &x, &gv_inverse);
  Ok((d, v, gv_inverse, poke2_proof, poe_proof))
}

/// Verifies the PoKE2 and PoE returned by `prove_nonmembership`.
pub fn verify_nonmembership<G: InvertibleGroup + CyclicGroup + Serialize>(
  acc: &G,
  elems: &[BigUint],
  d: &G,
  v: &G,
  gv_inverse: &G,
  poke2_proof: &PoKE2<G>,
  poe_proof: &PoE<G>,
) -> bool {
  let x = product(elems);
  poke2::verify_poke2(acc, v, poke2_proof) && poe::verify_poe(d, &x, gv_inverse, poe_proof)
}

/// Returns (a, b, GCD(x, y)).
fn bezout(x: &BigUint, y: &BigUint) -> (BigInt, BigInt, BigInt) {
  let (mut s, mut old_s): (BigInt, BigInt) = (num::zero(), num::one());
  let (mut t, mut old_t): (BigInt, BigInt) = (num::one(), num::zero());
  let (mut r, mut old_r) = (
    BigInt::from_biguint(Plus, y.clone()),
    BigInt::from_biguint(Plus, x.clone()),
  );

  while !r.is_zero() {
    let quotient = &old_r / &r;
    let (temp_r, temp_s, temp_t) = (r, s, t);

    r = &old_r - &quotient * &temp_r;
    s = &old_s - &quotient * &temp_s;
    t = &old_t - &quotient * &temp_t;

    old_r = temp_r;
    old_s = temp_s;
    old_t = temp_t;
  }

  (old_s, old_t, old_r)
}

fn product(elems: &[BigUint]) -> BigUint {
  elems.iter().fold(num::one(), |a, b| a * b)
}

/// Computes the (xy)th root of g given the xth and yth roots of g and (x, y) coprime).
fn shamir_trick<G: InvertibleGroup>(
  xth_root: &G,
  yth_root: &G,
  x: &BigUint,
  y: &BigUint,
) -> Option<G> {
  if xth_root.exp(&x) != yth_root.exp(&y) {
    return None;
  }

  let (a, b, gcd) = bezout(&x, &y);

  if !gcd.is_one() {
    return None;
  }

  Some(xth_root.exp_signed(&b).op(&yth_root.exp_signed(&a)))
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_bezout_simple() {
    let x = BigUint::from(7 as u16);
    let y = BigUint::from(165 as u16);
    let (a, b, gcd) = bezout(&x, &y);
    assert!(gcd.is_one());
    assert!(a == BigInt::from(-47 as i16));
    assert!(b == BigInt::from(2 as i16));
  }
}
