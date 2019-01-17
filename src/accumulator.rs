// NOTE: Remove underscores on variable names when implementing!
use super::group::{Generator, Inverse, Pow};
use super::proof::{poe, poke2, PoE, PoKE2};
use alga::general::{AbstractGroup, Operator};
use num;
use num::BigInt;
use num::BigUint;
use num_bigint::Sign::Plus;
use num_traits::identities::One;
use num_traits::identities::Zero;
use serde::ser::Serialize;

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
fn shamir_trick<O, G: AbstractGroup<O> + Inverse<O>>(
  xth_root: &G,
  yth_root: &G,
  x: &BigUint,
  y: &BigUint,
) -> Option<G>
where
  O: Operator,
{
  if xth_root.pow(&x) != yth_root.pow(&y) {
    return None;
  }

  let (a, b, gcd) = bezout(&x, &y);

  if !gcd.is_one() {
    return None;
  }

  Some(xth_root.pow_signed(&b).operate(&yth_root.pow_signed(&a)))
}

pub fn setup<O, G: AbstractGroup<O> + Generator<O>>() -> G
where
  O: Operator,
{
  G::generator()
}

pub fn add<O, G: AbstractGroup<O> + Pow<O> + Serialize>(acc: &G, elems: &[BigUint]) -> (G, PoE<G>)
where
  O: Operator,
{
  let x = product(elems);
  let new_acc = acc.pow(&x);
  let poe_proof = poe::compute_poe(acc, &x, &new_acc);
  (new_acc, poe_proof)
}

pub fn delete<O, G: AbstractGroup<O> + Inverse<O> + Serialize>(
  acc: &G,
  elem_witnesses: &[(BigUint, G)],
) -> Option<(G, PoE<G>)>
where
  O: Operator,
{
  if elem_witnesses.is_empty() {
    return None;
  }

  let (mut elem_aggregate, mut acc_next) = (&elem_witnesses[0]).clone();

  for (elem, witness) in elem_witnesses
    .split_first() // Chop off first entry.
    .expect("unexpected witnesses")
    .1
  {
    if &witness.pow(elem) != acc {
      return None;
    }

    let acc_next_option = shamir_trick(&acc_next, witness, &elem_aggregate, elem);
    match acc_next_option {
      Some(acc_next_value) => acc_next = acc_next_value,
      None => return None,
    };

    elem_aggregate = elem_aggregate * elem;
  }

  let poe_proof = poe::compute_poe(&acc_next, &elem_aggregate, acc);
  Some((acc_next, poe_proof))
}

pub fn prove_membership<O, G: AbstractGroup<O> + Inverse<O> + Serialize>(
  acc: &G,
  elem_witnesses: &[(BigUint, G)],
) -> Option<(G, PoE<G>)>
where
  O: Operator,
{
  delete(acc, elem_witnesses)
}

pub fn verify_membership<O, G: AbstractGroup<O> + Pow<O> + Serialize>(
  witness: &G,
  elems: &[BigUint],
  result: &G,
  proof: &PoE<G>,
) -> bool
where
  O: Operator,
{
  let exp = product(elems);
  poe::verify_poe(witness, &exp, result, proof)
}

pub fn prove_nonmembership<O, G: AbstractGroup<O> + Generator<O> + Inverse<O> + Serialize>(
  acc: &G,
  acc_set: &[BigUint],
  elems: &[BigUint],
) -> Option<(G, G, G, PoKE2<G>, PoE<G>)>
where
  O: Operator,
{
  let x = product(elems);
  let s = product(acc_set);
  let (a, b, gcd) = bezout(&x, &s);

  if !gcd.is_one() {
    return None;
  }

  let g = G::generator();
  let d = g.pow_signed(&a);
  let v = acc.pow_signed(&b);
  let gv_inverse = g.operate(&v.efficient_inverse(&num::one()));

  // let poke2_proof = poke2::compute_poke2(acc, &b, &v);
  let poe_proof = poe::compute_poe(&d, &x, &gv_inverse);
  unimplemented!()
}

pub fn verify_nonmembership<O, G: AbstractGroup<O> + Pow<O>>(
  _acc: &G,
  _elems: &[BigUint],
  _d: &G,
  _v: &G,
  _poke_proof: &PoKE2<G>,
  _poe_proof: &PoE<G>,
) -> bool
where
  O: Operator,
{
  unimplemented!()
}
