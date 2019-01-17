// NOTE: Remove underscores on variable names when implementing!

use super::group::Generator;
use super::group::Pow;
use super::proof::poe;
use super::proof::PoE;
use super::proof::PoKE2;
use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigUint;
use serde::ser::Serialize;

fn product(elems: &[BigUint]) -> BigUint {
  elems.iter().fold(BigUint::new(vec![1]), |a, b| a * b)
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

pub fn delete<O, G: AbstractGroup<O> + Pow<O>>(
  _acc: &G,
  _elem_witnesses: &[(BigUint, G)],
) -> Option<(G, PoE<G>)>
where
  O: Operator,
{
  unimplemented!()
}

pub fn prove_membership<O, G: AbstractGroup<O> + Pow<O>>(
  acc: &G,
  elem_witnesses: &[(BigUint, G)],
) -> Option<(G, PoE<G>)>
where
  O: Operator,
{
  delete(acc, elem_witnesses)
}

pub fn verify_membership<O, G: AbstractGroup<O> + Pow<O>>(
  _witness: &G,
  _elems: &[BigUint],
  _result: &G,
  _proof: &PoE<G>,
) -> Option<(G, G, PoKE2<G>, PoE<G>)>
where
  O: Operator,
{
  unimplemented!()
}

pub fn prove_nonmembership<O, G: AbstractGroup<O> + Pow<O>>(_acc: &G, _elems: &[BigUint]) -> bool
where
  O: Operator,
{
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
