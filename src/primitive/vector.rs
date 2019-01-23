// TODO (Also, how to aggregate?)
use super::accumulator::AccError;
use crate::group::Group;
use crate::proof::{poe::PoE, poke2::PoKE2};
use bitvec::BitVec;
use num::BigUint;

pub enum UpdateResult<G: Group> {
  NoChange(G::Elem),
  Update(G::Elem, PoE<G>),
}

// TODO: Enumerate error types.
pub enum OpenError {
  Error,
}

pub enum Proof {
  MembershipProof,
  NonMembershipProof,
}

#[allow(dead_code)]
struct MembershipProof<G: Group> {
  acc: G::Elem,
  proof: PoE<G>,
}

#[allow(dead_code)]
struct NonMembershipProof<G: Group> {
  d: G::Elem,
  v: G::Elem,
  gv_inverse: G::Elem,
  poke_proof: PoKE2<G>,
  poe_proof: PoE<G>,
}

pub fn setup<G: Group>() -> G::Elem {
  G::base_elem()
}

pub fn commit<G: Group>(_m: BitVec) -> G::Elem {
  G::base_elem()
}

pub fn update<G: Group>(_bits: &[bool], _indices: &[BigUint]) -> Result<UpdateResult<G>, AccError> {
  unimplemented!();
}

pub fn open<G: Group>(_bits: &[bool], _indices: &[BigUint]) -> Result<Vec<Proof>, OpenError> {
  unimplemented!();
}

pub fn verify<G: Group>(
  _acc: &G::Elem,
  _bits: &[bool],
  _indices: &[BigUint],
  _proofs: Vec<Proof>,
) -> bool {
  unimplemented!();
}
