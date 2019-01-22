// TODO
use super::accumulator::AccError;
use super::group::Group;
use super::proof::{poe::PoE, poke2::PoKE2};
use bitvec::BitVec;
use num::BigUint;

pub enum UpdateResult<G: Group> {
  NoChange(G::Elem),
  Update(G::Elem, PoE<G>),
}

// TODO enumerate error types
pub enum OpenError {
  Error,
}

pub enum Proof {
  MembershipProof,
  NonMembershipProof,
}

struct MembershipProof<G: Group> {
  acc: G::Elem,
  proof: PoE<G>,
}

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

pub fn commit<G: Group>(m: BitVec) -> G::Elem {
  G::base_elem()
}

pub fn update<G: Group>(bits: &[bool], indices: &[BigUint]) -> Result<UpdateResult<G>, AccError> {
  unimplemented!();
}

pub fn open<G: Group>(bits: &[bool], indices: &[BigUint]) -> Result<Vec<Proof>, OpenError> {
  unimplemented!();
}

pub fn verify<G: Group>(
  acc: &G::Elem,
  bits: &[bool],
  indices: &[BigUint],
  proofs: Vec<Proof>,
) -> bool {
  unimplemented!();
}
