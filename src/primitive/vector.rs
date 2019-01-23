// TODO (Also, how to aggregate?)
use super::super::hash::{hash_to_prime, Blake2b};
use super::accumulator;
use super::accumulator::AccError;
use crate::group::UnknownOrderGroup;
use crate::proof::{PoE, PoKE2};
use bitvec::BitVec;
use num::BigUint;

#[derive(Debug)]
pub enum UpdateResult<G: UnknownOrderGroup> {
  NoChange(G::Elem),
  Update((MembershipProof<G>, MembershipProof<G>)),
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
#[derive(Debug)]
pub struct MembershipProof<G: UnknownOrderGroup> {
  acc: G::Elem,
  proof: PoE<G>,
}

#[allow(dead_code)]
struct NonMembershipProof<G: UnknownOrderGroup> {
  d: G::Elem,
  v: G::Elem,
  gv_inverse: G::Elem,
  poke_proof: PoKE2<G>,
  poe_proof: PoE<G>,
}

pub fn setup<G: UnknownOrderGroup>() -> G::Elem {
  G::unknown_order_elem()
}

// TODO 1. Somehow check if element is already in accumulator
// Find a better way to pass a reference to vector of references, or remove this in accumulator.rs
// TODO Option type to allow for only deletion or only addition
pub fn update<G: UnknownOrderGroup>(
  acc: &G::Elem,
  bits: &BitVec,
  indices: &[&BigUint],
  del_witnesses: &[&G::Elem],
) -> Result<UpdateResult<G>, AccError> {
  // Must hold hash commitments in vec in order to pass by reference to accumulator fns
  let mut add_commitments = vec![];
  let mut del_commitments = vec![];
  for i in 0..bits.len() {
    if bits[i] {
      add_commitments.push(hash_to_prime(&Blake2b::default, indices[i]));
    } else {
      del_commitments.push(hash_to_prime(&Blake2b::default, indices[i]));
    }
  }
  let mut add_commitments_ref = vec![];
  for item in &add_commitments {
    add_commitments_ref.push(item);
  }
  let mut del_commitments_ref = vec![];
  for i in 0..del_commitments.len() {
    del_commitments_ref.push((&del_commitments[i], del_witnesses[i]));
  }
  let (added_acc, add_poe) = accumulator::add::<G>(acc, &add_commitments_ref);
  let add_proof = MembershipProof {
    acc: added_acc,
    proof: add_poe,
  };
  match accumulator::delete::<G>(acc, &del_commitments_ref) {
    Ok((del_acc, del_poe)) => {
      let del_proof = MembershipProof {
        acc: del_acc,
        proof: del_poe,
      };
      Ok(UpdateResult::Update((add_proof, del_proof)))
    }
    Err(n) => Err(n),
  }
}

pub fn open<G: UnknownOrderGroup>(
  _bits: &BitVec,
  _indices: &[BigUint],
) -> Result<Vec<Proof>, OpenError> {
  // let mut one_commitment = BigUint::One;
  // let mut zero_commitment;
  // for i in 0..bits.len() {
  //   if bits[i] {
  //     one_commitment *= h_prime(&blake2, indices[i].to_str_radix(16).as_bytes());
  //   } else {
  //     zero_commitment *= h_prime(&blake2, indices[i].to_str_radix(16).as_bytes());
  //   }
  // }
  // inclusion_witnesses.iter().map(|(index, witness)| { index = h_prime(&blake2, index.to_str_radix(16).as_bytes(); (index, witness) });
  // match accumulator::prove_membership(acc, inclusion_witnesses) {
  //   Ok(inclusion_proof) => {
  //     match accumulator::prove_nonmembership(acc, acc_set, )
  //   }
  //   Err(n) => Err(n),
  // }
  unimplemented!();
}

pub fn verify<G: UnknownOrderGroup>(
  _acc: &G::Elem,
  _bits: &[bool],
  _indices: &[BigUint],
  _proofs: Vec<Proof>,
) -> bool {
  unimplemented!();
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::DummyRSA;
  use crate::util::bu;

  #[test]
  fn test_update() {
    let vc = setup::<DummyRSA>();
    let mut bv: BitVec = BitVec::new();
    bv.push(true);
    let proofs = update::<DummyRSA>(&vc, &bv, &[&bu(2u8)], &[&DummyRSA::unknown_order_elem()]);
    println!("{:?}", proofs);
    //    G: UnknownOrderGroup>(
    //   acc: &G::Elem,
    //   bits: &BitVec,
    //   indices: &[&BigUint],
    //   del_witnesses: &[&G::Elem],
    // )
  }
}
