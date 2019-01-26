// TODO: Reconsider when/where to group elems by bit.
use super::accumulator;
use super::accumulator::{MembershipProof, NonmembershipProof};
use crate::group::UnknownOrderGroup;
use crate::hash::{hash_to_prime, Blake2b};
use rug::Integer;
use std::collections::HashSet;

#[derive(Debug)]
pub enum VCError {
  ConflictingIndicesError,
  InvalidOpenError,
}

pub struct VectorProof<G: UnknownOrderGroup> {
  membership_proof: MembershipProof<G>,
  nonmembership_proof: NonmembershipProof<G>,
}

pub fn setup<G: UnknownOrderGroup>() -> G::Elem {
  accumulator::setup::<G>()
}

fn group_elems_by_bit(bits: &[(bool, Integer)]) -> Result<(Vec<Integer>, Vec<Integer>), VCError> {
  let mut elems_with_one = vec![];
  let mut elems_with_zero = vec![];
  let mut seen_indices = HashSet::new();
  for (bit, i) in bits {
    if !seen_indices.insert(i) {
      return Err(VCError::ConflictingIndicesError);
    }
    if *bit {
      elems_with_one.push(hash_to_prime(&Blake2b::default, &i));
    } else {
      elems_with_zero.push(hash_to_prime(&Blake2b::default, &i));
    }
  }
  Ok((elems_with_zero, elems_with_one))
}

pub fn update<G: UnknownOrderGroup>(
  acc: G::Elem,
  acc_set: &[Integer],
  bits: &[(bool, Integer)],
) -> Result<(G::Elem, VectorProof<G>), VCError> {
  // Must hold hash commitments in vec in order to pass by reference to accumulator fns.
  let group_result = group_elems_by_bit(&bits);
  if let Err(e) = group_result {
    return Err(e);
  }
  let (elems_with_zero, elems_with_one) = group_result.unwrap();
  let (new_acc, membership_proof) = accumulator::add::<G>(acc, &elems_with_one);
  let nonmembership_proof =
    accumulator::prove_nonmembership(&new_acc, acc_set, &elems_with_zero).unwrap();
  Ok((
    new_acc,
    VectorProof {
      membership_proof,
      nonmembership_proof,
    },
  ))
}

pub fn open<G: UnknownOrderGroup>(
  acc: &G::Elem,
  acc_set: &[Integer],
  zero_bits: &[Integer],
  one_bit_witnesses: &[(Integer, G::Elem)],
) -> Result<VectorProof<G>, VCError> {
  let elems_with_zero: Vec<Integer> = zero_bits
    .iter()
    .map(|i| hash_to_prime(&Blake2b::default, &i))
    .collect();
  let elem_witnesses_with_one: Vec<(Integer, G::Elem)> = one_bit_witnesses
    .iter()
    .map(|(i, witness)| (hash_to_prime(&Blake2b::default, &i), witness.clone()))
    .collect();
  let membership_proof = accumulator::prove_membership::<G>(acc, &elem_witnesses_with_one);
  if membership_proof.is_err() {
    return Err(VCError::InvalidOpenError);
  }
  let nonmembership_proof = accumulator::prove_nonmembership::<G>(acc, acc_set, &elems_with_zero);
  if nonmembership_proof.is_err() {
    return Err(VCError::InvalidOpenError);
  }
  Ok(VectorProof {
    membership_proof: membership_proof.unwrap(),
    nonmembership_proof: nonmembership_proof.unwrap(),
  })
}

pub fn verify<G: UnknownOrderGroup>(
  acc: &G::Elem,
  bits: &[(bool, Integer)],
  VectorProof {
    membership_proof,
    nonmembership_proof,
  }: &VectorProof<G>,
) -> bool {
  let group_result = group_elems_by_bit(&bits);
  if group_result.is_err() {
    return false;
  }
  let (elems_with_zero, elems_with_one) = group_result.unwrap();
  let verified_membership = accumulator::verify_membership(acc, &elems_with_one, membership_proof);
  let verified_nonmembership =
    accumulator::verify_nonmembership(acc, &elems_with_zero, nonmembership_proof);
  verified_membership && verified_nonmembership
}

#[cfg(test)]
mod tests {}
