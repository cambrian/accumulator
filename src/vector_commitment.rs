//! Vector commitment library, built on a generic group interface. WIP.
use super::accumulator::{Accumulator, MembershipProof, NonmembershipProof};
use crate::group::UnknownOrderGroup;
use rug::Integer;
use std::collections::HashSet;

#[derive(Debug)]
pub enum VCError {
  ConflictingIndices,
  InvalidOpen,
  UnexpectedState,
}

// Represented internally by an accumulator of indices where the corresponding bit == 1.
#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct VectorCommitment<G: UnknownOrderGroup>(Accumulator<G, Integer>);

#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct VectorProof<G: UnknownOrderGroup> {
  membership_proof: MembershipProof<G, Integer>,
  nonmembership_proof: NonmembershipProof<G, Integer>,
}

fn group_elems_by_bit(bits: &[(bool, Integer)]) -> Result<(Vec<Integer>, Vec<Integer>), VCError> {
  let mut elems_with_one = vec![];
  let mut elems_with_zero = vec![];
  let mut seen_indices = HashSet::new();
  for (bit, i) in bits {
    if !seen_indices.insert(i) {
      return Err(VCError::ConflictingIndices);
    }
    if *bit {
      elems_with_one.push(i.clone());
    } else {
      elems_with_zero.push(i.clone());
    }
  }
  Ok((elems_with_zero, elems_with_one))
}

impl<G: UnknownOrderGroup> VectorCommitment<G> {
  pub fn empty() -> Self {
    VectorCommitment(Accumulator::<G, Integer>::empty())
  }

  // Uses a move instead of a `&self` reference to prevent accidental use of the old vc state.
  pub fn update(
    vc: Self,
    vc_acc_set: &[Integer],
    bits: &[(bool, Integer)],
  ) -> Result<(Self, VectorProof<G>), VCError> {
    // Must hold hash commitments in vec in order to pass by reference to accumulator fns.
    let (elems_with_zero, elems_with_one) = group_elems_by_bit(&bits)?;
    let (new_acc, membership_proof) = vc.0.add(&elems_with_one);
    let nonmembership_proof = new_acc
      .prove_nonmembership(vc_acc_set, &elems_with_zero)
      .map_err(|_| VCError::UnexpectedState)?;
    Ok((
      VectorCommitment(new_acc),
      VectorProof {
        membership_proof,
        nonmembership_proof,
      },
    ))
  }

  pub fn open(
    vc: &Self,
    vc_acc_set: &[Integer],
    zero_bits: &[Integer],
    one_bit_witnesses: &[(Integer, Accumulator<G, Integer>)],
  ) -> Result<VectorProof<G>, VCError> {
    let membership_proof = vc
      .0
      .prove_membership(one_bit_witnesses)
      .map_err(|_| VCError::InvalidOpen)?;
    let nonmembership_proof = vc
      .0
      .prove_nonmembership(vc_acc_set, zero_bits)
      .map_err(|_| VCError::InvalidOpen)?;
    Ok(VectorProof {
      membership_proof,
      nonmembership_proof,
    })
  }

  pub fn verify(
    vc: &Self,
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
    let verified_membership = vc.0.verify_membership(&elems_with_one, membership_proof);
    let verified_nonmembership = vc
      .0
      .verify_nonmembership(&elems_with_zero, nonmembership_proof);
    verified_membership && verified_nonmembership
  }
}

#[cfg(test)]
mod tests {}
