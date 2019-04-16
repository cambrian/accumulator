//! Vector commitment library, built on a generic group interface. **Very much a WIP.**
use super::accumulator::{Accumulator, MembershipProof, NonmembershipProof, Witness};
use crate::group::UnknownOrderGroup;
use rug::Integer;
use std::collections::HashSet;

#[derive(Debug)]
/// The different types of vector commitment errors.
pub enum VCError {
  /// When there are conflicting indices in the vector commitment.
  ConflictingIndices,
  /// When an opening fails.
  InvalidOpen,
  /// Unexpected state during an update.
  UnexpectedState,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
/// A vector commitment, wrapping an underlying accumulator. The accumulator contains indices of an
/// abstract vector where the corresponding bit is True.
pub struct VectorCommitment<G: UnknownOrderGroup>(Accumulator<G, Integer>);

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
/// A vector commitment proof.
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
  /// Initializes a new vector commitment (VC).
  pub fn empty() -> Self {
    Self(Accumulator::<G, Integer>::empty())
  }

  /// Updates a VC with a list of values and indices.
  ///
  /// # Arguments
  ///
  /// * `vc_acc_set` - All indices that are set (True).
  /// * `bits` - Tuples (truth value, bit index) to set.
  ///
  /// Uses a move instead of a `&self` reference to prevent accidental use of the old VC state.
  pub fn update(
    vc: Self,
    vc_acc_set: &[Integer],
    bits: &[(bool, Integer)],
  ) -> Result<(Self, VectorProof<G>), VCError> {
    let (elems_with_zero, elems_with_one) = group_elems_by_bit(&bits)?;
    let (new_acc, membership_proof) = vc.0.add_with_proof(&elems_with_one);
    let nonmembership_proof = new_acc
      .prove_nonmembership(vc_acc_set, &elems_with_zero)
      .map_err(|_| VCError::UnexpectedState)?;
    Ok((
      Self(new_acc),
      VectorProof {
        membership_proof,
        nonmembership_proof,
      },
    ))
  }

  /// Opens/generates a commitment to indices in the VC.
  ///
  /// # Arguments
  /// * `vc_acc_set` - All indices that are set (True).
  /// * `zero_bits` - Indices you want to prove are unset (False).
  /// * `one_bit_witnesses` - Indices you want to prove are set (True) and their witnesses.
  pub fn open(
    vc: &Self,
    vc_acc_set: &[Integer],
    zero_bits: &[Integer],
    one_bit_witnesses: &[(Integer, Witness<G, Integer>)],
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

  /// Verifies a commitment to indices in the VC.
  ///
  /// # Arguments
  ///
  /// * `bits` - Tuples (truth value, bit index) to verify.
  /// * `VectorProof` - A `VectorProof` to verify against.
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
    let verified_membership = vc
      .0
      .verify_membership_batch(&elems_with_one, membership_proof);
    let verified_nonmembership = vc
      .0
      .verify_nonmembership(&elems_with_zero, nonmembership_proof);
    verified_membership && verified_nonmembership
  }
}

// TODO: Write tests.
#[cfg(test)]
mod tests {}
