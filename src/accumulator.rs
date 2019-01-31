//! Accumulator library, built on a generic group interface.
//!
//! Operations that "mutate" the accumulator (add, delete) use moves instead of references so that
//! you don't accidentally use the old accumulator state.
use crate::group::UnknownOrderGroup;
use crate::proof::{PoE, PoKE2};
use crate::util::{int, shamir_trick};
use rug::Integer;

#[derive(Debug)]
pub enum AccError {
  BadWitness,
  InputsNotCoprime,
}

#[derive(PartialEq, Eq, Clone, Debug, Default)]
pub struct Accumulator<G: UnknownOrderGroup>(G::Elem);

pub struct MembershipProof<G: UnknownOrderGroup> {
  witness: Accumulator<G>,
  proof: PoE<G>,
}

pub struct NonmembershipProof<G: UnknownOrderGroup> {
  d: G::Elem,
  v: G::Elem,
  gv_inv: G::Elem,
  poke2_proof: PoKE2<G>,
  poe_proof: PoE<G>,
}

impl<G: UnknownOrderGroup> Accumulator<G> {
  /// Initializes the accumulator to a group element.
  pub fn new() -> Self {
    Accumulator(G::unknown_order_elem())
  }

  // The conciseness of accumulator.add() and low probability of confusion with implementations of
  // the Add trait probably justify this...
  #[allow(clippy::should_implement_trait)]
  /// Adds `elems` to the accumulator `acc`. Cannot check whether the elements are coprime with the
  /// accumulator, but it is up to clients to either ensure uniqueness or treat this as multiset.
  pub fn add(self, elems: &[Integer]) -> (Self, MembershipProof<G>) {
    let x = elems.iter().product();
    let new_acc = G::exp(&self.0, &x);
    let poe_proof = PoE::<G>::prove(&self.0, &x, &new_acc);
    (
      Accumulator(new_acc),
      MembershipProof {
        witness: self,
        proof: poe_proof,
      },
    )
  }

  /// Removes the elements in `elem_witnesses` from the accumulator `acc`.
  pub fn delete(
    self,
    elem_witnesses: &[(Integer, Self)],
  ) -> Result<(Self, MembershipProof<G>), AccError> {
    let mut elem_aggregate = int(1);
    let mut acc_next = self.0.clone();

    for (elem, witness) in elem_witnesses {
      if G::exp(&witness.0, elem) != self.0 {
        return Err(AccError::BadWitness);
      }

      acc_next = shamir_trick::<G>(&acc_next, &witness.0, &elem_aggregate, elem)
        .ok_or(AccError::InputsNotCoprime)?;
      elem_aggregate *= elem;
    }

    let poe_proof = PoE::<G>::prove(&acc_next, &elem_aggregate, &self.0);
    Ok((
      Accumulator(acc_next.clone()),
      MembershipProof {
        witness: Accumulator(acc_next),
        proof: poe_proof,
      },
    ))
  }

  /// Returns a proof (and associated variables) that `elem_witnesses` are aggregated in `acc`.
  pub fn prove_membership(
    &self,
    elem_witnesses: &[(Integer, Self)],
  ) -> Result<MembershipProof<G>, AccError> {
    Ok(self.clone().delete(elem_witnesses)?.1)
  }

  /// Verifies the PoE returned by `prove_membership`.
  pub fn verify_membership(
    &self,
    elems: &[Integer],
    MembershipProof { witness, proof }: &MembershipProof<G>,
  ) -> bool {
    let exp = elems.iter().product();
    PoE::verify(&witness.0, &exp, &self.0, proof)
  }

  /// Returns a proof (and associated variables) that `elems` are not in `acc_set`.
  pub fn prove_nonmembership(
    &self,
    acc_set: &[Integer],
    elems: &[Integer],
  ) -> Result<NonmembershipProof<G>, AccError> {
    let x: Integer = elems.iter().product();
    let s = acc_set.iter().product();
    let (gcd, a, b) = <(Integer, Integer, Integer)>::from(x.gcd_cofactors_ref(&s));

    if gcd != int(1) {
      return Err(AccError::InputsNotCoprime);
    }

    let g = G::unknown_order_elem();
    let d = G::exp(&g, &a);
    let v = G::exp(&self.0, &b);
    let gv_inv = G::op(&g, &G::inv(&v));

    let poke2_proof = PoKE2::prove(&self.0, &b, &v);
    let poe_proof = PoE::prove(&d, &x, &gv_inv);
    Ok(NonmembershipProof {
      d,
      v,
      gv_inv,
      poke2_proof,
      poe_proof,
    })
  }

  /// Verifies the PoKE2 and PoE returned by `prove_nonmembership`.
  pub fn verify_nonmembership(
    &self,
    elems: &[Integer],
    NonmembershipProof {
      d,
      v,
      gv_inv,
      poke2_proof,
      poe_proof,
    }: &NonmembershipProof<G>,
  ) -> bool {
    let x = elems.iter().product();
    PoKE2::verify(&self.0, v, poke2_proof) && PoE::verify(d, &x, gv_inv, poe_proof)
  }
}

// TODO: Add test for `prove_membership`.
#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, RSA2048};
  use crate::util::int;

  fn init_acc<G: UnknownOrderGroup>() -> Accumulator<G> {
    Accumulator::<G>::new().add(&[int(41), int(67), int(89)]).0
  }

  #[test]
  fn test_add() {
    let acc = init_acc::<RSA2048>();
    let new_elems = [int(5), int(7), int(11)];
    let (new_acc, proof) = acc.add(&new_elems);
    let expected_acc = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(94_125_955));
    assert!(new_acc.0 == expected_acc);
    assert!(new_acc.verify_membership(&new_elems, &proof));
  }

  #[test]
  fn test_delete() {
    let acc = init_acc::<RSA2048>();
    let y_witness = Accumulator::<RSA2048>::new().add(&[int(3649)]).0;
    let z_witness = Accumulator::<RSA2048>::new().add(&[int(2747)]).0;
    let (new_acc, proof) = acc
      .clone()
      .delete(&[(int(67), y_witness), (int(89), z_witness)])
      .expect("valid delete expected");
    let expected_acc = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(41));
    assert!(new_acc.0 == expected_acc);
    assert!(acc.verify_membership(&[int(67), int(89)], &proof));
  }

  #[test]
  fn test_delete_empty() {
    let acc = init_acc::<RSA2048>();
    let (new_acc, proof) = acc.clone().delete(&[]).expect("valid delete expected");
    assert!(new_acc == acc);
    assert!(acc.verify_membership(&[], &proof));
  }

  #[should_panic(expected = "BadWitness")]
  #[test]
  fn test_delete_bad_witness() {
    let acc = init_acc::<RSA2048>();
    let y_witness = Accumulator::<RSA2048>::new().add(&[int(3648)]).0;
    let z_witness = Accumulator::<RSA2048>::new().add(&[int(2746)]).0;
    acc
      .delete(&[(int(67), y_witness), (int(89), z_witness)])
      .unwrap();
  }

  #[test]
  fn test_prove_nonmembership() {
    let acc = init_acc::<RSA2048>();
    let acc_set = [int(41), int(67), int(89)];
    let elems = [int(5), int(7), int(11)];
    let proof = acc
      .prove_nonmembership(&acc_set, &elems)
      .expect("valid proof expected");
    assert!(acc.verify_nonmembership(&elems, &proof));
  }

  #[should_panic(expected = "InputsNotCoprime")]
  #[test]
  fn test_prove_nonmembership_failure() {
    let acc = init_acc::<RSA2048>();
    let acc_set = [int(41), int(67), int(89)];
    let elems = [int(41), int(7), int(11)];
    acc.prove_nonmembership(&acc_set, &elems).unwrap();
  }
}
