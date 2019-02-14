//! Accumulator library, built on a generic group interface.
//!
//! Operations that "mutate" the accumulator (add, delete) use moves instead of references so that
//! you don't accidentally use the old accumulator state.
use crate::group::UnknownOrderGroup;
use crate::proof::{Poe, Poke2};
use crate::util::{int, shamir_trick};
use rug::Integer;

#[derive(Debug)]
pub enum AccError {
  BadWitness,
  FailedDivision,
  InputsNotCoprime,
}

#[derive(PartialEq, Eq, Clone, Debug, Default, Hash)]
pub struct Accumulator<G: UnknownOrderGroup>(G::Elem);

#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct MembershipProof<G: UnknownOrderGroup> {
  pub witness: Accumulator<G>,
  proof: Poe<G>,
}

#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct NonmembershipProof<G: UnknownOrderGroup> {
  d: G::Elem,
  v: G::Elem,
  gv_inv: G::Elem,
  poke2_proof: Poke2<G>,
  poe_proof: Poe<G>,
}

impl<G: UnknownOrderGroup> Accumulator<G> {
  /// Initializes the accumulator to a group element.
  pub fn new() -> Self {
    Accumulator(G::unknown_order_elem())
  }

  // Computes `self ^ (numerator / denominator)`.
  pub fn exp_quotient(self, numerator: Integer, denominator: Integer) -> Result<Self, AccError> {
    if denominator == int(0) {
      return Err(AccError::FailedDivision);
    }

    let (quotient, remainder) = numerator.div_rem(denominator);

    if remainder != int(0) {
      return Err(AccError::FailedDivision);
    }

    Ok(Accumulator(G::exp(&self.0, &quotient)))
  }

  // The conciseness of accumulator.add() and low probability of confusion with implementations of
  // the Add trait probably justify this...
  #[allow(clippy::should_implement_trait)]
  /// Adds `elems` to the accumulator `acc`. Cannot check whether the elements are coprime with the
  /// accumulator, but it is up to clients to either ensure uniqueness or treat this as multiset.
  pub fn add(self, elems: &[Integer]) -> (Self, MembershipProof<G>) {
    let x = elems.iter().product();
    let acc_new = G::exp(&self.0, &x);
    let poe_proof = Poe::<G>::prove(&self.0, &x, &acc_new);
    (
      Accumulator(acc_new),
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

    let poe_proof = Poe::<G>::prove(&acc_next, &elem_aggregate, &self.0);
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
    Poe::verify(&witness.0, &exp, &self.0, proof)
  }

  /// Updates a membership witness for some set. See Section 4.2 in the Li, Li, Xue paper.
  pub fn update_membership_witness(
    self,
    acc_new: Self,
    witness_set: &[Integer],
    untracked_additions: &[Integer],
    untracked_deletions: &[Integer],
  ) -> Result<Self, AccError> {
    let witness_set_product: Integer = witness_set.iter().product();
    let untracked_delete_product: Integer = untracked_deletions.iter().product();

    let (gcd, a, b) = <(Integer, Integer, Integer)>::from(
      witness_set_product.gcd_cofactors_ref(&untracked_delete_product),
    );

    if gcd != int(1) {
      return Err(AccError::InputsNotCoprime);
    }

    let witness_post_add = self.add(untracked_additions);
    let w_to_b = G::exp(&(witness_post_add.0).0, &b);
    let acc_new_to_a = G::exp(&acc_new.0, &a);
    Ok(Accumulator(G::op(&w_to_b, &acc_new_to_a)))
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

    let poke2_proof = Poke2::prove(&self.0, &b, &v);
    let poe_proof = Poe::prove(&d, &x, &gv_inv);
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
    Poke2::verify(&self.0, v, poke2_proof) && Poe::verify(d, &x, gv_inv, poe_proof)
  }

  #[allow(non_snake_case)]
  /// For accumulator with value `g` and elems `[x_1, ..., x_n]`, computes a membership witness for
  /// each `x_i` in accumulator `g^{x_1 * ... * x_n}`, namely `g^{x_1 * ... * x_n / x_i}`, in O(N
  /// log N) time.
  pub fn root_factor(&self, elems: &[Integer]) -> Vec<Accumulator<G>> {
    if elems.len() == 1 {
      return vec![self.clone()];
    }
    let half_n = elems.len() / 2;
    let g_l = elems[..half_n]
      .iter()
      .fold(self.clone(), |sum, x| Accumulator(G::exp(&sum.0, x)));
    let g_r = elems[half_n..]
      .iter()
      .fold(self.clone(), |sum, x| Accumulator(G::exp(&sum.0, x)));
    let mut L = g_r.root_factor(&Vec::from(&elems[..half_n]));
    let mut R = g_l.root_factor(&Vec::from(&elems[half_n..]));
    L.append(&mut R);
    L
  }
}

// TODO: Add test for `prove_membership`.
#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, Rsa2048};
  use crate::util::int;

  fn init_acc<G: UnknownOrderGroup>() -> Accumulator<G> {
    Accumulator::<G>::new().add(&[int(41), int(67), int(89)]).0
  }

  #[test]
  fn test_exp_quotient() {
    let empty_acc = Accumulator::<Rsa2048>::new();
    let exp_quotient_result = empty_acc
      .exp_quotient(int(17 * 41 * 67 * 89), int(17 * 89))
      .unwrap();
    let exp_quotient_expected =
      Accumulator(Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(41 * 67)));
    assert!(exp_quotient_result == exp_quotient_expected);
  }

  #[test]
  #[should_panic(expected = "FailedDivision")]
  fn test_exp_quotient_zero() {
    Accumulator::<Rsa2048>::new()
      .exp_quotient(int(17 * 41 * 67 * 89), int(0))
      .unwrap();
  }

  #[test]
  #[should_panic(expected = "FailedDivision")]
  fn test_exp_quotient_remainder() {
    Accumulator::<Rsa2048>::new()
      .exp_quotient(int(17 * 41 * 67 * 89), int(5))
      .unwrap();
  }

  #[test]
  fn test_add() {
    let acc = init_acc::<Rsa2048>();
    let new_elems = [int(5), int(7), int(11)];
    let (acc_new, proof) = acc.add(&new_elems);
    let acc_expected = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(94_125_955));
    assert!(acc_new.0 == acc_expected);
    assert!(acc_new.verify_membership(&new_elems, &proof));
  }

  #[test]
  fn test_delete() {
    let acc = init_acc::<Rsa2048>();
    let y_witness = Accumulator::<Rsa2048>::new().add(&[int(3649)]).0;
    let z_witness = Accumulator::<Rsa2048>::new().add(&[int(2747)]).0;
    let (acc_new, proof) = acc
      .clone()
      .delete(&[(int(67), y_witness), (int(89), z_witness)])
      .expect("valid delete expected");
    let acc_expected = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(41));
    assert!(acc_new.0 == acc_expected);
    assert!(acc.verify_membership(&[int(67), int(89)], &proof));
  }

  #[test]
  fn test_delete_empty() {
    let acc = init_acc::<Rsa2048>();
    let (acc_new, proof) = acc.clone().delete(&[]).expect("valid delete expected");
    assert!(acc_new == acc);
    assert!(acc.verify_membership(&[], &proof));
  }

  #[should_panic(expected = "BadWitness")]
  #[test]
  fn test_delete_bad_witness() {
    let acc = init_acc::<Rsa2048>();
    let y_witness = Accumulator::<Rsa2048>::new().add(&[int(3648)]).0;
    let z_witness = Accumulator::<Rsa2048>::new().add(&[int(2746)]).0;
    acc
      .delete(&[(int(67), y_witness), (int(89), z_witness)])
      .unwrap();
  }

  #[test]
  fn test_update_membership_witness() {
    // Original accumulator has [3, 5, 11, 13].
    // Witness is tracking elements [3, 5] and eventually [7].
    let acc_new = Accumulator::<Rsa2048>::new()
      .add(&[int(3), int(7), int(11), int(17)])
      .0;
    let witness = Accumulator::<Rsa2048>::new().add(&[int(11), int(13)]).0;
    let witness_new = witness
      .update_membership_witness(acc_new.clone(), &[int(3), int(7)], &[int(17)], &[int(13)])
      .unwrap();
    assert!(witness_new.add(&[int(3), int(7)]).0 == acc_new);
  }

  #[test]
  fn test_prove_nonmembership() {
    let acc = init_acc::<Rsa2048>();
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
    let acc = init_acc::<Rsa2048>();
    let acc_set = [int(41), int(67), int(89)];
    let elems = [int(41), int(7), int(11)];
    acc.prove_nonmembership(&acc_set, &elems).unwrap();
  }

  #[test]
  fn test_root() {
    let acc = Accumulator::<Rsa2048>::new();
    let (acc, _) = acc.add(&[int(41), int(67), int(89)]);
    let factors = [int(97), int(101), int(103), int(107), int(109)];
    let witnesses = acc.root_factor(&factors);
    for (i, witness) in witnesses.iter().enumerate() {
      let partial_product = factors.iter().product::<Integer>() / factors[i].clone();
      let expected = acc.clone().add(&[partial_product]).0;
      assert_eq!(*witness, expected);
    }
  }
}
