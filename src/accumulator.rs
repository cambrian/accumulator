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

/// Initializes the accumulator to a group element.
pub fn setup<G: UnknownOrderGroup>() -> G::Elem {
  G::unknown_order_elem()
}

pub struct MembershipProof<G: UnknownOrderGroup> {
  witness: G::Elem,
  proof: PoE<G>,
}

/// Adds `elems` to the accumulator `acc`. Cannot check whether the elements are co-prime with the
/// accumulator, but it is up to clients to either ensure uniqueness or treat this as multi-set.
pub fn add<G: UnknownOrderGroup>(acc: G::Elem, elems: &[Integer]) -> (G::Elem, MembershipProof<G>) {
  let x = elems.iter().product();
  let new_acc = G::exp(&acc, &x);
  let poe_proof = PoE::<G>::prove(&acc, &x, &new_acc);
  (
    new_acc,
    MembershipProof {
      witness: acc.clone(),
      proof: poe_proof,
    },
  )
}

/// Removes the elements in `elem_witnesses` from the accumulator `acc`.
pub fn delete<G: UnknownOrderGroup>(
  acc: G::Elem,
  elem_witnesses: &[(Integer, G::Elem)],
) -> Result<(G::Elem, MembershipProof<G>), AccError> {
  let mut elem_aggregate = int(1);
  let mut acc_next = acc.clone();

  for (elem, witness) in elem_witnesses {
    if G::exp(witness, elem) != acc {
      return Err(AccError::BadWitness);
    }

    let acc_next_option = shamir_trick::<G>(&acc_next, witness, &elem_aggregate, elem);
    match acc_next_option {
      Some(acc_next_value) => acc_next = acc_next_value,
      None => return Err(AccError::InputsNotCoprime),
    };

    elem_aggregate *= elem;
  }

  let poe_proof = PoE::<G>::prove(&acc_next, &elem_aggregate, &acc);
  Ok((
    acc_next.clone(),
    MembershipProof {
      witness: acc_next,
      proof: poe_proof,
    },
  ))
}

/// Returns a proof (and associated variables) that `elem_witnesses` are aggregated in `acc`.
pub fn prove_membership<G: UnknownOrderGroup>(
  acc: &G::Elem,
  elem_witnesses: &[(Integer, G::Elem)],
) -> Result<MembershipProof<G>, AccError> {
  let delete_result = delete::<G>(acc.clone(), elem_witnesses);
  match delete_result {
    Ok((_, membership_proof)) => Ok(membership_proof),
    Err(e) => Err(e),
  }
}

/// Verifies the PoE returned by `prove_membership`.
pub fn verify_membership<G: UnknownOrderGroup>(
  acc: &G::Elem,
  elems: &[Integer],
  MembershipProof { witness, proof }: &MembershipProof<G>,
) -> bool {
  let exp = elems.iter().product();
  PoE::verify(witness, &exp, acc, proof)
}

pub struct NonmembershipProof<G: UnknownOrderGroup> {
  d: G::Elem,
  v: G::Elem,
  gv_inv: G::Elem,
  poke2_proof: PoKE2<G>,
  poe_proof: PoE<G>,
}

/// Returns a proof (and associated variables) that `elems` are not in `acc_set`.
pub fn prove_nonmembership<G: UnknownOrderGroup>(
  acc: &G::Elem,
  acc_set: &[Integer],
  elems: &[Integer],
) -> Result<NonmembershipProof<G>, AccError> {
  let x: Integer = elems.iter().product();
  let s = acc_set.iter().product();
  let (gcd, a, b) = x.clone().gcd_cofactors(s, Integer::new());

  if gcd != int(1) {
    return Err(AccError::InputsNotCoprime);
  }

  let g = G::unknown_order_elem();
  let d = G::exp(&g, &a);
  let v = G::exp(acc, &b);
  let gv_inv = G::op(&g, &G::inv(&v));

  let poke2_proof = PoKE2::prove(acc, &b, &v);
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
pub fn verify_nonmembership<G: UnknownOrderGroup>(
  acc: &G::Elem,
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
  PoKE2::verify(acc, v, poke2_proof) && PoE::verify(d, &x, gv_inv, poe_proof)
}

// TODO: Add test for `prove_membership`.
#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{ClassGroup, Group, RSA2048};
  use crate::util::int;

  fn init_acc<G: UnknownOrderGroup>() -> G::Elem {
    G::exp(&setup::<G>(), &(int(41) * &int(67) * &int(89)))
  }

  #[test]
  fn test_add() {
    fn test_add_generic<G: UnknownOrderGroup>() {
      let acc = init_acc::<G>();
      let new_elems = [int(5), int(7), int(11)];
      let (new_acc, proof) = add::<G>(acc.clone(), &new_elems);
      let expected_acc = G::exp(&G::unknown_order_elem(), &int(94_125_955));
      assert!(new_acc == expected_acc);
      assert!(verify_membership(&new_acc, &new_elems, &proof));
    }
    test_add_generic::<RSA2048>();
    test_add_generic::<ClassGroup>();
  }

  #[test]
  fn test_delete() {
    fn test_delete_generic<G: UnknownOrderGroup>() {
      let acc = init_acc::<G>();
      let y_witness = G::exp(&G::unknown_order_elem(), &int(3649));
      let z_witness = G::exp(&G::unknown_order_elem(), &int(2747));
      let (new_acc, proof) =
        delete::<G>(acc.clone(), &[(int(67), y_witness), (int(89), z_witness)])
          .expect("valid delete expected");
      let expected_acc = G::exp(&G::unknown_order_elem(), &int(41));
      assert!(new_acc == expected_acc);
      assert!(verify_membership(&acc, &[int(67), int(89)], &proof));
    }
    test_delete_generic::<RSA2048>();
    test_delete_generic::<ClassGroup>();
  }

  #[test]
  fn test_delete_empty() {
    fn test_delete_empty_generic<G: UnknownOrderGroup>() {
      let acc = init_acc::<G>();
      let (new_acc, proof) = delete::<G>(acc.clone(), &[]).expect("valid delete expected");
      assert!(new_acc == acc);
      assert!(verify_membership(&acc, &[], &proof));
    }
    test_delete_empty_generic::<RSA2048>();
    test_delete_empty_generic::<ClassGroup>();
  }

  #[should_panic(expected = "BadWitness")]
  #[test]
  fn test_delete_bad_witness_rsa() {
    let acc = init_acc::<RSA2048>();
    let y_witness = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(3648));
    let z_witness = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(2746));
    delete::<RSA2048>(acc, &[(int(67), y_witness), (int(89), z_witness)]).unwrap();
  }

  #[should_panic(expected = "BadWitness")]
  #[test]
  fn test_delete_bad_witness_class() {
    // TODO: figure out how to merge this and previous test into one single test.
    let acc = init_acc::<ClassGroup>();
    let y_witness = ClassGroup::exp(&ClassGroup::unknown_order_elem(), &int(3648));
    let z_witness = ClassGroup::exp(&ClassGroup::unknown_order_elem(), &int(2746));
    delete::<ClassGroup>(acc, &[(int(67), y_witness), (int(89), z_witness)]).unwrap();
  }

  #[test]
  fn test_prove_nonmembership() {
    fn test_prove_nonmembership_generic<G: UnknownOrderGroup>() {
      let acc = init_acc::<G>();
      let acc_set = [int(41), int(67), int(89)];
      let elems = [int(5), int(7), int(11)];
      let proof = prove_nonmembership::<G>(&acc, &acc_set, &elems).expect("valid proof expected");
      assert!(verify_nonmembership::<G>(&acc, &elems, &proof));
    }
    test_prove_nonmembership_generic::<RSA2048>();
    test_prove_nonmembership_generic::<ClassGroup>();
  }

  #[should_panic(expected = "InputsNotCoprime")]
  #[test]
  fn test_prove_nonmembership_failure_rsa() {
    let acc = init_acc::<RSA2048>();
    let acc_set = [int(41), int(67), int(89)];
    let elems = [int(41), int(7), int(11)];
    prove_nonmembership::<RSA2048>(&acc, &acc_set, &elems).unwrap();
  }

  #[should_panic(expected = "InputsNotCoprime")]
  #[test]
  fn test_prove_nonmembership_failure_class() {
    // TODO: figure out how to merge this and previous test into one single test.
    let acc = init_acc::<ClassGroup>();
    let acc_set = [int(41), int(67), int(89)];
    let elems = [int(41), int(7), int(11)];
    prove_nonmembership::<ClassGroup>(&acc, &acc_set, &elems).unwrap();
  }
}
