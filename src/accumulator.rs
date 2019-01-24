//! Accumulator library, built on a generic group interface.
//!
//! Operations that "mutate" the accumulator (add, delete) use moves instead of references so that
//! you don't accidentally use the old accumulator state.
use crate::group::UnknownOrderGroup;
use crate::proof::{PoE, PoKE2};
use crate::util::{bezout, int, product, shamir_trick};
use rug::Integer;
use num_traits::identities::One;

#[derive(Debug)]
pub enum AccError {
  BadWitness,
  InputsNotCoPrime,
}

/// Initializes the accumulator to a group element.
pub fn setup<G: UnknownOrderGroup>() -> G::Elem {
  G::unknown_order_elem()
}

/// Adds `elems` to the accumulator `acc`. Cannot check whether the elements are co-prime with the
/// accumulator, but it is up to clients to either ensure uniqueness or treat this as multi-set.
pub fn add<G: UnknownOrderGroup>(acc: G::Elem, elems: &[Integer]) -> (G::Elem, PoE<G>) {
  let x = product(elems);
  let new_acc = G::exp(&acc, &x);
  let poe_proof = PoE::<G>::prove(&acc, &x, &new_acc);
  (new_acc, poe_proof)
}

/// Removes the elements in `elem_witnesses` from the accumulator `acc`.
pub fn delete<G: UnknownOrderGroup>(
  acc: G::Elem,
  elem_witnesses: &[(Integer, G::Elem)],
) -> Result<(G::Elem, PoE<G>), AccError> {
  let mut elem_aggregate = int(1);
  let mut acc_next = acc.clone();

  for (elem, witness) in elem_witnesses {
    if G::exp(witness, elem) != acc {
      return Err(AccError::BadWitness);
    }

    let acc_next_option = shamir_trick::<G>(&acc_next, witness, &elem_aggregate, elem);
    match acc_next_option {
      Some(acc_next_value) => acc_next = acc_next_value,
      None => return Err(AccError::InputsNotCoPrime),
    };

    elem_aggregate *= elem;
  }

  let poe_proof = PoE::<G>::prove(&acc_next, &elem_aggregate, &acc);
  Ok((acc_next, poe_proof))
}

/// See `delete`.
pub fn prove_membership<G: UnknownOrderGroup>(
  acc: &G::Elem,
  elem_witnesses: &[(Integer, G::Elem)],
) -> Result<(G::Elem, PoE<G>), AccError> {
  delete::<G>(acc.clone(), elem_witnesses)
}

/// Verifies the PoE returned by `prove_membership` s.t. `witness` ^ `elems` = `result`.
pub fn verify_membership<G: UnknownOrderGroup>(
  witness: &G::Elem,
  elems: &[Integer],
  result: &G::Elem,
  proof: &PoE<G>,
) -> bool {
  let exp = product(elems);
  PoE::verify(witness, &exp, result, proof)
}

pub struct NonMembershipProof<G: UnknownOrderGroup> {
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
) -> Result<NonMembershipProof<G>, AccError> {
  let x = product(elems);
  let s = product(acc_set);
  let (a, b, gcd) = bezout(&x, &s);

  if !gcd.is_one() {
    return Err(AccError::InputsNotCoPrime);
  }

  let g = G::unknown_order_elem();
  let d = G::exp_signed(&g, &a);
  let v = G::exp_signed(acc, &b);
  let gv_inv = G::op(&g, &G::inv(&v));

  let poke2_proof = PoKE2::prove(acc, &b, &v);
  let poe_proof = PoE::prove(&d, &x, &gv_inv);
  Ok(NonMembershipProof {
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
  NonMembershipProof {
    d,
    v,
    gv_inv,
    poke2_proof,
    poe_proof,
  }: &NonMembershipProof<G>,
) -> bool {
  let x = product(elems);
  PoKE2::verify(acc, v, poke2_proof) && PoE::verify(d, &x, gv_inv, poe_proof)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, RSA2048};
  use crate::util::int;

  fn init_acc<G: UnknownOrderGroup>() -> G::Elem {
    G::exp(&setup::<G>(), &(int(41) * &int(67) * &int(89)))
  }

  #[test]
  fn test_shamir_trick() {
    let (x, y, z) = (&int(13), &int(17), &int(19));
    let xth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &(y * z));
    let yth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &(x * z));
    let xyth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), z);
    assert!(shamir_trick::<RSA2048>(&xth_root, &yth_root, x, y) == Some(xyth_root));
  }

  #[test]
  fn test_shamir_trick_failure() {
    let (x, y, z) = (&int(7), &int(14), &int(19)); // Inputs not co-prime.
    let xth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &(y * z));
    let yth_root = RSA2048::exp(&RSA2048::unknown_order_elem(), &(x * z));
    assert!(shamir_trick::<RSA2048>(&xth_root, &yth_root, x, y) == None);
  }

  #[test]
  fn test_add() {
    let acc = init_acc::<RSA2048>();
    let new_elems = [int(5), int(7), int(11)];
    let (new_acc, poe) = add::<RSA2048>(acc.clone(), &new_elems);
    let expected_acc = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(94_125_955u32));
    assert!(new_acc == expected_acc);
    assert!(PoE::verify(&acc, &int(385u16), &new_acc, &poe));
  }

  #[test]
  fn test_delete() {
    let acc = init_acc::<RSA2048>();
    let y_witness = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(3649u16));
    let z_witness = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(2747u16));
    let (new_acc, poe) =
      delete::<RSA2048>(acc.clone(), &[(int(67), y_witness), (int(89), z_witness)])
        .expect("valid delete expected");
    let expected_acc = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(41));
    assert!(new_acc == expected_acc);
    assert!(PoE::verify(&new_acc, &int(5963u16), &acc, &poe));
  }

  #[test]
  fn test_delete_empty() {
    let acc = init_acc::<RSA2048>();
    let (new_acc, poe) = delete::<RSA2048>(acc.clone(), &[]).expect("valid delete expected");
    assert!(new_acc == acc);
    assert!(PoE::verify(&new_acc, &int(1), &acc, &poe));
  }

  #[should_panic(expected = "BadWitness")]
  #[test]
  fn test_delete_bad_witness() {
    let acc = init_acc::<RSA2048>();
    let y_witness = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(3648u16));
    let z_witness = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(2746u16));
    delete::<RSA2048>(acc, &[(int(67), y_witness), (int(89), z_witness)]).unwrap();
  }

  #[test]
  fn test_prove_nonmembership() {
    let acc = init_acc::<RSA2048>();
    let acc_set = [int(41), int(67), int(89)];
    let elems = [int(5), int(7), int(11)];
    let proof =
      prove_nonmembership::<RSA2048>(&acc, &acc_set, &elems).expect("valid proof expected");
    assert!(verify_nonmembership::<RSA2048>(&acc, &elems, &proof));
  }

  #[should_panic(expected = "InputsNotCoPrime")]
  #[test]
  fn test_prove_nonmembership_failure() {
    let acc = init_acc::<RSA2048>();
    let acc_set = [int(41), int(67), int(89)];
    let elems = [int(41), int(7), int(11)];
    prove_nonmembership::<RSA2048>(&acc, &acc_set, &elems).unwrap();
  }
}
