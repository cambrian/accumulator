use super::group::{Group, InvertibleGroup};
use super::proof::{poe::PoE, poke2::PoKE2};
use super::util{product, bezout, shamir_trick};
use num;
use num::BigUint;
use num_traits::identities::{One, Zero};

#[derive(Debug)]
pub enum AccError {
  BadWitness,
  InputsNotCoPrime,
}

/// Initializes the accumulator to a group element.
pub fn setup<G: Group>() -> G::Elem {
  G::base_elem()
}

/// Adds `elems` to the accumulator `acc`.
pub fn add<G: Group>(acc: G::Elem, elems: &[&BigUint]) -> (G::Elem, PoE<G>) {
  let x = product(elems);
  let new_acc = G::exp(&acc, &x);
  let poe_proof = PoE::prove(&acc, &x, &new_acc);
  (new_acc, poe_proof)
}

/// Removes the elements in `elem_witnesses` from the accumulator `acc`.
pub fn delete<G: InvertibleGroup>(
  acc: G::Elem,
  elem_witnesses: &[&(BigUint, G::Elem)],
) -> Result<(G::Elem, PoE<G>), AccError> {
  // REVIEW: It should be possible to restructure the loop in such a way that this check is
  // unnecessary
  if elem_witnesses.is_empty() {
    let poe_proof = PoE::prove(&acc, &BigUint::zero(), &acc);
    return Ok((acc.clone(), poe_proof));
  }

  let mut elem_aggregate = elem_witnesses[0].0.clone();
  let mut acc_next = elem_witnesses[0].1.clone();

  for (elem, witness) in elem_witnesses
    .split_first() // Chop off first entry.
    .expect("unexpected witnesses")
    .1
  {
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

  let poe_proof = PoE::prove(&acc_next, &elem_aggregate, &acc);
  Ok((acc_next, poe_proof))
}

/// See `delete`.
pub fn prove_membership<G: InvertibleGroup>(
  acc: &G::Elem,
  elem_witnesses: &[&(BigUint, G::Elem)],
) -> Result<(G::Elem, PoE<G>), AccError> {
  delete::<G>(acc.clone(), elem_witnesses)
}

/// Verifies the PoE returned by `prove_membership` s.t. `witness` ^ `elems` = `result`.
pub fn verify_membership<G: Group>(
  witness: &G::Elem,
  elems: &[&BigUint],
  result: &G::Elem,
  proof: &PoE<G>,
) -> bool {
  let exp = product(elems);
  PoE::verify(witness, &exp, result, proof)
}

pub struct NonMembershipProof<G: Group> {
  d: G::Elem,
  v: G::Elem,
  gv_inv: G::Elem,
  poke2_proof: PoKE2<G>,
  poe_proof: PoE<G>,
}

/// Returns a proof (and associated variables) that `elems` are not in `acc_set`.
#[allow(clippy::type_complexity)]
pub fn prove_nonmembership<G: InvertibleGroup>(
  acc: &G::Elem,
  acc_set: &[&BigUint],
  elems: &[&BigUint],
) -> Result<NonMembershipProof<G>, AccError> {
  let x = product(elems);
  let s = product(acc_set);
  let (a, b, gcd) = bezout(&x, &s);

  if !gcd.is_one() {
    return Err(AccError::InputsNotCoPrime);
  }

  let g = G::base_elem();
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
pub fn verify_nonmembership<G: Group>(
  acc: &G::Elem,
  elems: &[&BigUint],
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
  use super::super::group::dummy::DummyRSA;
  use super::*;

  fn big(val: u64) -> BigUint {
    BigUint::from(val)
  }

  fn init_acc<G: Group>() -> G::Elem {
    G::exp(&setup::<G>(), &(big(41) * &big(67) * &big(89)))
  }

  #[test]
  fn test_shamir_trick() {
    let (x, y, z) = (&big(13), &big(17), &big(19));
    let xth_root = DummyRSA::exp(&DummyRSA::base_elem(), &(y * z));
    let yth_root = DummyRSA::exp(&DummyRSA::base_elem(), &(x * z));
    let xyth_root = DummyRSA::exp(&DummyRSA::base_elem(), z);
    assert!(shamir_trick::<DummyRSA>(&xth_root, &yth_root, x, y) == Some(xyth_root));
  }

  #[test]
  fn test_shamir_trick_failure() {
    let (x, y, z) = (&big(7), &big(14), &big(19)); // Inputs not co-prime.
    let xth_root = DummyRSA::exp(&DummyRSA::base_elem(), &(y * z));
    let yth_root = DummyRSA::exp(&DummyRSA::base_elem(), &(x * z));
    assert!(shamir_trick::<DummyRSA>(&xth_root, &yth_root, x, y) == None);
  }

  #[test]
  fn test_add() {
    let acc = init_acc::<DummyRSA>();
    let new_elems = [&big(5), &big(7), &big(11)];
    let (new_acc, poe) = add::<DummyRSA>(acc.clone(), &new_elems);
    let expected_acc = DummyRSA::exp(&DummyRSA::base_elem(), &big(94_125_955));
    assert!(new_acc == expected_acc);
    assert!(PoE::verify(&acc, &big(385), &new_acc, &poe));
  }

  #[test]
  fn test_delete() {
    let acc = init_acc::<DummyRSA>();
    let y_witness = DummyRSA::exp(&DummyRSA::base_elem(), &big(3649));
    let z_witness = DummyRSA::exp(&DummyRSA::base_elem(), &big(2747));
    let (new_acc, poe) =
      delete::<DummyRSA>(acc.clone(), &[&(big(67), y_witness), &(big(89), z_witness)])
        .expect("valid delete expected");
    let expected_acc = DummyRSA::exp(&DummyRSA::base_elem(), &big(41));
    assert!(new_acc == expected_acc);
    assert!(PoE::verify(&new_acc, &big(5963), &acc, &poe));
  }

  #[should_panic(expected = "BadWitness")]
  #[test]
  fn test_delete_bad_witness() {
    let acc = init_acc::<DummyRSA>();
    let y_witness = DummyRSA::exp(&DummyRSA::base_elem(), &big(3648));
    let z_witness = DummyRSA::exp(&DummyRSA::base_elem(), &big(2746));
    delete::<DummyRSA>(acc, &[&(big(67), y_witness), &(big(89), z_witness)]).unwrap();
  }

  #[test]
  fn test_prove_nonmembership() {
    let acc = init_acc::<DummyRSA>();
    let acc_set = [&big(41), &big(67), &big(89)];
    let elems = [&big(5), &big(7), &big(11)];
    let proof =
      prove_nonmembership::<DummyRSA>(&acc, &acc_set, &elems).expect("valid proof expected");
    assert!(verify_nonmembership::<DummyRSA>(&acc, &elems, &proof));
  }

  #[should_panic(expected = "InputsNotCoPrime")]
  #[test]
  fn test_prove_nonmembership_failure() {
    let acc = init_acc::<DummyRSA>();
    let acc_set = [&big(41), &big(67), &big(89)];
    let elems = [&big(41), &big(7), &big(11)];
    prove_nonmembership::<DummyRSA>(&acc, &acc_set, &elems).unwrap();
  }
}
