//! Accumulator library, built on a generic group interface.
use crate::group::UnknownOrderGroup;
use crate::hash::hash_to_prime;
use crate::proof::{Poe, Poke2};
use crate::util::{divide_and_conquer, int, prime_hash_product, shamir_trick};
use rug::Integer;
use std::hash::Hash;
use std::marker::PhantomData;

#[derive(Debug)]
pub enum AccError {
  BadWitness,
  BadWitnessUpdate,
  DivisionByZero,
  InexactDivision,
  InputsNotCoprime,
}

// See https://doc.rust-lang.org/std/marker/struct.PhantomData.html#ownership-and-the-drop-check
// for recommendations re phantom types.
#[derive(PartialEq, Eq, Debug, Hash)]
pub struct Accumulator<G: UnknownOrderGroup, T: Hash + Eq> {
  phantom: PhantomData<*const T>,
  value: G::Elem,
}

// Manual clone impl required because Rust's type inference is not good. See
// https://github.com/rust-lang/rust/issues/26925
impl<G: UnknownOrderGroup, T: Hash + Eq> Clone for Accumulator<G, T> {
  fn clone(&self) -> Self {
    Accumulator {
      phantom: PhantomData,
      value: self.value.clone(),
    }
  }
}

#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct Witness<G: UnknownOrderGroup, T: Hash + Eq>(Accumulator<G, T>);

// // Manual clone impl required because Rust's type inference is not good. See
// // https://github.com/rust-lang/rust/issues/26925
// impl<G: UnknownOrderGroup, T: Hash + Eq> Clone for Witness<G, T> {
//   fn clone(&self) -> Self {
//     Witness(self.0.clone())
//   }
// }

#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct MembershipProof<G: UnknownOrderGroup, T: Hash + Eq> {
  pub witness: Witness<G, T>,
  proof: Poe<G>,
}

#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct NonmembershipProof<G: UnknownOrderGroup, T: Hash + Eq> {
  phantom: PhantomData<*const T>,
  d: G::Elem,
  v: G::Elem,
  gv_inv: G::Elem,
  poke2_proof: Poke2<G>,
  poe_proof: Poe<G>,
}

impl<G: UnknownOrderGroup, T: Hash + Eq> Accumulator<G, T> {
  /// Create a new, empty accumulator
  pub fn empty() -> Self {
    Accumulator {
      phantom: PhantomData,
      value: G::unknown_order_elem(),
    }
  }

  /// Given a witness for old_witness_set, returns a witness for new_witness_set
  /// new_witness_set must be a subset of old_witness_set.
  /// REVIEW: explicitly check that new_witness_set is a subset of old_witness_set, and improve
  /// the error type
  pub fn compute_sub_witness(
    self,
    old_witness_set: &[T],
    new_witness_set: &[T],
  ) -> Result<Self, AccError> {
    let numerator = prime_hash_product(old_witness_set);
    let denominator = prime_hash_product(new_witness_set);

    let (quotient, remainder) = numerator.div_rem(denominator);

    if remainder != int(0) {
      return Err(AccError::InexactDivision);
    }

    Ok(Accumulator {
      phantom: PhantomData,
      value: G::exp(&self.value, &quotient),
    })
  }

  // Internal fn that also returns the prime hash product of the added elements, to enable
  // efficient add_with_proof.
  // Uses a move instead of a `&self` reference to prevent accidental use of the old accumulator
  // state.
  fn add_(self, elems: &[T]) -> (Self, Integer, Witness<G, T>) {
    let x = prime_hash_product(elems);
    let acc_elem = G::exp(&self.value, &x);
    (
      Accumulator {
        phantom: PhantomData,
        value: acc_elem,
      },
      x,
      Witness(self),
    )
  }

  // The conciseness of `accumulator.add()` and low probability of confusion with implementations of
  // the `Add` trait probably justify this...
  #[allow(clippy::should_implement_trait)]
  /// Adds `elems` to the accumulator `acc`. Cannot check whether the elements have not already been
  /// added. It is up to clients to either ensure uniqueness or treat this as multiset.
  pub fn add(self, elems: &[T]) -> (Self, Witness<G, T>) {
    let (acc, _, witness) = self.add_(elems);
    (acc, witness)
  }

  /// Adds `elems` to the accumulator `acc`. Cannot check whether the elements have not already been
  /// added. It is up to clients to either ensure uniqueness or treat this as multiset.
  /// Also returns a batch membership proof for elems in the new accumulator.
  pub fn add_with_proof(self, elems: &[T]) -> (Self, MembershipProof<G, T>) {
    let (acc, x, witness) = self.add_(elems);
    let proof = Poe::<G>::prove(&witness.0.value, &x, &acc.value);
    (acc, MembershipProof { witness, proof })
  }

  // Internal fn that also returns the prime hash product of the deleted elements, to enable
  // efficient delete_with_proof.
  // Uses a divide-and-conquer approach to running the ShamirTrick, which keeps the average input
  // smaller: For `[a, b, c, d]` do `S(S(a, b), S(c, d))` instead of `S(S(S(a, b), c), d)`.
  // Uses a move instead of a `&self` reference to prevent accidental use of the old accumulator
  // state.
  pub fn delete_(self, elem_witnesses: &[(T, Witness<G, T>)]) -> Result<(Self, Integer), AccError> {
    let prime_witnesses = elem_witnesses
      .iter()
      .map(|(elem, witness)| (hash_to_prime(elem), witness.0.value.clone()))
      .collect::<Vec<_>>(); // doesn't cooperate when we try to collect to a &[_]

    for (p, witness_elem) in &prime_witnesses {
      if G::exp(&witness_elem, &p) != self.value {
        return Err(AccError::BadWitness);
      }
    }

    let (prime_product, acc_elem) = divide_and_conquer(
      |(p1, v1), (p2, v2)| Ok((int(p1 * p2), shamir_trick::<G>(&v1, &v2, p1, p2).unwrap())),
      (int(1), self.value),
      &prime_witnesses[..],
    )?;

    Ok((
      Accumulator {
        phantom: PhantomData,
        value: acc_elem.clone(),
      },
      prime_product,
    ))
  }

  /// Removes the elements in `elem_witnesses` from the accumulator.
  pub fn delete(self, elem_witnesses: &[(T, Witness<G, T>)]) -> Result<Self, AccError> {
    Ok(self.delete_(elem_witnesses)?.0)
  }

  /// Removes the elements in `elem_witnesses` from the accumulator. Also returns inclusion proof
  /// for the deleted elements in the original accumulator, using the new accumulator as witness.
  pub fn delete_with_proof(
    self,
    elem_witnesses: &[(T, Witness<G, T>)],
  ) -> Result<(Self, MembershipProof<G, T>), AccError> {
    let (acc, prime_product) = self.clone().delete_(elem_witnesses)?;
    let proof = Poe::<G>::prove(&acc.value, &prime_product, &self.value);
    Ok((
      acc.clone(),
      MembershipProof {
        witness: Witness(acc),
        proof,
      },
    ))
  }

  /// Compute the batch MembershipProof for `elem_witnesses`
  pub fn prove_membership(
    &self,
    elem_witnesses: &[(T, Witness<G, T>)],
  ) -> Result<MembershipProof<G, T>, AccError> {
    let witness_accum = self.clone().delete(elem_witnesses)?;
    let prod = elem_witnesses
      .iter()
      .map(|(t, _)| hash_to_prime(t))
      .product();
    let proof = Poe::<G>::prove(&witness_accum.value, &prod, &self.value);
    Ok(MembershipProof {
      witness: Witness(witness_accum),
      proof,
    })
  }

  pub fn verify_membership(
    &self,
    t: &T,
    MembershipProof { witness, proof }: &MembershipProof<G, T>,
  ) -> bool {
    let exp = hash_to_prime(t);
    Poe::verify(&witness.0.value, &exp, &self.value, proof)
  }

  pub fn verify_aggregate_membership(
    &self,
    elems: &[T],
    MembershipProof { witness, proof }: &MembershipProof<G, T>,
  ) -> bool {
    let exp = prime_hash_product(elems);
    Poe::verify(&witness.0.value, &exp, &self.value, proof)
  }

  /// See Section 4.2 in the Li, Li, Xue paper.
  pub fn update_membership_witness(
    self,
    acc_new: &Self,
    tracked_elems: &[T],
    untracked_additions: &[T],
    untracked_deletions: &[T],
  ) -> Result<Self, AccError> {
    let x = prime_hash_product(tracked_elems);
    let x_hat = prime_hash_product(untracked_deletions);

    for elem in tracked_elems {
      if untracked_additions.contains(elem) || untracked_deletions.contains(elem) {
        return Err(AccError::BadWitnessUpdate);
      }
    }

    let (gcd, a, b) = <(Integer, Integer, Integer)>::from(x.gcd_cofactors_ref(&x_hat));
    assert!(gcd == int(1));

    let w = self.add(untracked_additions).0;
    let w_to_b = G::exp(&w.value, &b);
    let acc_new_to_a = G::exp(&acc_new.value, &a);
    Ok(Accumulator {
      phantom: PhantomData,
      value: G::op(&w_to_b, &acc_new_to_a),
    })
  }

  /// Returns a proof (and associated variables) that `elems` are not in `acc_set`.
  pub fn prove_nonmembership(
    &self,
    acc_set: &[T],
    elems: &[T],
  ) -> Result<NonmembershipProof<G, T>, AccError> {
    let x: Integer = elems.iter().map(hash_to_prime).product();
    let s = acc_set.iter().map(hash_to_prime).product();
    let (gcd, a, b) = <(Integer, Integer, Integer)>::from(x.gcd_cofactors_ref(&s));

    if gcd != int(1) {
      return Err(AccError::InputsNotCoprime);
    }

    let g = G::unknown_order_elem();
    let d = G::exp(&g, &a);
    let v = G::exp(&self.value, &b);
    let gv_inv = G::op(&g, &G::inv(&v));

    let poke2_proof = Poke2::prove(&self.value, &b, &v);
    let poe_proof = Poe::prove(&d, &x, &gv_inv);
    Ok(NonmembershipProof {
      phantom: PhantomData,
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
    elems: &[T],
    NonmembershipProof {
      d,
      v,
      gv_inv,
      poke2_proof,
      poe_proof,
      ..
    }: &NonmembershipProof<G, T>,
  ) -> bool {
    let x = elems.iter().map(hash_to_prime).product();
    Poke2::verify(&self.value, v, poke2_proof) && Poe::verify(d, &x, gv_inv, poe_proof)
  }

  /// For accumulator with elems `[x_1, ..., x_n]`, computes a membership witness for each `x_i` in
  /// accumulator `g^{x_1 * ... * x_n}`, namely `g^{x_1 * ... * x_n / x_i}`, in O(N
  /// log N) time using the root factor algorithm.
  pub fn compute_individual_witnesses<'a>(elems: &'a [T]) -> Vec<(&'a T, Witness<G, T>)> {
    let primes = elems.iter().map(hash_to_prime).collect::<Vec<_>>();
    let witnesses = Self::root_factor(&G::unknown_order_elem(), &primes);
    // why is it necessary to split this calculation into 2 lines??
    let witnesses = witnesses.iter().map(|value| {
      Witness(Accumulator {
        phantom: PhantomData,
        value: value.clone(),
      })
    });
    elems.iter().zip(witnesses).collect::<Vec<_>>()
  }

  #[allow(non_snake_case)]
  fn root_factor(g: &G::Elem, primes: &[Integer]) -> Vec<G::Elem> {
    dbg!((&g, &primes));
    if primes.len() == 1 {
      return vec![g.clone()];
    }
    let half_n = primes.len() / 2;
    let g_l = primes[..half_n]
      .iter()
      .fold(g.clone(), |g_, x| G::exp(&g_, x));
    let g_r = primes[half_n..]
      .iter()
      .fold(g.clone(), |g_, x| G::exp(&g_, x));
    let mut L = Self::root_factor(&g_r, &primes[..half_n]);
    let mut R = Self::root_factor(&g_l, &primes[half_n..]);
    L.append(&mut R);
    L
  }
}

impl<G: UnknownOrderGroup, T: Hash + Eq> From<&[T]> for Accumulator<G, T> {
  fn from(ts: &[T]) -> Self {
    Accumulator::<G, T>::empty().add(ts).0
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{ClassGroup, ElemFrom, Rsa2048};

  // For some reason this doesn't work with `from`
  fn new_acc<G: UnknownOrderGroup, T: Hash + Eq>(data: &[T]) -> Accumulator<G, T> {
    Accumulator::<G, T>::empty().add(data).0
  }

  macro_rules! test_all_groups {
    ($test_func:ident, $func_name_rsa:ident, $func_name_class:ident, $($attr:meta)*) => {
      #[test]
      $(
        #[$attr]
      )*
      fn $func_name_rsa() {
        $test_func::<Rsa2048>();
      }

      #[test]
      $(
        #[$attr]
      )*
      fn $func_name_class() {
        $test_func::<ClassGroup>();
      }
    };
  }

  test_all_groups!(
    test_compute_sub_witness,
    test_compute_sub_witness_rsa2048,
    test_compute_sub_witness_class,
  );
  fn test_compute_sub_witness<G: UnknownOrderGroup>() {
    let empty_acc = Accumulator::<G, &'static str>::empty();
    let sub_witness = empty_acc.compute_sub_witness(&["a", "b"], &["a"]).unwrap();
    let exp_quotient_expected = new_acc::<G, &'static str>(&["b"]);
    assert!(sub_witness == exp_quotient_expected);
  }

  test_all_groups!(test_add, test_add_rsa2048, test_add_class,);
  fn test_add<G: UnknownOrderGroup>() {
    let acc = new_acc::<G, &'static str>(&["a", "b"]);
    let new_elems = ["c", "d"];
    let (acc_new, proof) = acc.add_with_proof(&new_elems);
    let acc_expected = G::exp(
      &G::unknown_order_elem(),
      &prime_hash_product(&["a", "b", "c", "d"]),
    );
    assert!(acc_new.value == acc_expected);
    assert!(acc_new.verify_aggregate_membership(&new_elems, &proof));
  }

  test_all_groups!(test_delete, test_delete_rsa2048, test_delete_class,);
  fn test_delete<G: UnknownOrderGroup>() {
    let acc_0 = new_acc::<G, &'static str>(&["a", "b"]);
    let (acc_1, c_witness) = acc_0.clone().add(&["c"]);
    let (acc_2, proof) = acc_1
      .clone()
      .delete_with_proof(&[("c", c_witness)])
      .expect("valid delete expected");
    assert!(acc_2 == acc_0);
    assert!(acc_1.verify_membership(&"c", &proof));
  }

  test_all_groups!(
    test_delete_empty,
    test_delete_empty_rsa2048,
    test_delete_empty_class,
  );
  fn test_delete_empty<G: UnknownOrderGroup>() {
    let acc = new_acc::<G, &'static str>(&["a", "b"]);
    let (acc_new, proof) = acc
      .clone()
      .delete_with_proof(&[])
      .expect("valid delete expected");
    assert!(acc_new == acc);
    assert!(acc.verify_aggregate_membership(&[], &proof));
  }

  // test_all_groups!(
  //   test_delete_bad_witness,
  //   test_delete_bad_witness_rsa2048,
  //   test_delete_bad_witness_class,
  //   should_panic(expected = "BadWitness")
  // );
  // fn test_delete_bad_witness<G: UnknownOrderGroup>() {
  //   let acc = Accumulator::<G, &'static str>::empty();
  //   let a_witness = new_acc::<G, &'static str>(&["b", "c"]);
  //   let b_witness = new_acc::<G, &'static str>(&["a", "c"]);
  //   acc.delete(&[("a", a_witness), ("b", b_witness)]).unwrap();
  // }

  test_all_groups!(
    test_update_membership_witness,
    test_update_membership_witness_rsa2048,
    test_update_membership_witness_class,
  );
  fn test_update_membership_witness<G: UnknownOrderGroup>() {
    let acc = new_acc::<G, &'static str>(&["a", "b", "c"]);
    let witness = new_acc::<G, &'static str>(&["c", "d"]);
    let witness_new = witness
      .update_membership_witness(&acc, &["a"], &["b"], &["d"])
      .unwrap();
    assert!(witness_new.add(&["a"]).0 == acc);
  }

  test_all_groups!(
    test_update_membership_witness_failure,
    test_update_membership_witness_failure_rsa2048,
    test_update_membership_witness_failure_class,
    should_panic(expected = "BadWitnessUpdate")
  );
  fn test_update_membership_witness_failure<G: UnknownOrderGroup>() {
    let acc = new_acc::<G, &'static str>(&["a", "b", "c"]);
    let witness = new_acc::<G, &'static str>(&["c", "d"]);
    witness
      .update_membership_witness(&acc, &["a"], &["b"], &["a"])
      .unwrap();
  }

  test_all_groups!(
    test_prove_nonmembership,
    test_prove_nonmembership_rsa2048,
    test_prove_nonmembership_class,
  );
  fn test_prove_nonmembership<G: UnknownOrderGroup>() {
    let acc_set = ["a", "b"];
    let acc = new_acc::<G, &'static str>(&acc_set);
    let non_members = ["c", "d"];
    let proof = acc
      .prove_nonmembership(&acc_set, &non_members)
      .expect("valid proof expected");
    assert!(acc.verify_nonmembership(&non_members, &proof));
  }

  fn test_compute_individual_witnesses<G: UnknownOrderGroup + ElemFrom<u32>>() {
    let elems = ["a", "b", "c"];
    let acc = new_acc::<G, &'static str>(&elems);
    let witnesses = Accumulator::<G, &'static str>::compute_individual_witnesses(&elems);
    for (elem, witness) in witnesses {
      assert_eq!(acc.value, G::exp(&witness.0.value, &hash_to_prime(*elem)));
    }
  }

  #[test]
  fn test_compute_individual_witnesses_rsa2048() {
    // Class version takes too long for a unit test.
    test_compute_individual_witnesses::<Rsa2048>();
  }
}
