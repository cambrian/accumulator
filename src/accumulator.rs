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
#[derive(PartialEq, Eq, Debug, Default, Hash)]
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
pub struct MembershipProof<G: UnknownOrderGroup, T: Hash + Eq> {
  pub witness: Accumulator<G, T>,
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
    new_witness_set: &[T],
    old_witness_set: &[T],
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

  // The conciseness of `accumulator.add()` and low probability of confusion with implementations of
  // the `Add` trait probably justify this...
  #[allow(clippy::should_implement_trait)]
  /// Adds `elems` to the accumulator `acc`. Cannot check whether the elements are coprime with the
  /// accumulator, but it is up to clients to either ensure uniqueness or treat this as multiset.
  // Uses a move instead of a `&self` reference to prevent accidental use of the old accumulator
  // state.
  pub fn add(self, elems: &[T]) -> (Self, MembershipProof<G, T>) {
    let x = prime_hash_product(elems);
    let acc_new = G::exp(&self.value, &x);
    let poe_proof = Poe::<G>::prove(&self.value, &x, &acc_new);
    (
      Accumulator {
        phantom: PhantomData,
        value: acc_new,
      },
      MembershipProof {
        witness: self,
        proof: poe_proof,
      },
    )
  }

  /// Removes the elements in `elem_witnesses` from the accumulator `acc`.
  /// Uses a divide-and-conquer approach to running the ShamirTrick, which keeps the average input
  /// smaller: For `[a, b, c, d]` do `S(S(a, b), S(c, d))` instead of `S(S(S(a, b), c), d)`.
  // Uses a move instead of a `&self` reference to prevent accidental use of the old accumulator
  // state.
  pub fn delete(
    self,
    elem_witnesses: &[(T, Self)],
  ) -> Result<(Self, MembershipProof<G, T>), AccError> {
    let prime_witnesses = elem_witnesses
      .iter()
      .map(|(elem, witness)| (hash_to_prime(elem), witness.clone()))
      .collect::<Vec<_>>(); // doesn't cooperate when we try to collect to a &[_]
    for (p, witness) in &prime_witnesses {
      if G::exp(&witness.value, &p) != self.value {
        return Err(AccError::BadWitness);
      }
    }

    let (elem_aggregate, acc_next) = divide_and_conquer(
      |(p1, v1), (p2, v2)| {
        Ok((
          int(p1 * p2),
          Accumulator {
            phantom: PhantomData,
            value: shamir_trick::<G>(&v1.value, &v2.value, p1, p2)
              .ok_or(AccError::InputsNotCoprime)?,
          },
        ))
      },
      (int(1), self.clone()),
      &prime_witnesses[..],
    )?;

    let poe_proof = Poe::<G>::prove(&acc_next.value, &elem_aggregate, &self.value);
    Ok((
      acc_next.clone(),
      MembershipProof {
        witness: acc_next,
        proof: poe_proof,
      },
    ))
  }

  /// Returns a proof (and associated variables) that `elem_witnesses` are aggregated in `acc`.
  pub fn prove_membership(
    &self,
    elem_witnesses: &[(T, Self)],
  ) -> Result<MembershipProof<G, T>, AccError> {
    Ok(self.clone().delete(elem_witnesses)?.1)
  }

  /// Verifies the PoE returned by `prove_membership`.
  pub fn verify_membership(
    &self,
    elems: &[T],
    MembershipProof { witness, proof }: &MembershipProof<G, T>,
  ) -> bool {
    let exp = prime_hash_product(elems);
    Poe::verify(&witness.value, &exp, &self.value, proof)
  }

  /// Updates a membership witness for some set. See Section 4.2 in the Li, Li, Xue paper.
  pub fn update_membership_witness(
    self,
    acc_new: &Self,
    witness_set: &[T],
    untracked_additions: &[T],
    untracked_deletions: &[T],
  ) -> Result<Self, AccError> {
    let x = prime_hash_product(witness_set);
    let x_hat = prime_hash_product(untracked_deletions);

    for witness_elem in witness_set {
      if untracked_additions.contains(witness_elem) || untracked_deletions.contains(witness_elem) {
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

  /// For accumulator with value `g` and elems `[x_1, ..., x_n]`, computes a membership witness for
  /// each `x_i` in accumulator `g^{x_1 * ... * x_n}`, namely `g^{x_1 * ... * x_n / x_i}`, in O(N
  /// log N) time using the root factor algorithm.
  pub fn compute_all_witnesses<'a>(&self, elems: &'a [T]) -> Vec<(&'a T, Self)> {
    let primes = elems.iter().map(hash_to_prime).collect::<Vec<_>>();
    let witnesses = Self::root_factor(&self.value, &primes); // why is this necessary??
    let witnesses = witnesses.iter().map(|value| Accumulator {
      phantom: PhantomData,
      value: value.clone(),
    });
    elems.iter().zip(witnesses).collect::<Vec<_>>()
  }

  #[allow(non_snake_case)]
  fn root_factor(g: &G::Elem, primes: &[Integer]) -> Vec<G::Elem> {
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
  use crate::group::{ClassGroup, Rsa2048};
  use crate::hash;
  use crate::util::int;

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

  // test_all_groups!(
  //   test_exp_quotient,
  //   test_exp_quotient_rsa2048,
  //   test_exp_quotient_class,
  // );
  // fn test_exp_quotient<G: UnknownOrderGroup>() {
  //   let empty_acc = Accumulator::<G>::new();
  //   let exp_quotient_result = empty_acc
  //     .exp_quotient(int(17 * 41 * 67 * 89), int(17 * 89))
  //     .unwrap();
  //   let exp_quotient_expected = Accumulator(G::exp(&G::unknown_order_elem(), &int(41 * 67)));
  //   assert!(exp_quotient_result == exp_quotient_expected);
  // }

  // test_all_groups!(
  //   test_exp_quotient_zero,
  //   test_exp_quotient_zero_rsa2048,
  //   test_exp_quotient_zero_class,
  //   should_panic(expected = "DivisionByZero")
  // );
  // fn test_exp_quotient_zero<G: UnknownOrderGroup>() {
  //   Accumulator::<G>::new()
  //     .exp_quotient(int(17 * 41 * 67 * 89), int(0))
  //     .unwrap();
  // }

  // test_all_groups!(
  //   test_exp_quotient_inexact,
  //   test_exp_quotient_inexact_rsa2048,
  //   test_exp_quotient_inexact_class,
  //   should_panic(expected = "InexactDivision")
  // );
  // fn test_exp_quotient_inexact<G: UnknownOrderGroup>() {
  //   Accumulator::<G>::new()
  //     .exp_quotient(int(17 * 41 * 67 * 89), int(5))
  //     .unwrap();
  // }

  test_all_groups!(test_add, test_add_rsa2048, test_add_class,);
  fn test_add<G: UnknownOrderGroup>() {
    let acc = new_acc::<G, &'static str>(&["a", "b"]);
    let new_elems = ["c", "d"];
    let (acc_new, proof) = acc.add(&new_elems);
    let acc_expected = G::exp(
      &G::unknown_order_elem(),
      &prime_hash_product(&["a", "b", "c", "d"]),
    );
    assert!(acc_new.value == acc_expected);
    assert!(acc_new.verify_membership(&new_elems, &proof));
  }

  // test_all_groups!(test_delete, test_delete_rsa2048, test_delete_class,);
  // fn test_delete<G: UnknownOrderGroup>() {
  //   let acc = init_acc::<G, &'static str>(&["a", "b"]);
  //   let y_witness = init_acc::<G, &'static str>(&["a", "b"]);
  //   let z_witness = Accumulator::<G>::new().add(&[int(2747)]).0;
  //   let (acc_new, proof) = acc
  //     .clone()
  //     .delete(&[(int(67), y_witness), (int(89), z_witness)])
  //     .expect("valid delete expected");
  //   let acc_expected = G::exp(&G::unknown_order_elem(), &int(41));
  //   assert!(acc_new.0 == acc_expected);
  //   assert!(acc.verify_membership(&[int(67), int(89)], &proof));
  // }

  test_all_groups!(
    test_delete_empty,
    test_delete_empty_rsa2048,
    test_delete_empty_class,
  );
  fn test_delete_empty<G: UnknownOrderGroup>() {
    let acc = new_acc::<G, &'static str>(&["a", "b"]);
    let (acc_new, proof) = acc.clone().delete(&[]).expect("valid delete expected");
    assert!(acc_new == acc);
    assert!(acc.verify_membership(&[], &proof));
  }

  // test_all_groups!(
  //   test_delete_bad_witness,
  //   test_delete_bad_witness_rsa2048,
  //   test_delete_bad_witness_class,
  //   should_panic(expected = "BadWitness")
  // );
  // fn test_delete_bad_witness<G: UnknownOrderGroup>() {
  //   let acc = init_acc::<G>();
  //   let y_witness = Accumulator::<G>::new().add(&[int(3648)]).0;
  //   let z_witness = Accumulator::<G>::new().add(&[int(2746)]).0;
  //   acc
  //     .delete(&[(int(67), y_witness), (int(89), z_witness)])
  //     .unwrap();
  // }

  // test_all_groups!(
  //   test_update_membership_witness,
  //   test_update_membership_witness_rsa2048,
  //   test_update_membership_witness_class,
  // );
  // fn test_update_membership_witness<G: UnknownOrderGroup>() {
  //   // Original accumulator has [3, 5, 11, 13].
  //   // Witness is tracking elements [3, 5] and eventually [7].
  //   let acc_new = Accumulator::<G>::new()
  //     .add(&[int(3), int(7), int(11), int(17)])
  //     .0;
  //   let witness = Accumulator::<G>::new().add(&[int(11), int(13)]).0;
  //   let witness_new = witness
  //     .update_membership_witness(&acc_new, &[int(3), int(7)], &[int(17)], &[int(13)])
  //     .unwrap();
  //   assert!(witness_new.add(&[int(3), int(7)]).0 == acc_new);
  // }

  // test_all_groups!(
  //   test_update_membership_witness_failure,
  //   test_update_membership_witness_failure_rsa2048,
  //   test_update_membership_witness_failure_class,
  //   should_panic(expected = "BadWitnessUpdate")
  // );
  // fn test_update_membership_witness_failure<G: UnknownOrderGroup>() {
  //   let acc_new = Accumulator::<G>::new()
  //     .add(&[int(3), int(7), int(11), int(17)])
  //     .0;
  //   let witness = Accumulator::<G>::new().add(&[int(11), int(13)]).0;
  //   witness
  //     .update_membership_witness(&acc_new, &[int(3), int(7)], &[int(3)], &[int(13)])
  //     .unwrap();
  // }

  // test_all_groups!(
  //   test_prove_nonmembership,
  //   test_prove_nonmembership_rsa2048,
  //   test_prove_nonmembership_class,
  // );
  // fn test_prove_nonmembership<G: UnknownOrderGroup>() {
  //   let acc = init_acc::<G>();
  //   let acc_set = [int(41), int(67), int(89)];
  //   let elems = [int(5), int(7), int(11)];
  //   let proof = acc
  //     .prove_nonmembership(&acc_set, &elems)
  //     .expect("valid proof expected");
  //   assert!(acc.verify_nonmembership(&elems, &proof));
  // }

  // test_all_groups!(
  //   test_prove_nonmembership_failure,
  //   test_prove_nonmembership_failure_rsa2048,
  //   test_prove_nonmembership_failure_class,
  //   should_panic(expected = "InputsNotCoprime")
  // );
  // fn test_prove_nonmembership_failure<G: UnknownOrderGroup>() {
  //   let acc = init_acc::<G>();
  //   let acc_set = [int(41), int(67), int(89)];
  //   let elems = [int(41), int(7), int(11)];
  //   acc.prove_nonmembership(&acc_set, &elems).unwrap();
  // }

  // fn test_root_factor<G: UnknownOrderGroup>() {
  //   let acc = init_acc::<G>();
  //   let orig = [
  //     97 as usize,
  //     101 as usize,
  //     103 as usize,
  //     107 as usize,
  //     109 as usize,
  //   ];
  //   let factors: Vec<Integer> = orig.iter().map(|x| hash::hash_to_prime(x)).collect();
  //   let witnesses = acc.root_factor(&factors, &orig);
  //   for (i, (_, witness)) in witnesses.iter().enumerate() {
  //     let partial_product = factors.iter().product::<Integer>() / factors[i].clone();
  //     let expected = acc.clone().add(&[partial_product]).0;
  //     assert_eq!(*witness, expected);
  //   }
  // }

  // #[test]
  // fn test_root_factor_rsa2048() {
  //   // Class version takes too long for a unit test.
  //   test_root_factor::<Rsa2048>();
  // }
}
