use crypto::accumulator::{
  add, delete, prove_membership, prove_nonmembership, setup, verify_membership,
  verify_nonmembership,
};
use crypto::group::{Group, UnknownOrderGroup, RSA2048};
use crypto::hash::hash_to_prime;
use crypto::util::int;
use rand::Rng;

fn init_acc<G: UnknownOrderGroup>() -> G::Elem {
  G::exp(&setup::<G>(), &(int(41) * &int(67) * &int(89)))
}
#[test]
fn test_add() {
  let acc = init_acc::<RSA2048>();
  let new_elems = [int(5), int(7), int(11)];
  let (new_acc, proof) = add::<RSA2048>(acc.clone(), &new_elems);
  let expected_acc = RSA2048::exp(&RSA2048::unknown_order_elem(), &int(94_125_955));
  assert!(new_acc == expected_acc);
  assert!(verify_membership(&new_acc, &new_elems, &proof));
}

/// Adds 10,000 random primes to accumulator (unverified), then tests 100 more random additions
/// (with verification) and 100 random elements are verified to be nonmembers.
/// Takes about 5 minutes
#[test]
fn stress_test() {
  let mut acc_set = Vec::new();
  let mut acc = setup::<RSA2048>();
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&random_bytes);
    acc_set.push(prime);
  }
  println!("Starting add");
  let (holder, _) = add::<RSA2048>(acc.clone(), &acc_set);
  acc = holder;
  println!("{}", acc_set.len());
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&random_bytes);
    assert!(!acc_set.contains(&prime));
    let (new_acc, add_proof) = add::<RSA2048>(acc.clone(), &[prime.clone()]);
    assert!(verify_membership(&new_acc, &[prime.clone()], &add_proof));
    let (_, del_proof) =
      delete::<RSA2048>(new_acc.clone(), &[(prime.clone(), acc.clone())]).unwrap();
    assert!(verify_membership(&new_acc, &[prime.clone()], &del_proof));
    let nonmem_proof =
      prove_nonmembership::<RSA2048>(&acc, &acc_set, &[prime.clone()]).expect("It works");
    assert!(verify_nonmembership::<RSA2048>(
      &acc,
      &[prime.clone()],
      &nonmem_proof
    ));
    for _ in 0..100 {
      let random_bytes_2 = rand::thread_rng().gen::<[u8; 32]>();
      let random_exp = hash_to_prime(&random_bytes_2);
      let false_witness = RSA2048::exp(&RSA2048::unknown_order_elem(), &random_exp);
      assert!(prove_membership::<RSA2048>(&acc, &[(prime.clone(), false_witness)]).is_err());
    }
  }
}
