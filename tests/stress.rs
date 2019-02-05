use crypto::accumulator::Accumulator;
use crypto::group::RSA2048;
use crypto::hash::hash_to_prime;
use rand::Rng;

/// Adds 10,000 random primes to accumulator (unverified), then tests 100 more random additions
/// (with verification) and 100 random elements are verified to be nonmembers.
/// Takes about 5 minutes. (TODO: Disable by default since it takes so long?)
// #[test]
#[allow(dead_code)]
fn stress_test() {
  let mut acc_set = Vec::new();
  let mut acc = Accumulator::<RSA2048>::new();
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&random_bytes);
    acc_set.push(prime);
  }
  println!("Starting add");
  let (holder, _) = acc.clone().add(&acc_set);
  acc = holder;
  println!("{}", acc_set.len());
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&random_bytes);
    assert!(!acc_set.contains(&prime));
    let (new_acc, add_proof) = acc.clone().add(&[prime.clone()]);
    assert!(new_acc.verify_membership(&[prime.clone()], &add_proof));
    let (_, del_proof) = new_acc
      .clone()
      .delete(&[(prime.clone(), acc.clone())])
      .unwrap();
    assert!(new_acc.verify_membership(&[prime.clone()], &del_proof));
    let nonmem_proof = acc
      .prove_nonmembership(&acc_set, &[prime.clone()])
      .expect("It works");
    assert!(acc.verify_nonmembership(&[prime.clone()], &nonmem_proof));
    for _ in 0..100 {
      let random_bytes_2 = rand::thread_rng().gen::<[u8; 32]>();
      let random_exp = hash_to_prime(&random_bytes_2);
      let false_witness = Accumulator::<RSA2048>::new().add(&[random_exp]).0;
      assert!(acc
        .prove_membership(&[(prime.clone(), false_witness)])
        .is_err());
    }
  }
}
