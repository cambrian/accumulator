use accumulator::group::Rsa2048;
use accumulator::Accumulator;
use rand::Rng;

/// Adds 10,000 random primes to accumulator (unverified), then tests 100 more random additions
/// (with verification) and 100 random elements are verified to be nonmembers.
/// Takes about 5 minutes.
#[test]
#[ignore]
fn stress_test() {
  let mut acc_set = Vec::new();
  let mut acc = Accumulator::<Rsa2048, [u8; 32]>::empty();
  for _ in 0..100 {
    let random_elem = rand::thread_rng().gen::<[u8; 32]>();
    acc_set.push(random_elem);
  }
  println!("Starting add");
  let (holder, _) = acc.clone().add(&acc_set);
  acc = holder;
  println!("{}", acc_set.len());
  for _ in 0..100 {
    let new_elem = rand::thread_rng().gen::<[u8; 32]>();
    assert!(!acc_set.contains(&new_elem));
    let (new_acc, add_proof) = acc.clone().add(&[new_elem]);
    assert!(new_acc.verify_membership(&[new_elem], &add_proof));
    let (_, del_proof) = new_acc.clone().delete(&[(new_elem, acc.clone())]).unwrap();
    assert!(new_acc.verify_membership(&[new_elem], &del_proof));
    let nonmem_proof = acc
      .prove_nonmembership(&acc_set, &[new_elem])
      .expect("It works");
    assert!(acc.verify_nonmembership(&[new_elem], &nonmem_proof));
    for _ in 0..100 {
      let new_elem_2 = rand::thread_rng().gen::<[u8; 32]>();
      let false_witness = Accumulator::<Rsa2048, [u8; 32]>::empty()
        .add(&[new_elem_2])
        .0;
      assert!(acc.prove_membership(&[(new_elem, false_witness)]).is_err());
    }
  }
}
