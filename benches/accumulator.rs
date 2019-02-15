/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::group::{Rsa2048, UnknownOrderGroup};
use accumulator::hash::hash_to_prime;
use accumulator::util::int;
use accumulator::{Accumulator, MembershipProof};
use criterion::Criterion;
use rand::Rng;
use rug::Integer;

fn bench_delete<G: UnknownOrderGroup>(acc: Accumulator<G>, witness: &[(Integer, Accumulator<G>)]) {
  acc.delete(witness).expect("valid delete expected");
}

#[allow(dead_code)]
fn bench_add(elems: &[Integer]) {
  let acc = Accumulator::<Rsa2048>::new();
  acc.add(elems);
}

#[allow(dead_code)]
fn bench_verify<G: UnknownOrderGroup>(
  acc: &Accumulator<G>,
  elems: &[Integer],
  proof: &MembershipProof<G>,
) {
  assert!(acc.verify_membership(elems, proof));
}

#[allow(dead_code)]
fn bench_iterative_add(elems: &[Integer]) {
  let mut acc = Accumulator::<Rsa2048>::new();
  for elem in elems.chunks(1) {
    acc = acc.add(elem).0;
  }
}

fn init_acc<G: UnknownOrderGroup>() -> (Accumulator<G>, MembershipProof<G>, Vec<Integer>) {
  let mut elems = Vec::new();
  for _ in 0..50 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&random_bytes);
    elems.push(prime);
  }
  let acc = Accumulator::<G>::new();
  let (acc, proof) = acc.clone().add(&elems);
  // for _ in 0..100 {
  // elems = vec![];
  // for _ in 0..100 {
  //   let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  //   let prime = hash_to_prime(&random_bytes);
  //   elems.push(prime);
  // }
  // let (curr_acc, curr_proof) = acc.add(&elems);
  // acc = curr_acc;
  // proof = curr_proof;
  // }
  (acc, proof, elems)
}

fn criterion_benchmark(c: &mut Criterion) {
  // Test verification on lots of elements. Added in batches to not go crazy with exponent size.
  let (acc, proof, elems) = init_acc::<Rsa2048>();
  let elems_1 = elems.clone();
  let elems_2 = elems.clone();
  let elems_3 = elems.clone();
  let elems_4 = elems.clone();

  // c.bench_function("add_1", move |b| b.iter(|| bench_add(&elems_1[0..1])));
  // c.bench_function("add_10", move |b| b.iter(|| bench_add(&elems_2[0..10])));
  // c.bench_function("add_100", move |b| b.iter(|| bench_add(&elems_3)));
  // c.bench_function("verify", move |b| {
  //   b.iter(|| bench_verify(&acc, &elems, &proof))
  // });
  let mut delete_arg = Vec::new();
  for elem in elems_4 {
    let exp: Integer = elems.iter().product::<Integer>() / elem.clone();
    let witness = Accumulator::<Rsa2048>::new();
    let witness = witness.exp_quotient(exp, int(1)).unwrap();
    delete_arg.push((elem, witness));
  }
  let acc_2 = acc.clone();
  c.bench_function("delete_50", move |b| {
    b.iter(|| bench_delete(acc_2.clone(), &delete_arg[..]))
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
