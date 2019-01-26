/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::accumulator;
use crypto::accumulator::MembershipProof;
use crypto::group::{UnknownOrderGroup, RSA2048};
use crypto::hash::{hash_to_prime, Blake2b};
use rand::Rng;
use rug::Integer;

#[allow(dead_code)]
fn bench_delete<G: UnknownOrderGroup>(acc: G::Elem, witness: &[(Integer, G::Elem)]) {
  accumulator::delete::<G>(acc, witness).expect("valid delete expected");
}

fn bench_add(elems: &[Integer]) {
  let acc = accumulator::setup::<RSA2048>();
  accumulator::add::<RSA2048>(acc, elems);
}

fn bench_verify<G: UnknownOrderGroup>(
  acc: &G::Elem,
  elems: &[Integer],
  proof: &MembershipProof<G>,
) {
  assert!(accumulator::verify_membership::<G>(acc, elems, proof));
}

#[allow(dead_code)]
fn bench_iterative_add(elems: &[Integer]) {
  let mut acc = accumulator::setup::<RSA2048>();
  for elem in elems.chunks(1) {
    acc = accumulator::add::<RSA2048>(acc, elem).0;
  }
}

fn criterion_benchmark(c: &mut Criterion) {
  // Setup (not included in the benchmark time).
  // REVIEW: Consider breaking up setup into functions (if it's not annoying)?
  let mut elems = Vec::new();
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&Blake2b::default, &random_bytes);
    elems.push(prime);
  }
  let elems_1 = [elems[0].clone()];
  let elems_2 = elems.clone();
  let elems_3 = elems.clone();
  let acc = accumulator::setup::<RSA2048>();
  let mut new_acc;
  let mut proof;
  let (holder, proof_holder) = accumulator::add::<RSA2048>(acc.clone(), &elems.clone());
  new_acc = holder;
  proof = proof_holder;
  // Test verification on lots of elements. Added in batches to not go crazy with exponent size.
  for _ in 0..100 {
    elems = vec![];
    for _ in 0..100 {
      let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
      let prime = hash_to_prime(&Blake2b::default, &random_bytes);
      elems.push(prime);
    }
    let (curr_acc, curr_proof) = accumulator::add::<RSA2048>(new_acc.clone(), &elems.clone());
    new_acc = curr_acc;
    proof = curr_proof;
  }

  c.bench_function("add_1", move |b| b.iter(|| bench_add(&elems_1)));
  c.bench_function("add_10", move |b| b.iter(|| bench_add(&elems_2[0..10])));
  c.bench_function("add_100", move |b| b.iter(|| bench_add(&elems_3)));
  c.bench_function("verify", move |b| {
    b.iter(|| bench_verify(&new_acc, &elems, &proof))
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
