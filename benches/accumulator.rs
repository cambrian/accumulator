/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::accumulator::{add, setup};
use crypto::group::DummyRSA;
use crypto::hash::{hash_to_prime, Blake2b};
use num::bigint::BigUint;
use rand::Rng;

fn bench_add(elems: &[BigUint]) {
  let acc = setup::<DummyRSA>();
  add::<DummyRSA>(acc, elems);
}

fn bench_iterative_add(elems: &[BigUint]) {
  let mut acc = setup::<DummyRSA>();
  for elem in elems.chunks(1) {
    acc = add::<DummyRSA>(acc, elem).0;
  }
}

fn criterion_benchmark(c: &mut Criterion) {
  let mut elems = Vec::new();
  for _ in 0..10 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&Blake2b::default, &random_bytes);
    elems.push(prime);
  }
  let elems_2 = elems.clone();
  c.bench_function("add", move |b| b.iter(|| bench_add(&elems)));
  c.bench_function("iterative_add", move |b| {
    b.iter(|| bench_iterative_add(&elems_2))
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
