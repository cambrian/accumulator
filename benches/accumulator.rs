/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::accumulator::{add, setup};
use crypto::group::RSA2048;
use crypto::hash::{hash_to_prime, Blake2b};
use num::bigint::BigUint;
use rand::Rng;

fn bench_add(elems: &[BigUint]) {
  let acc = setup::<RSA2048>();
  add::<RSA2048>(acc, elems);
}

fn criterion_benchmark(c: &mut Criterion) {
  let mut elems = Vec::new();
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&Blake2b::default, &random_bytes);
    elems.push(prime);
  }
  let elems_1 = [elems[0].clone()];
  let elems_2 = elems.clone();
  let elems_3 = elems.clone();
  c.bench_function("add_1", move |b| b.iter(|| bench_add(&elems_1)));
  c.bench_function("add_10", move |b| b.iter(|| bench_add(&elems_2[0..10])));
  c.bench_function("add_100", move |b| b.iter(|| bench_add(&elems_3)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
