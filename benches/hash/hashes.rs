/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use accumulator::hash::{blake2b, hash_to_prime};
use rand::Rng;

fn bench_blake2() {
  blake2b("werg");
}

fn bench_hash_to_prime() {
  let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  hash_to_prime(&random_bytes);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("blake2", |b| b.iter(bench_blake2));
  c.bench_function("hash_to_prime", |b| b.iter(bench_hash_to_prime));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
