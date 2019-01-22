/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::hash::{blake2, sha256, h_prime};
use rand::Rng;

fn bench_blake2() {
  blake2(b"test", None);
}

fn bench_sha256() {
  sha256(b"test", None);
}

fn bench_h_prime() {
  let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  h_prime(&blake2, &random_bytes);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("blake2", |b| b.iter(|| bench_blake2()));
  c.bench_function("sha256", |b| b.iter(|| bench_sha256()));
  c.bench_function("h_prime", |b| b.iter(|| bench_h_prime()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
