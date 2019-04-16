/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::hash::blake2b;
use accumulator::uint::U256;
use criterion::{black_box, Criterion};
use rug::integer::Order;
use std::ops::Mul;

fn bench_mul<T: Mul>(a: T, b: T) {
  black_box(a * b);
}

fn criterion_benchmark(c: &mut Criterion) {
  let int = blake2b("data");
  let mut bytes = [0; 4];
  int.write_digits(&mut bytes, Order::LsfBe);
  let u256 = U256::from(bytes);
  c.bench_function("mul_rug", move |b| b.iter(|| bench_mul(&int, &int)));
  c.bench_function("mul_u256", move |b| b.iter(|| bench_mul(u256, u256)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
