/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::hash::primality::{passes_lucas, passes_miller_rabin_base_2};
use accumulator::uint::u256;
use criterion::Criterion;
use rand::Rng;
use rug::integer::Order;
use rug::Integer;

fn bench_mr2(bytes: &[u8; 32]) {
  passes_miller_rabin_base_2(&u256(bytes));
}

fn bench_mr2_rug(bytes: &[u8; 32]) {
  let n = Integer::from_digits(bytes, Order::Lsf);
  // GMP does not let us demand a base-2 Fermat test so we just do one of random base.
  n.is_probably_prime(1);
}

fn bench_lucas(bytes: &[u8; 32]) {
  passes_lucas(&u256(bytes));
}

fn criterion_benchmark(c: &mut Criterion) {
  let mut random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  random_bytes[0] |= 1;
  c.bench_function("mr2", move |b| b.iter(|| bench_mr2(&random_bytes)));
  c.bench_function("mr2_rug", move |b| b.iter(|| bench_mr2_rug(&random_bytes)));
  c.bench_function("lucas", move |b| b.iter(|| bench_lucas(&random_bytes)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
