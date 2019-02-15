/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::hash::primality::{passes_miller_rabin_base_2, passes_miller_rabin_base_22};
use accumulator::{i256, I256};
use criterion::Criterion;
use rand::Rng;
use rug::integer::Order;
use rug::Integer;

const NUM_JACOBI_AS: u64 = 100;

fn bench_jacobi_rug() {
  let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  let n = Integer::from_digits(&random_bytes, Order::MsfLe);
  let mut a = Integer::from(5);
  for _ in 0..NUM_JACOBI_AS {
    a.jacobi(&n);
    if a < 0 {
      a = -a + 2;
    } else {
      a = -a - 2;
    }
  }
}

fn bench_mr2_pablo(bytes: &[u8; 32]) {
  let n = Integer::from_digits(bytes, Order::LsfBe);
  passes_miller_rabin_base_2(&n);
}

fn bench_mr2_rug(bytes: &[u8; 32]) {
  let n = Integer::from_digits(bytes, Order::LsfBe);
  // GMP does not let us demand a base-2 Fermat test so we just do 1 of random base
  n.is_probably_prime(1);
}

fn bench_mr2_zero(bytes: &[u8; 32]) {
  // GMP does not let us demand a base-2 Fermat test so we just do 1 of random base
  passes_miller_rabin_base_22(&i256(bytes));
}

fn criterion_benchmark(c: &mut Criterion) {
  let mut random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  random_bytes[0] |= 1;
  c.bench_function("jacobi_rug", |b| b.iter(bench_jacobi_rug));
  c.bench_function("mr2_pablo", move |b| {
    b.iter(|| bench_mr2_pablo(&random_bytes))
  });
  c.bench_function("mr2_rug", move |b| b.iter(|| bench_mr2_rug(&random_bytes)));
  c.bench_function("mr2_zero", move |b| {
    b.iter(|| bench_mr2_zero(&random_bytes))
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
