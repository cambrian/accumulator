/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::hash::primality::passes_miller_rabin_base_2;
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

fn bench_mr2_pablo() {
  let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  let n = Integer::from_digits(&random_bytes, Order::MsfBe);
  passes_miller_rabin_base_2(&n);
}

fn bench_mr2_rug() {
  let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
  let n = Integer::from_digits(&random_bytes, Order::MsfBe);
  // GMP does not let us demand a base-2 Fermat test so we just do 1 of random base
  n.is_probably_prime(1);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("jacobi_rug", |b| b.iter(bench_jacobi_rug));
  c.bench_function("mr2_pablo", |b| b.iter(bench_mr2_pablo));
  c.bench_function("mr2_rug", |b| b.iter(bench_mr2_rug));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
