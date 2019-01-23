/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::group::UnknownOrderGroup;
use crypto::group::{DummyRSA, RSA2048};
use crypto::proof::PoE;
use crypto::util::bu;

fn bench_poe_dummy() {
  let base = DummyRSA::elem_of(2);
  let exp = bu(20u8);
  let result = DummyRSA::elem_of(1_048_576);
  let proof = PoE::<DummyRSA>::prove(&base, &exp, &result);
  PoE::<DummyRSA>::verify(&base, &exp, &result, &proof);
}

fn bench_poe_rsa() {
  let base = RSA2048::unknown_order_elem();
  let exp = bu(20u8);
  let result = RSA2048::elem_of(1_048_576);
  let proof = PoE::<RSA2048>::prove(&base, &exp, &result);
  PoE::<RSA2048>::verify(&base, &exp, &result, &proof);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("poe_dummy", |b| b.iter(bench_poe_dummy));
  c.bench_function("poe_rsa", |b| b.iter(bench_poe_rsa));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
