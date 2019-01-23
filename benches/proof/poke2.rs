/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::group::UnknownOrderGroup;
use crypto::group::{DummyRSA, RSA2048};
use crypto::proof::PoKE2;
use crypto::util::bi;

fn bench_poke2_dummy() {
  let base = DummyRSA::unknown_order_elem();
  let exp = bi(20);
  let result = DummyRSA::elem_of(1_048_576);
  let proof = PoKE2::<DummyRSA>::prove(&base, &exp, &result);
  PoKE2::verify(&base, &result, &proof);
}

fn bench_poke2_rsa() {
  let base = RSA2048::unknown_order_elem();
  let exp = bi(20);
  let result = RSA2048::elem_of(1_048_576);
  let proof = PoKE2::<RSA2048>::prove(&base, &exp, &result);
  PoKE2::verify(&base, &result, &proof);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("poke2_dummy", |b| b.iter(bench_poke2_dummy));
  c.bench_function("poke2_rsa", |b| b.iter(bench_poke2_rsa));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
