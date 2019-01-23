/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::group::{DummyRSA, RSA2048};
use crypto::proof::PoKCR;
use crypto::util::bi;

fn bench_pokcr_dummy() {
  let witnesses = [&DummyRSA::elem_of(2), &DummyRSA::elem_of(3)];
  let x = [&bi(2), &bi(2)];
  let alphas = [&DummyRSA::elem_of(4), &DummyRSA::elem_of(9)];
  let proof = PoKCR::<DummyRSA>::prove(&witnesses);
  PoKCR::verify(&alphas, &x, &proof);
}

fn bench_pokcr_rsa() {
  let witnesses = [&RSA2048::elem_of(2), &RSA2048::elem_of(3)];
  let x = [&bi(2), &bi(2)];
  let alphas = [&RSA2048::elem_of(4), &RSA2048::elem_of(9)];
  let proof = PoKCR::<RSA2048>::prove(&witnesses);
  PoKCR::verify(&alphas, &x, &proof);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("pokcr_dummy", |b| b.iter(bench_pokcr_dummy));
  c.bench_function("pokcr_rsa", |b| b.iter(bench_pokcr_rsa));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
