/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::group::{ElemFrom, RSA2048};
use crypto::proof::PoKCR;
use crypto::util::int;

fn bench_pokcr_rsa() {
  let witnesses = [RSA2048::elem(2), RSA2048::elem(3)];
  let x = [int(2), int(2)];
  let alphas = [RSA2048::elem(4), RSA2048::elem(9)];
  let proof = PoKCR::<RSA2048>::prove(&witnesses);
  PoKCR::verify(&alphas, &x, &proof);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("pokcr_rsa", |b| b.iter(bench_pokcr_rsa));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
