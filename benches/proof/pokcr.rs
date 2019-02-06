/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::group::{ElemFrom, Rsa2048};
use accumulator::proof::Pokcr;
use accumulator::util::int;
use criterion::Criterion;

fn bench_pokcr_rsa() {
  let witnesses = [Rsa2048::elem(2), Rsa2048::elem(3)];
  let x = [int(2), int(2)];
  let alphas = [Rsa2048::elem(4), Rsa2048::elem(9)];
  let proof = Pokcr::<Rsa2048>::prove(&witnesses);
  Pokcr::verify(&alphas, &x, &proof);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("pokcr_rsa", |b| b.iter(bench_pokcr_rsa));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
