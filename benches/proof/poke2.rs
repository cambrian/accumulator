/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::internal::group::Rsa2048;
use accumulator::internal::group::{ElemFrom, UnknownOrderGroup};
use accumulator::internal::proof::Poke2;
use accumulator::internal::util::int;
use criterion::Criterion;

fn bench_poke2_rsa() {
  let base = Rsa2048::unknown_order_elem();
  let exp = int(20);
  let result = Rsa2048::elem(1_048_576);
  let proof = Poke2::<Rsa2048>::prove(&base, &exp, &result);
  Poke2::verify(&base, &result, &proof);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("poke2_rsa", |b| b.iter(bench_poke2_rsa));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
