/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::group::{ElemFrom, UnknownOrderGroup, RSA2048};
use crypto::proof::PoE;
use crypto::util::int;

fn bench_poe_rsa() {
  let base = RSA2048::unknown_order_elem();
  let exp = int(20);
  let result = RSA2048::elem(1_048_576);
  let proof = PoE::<RSA2048>::prove(&base, &exp, &result);
  PoE::<RSA2048>::verify(&base, &exp, &result, &proof);
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("poe_rsa", |b| b.iter(bench_poe_rsa));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
