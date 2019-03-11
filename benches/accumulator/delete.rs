/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::group::{ClassGroup, Rsa2048, UnknownOrderGroup};
use accumulator::{Accumulator, MembershipProof};
use criterion::Criterion;

fn bench_delete<G: UnknownOrderGroup>(
  acc: &Accumulator<G, &'static str>,
  c_proof: &MembershipProof<G, &'static str>,
) {
  acc
    .clone()
    .delete_with_proof(&[("c", c_proof.clone().witness)])
    .expect("Valid delete expected.");
}

macro_rules! benchmark_delete {
  ($group_type : ty, $criterion: ident) => {
    let group_type_str = String::from(stringify!($group_type)).to_lowercase();
    let acc_0 = Accumulator::<$group_type, &'static str>::empty().add(&["a", "b"]);
    let (acc_1, c_proof) = acc_0.clone().add_with_proof(&["c"]);
    $criterion.bench_function(format! {"{}_delete", group_type_str}.as_str(), move |b| {
      b.iter(|| bench_delete(&acc_1.clone(), &c_proof))
    });
  };
}

fn criterion_benchmark(c: &mut Criterion) {
  benchmark_delete! {Rsa2048, c};
  benchmark_delete! {ClassGroup, c};
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
