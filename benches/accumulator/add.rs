/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::group::{ClassGroup, Rsa2048, UnknownOrderGroup};
use accumulator::hash::hash_to_prime;
use accumulator::{Accumulator, MembershipProof};
use criterion::Criterion;
use rand::Rng;
use rug::Integer;

fn bench_add<G: UnknownOrderGroup>(elems: &[Integer]) {
  let acc = Accumulator::<G>::new();
  acc.add(elems);
}

fn bench_verify<G: UnknownOrderGroup>(
  acc: &Accumulator<G>,
  elems: &[Integer],
  proof: &MembershipProof<G>,
) {
  assert!(acc.verify_membership(elems, proof));
}

#[allow(dead_code)]
fn bench_iterative_add(elems: &[Integer]) {
  let mut acc = Accumulator::<Rsa2048>::new();
  for elem in elems.chunks(1) {
    acc = acc.add(elem).0;
  }
}

fn init_acc<G: UnknownOrderGroup>() -> (Accumulator<G>, MembershipProof<G>, Vec<Integer>) {
  let mut elems = Vec::new();
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&random_bytes);
    elems.push(prime);
  }
  let acc = Accumulator::<G>::new();
  let (mut acc, mut proof) = acc.clone().add(&elems);
  for _ in 0..100 {
    elems = vec![];
    for _ in 0..100 {
      let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
      let prime = hash_to_prime(&random_bytes);
      elems.push(prime);
    }
    let (curr_acc, curr_proof) = acc.add(&elems);
    acc = curr_acc;
    proof = curr_proof;
  }
  (acc, proof, elems)
}

macro_rules! benchmark_add {
  ($group_type : ty, $criterion: ident) => {
    let group_type_str = String::from(stringify!($group_type)).to_lowercase();
    let (acc, proof, elems) = init_acc::<$group_type>();
    let elems_1 = elems.clone();
    let elems_2 = elems.clone();
    let elems_3 = elems.clone();

    $criterion.bench_function(format!("{}_add_1", group_type_str).as_str(), move |b| {
      b.iter(|| bench_add::<$group_type>(&elems_1[0..1]))
    });
    $criterion.bench_function(format!("{}_add_10", group_type_str).as_str(), move |b| {
      b.iter(|| bench_add::<$group_type>(&elems_2[0..10]))
    });
    $criterion.bench_function(format!("{}_add_100", group_type_str).as_str(), move |b| {
      b.iter(|| bench_add::<$group_type>(&elems_3))
    });
    $criterion.bench_function(format!("{}_verify", group_type_str).as_str(), move |b| {
      b.iter(|| bench_verify(&acc, &elems, &proof))
    });
  };
}

fn criterion_benchmark(c: &mut Criterion) {
  benchmark_add! {Rsa2048, c};
  benchmark_add! {ClassGroup, c};
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
