/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::internal::group::{ClassGroup, Rsa2048, UnknownOrderGroup};
use accumulator::internal::hash::hash_to_prime;
use accumulator::internal::util::int;
use accumulator::{Accumulator, MembershipProof};
use criterion::Criterion;
use rand::Rng;
use rug::Integer;

fn bench_delete<G: UnknownOrderGroup>(acc: Accumulator<G>, witness: &[(Integer, Accumulator<G>)]) {
  acc.delete(witness).expect("valid delete expected");
}

fn init_acc<G: UnknownOrderGroup>() -> (Accumulator<G>, MembershipProof<G>, Vec<Integer>) {
  let mut elems = Vec::new();
  for _ in 0..50 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&random_bytes);
    elems.push(prime);
  }
  let acc = Accumulator::<G>::new();
  let (acc, proof) = acc.clone().add(&elems);
  (acc, proof, elems)
}

macro_rules! benchmark_delete {
  ($group_type : ty, $criterion: ident) => {
    let group_type_str = String::from(stringify!($group_type)).to_lowercase();
    // Test verification on lots of elements. Added in batches to not go crazy with exponent size.
    let (acc, _, elems) = init_acc::<$group_type>();
    let mut delete_arg = Vec::new();
    for elem in elems.clone() {
      let exp: Integer = elems.iter().product::<Integer>() / elem.clone();
      let witness = Accumulator::<$group_type>::new();
      let witness = witness.exp_quotient(exp, int(1)).unwrap();
      delete_arg.push((elem, witness));
    }
    let acc_2 = acc.clone();
    $criterion.bench_function(
      format! {"{}_delete_50", group_type_str}.as_str(),
      move |b| b.iter(|| bench_delete(acc_2.clone(), &delete_arg[..])),
    );
  };
}

fn criterion_benchmark(c: &mut Criterion) {
  benchmark_delete! {Rsa2048, c};
  benchmark_delete! {ClassGroup, c};
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
