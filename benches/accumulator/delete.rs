/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use accumulator::group::{Rsa2048, UnknownOrderGroup};
use accumulator::hash::hash_to_prime;
use accumulator::util::int;
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

fn criterion_benchmark(c: &mut Criterion) {
  // Test verification on lots of elements. Added in batches to not go crazy with exponent size.
  let (acc, _, elems) = init_acc::<Rsa2048>();

  let mut delete_arg = Vec::new();
  for elem in elems.clone() {
    let exp: Integer = elems.iter().product::<Integer>() / elem.clone();
    let witness = Accumulator::<Rsa2048>::new();
    let witness = witness.exp_quotient(exp, int(1)).unwrap();
    delete_arg.push((elem, witness));
  }
  let acc_2 = acc.clone();
  c.bench_function("delete_50", move |b| {
    b.iter(|| bench_delete(acc_2.clone(), &delete_arg[..]))
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
