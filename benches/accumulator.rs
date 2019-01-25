/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use criterion::Fun;
use crypto::accumulator::{add, delete, setup, verify_membership};
use crypto::group::{DummyRSA2048, UnknownOrderGroup, RSA2048};
use crypto::hash::{hash_to_prime, Blake2b};
use crypto::proof::PoE;
use num::bigint::BigUint;
use rand::Rng;

// fn bench_add_ring(elems: &[BigUint]) {
//   let acc = setup::<RSA2048>();
//   add::<RSA2048>(acc, elems);
// }

fn bench_delete<G: UnknownOrderGroup>(acc: G::Elem, witness: &[(BigUint, G::Elem)]) {
  delete::<G>(acc, witness);
}

fn bench_add_dummy(elems: &[BigUint]) {
  let acc = setup::<DummyRSA2048>();
  add::<DummyRSA2048>(acc, elems);
}

fn bench_verify<G: UnknownOrderGroup>(
  witness: &G::Elem,
  elems: &[BigUint],
  result: &G::Elem,
  proof: &PoE<G>,
) {
  assert!(verify_membership::<G>(witness, elems, result, proof));
}

fn criterion_benchmark(c: &mut Criterion) {
  let mut elems = Vec::new();
  for _ in 0..100 {
    let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
    let prime = hash_to_prime(&Blake2b::default, &random_bytes);
    elems.push(prime);
  }
  let elems_1 = [elems[0].clone()];
  let elems_2 = elems.clone();
  let elems_3 = elems.clone();
  let elems_4 = elems_1.clone();
  let elems_5 = elems.clone();
  let elems_6 = elems.clone();
  let acc_dummy = setup::<DummyRSA2048>();
  let (new_acc_dummy, poe_dummy) = add::<DummyRSA2048>(acc_dummy.clone(), &elems.clone()[1..]);
  let acc_ring = setup::<RSA2048>();
  let (new_acc_ring, poe_ring) = add::<RSA2048>(acc_ring.clone(), &elems.clone()[1..]);
  let (del_acc_dummy, _) = add::<DummyRSA2048>(new_acc_dummy, &[elems[0].clone()]);
  let (del_acc_ring, _) = add::<RSA2048>(new_acc_ring, &[elems[0].clone()]);

  // c.bench_function("add_ring_1", move |b| b.iter(|| bench_add_ring(&elems_1)));
  // c.bench_function("add_ring_10", move |b| {
  //   b.iter(|| bench_add_ring(&elems_2[0..10]))
  // });
  // c.bench_function("add_ring_100", move |b| b.iter(|| bench_add_ring(&elems_3)));
  // c.bench_function("add_dummy_1", move |b| b.iter(|| bench_add_dummy(&elems_4)));
  // c.bench_function("add_dummy_10", move |b| {
  //   b.iter(|| bench_add_dummy(&elems_5[0..10]))
  // });
  // c.bench_function("add_dummy_100", move |b| {
  //   b.iter(|| bench_add_dummy(&elems_6))
  // });
  // c.bench_function("verify_dummy", move |b| {
  //   b.iter(|| bench_verify(&acc_dummy, &elems, &new_acc_dummy, &poe_dummy))
  // });
  // let delete_fun = Fun::new("Delete_Ring", |b, (new_acc_ring, elem, del_acc_ring)| {
  //   b.iter(|| bench_delete::<RSA2048>(*new_acc_ring, &[(*elem, *del_acc_ring)]))
  // });
  // let del_funs = vec![delete_fun];
  // c.bench_functions(
  //   "delete",
  //   del_funs,
  //   (new_acc_ring.clone(), elems[0].clone(), del_acc_ring),
  // );
  // c.bench_function("verify_ring", move |b| {
  //   b.iter(|| bench_verify(&acc_ring, &elems_2, &new_acc_ring, &poe_ring))
  // });
  // c.bench_function("delete_ring", move |b| {
  //   b.iter(|| bench_delete::<RSA2048>(new_acc_ring, &[(elems_1[0], del_acc_ring)]))
  // });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
