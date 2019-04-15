/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

#[macro_use]
mod group_bench_util;

use accumulator::group::{ClassGroup, ElemFrom, Group, UnknownOrderGroup};
use accumulator::num::mpz::Mpz;
use criterion::Criterion;
use rug::Integer;
use std::str::FromStr;

#[derive(Clone)]
struct ClassBenchEnv {
  op_l: <ClassGroup as Group>::Elem,
  op_r: <ClassGroup as Group>::Elem,
  exp_base: <ClassGroup as Group>::Elem,
  exp: Integer,
  elem_to_square: <ClassGroup as Group>::Elem,
  elem_to_inv: <ClassGroup as Group>::Elem,
  elem_to_reduce: (Mpz, Mpz, Mpz),
  elem_to_normalize: (Mpz, Mpz, Mpz),
}

// Initialize all the elements we need here so that initialization logic
// does not pollute the benchmarks.
fn init_env() -> ClassBenchEnv {
  let left = ClassGroup::elem((
    Mpz::from_str("16").unwrap(),
    Mpz::from_str("9").unwrap(),
    Mpz::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036",
    )
    .unwrap(),
  ));
  let right = left.clone();

  let base = ClassGroup::unknown_order_elem();
  let exp = Integer::from_str(
    "6531513683389606180955725446695124007119189061243576857500117325602044754680002922154438028",
  )
  .unwrap();

  let aa = Mpz::from_str("16").unwrap();
  let bb = Mpz::from_str("105").unwrap();
  let cc = Mpz::from_str(
    "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
     672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
     531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
     198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
     494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
     038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
     9057462766047140854869124473221137588347335081555186814207",
  )
  .unwrap();

  // Element which requires one iteration to reduce, represented as a tuple
  // here, since only reduced representations of ClassElem are allowed.
  let g_red = (cc.clone(), bb.clone(), aa.clone());
  let g_norm = (aa, bb, cc);

  ClassBenchEnv {
    op_l: left,
    op_r: right,
    exp_base: base.clone(),
    exp: exp,
    elem_to_inv: base.clone(),
    elem_to_square: ClassGroup::unknown_order_elem(),
    elem_to_reduce: g_red,
    elem_to_normalize: g_norm,
  }
}

fn criterion_benchmark(c: &mut Criterion) {
  let env = init_env();
  c.bench_function(
    "group_class_op",
    enclose!(
      (env) move |b| {
        b.iter(|| ClassGroup::op(&env.op_l, &env.op_r))
      }
    ),
  );

  c.bench_function(
    "group_class_exp",
    enclose!(
      (env) move |b| {
        b.iter(|| ClassGroup::exp(&env.exp_base, &env.exp))
      }
    ),
  );

  c.bench_function(
    "group_class_inv",
    enclose!(
      (env) move |b| {
        b.iter(|| ClassGroup::inv(&env.elem_to_inv))
      }
    ),
  );

  c.bench_function(
    "group_class_normalize",
    enclose!(
      (env) move |b| {
        b.iter_with_setup(
          || env.elem_to_normalize.clone(),
          |g| ClassGroup::normalize(g.0, g.1, g.2)
        )
      }
    ),
  );

  c.bench_function(
    "group_class_reduce",
    enclose!(
      (env) move |b| {
        b.iter_with_setup(
          || env.elem_to_reduce.clone(),
          |g| ClassGroup::reduce(g.0, g.1, g.2)
        )
      }
    ),
  );

  c.bench_function(
    "group_class_square",
    enclose!(
      (env) move |b| {
        b.iter_with_setup(
          || env.elem_to_square.clone(),
          |mut g| ClassGroup::square(&mut g)
        )
      }
    ),
  );
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
