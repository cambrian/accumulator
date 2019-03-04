/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

#[macro_use]
mod group_bench_util;

use accumulator::group::{ElemFrom, Group, Rsa2048, UnknownOrderGroup};
use criterion::Criterion;
use rug::Integer;
use std::str::FromStr;

#[derive(Clone)]
struct RsaBenchEnv {
  op_l: <Rsa2048 as Group>::Elem,
  op_r: <Rsa2048 as Group>::Elem,
  op_large_l: <Rsa2048 as Group>::Elem,
  op_large_r: <Rsa2048 as Group>::Elem,
  exp_base: <Rsa2048 as Group>::Elem,
  exp: Integer,
  elem_to_square: <Rsa2048 as Group>::Elem,
  elem_to_inv: <Rsa2048 as Group>::Elem,
}

// Initialize all the elements we need here so that initialization logic
// does not pollute the benchmarks.
fn init_env() -> RsaBenchEnv {
  let op_l = Rsa2048::elem(
    Integer::from_str(
      "111066521363124532171649626395987136074128970245601106158251038731392583290069",
    )
    .unwrap(),
  );
  let op_r = Rsa2048::elem(
    Integer::from_str(
      "106610920435831899020588753249099054915951032185883121197718271189872278955399",
    )
    .unwrap(),
  );

  let op_large_l = Rsa2048::elem(
    Integer::from_str(
      "2172073899553954285893691587818692186975191598984015216589930386158248724081087849265975174\
       9672737203717627738047648700009977053044057502917091973287111671693426065546612150833232954\
       3615367099810550371217642707848747209719337160655740326150736137284544974770721296865388733\
       3057277396369601863707823088589609031265453680152037285312247125429494632830592984498231941\
       6384204134056551840145916685870951507887895129356414704422748714217113880489703934147612551\
       9380825017530552968018297030172607314398711102156189885095451290884843968486448057303474665\
       81515692959313583208325725034506693916571047785061884094866050395109710",
    )
    .unwrap(),
  );

  let op_large_r = Rsa2048::elem(
    Integer::from_str(
      "3172073899553954285893691587818692186975191598984015216589930386158248724081087849265975174\
       9672737203717627738047648700009977053044057502917091973287111671693426065546612150833232954\
       3615367099810550371217642707848747209719337160655740326150736137284544974770721296865388733\
       3057277396369601863707823088589609031265453680152037285312247125429494632830592984498231941\
       6384204134056551840145916685870951507887895129356414704422748714217113880489703934147612551\
       9380825017530552968018297030172607314398711102156189885095451290884843968486448057303474665\
       81515692959313583208325725034506693916571047785061884094866050395109710",
    )
    .unwrap(),
  );

  let exp_base = Rsa2048::unknown_order_elem();
  let exp = Integer::from_str(
    "6531513683389606180955725446695124007119189061243576857500117325602044754680002922154438028\
     8474666886816442984548106882909827295319824031764930714696522619672276938781971873901815262\
     4216545626917306691611266738335435709225561930968971212874444236961226918266618788498569915\
     09472508677693535083051665283493383",
  )
  .unwrap();

  RsaBenchEnv {
    op_l: op_l,
    op_r: op_r,
    op_large_l: op_large_l,
    op_large_r: op_large_r,
    exp_base: exp_base,
    exp: exp,
    elem_to_square: Rsa2048::unknown_order_elem(),
    elem_to_inv: Rsa2048::unknown_order_elem(),
  }
}

fn criterion_benchmark(c: &mut Criterion) {
  let env = init_env();
  c.bench_function(
    "group_rsa_op",
    enclose!((env) move |b| b.iter(|| Rsa2048::op(&env.op_l, &env.op_r))),
  );
  c.bench_function(
    "group_rsa_op_large",
    enclose!((env) move |b| b.iter(|| Rsa2048::op(&env.op_large_l, &env.op_large_r))),
  );
  c.bench_function(
    "group_rsa_exp",
    enclose!((env) move |b| b.iter(|| Rsa2048::exp(&env.exp_base, &env.exp))),
  );
  c.bench_function(
    "group_rsa_inv",
    enclose!((env) move |b| b.iter(|| Rsa2048::inv(&env.elem_to_inv))),
  );
  c.bench_function(
    "group_rsa_square",
    enclose!((env) move |b| b.iter(|| Rsa2048::op(&env.elem_to_square, &env.elem_to_square))),
  );
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
