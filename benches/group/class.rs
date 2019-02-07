/// See https://bheisler.github.io/criterion.rs/book/getting_started.html to add more benchmarks.
#[macro_use]
extern crate criterion;

use criterion::Criterion;
use crypto::group::{ClassGroup, ElemFrom, Group};
use rug::Integer;
use std::str::FromStr;

fn bench_op_class() {
  ClassGroup::op(
    &ClassGroup::elem(
      (Integer::from_str("16").unwrap(),
      Integer::from_str("105").unwrap(),
      Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",
    ).unwrap(),
    )),
    &ClassGroup::elem(
      (Integer::from_str("16").unwrap(),
      Integer::from_str("9").unwrap(),
      Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036",
    ).unwrap(),
    )),
  );
}

fn bench_exp_class() {
  // generator obtained from hard-coded discriminant
  ClassGroup::exp(
    &ClassGroup::elem(
      (Integer::from(2),
      Integer::from(1),
      Integer::from_str("3827008629350940493386707189540101901936609547020633487839623582225300004666489306027281448853763777368990198117880164809708227406024703459009725115772610407878810521392085902015295545562523958711866779371531088132889638114041946661849770572154226710985917599916457066302682148335909706585071959150959814546206265435103373673496943574788744935795178127732520127531075979159538289365466373182137158779382092647246679657171935507126728878971929489212668908199079072163111583975633638661816714659180109107951783005735418950482497851235754121794548776139119565032459702128377126838952995785769100706778680652441494512278").unwrap(),
    )),
    &Integer::from_str(
      "6531513683389606180955725446695124007119189061243576857500117325602044754680002922154438028\
       8474666886816442984548106882909827295319824031764930714696522619672276938781971873901815262\
       4216545626917306691611266738335435709225561930968971212874444236961226918266618788498569915\
       09472508677693535083051665283493383",
    )
    .unwrap(),
  );
}

fn bench_square_class() {
  let mut g = ClassGroup::elem(
    (Integer::from_str("16").unwrap(),
    Integer::from_str("105").unwrap(),
    Integer::from_str(
    "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
    672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
    531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
    198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
    494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
    038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
    9057462766047140854869124473221137588347335081555186814207",
  ).unwrap(),
  ));
  g.bench_square();
}

fn bench_normalize_class() {
  let mut g = ClassGroup::elem(
    (Integer::from_str("16").unwrap(),
    Integer::from_str("105").unwrap(),
    Integer::from_str(
    "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
    672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
    531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
    198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
    494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
    038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
    9057462766047140854869124473221137588347335081555186814207",
  ).unwrap(),
  ));
  g.bench_normalize();
}

fn bench_reduce_class() {
  let mut g = ClassGroup::elem(
    (Integer::from_str("16").unwrap(),
    Integer::from_str("105").unwrap(),
    Integer::from_str(
    "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
    672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
    531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
    198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
    494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
    038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
    9057462766047140854869124473221137588347335081555186814207",
  ).unwrap(),
  ));
  g.bench_reduce();
}

fn bench_inv_class() {
  // generator obtained from hard-coded discriminant
  ClassGroup::inv(&ClassGroup::elem(
    (Integer::from(2),
    Integer::from(1),
    Integer::from_str("3827008629350940493386707189540101901936609547020633487839623582225300004666489306027281448853763777368990198117880164809708227406024703459009725115772610407878810521392085902015295545562523958711866779371531088132889638114041946661849770572154226710985917599916457066302682148335909706585071959150959814546206265435103373673496943574788744935795178127732520127531075979159538289365466373182137158779382092647246679657171935507126728878971929489212668908199079072163111583975633638661816714659180109107951783005735418950482497851235754121794548776139119565032459702128377126838952995785769100706778680652441494512278").unwrap(),
    )));
}

fn criterion_benchmark(c: &mut Criterion) {
  c.bench_function("group_class_op", |b| b.iter(bench_op_class));
  c.bench_function("group_class_exp", |b| b.iter(bench_exp_class));
  c.bench_function("group_class_inv", |b| b.iter(bench_inv_class));
  c.bench_function("group_class_square", |b| b.iter(bench_square_class));
  c.bench_function("group_class_normalize", |b| b.iter(bench_normalize_class));
  c.bench_function("group_class_reduce", |b| b.iter(bench_reduce_class));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
