//! Class Group implementation
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::util;
use crate::util::{int, TypeRep};
use gmp_mpfr_sys::gmp::{
  mpz_add, mpz_cmp, mpz_fdiv_q, mpz_fdiv_q_ui, mpz_fdiv_qr, mpz_gcd, mpz_gcdext, mpz_init, mpz_mod,
  mpz_mul, mpz_neg, mpz_set, mpz_set_str, mpz_set_ui, mpz_sub, mpz_t,
};
use rug::{Assign, Integer};
use std::cell::RefCell;
use std::ffi::CString;
use std::hash::{Hash, Hasher};
use std::mem::uninitialized;
use std::str::FromStr;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ClassGroup {}

// REVIEW: It appears that this group implementation always expects its elements to be reduced.
// However, the type signature fn reduce(&mut self) implies that there is a valid representation
// for un-reduced elements. This is suboptimal. Can you restructure the code in such a way that
// ClassElem instances are only created for valid elements? Likewise for normalize.
// Hint: change reduce to look something like
// fn reduce(a: Integer, b: Integer, c: Integer) -> ClassElem { ... }

// 2048-bit prime, negated, congruent to 3 mod 4.  Generated using OpenSSL.

// According to "A Survey of IQ Cryptography" (Buchmann & Hamdy) Table 1,
// IQ-MPQS for computing discrete logarithms in class groups with a 2048-bit discriminant is
// comparable in complexity to GNFS for factoring a 4096-bit integer.
const DISCRIMINANT2048_DECIMAL: &str =
  "-30616069034807523947093657516320815215492876376165067902716988657802400037331914448218251590830\
  1102189519215849430413184776658192481976276720778009261808832630304841711366872161223643645001916\
  6969493423497224870506311710491233557329479816457723381368788734079933165653042145718668727765268\
  0575673207678516369650123480826989387975548598309959486361425021860161020248607833276306314923730\
  9854570972702350567411779734372573754840570138310317754359137013512655926325773048926718050691092\
  9453371727344087286361426404588335160385998280988603297435639020911295652025967761702701701471162\
  3966286152805654229445219531956098223";

//lazy_static! {
static CLASS_GROUP_DISCRIMINANT: mpz_t = {
  let d = new_mpz();
  let d_str = CString::new(DISCRIMINANT2048_DECIMAL).unwrap();
  unsafe {
    mpz_set_str(&mut d, d_str.as_ptr(), 10);
  }
  d
};
//}
static CTX: RefCell<Ctx> = Default::default();

//#[derive(Clone, Debug, Eq)]
pub struct ClassElem {
  a: mpz_t,
  b: mpz_t,
  c: mpz_t,
}

impl Default for ClassElem {
  fn default() -> Self {
    ClassElem {
      a: new_mpz(),
      b: new_mpz(),
      c: new_mpz(),
    }
  }
}

fn new_mpz() -> mpz_t {
  unsafe {
    let ret = uninitialized();
    mpz_init(&mut ret);
    ret
  }
}

pub struct Ctx {
  negative_a: mpz_t,
  r: mpz_t,
  denom: mpz_t,
  old_b: mpz_t,
  ra: mpz_t,
  s: mpz_t,
  x: mpz_t,
  old_a: mpz_t,
  g: mpz_t,
  d: mpz_t,
  e: mpz_t,
  q: mpz_t,
  w: mpz_t,
  u: mpz_t,
  a: mpz_t,
  b: mpz_t,
  m: mpz_t,
  k: mpz_t,
  mu: mpz_t,
  v: mpz_t,
  sigma: mpz_t,
  lambda: mpz_t,
  h: mpz_t,
  t: mpz_t,
  l: mpz_t,
  j: mpz_t,
}

impl Default for Ctx {
  fn default() -> Self {
    Ctx {
      negative_a: new_mpz(),
      r: new_mpz(),
      denom: new_mpz(),
      old_b: new_mpz(),
      ra: new_mpz(),
      s: new_mpz(),
      x: new_mpz(),
      old_a: new_mpz(),
      g: new_mpz(),
      d: new_mpz(),
      e: new_mpz(),
      q: new_mpz(),
      w: new_mpz(),
      u: new_mpz(),
      a: new_mpz(),
      b: new_mpz(),
      m: new_mpz(),
      k: new_mpz(),
      mu: new_mpz(),
      v: new_mpz(),
      sigma: new_mpz(),
      lambda: new_mpz(),
      h: new_mpz(),
      t: new_mpz(),
      l: new_mpz(),
      j: new_mpz(),
    }
  }
}

// TODO: Use seperate context object for linear congruences and move into util.rs
// TODO: Check for solution
pub fn solve_linear_congruence_mpz(
  ctx: RefCell<Ctx>,
  mu: &mut mpz_t,
  v: &mut mpz_t,
  a: &mpz_t,
  b: &mpz_t,
  m: &mpz_t,
) {
  let ctx = ctx.borrow_mut();
  // g = gcd(a, m) => da + em = g
  mpz_gcdext(&mut ctx.g, &mut ctx.d, &mut ctx.e, a, m);
  // q = floor_div(b, g)
  // r = b % g
  mpz_fdiv_qr(&mut ctx.q, &mut ctx.r, b, &ctx.g);
  // mu = (q * d) % m
  mpz_mul(mu, &ctx.q, &ctx.d);
  mpz_mod(mu, mu, m);
  // v = m / g
  mpz_fdiv_q(v, m, &ctx.g);
}

// ClassElem and ClassGroup ops based on Chia's fantastic doc explaining applied class groups:
//  https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf.

impl ClassElem {
  fn normalize(&mut self) {
    if self.is_normalized() {
      return;
    }
    // r = floor_div((a - b), 2a)
    // (a, b, c) = (a, b + 2ra, ar^2 + br + c)
    let (r, _) = Integer::from(&self.a - &self.b).div_rem_floor(Integer::from(2 * &self.a));
    let new_b = &self.b + 2 * Integer::from(&r * &self.a);
    let new_c = &self.c + Integer::from(&self.b * &r) + &self.a * r.square();
    self.b = new_b;
    self.c = new_c;
  }

  fn reduce(&mut self) {
    self.normalize();
    while !self.is_reduced() {
      // s = floor_div(c + b, 2c)
      let (s, _) = Integer::from(&self.c + &self.b).div_rem_floor(Integer::from(2 * &self.c));

      // (a, b, c) = (c, −b + 2sc, cs^2 − bs + a)
      let new_a = self.c.clone();
      let new_b = Integer::from(-&self.b) + 2 * Integer::from(&s * &new_a);
      let new_c = -Integer::from(&self.b * &s) + &self.a + &self.c * s.square();
      self.a = new_a;
      self.b = new_b;
      self.c = new_c;
    }
    self.normalize();
  }

  #[allow(non_snake_case)]
  fn square(&mut self) {
    // Solve `bk = c mod a` for k, represented by mu, v and any integer n s.t. k = mu + v * n
    //
    let (mu, _) = util::solve_linear_congruence(&self.b, &self.c, &self.a).unwrap();

    // A = a^2
    // B = b - 2a * mu
    // tmp = (b * mu) / a
    // C = mu^2 - tmp
    let A = Integer::from(self.a.square_ref());
    let B = &self.b - Integer::from(2 * &self.a) * &mu;
    let (tmp, _) = <(Integer, Integer)>::from(
      Integer::from((&self.b * &mu) - &self.c).div_rem_floor_ref(&self.a),
    );
    let C = mu.square() - tmp;
    self.a = A;
    self.b = B;
    self.c = C;
    self.reduce();
  }

  fn discriminant(&self) -> Integer {
    Integer::from(self.b.square_ref()) - Integer::from(4) * &self.a * &self.c
  }

  fn validate(&self) -> bool {
    &self.discriminant() == ClassGroup::rep()
  }

  fn is_reduced(&self) -> bool {
    self.is_normalized() && (self.a <= self.c && !(self.a == self.c && self.b < 0))
  }

  fn is_normalized(&self) -> bool {
    -Integer::from(&self.a) < self.b && self.b <= self.a
  }

  #[inline]
  #[cfg(feature = "benchmark")]
  pub fn bench_square(&mut self) {
    self.square()
  }

  #[inline]
  #[cfg(feature = "benchmark")]
  pub fn bench_normalize(&mut self) {
    self.normalize()
  }

  #[inline]
  #[cfg(feature = "benchmark")]
  pub fn bench_reduce(&mut self) {
    self.reduce()
  }

  #[cfg(feature = "benchmark")]
  pub fn assign(&mut self, a: &Integer, b: &Integer, c: &Integer) {
    self.a = a.clone();
    self.b = b.clone();
    self.c = c.clone();
  }
}

pub struct ClassRep {
  discriminant: mpz_t,
  ctx: RefCell<Ctx>,
}
impl TypeRep for ClassGroup {
  type Rep = ClassRep;
  fn rep() -> &'static Self::Rep {
    &ClassRep {
      discriminant: CLASS_GROUP_DISCRIMINANT,
      ctx: CTX,
    }
  }
}

impl Group for ClassGroup {
  type Elem = ClassElem;

  #[allow(non_snake_case)]
  fn op_(rep: &ClassRep, x: &ClassElem, y: &ClassElem) -> ClassElem {
    // TODO: Assert discriminants equal
    let ctx = rep.ctx.borrow_mut();

    // g = (b1 + b2) / 2
    mpz_add(&mut ctx.g, &x.b, &y.b);
    mpz_fdiv_q_ui(&mut ctx.g, &ctx.g, 2);
    // h = (b2 - b1) / 2
    mpz_sub(&mut ctx.h, &y.b, &x.b);
    // w = gcd(a1, a2, g)
    mpz_gcd(&mut ctx.w, &x.a, &y.a);
    mpz_gcd(&mut ctx.w, &ctx.w, &ctx.g);

    // j = w
    mpz_set(&mut ctx.j, &ctx.w);
    // r = 0
    mpz_set_ui(&mut ctx.r, 0);
    // s = a1 / w
    mpz_fdiv_q(&mut ctx.s, &x.a, &ctx.w);
    // t = a2 / w
    mpz_fdiv_q(&mut ctx.t, &y.a, &ctx.w);
    // u = g / w
    mpz_fdiv_q(&mut ctx.u, &ctx.g, &ctx.w);

    // a = tu
    mpz_mul(&mut ctx.a, &ctx.t, &ctx.u);

    // b = hu + sc
    mpz_mul(&mut ctx.b, &ctx.h, &ctx.u);
    mpz_mul(&mut ctx.m, &ctx.s, &x.c);
    mpz_add(&mut ctx.b, &ctx.b, &ctx.m);
    // m = st
    mpz_mul(&mut ctx.m, &ctx.s, &ctx.t);
    // Solve linear congruence `(tu)k = hu + sc mod st` or `ak = b mod m` for solutions k.
    solve_linear_congruence_mpz(rep.ctx, &mut ctx.mu, &mut ctx.v, &ctx.a, &ctx.b, &ctx.m);

    // a = tv
    mpz_mul(&mut ctx.m, &ctx.t, &ctx.v);
    // b = h - t * mu
    mpz_mul(&mut ctx.m, &ctx.t, &ctx.mu);
    mpz_sub(&mut ctx.b, &ctx.h, &ctx.m);
    // m = s
    mpz_set(&mut ctx.m, &ctx.s);
    // Solve linear congruence `(tv)k = h - t * mu mod s` or `ak = b mod m` for solutions k
    solve_linear_congruence_mpz(
      rep.ctx,
      &mut ctx.lambda,
      &mut ctx.sigma,
      &ctx.a,
      &ctx.b,
      &ctx.m,
    );

    // k = mu + v * lambda
    mpz_mul(&mut ctx.a, &ctx.v, &ctx.lambda);
    mpz_add(&mut ctx.k, &ctx.mu, &ctx.a);
    // l = (k * t - h) / s
    mpz_mul(&mut ctx.l, &ctx.k, &ctx.t);
    mpz_sub(&mut ctx.l, &ctx.l, &ctx.h);
    mpz_fdiv_q(&mut ctx.l, &ctx.l, &ctx.s);
    // m = (tuk - hu - cs) / st
    mpz_mul(&mut ctx.m, &ctx.t, &ctx.u);
    mpz_mul(&mut ctx.m, &ctx.m, &ctx.k);
    mpz_mul(&mut ctx.a, &ctx.h, &ctx.u);
    mpz_sub(&mut ctx.m, &ctx.m, &ctx.a);
    mpz_mul(&mut ctx.a, &x.c, &ctx.s);
    mpz_sub(&mut ctx.m, &ctx.m, &ctx.a);
    mpz_mul(&mut ctx.a, &ctx.s, &ctx.t);
    mpz_fdiv_q(&mut ctx.m, &ctx.m, &ctx.a);

    let ret = ClassElem::default();
    // A = st - ru
    mpz_mul(&mut ret.a, &ctx.s, &ctx.t);
    mpz_mul(&mut ctx.a, &ctx.r, &ctx.u);
    mpz_sub(&mut ret.a, &ret.a, &ctx.a);

    // B = ju - kt + ls
    mpz_mul(&mut ret.b, &ctx.j, &ctx.u);
    mpz_mul(&mut ret.a, &ctx.m, &ctx.r);
    mpz_add(&mut ret.b, &ret.b, &ctx.a);
    mpz_mul(&mut ctx.a, &ctx.k, &ctx.t);
    mpz_sub(&mut ret.b, &ret.b, &ctx.a);
    mpz_mul(&mut ctx.a, &ctx.l, &ctx.s);
    mpz_sub(&mut ret.b, &ret.b, &ctx.a);

    // C = kl - jm
    mpz_mul(&mut ret.c, &ctx.k, &ctx.l);
    mpz_mul(&mut ctx.a, &ctx.j, &ctx.m);
    mpz_sub(&mut ret.c, &ret.c, &ctx.a);

    ret.reduce();
    ret
  }

  fn id_(rep: &ClassRep) -> ClassElem {
    let ctx = rep.ctx.borrow_mut();
    let d = rep.discriminant;
    let ret = ClassElem::default();
    mpz_set_ui(&mut ret.a, 1);
    mpz_set_ui(&mut ret.b, 1);
    // c = (b * b - d) / 4a
    mpz_sub(&mut ctx.a, &ret.b, &d); // b == b*b
    mpz_fdiv_q_ui(&mut ret.c, &ctx.a, 4);
    ret
  }

  fn inv_(_: &ClassRep, x: &ClassElem) -> ClassElem {
    let ret = ClassElem::default();
    mpz_set(&mut ret.a, &x.a);
    mpz_neg(&mut ret.b, &x.b);
    mpz_set(&mut ret.c, &x.c);
    ret
  }

  fn exp_(_: &Integer, a: &ClassElem, n: &Integer) -> ClassElem {
    let (mut val, mut a, mut n) = {
      if *n < int(0) {
        (Self::id(), Self::inv(a), int(-n))
      } else {
        (Self::id(), a.clone(), n.clone())
      }
    };
    loop {
      if n == int(0) {
        return val;
      }
      if n.is_odd() {
        val = Self::op(&val, &a);
      }
      a.square();
      n >>= 1;
    }
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Integer) -> ClassElem {
    // a = 2
    // b = 1
    // c = (b * b - d) / 4a
    let a = Integer::from(2);
    let b = Integer::from(1);
    let c = Integer::from(1 - d) / Integer::from(8);
    let mut ret = ClassElem { a, b, c };
    ret.reduce();
    ret
  }
}

impl Hash for ClassElem {
  fn hash<H: Hasher>(&self, state: &mut H) {
    let mut self_reduced = self.clone();
    self_reduced.reduce();
    self_reduced.a.hash(state);
    self_reduced.b.hash(state);
    self_reduced.c.hash(state);
  }
}

impl PartialEq for ClassElem {
  fn eq(&self, other: &ClassElem) -> bool {
    let mut self_reduced = self.clone();
    self_reduced.reduce();
    let mut other_reduced = other.clone();
    other_reduced.reduce();
    self_reduced.a == other_reduced.a
      && self_reduced.b == other_reduced.b
      && self_reduced.c == other_reduced.c
  }
}

impl<T> ElemFrom<(T, T, T)> for ClassGroup
where
  Integer: From<T>,
{
  fn elem(t: (T, T, T)) -> ClassElem {
    let mut class_elem = ClassElem {
      a: Integer::from(t.0),
      b: Integer::from(t.1),
      c: Integer::from(t.2),
    };

    // Ideally, this should return an error and the
    // return type of ElemFrom should be Result<Self::Elem, Self:err>,
    // but this would require a lot of ugly "unwraps" in the accumulator
    // library. Besides, users should not need to create new class group
    // elements, so an invalid ElemFrom here should signal a severe internal error.
    assert!(class_elem.validate());

    class_elem.reduce();
    class_elem
  }
}

// Caveat: tests that use "ground truth" use outputs from
//  Chia's sample implementation in python:
//    https://github.com/Chia-Network/vdf-competition/blob/master/inkfish/classgroup.py.
#[cfg(test)]
mod tests {
  use super::*;
  use std::collections::hash_map::DefaultHasher;

  // Makes a class elem tuple but does not reduce.
  fn construct_raw_elem_from_strings(a: &str, b: &str, c: &str) -> ClassElem {
    ClassElem {
      a: Integer::from_str(a).unwrap(),
      b: Integer::from_str(b).unwrap(),
      c: Integer::from_str(c).unwrap(),
    }
  }

  #[should_panic]
  #[test]
  fn test_bad_elem() {
    let _ = ClassGroup::elem((1, 2, 3));
  }

  #[test]
  fn test_elem_from() {
    let a1 = Integer::from_str("16").unwrap();
    let b1 = Integer::from_str("105").unwrap();
    let c1 = Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",
    )
    .unwrap();

    let a2 = Integer::from_str("16").unwrap();
    let b2 = Integer::from_str("9").unwrap();
    let c2 = Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036",
    )
    .unwrap();

    let reduced_elem = ClassGroup::elem((a1, b1, c1));
    let also_reduced_elem = ClassGroup::elem((a2, b2, c2));
    assert_eq!(reduced_elem, also_reduced_elem);
  }

  #[test]
  fn test_equality() {
    let not_reduced = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    let diff_elem = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    assert!(not_reduced == reduced_ground_truth);
    assert!(not_reduced == not_reduced.clone());
    assert!(reduced_ground_truth == reduced_ground_truth.clone());
    assert!(not_reduced != diff_elem);
    assert!(reduced_ground_truth != diff_elem);
  }

  #[test]
  fn test_hash() {
    let not_reduced = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    let diff_elem = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    let mut hasher_lh = DefaultHasher::new();
    let mut hasher_rh = DefaultHasher::new();
    not_reduced.hash(&mut hasher_lh);
    reduced_ground_truth.hash(&mut hasher_rh);
    assert!(hasher_lh.finish() == hasher_rh.finish());
    assert!(hasher_lh.finish() == hasher_lh.finish());
    assert!(hasher_rh.finish() == hasher_rh.finish());

    hasher_lh = DefaultHasher::new();
    hasher_rh = DefaultHasher::new();
    not_reduced.hash(&mut hasher_lh);
    diff_elem.hash(&mut hasher_rh);
    assert!(hasher_lh.finish() != hasher_rh.finish());
  }

  #[test]
  fn test_reduce_basic() {
    let mut to_reduce = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    to_reduce.reduce();
    assert_eq!(to_reduce, reduced_ground_truth);

    let mut already_reduced = reduced_ground_truth.clone();
    already_reduced.reduce();
    assert_eq!(already_reduced, reduced_ground_truth);
  }

  #[test]
  fn test_normalize_basic() {
    let mut unnormalized = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",

    );

    let normalized_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "4783760786688675616733383986925127377420761933775791859799529477781625005833111632534101811\
       06720472171123774764735020601213528425753087932376215639471576300984851315174010737751911943\
       19531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526\
       85419887133231339948938699768182757831793879217091871179468485931169743972659665650159413844\
       97394942286170683296647767144847422761580905834957146491938390841109871491186151583613524884\
       88402038894799695420483272708933239751363849397287571692736881031223140446926522431859701738\
       9945629057462766047140854869124473221137588347335081555186814036",
    );

    unnormalized.normalize();
    assert_eq!(normalized_ground_truth, unnormalized);
  }

  #[test]
  fn test_discriminant_basic() {
    let g = ClassGroup::unknown_order_elem();
    assert_eq!(&g.discriminant(), ClassGroup::rep());
  }

  #[test]
  fn test_discriminant_across_ops() {
    let id = ClassGroup::id();
    let g1 = ClassGroup::unknown_order_elem();
    let g2 = ClassGroup::op(&g1, &g1);
    let g3 = ClassGroup::op(&id, &g2);
    let g3_inv = ClassGroup::inv(&g3);

    assert!(id.validate());
    assert!(g1.validate());
    assert!(g2.validate());
    assert!(g3.validate());
    assert!(g3_inv.validate());
  }

  #[test]
  fn test_op_single() {
    let a = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    let b = construct_raw_elem_from_strings(
      "16",
      "41",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814061"
    );

    let ground_truth = construct_raw_elem_from_strings(
      "64",
      "9",
      "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      50972369992385512081817723330993784096234932189292318422025780578511173163060796492543474864\
      07264365691511785213717281118305284397086833770388796703509"
    );

    assert_eq!(ClassGroup::op(&a, &b), ground_truth);
  }

  #[test]
  fn test_op_alternating() {
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();
    let mut g_star = ClassGroup::id();

    // g
    g = ClassGroup::op(&g_anchor, &g);

    // g^2, g^* = g^2
    g = ClassGroup::op(&g_anchor, &g);
    g_star = ClassGroup::op(&g, &g_star);

    // g^3
    g = ClassGroup::op(&g_anchor, &g);

    // g^4, g^* = g^2 * g^4 = g^6
    g = ClassGroup::op(&g_anchor, &g);
    g_star = ClassGroup::op(&g, &g_star);

    let ground_truth = construct_raw_elem_from_strings(
      "64",
      "9",
      "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      509723699923855120818177233309937840962349321892923184220257805785111731630607964925434748640\
      7264365691511785213717281118305284397086833770388796703509"
    );

    assert_eq!(ground_truth, g_star);
  }

  #[test]
  fn test_op_complex() {
    // 1. Take g^100, g^200, ..., g^1000.
    // 2. Compute g^* = g^100 * ... * g^1000
    // 3. For each of g^100, g^200, ..., g^1000
    //    compute the inverse of that element and
    //    assert that g^* * current_inverse = product of g^100, g^200, ..., g^1000
    //    without the inversed-out element.
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    let mut gs = vec![];
    let mut gs_invs = vec![];

    let mut g_star = ClassGroup::id();
    for i in 1..=1000 {
      g = ClassGroup::op(&g_anchor, &g);
      assert!(g.validate());
      if i % 100 == 0 {
        gs.push(g.clone());
        gs_invs.push(ClassGroup::inv(&g));
        g_star = ClassGroup::op(&g, &g_star);
        assert!(g_star.validate());
      }
    }

    let elems_n_invs = gs.iter().zip(gs_invs.iter());
    for (g_elem, g_inv) in elems_n_invs {
      assert!(g_elem.validate());
      assert!(g_inv.validate());
      let mut curr_prod = ClassGroup::id();
      for elem in &gs {
        if elem != g_elem {
          curr_prod = ClassGroup::op(&curr_prod, &elem);
          assert!(curr_prod.validate());
        }
      }
      assert_eq!(ClassGroup::id(), ClassGroup::op(&g_inv, &g_elem));
      assert_eq!(curr_prod, ClassGroup::op(&g_inv, &g_star));
    }
  }

  #[test]
  fn test_id_basic() {
    let g = ClassGroup::unknown_order_elem();
    let id = ClassGroup::id();
    assert_eq!(g, ClassGroup::op(&g, &id));
    assert_eq!(g, ClassGroup::op(&id, &g));
    assert_eq!(id, ClassGroup::op(&id, &id));
  }

  #[test]
  fn test_id_repeated() {
    let mut id = ClassGroup::id();
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();
    for _ in 0..1000 {
      id = ClassGroup::op(&id, &id);
      assert_eq!(id, ClassGroup::id());
      g = ClassGroup::op(&g, &ClassGroup::id());
      assert_eq!(g, g_anchor);
    }
  }

  #[test]
  fn test_inv() {
    let id = ClassGroup::id();
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();

    for _ in 0..1000 {
      g = ClassGroup::op(&g, &g_anchor);
      let g_inv = ClassGroup::inv(&g);
      assert_eq!(id, ClassGroup::op(&g_inv, &g));
      assert_eq!(id, ClassGroup::op(&g, &g_inv));
      assert_eq!(g, ClassGroup::inv(&g_inv));
    }
  }

  #[test]
  fn test_exp_basic() {
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    for i in 1..=1000 {
      g = ClassGroup::op(&g, &g_anchor);
      assert_eq!(&g, &ClassGroup::exp(&g_anchor, &Integer::from(i)));
    }
  }

  #[test]
  fn test_square_basic() {
    let g = ClassGroup::unknown_order_elem();
    let mut g4 = ClassGroup::id();

    // g^4
    for _ in 0..4 {
      g4 = ClassGroup::op(&g, &g4);
    }

    // g^2
    let mut g2 = ClassGroup::op(&g, &g);

    // g^4
    g2.square();

    assert_eq!(&g2, &g4);
  }
}
