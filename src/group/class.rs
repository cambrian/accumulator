//! Class Group implementation
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::util;
use crate::util::{int, TypeRep};
use gmp_mpfr_sys::gmp::{
  mpz_add, mpz_cmp, mpz_cmp_si, mpz_cmp_ui, mpz_fdiv_q, mpz_fdiv_q_ui, mpz_fdiv_qr, mpz_gcd,
  mpz_gcdext, mpz_get_str, mpz_init, mpz_mod, mpz_mul, mpz_mul_ui, mpz_neg, mpz_set, mpz_set_str,
  mpz_set_ui, mpz_sub, mpz_t,
};
use rug::Integer;
use std::cell::RefCell;
use std::ffi::CStr;
use std::ffi::CString;
use std::hash::{Hash, Hasher};
use std::mem::uninitialized;
use std::slice;
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

// need wrapper to make mpz_t thread-safe
#[cfg_attr(repr_transparent, repr(transparent))]
pub struct Discriminant {
  inner: mpz_t,
}

unsafe impl Send for Discriminant {}
unsafe impl Sync for Discriminant {}

impl Default for Discriminant {
  fn default() -> Self {
    Discriminant { inner: new_mpz() }
  }
}

lazy_static! {
  pub static ref CLASS_GROUP_DISCRIMINANT: Discriminant = {
    let mut d = Discriminant::default();
    let d_str = CString::new(DISCRIMINANT2048_DECIMAL).unwrap();
    unsafe {
      mpz_set_str(&mut d.inner, d_str.as_ptr(), 10);
    }
    d
  };
}

thread_local! {
  static CTX: RefCell<Ctx> = Default::default();
}

pub fn with_context<T, U>(f: T) -> U
where
  T: FnOnce(&mut Ctx) -> U,
{
  let mut opt = None;
  CTX.with(|x| opt = Some(f(&mut x.borrow_mut())));
  opt.unwrap()
}

#[derive(Debug)]
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

impl Clone for ClassElem {
  fn clone(&self) -> Self {
    let mut ret = ClassElem::default();
    unsafe {
      mpz_set(&mut ret.a, &self.a);
      mpz_set(&mut ret.b, &self.b);
      mpz_set(&mut ret.c, &self.c);
    };
    ret
  }
}

impl PartialEq for ClassElem {
  fn eq(&self, other: &ClassElem) -> bool {
    let mut r_self = self.clone();
    let mut r_other = other.clone();
    r_self.reduce();
    r_other.reduce();
    unsafe {
      mpz_cmp(&r_self.a, &r_other.a) == 0
        && mpz_cmp(&r_self.b, &r_other.b) == 0
        && mpz_cmp(&r_self.c, &r_other.c) == 0
    }
  }
}

impl Eq for ClassElem {}

unsafe impl Send for ClassElem {}
unsafe impl Sync for ClassElem {}

fn new_mpz() -> mpz_t {
  unsafe {
    let mut ret = uninitialized();
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
pub fn solve_linear_congruence_mpz(ctx: &mut Ctx) {
  unsafe {
    // g = gcd(a, m) => da + em = g
    mpz_gcdext(&mut ctx.g, &mut ctx.d, &mut ctx.e, &ctx.a, &ctx.m);
    // q = floor_div(b, g)
    // r = b % g
    mpz_fdiv_qr(&mut ctx.q, &mut ctx.r, &ctx.b, &ctx.g);
    // mu = (q * d) % m
    mpz_mul(&mut ctx.mu, &ctx.q, &ctx.d);
    mpz_mod(&mut ctx.mu, &ctx.mu, &ctx.m);
    // v = m / g
    mpz_fdiv_q(&mut ctx.v, &ctx.m, &ctx.g);
  }
}

// ClassElem and ClassGroup ops based on Chia's fantastic doc explaining applied class groups:
//  https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf.

impl ClassElem {
  fn normalize(&mut self) {
    with_context(|x| self.normalize_ctx(x))
  }

  fn normalize_ctx(&mut self, ctx: &mut Ctx) {
    if self.is_normalized(ctx) {
      return;
    }
    // r = floor_div((a - b), 2a)
    // (a, b, c) = (a, b + 2ra, ar^2 + br + c)
    unsafe {
      mpz_sub(&mut ctx.r, &self.a, &self.b);
      mpz_mul_ui(&mut ctx.denom, &self.a, 2);

      mpz_fdiv_q(&mut ctx.r, &ctx.r, &ctx.denom);

      mpz_set(&mut ctx.old_b, &self.b);

      mpz_mul(&mut ctx.ra, &ctx.r, &self.a);
      mpz_add(&mut self.b, &self.b, &ctx.ra);
      mpz_add(&mut self.b, &self.b, &ctx.ra);

      mpz_mul(&mut ctx.ra, &ctx.ra, &ctx.r);
      mpz_add(&mut self.c, &self.c, &ctx.ra);

      mpz_set(&mut ctx.ra, &ctx.r);
      mpz_mul(&mut ctx.ra, &ctx.ra, &ctx.old_b);
      mpz_add(&mut self.c, &self.c, &ctx.ra);
    }
  }

  fn reduce(&mut self) {
    with_context(|x| self.reduce_ctx(x))
  }

  // TODO: Ensure elements only valid if reduced
  fn reduce_ctx(&mut self, ctx: &mut Ctx) {
    self.normalize_ctx(ctx);
    let mut count = 0;
    while !self.is_reduced() {
      unsafe {
        println!("reduce count: {}", count);
        count += 1;
        // s = floor_div(c + b, 2c)
        mpz_add(&mut ctx.s, &self.c, &self.b);
        // x = 2c
        mpz_mul_ui(&mut ctx.x, &self.c, 2);
        mpz_fdiv_q(&mut ctx.s, &ctx.s, &ctx.x);
        mpz_set(&mut ctx.old_a, &self.a);
        mpz_set(&mut ctx.old_b, &self.b);
        mpz_set(&mut self.a, &self.c);
        mpz_neg(&mut self.b, &self.b);
        // x = 2sc
        mpz_mul(&mut ctx.x, &ctx.s, &self.c);
        mpz_mul_ui(&mut ctx.x, &ctx.x, 2);
        // b = 2sc - b
        mpz_add(&mut self.b, &self.b, &ctx.x);

        // c = cs^2
        mpz_mul(&mut self.c, &self.c, &ctx.s);
        mpz_mul(&mut self.c, &self.c, &ctx.s);
        // x = bs
        mpz_mul(&mut ctx.x, &ctx.old_b, &ctx.s);
        // c -= bs
        mpz_sub(&mut self.c, &self.c, &ctx.x);
        // c += a
        mpz_add(&mut self.c, &self.c, &ctx.old_a);
      }
    }
    self.normalize_ctx(ctx);
  }

  fn square(&mut self) {
    with_context(|x| self.square_ctx(x))
  }

  fn square_ctx(&mut self, ctx: &mut Ctx) {
    unsafe {
      // Solve `bk = c mod a` for k, represented by mu, v and any integer n s.t. k = mu + v * n
      solve_linear_congruence_mpz(ctx);

      // tmp = (b * mu) / a
      mpz_mul(&mut ctx.m, &self.b, &ctx.mu);
      mpz_sub(&mut ctx.m, &ctx.m, &self.c);
      mpz_fdiv_q(&mut ctx.m, &ctx.m, &self.a);

      // A = a^2
      mpz_set(&mut ctx.old_a, &self.a);
      mpz_mul(&mut self.a, &self.a, &self.a);

      // B = b - 2a * mu
      mpz_mul(&mut ctx.a, &ctx.mu, &ctx.old_a);
      mpz_mul_ui(&mut ctx.a, &ctx.a, 2);
      mpz_sub(&mut self.b, &self.b, &ctx.a);

      // C = mu^2 - tmp
      mpz_mul(&mut self.c, &ctx.mu, &ctx.mu);
      mpz_sub(&mut self.c, &self.c, &ctx.m);
    }

    self.reduce_ctx(ctx);
  }

  fn discriminant(&self, d: &mut mpz_t) -> mpz_t {
    let mut tmp = new_mpz();
    unsafe {
      mpz_mul(d, &self.b, &self.b);
      mpz_mul(&mut tmp, &self.a, &self.c);
      mpz_mul_ui(&mut tmp, &tmp, 4);
      mpz_sub(d, d, &tmp);
      *d
    }
  }

  fn validate(&self) -> bool {
    let mut d = new_mpz();
    unsafe {
      &self.discriminant(&mut d);
      mpz_cmp(&d, &ClassGroup::rep().inner) == 0
    }
  }

  // expects normalized element
  fn is_reduced(&self) -> bool {
    unsafe {
      !((mpz_cmp(&self.a, &self.c) > 0)
        || (mpz_cmp(&self.a, &self.c) == 0 && mpz_cmp_si(&self.b, 0) < 0))
    }
  }

  fn is_normalized(&self, ctx: &mut Ctx) -> bool {
    unsafe {
      mpz_neg(&mut ctx.negative_a, &self.a);
      mpz_cmp(&self.b, &ctx.negative_a) > 0 && mpz_cmp(&self.b, &self.a) <= 0
    }
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

impl TypeRep for ClassGroup {
  type Rep = Discriminant;
  fn rep() -> &'static Self::Rep {
    &CLASS_GROUP_DISCRIMINANT
  }
}

impl Group for ClassGroup {
  type Elem = ClassElem;

  fn op_(_: &Discriminant, x: &ClassElem, y: &ClassElem) -> ClassElem {
    with_context(|ctx| {
      unsafe {
        // TODO: Assert discriminants equal
        // g = (b1 + b2) / 2
        mpz_add(&mut ctx.g, &x.b, &y.b);
        mpz_fdiv_q_ui(&mut ctx.g, &ctx.g, 2);
        // h = (b2 - b1) / 2
        mpz_sub(&mut ctx.h, &y.b, &x.b);
        mpz_fdiv_q_ui(&mut ctx.h, &ctx.h, 2);
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
        solve_linear_congruence_mpz(ctx);

        // a = tv
        mpz_mul(&mut ctx.m, &ctx.t, &ctx.v);
        // b = h - t * mu
        mpz_mul(&mut ctx.m, &ctx.t, &ctx.mu);
        mpz_sub(&mut ctx.b, &ctx.h, &ctx.m);
        // m = s
        mpz_set(&mut ctx.m, &ctx.s);
        // Solve linear congruence `(tv)k = h - t * mu mod s` or `ak = b mod m` for solutions k
        solve_linear_congruence_mpz(ctx);

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

        let mut ret = ClassElem::default();
        // A = st - ru
        mpz_mul(&mut ret.a, &ctx.s, &ctx.t);
        mpz_mul(&mut ctx.a, &ctx.r, &ctx.u);
        mpz_sub(&mut ret.a, &ret.a, &ctx.a);

        // B = ju - kt + ls
        mpz_mul(&mut ret.b, &ctx.j, &ctx.u);
        mpz_mul(&mut ctx.a, &ctx.m, &ctx.r);
        mpz_add(&mut ret.b, &ret.b, &ctx.a);
        mpz_mul(&mut ctx.a, &ctx.k, &ctx.t);
        mpz_sub(&mut ret.b, &ret.b, &ctx.a);
        mpz_mul(&mut ctx.a, &ctx.l, &ctx.s);
        mpz_sub(&mut ret.b, &ret.b, &ctx.a);

        // C = kl - jm
        mpz_mul(&mut ret.c, &ctx.k, &ctx.l);
        mpz_mul(&mut ctx.a, &ctx.j, &ctx.m);
        mpz_sub(&mut ret.c, &ret.c, &ctx.a);
        ret.reduce_ctx(ctx);
        ret
      }
    })
  }

  fn id_(d: &Discriminant) -> ClassElem {
    with_context(|ctx| {
      let mut ret = ClassElem::default();
      unsafe {
        mpz_set_ui(&mut ret.a, 1);
        mpz_set_ui(&mut ret.b, 1);
        // c = (b * b - d) / 4a
        mpz_sub(&mut ctx.a, &ret.b, &d.inner); // b == b*b
        mpz_fdiv_q_ui(&mut ret.c, &ctx.a, 4);
      }
      ret
    })
  }

  fn inv_(_: &Discriminant, x: &ClassElem) -> ClassElem {
    let mut ret = ClassElem::default();
    unsafe {
      mpz_set(&mut ret.a, &x.a);
      mpz_neg(&mut ret.b, &x.b);
      mpz_set(&mut ret.c, &x.c);
    }
    ret
  }

  fn exp_(_: &Discriminant, a: &ClassElem, n: &Integer) -> ClassElem {
    with_context(|ctx| {
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
        a.square_ctx(ctx);
        n >>= 1;
      }
    })
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Discriminant) -> ClassElem {
    // a = 2
    // b = 1
    // c = (b * b - d) / 4a
    let mut ret = ClassElem {
      a: new_mpz(),
      b: new_mpz(),
      c: new_mpz(),
    };
    unsafe {
      mpz_set_ui(&mut ret.a, 2);
      mpz_set_ui(&mut ret.b, 1);

      mpz_set_ui(&mut ret.c, 1);
      mpz_sub(&mut ret.c, &ret.c, &d.inner);
      mpz_fdiv_q_ui(&mut ret.c, &ret.c, 8);
    }

    ret.reduce();
    ret
  }
}

impl Hash for ClassElem {
  fn hash<H: Hasher>(&self, state: &mut H) {
    let size = self.a.size;
    size.hash(state);
    if size != 0 {
      let limbs = size.checked_abs().expect("overflow") as usize;
      let slice = unsafe { slice::from_raw_parts(self.a.d, limbs) };
      slice.hash(state);
    }
  }
}

// TODO: This just returns (0,0,0) right now
impl<T> ElemFrom<(T, T, T)> for ClassGroup
//where
//  Integer: From<T>,
{
  fn elem(t: (T, T, T)) -> ClassElem {
    let mut class_elem = ClassElem {
      a: new_mpz(),
      b: new_mpz(),
      c: new_mpz(),
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
    let mut ret = ClassElem {
      a: new_mpz(),
      b: new_mpz(),
      c: new_mpz(),
    };

    let a_str = CString::new(a).unwrap();
    let b_str = CString::new(b).unwrap();
    let c_str = CString::new(c).unwrap();
    unsafe {
      mpz_set_str(&mut ret.a, a_str.as_ptr(), 10);
      mpz_set_str(&mut ret.b, b_str.as_ptr(), 10);
      mpz_set_str(&mut ret.c, c_str.as_ptr(), 10);
    }

    ret
  }

  #[should_panic]
  #[test]
  fn test_bad_elem() {
    let _ = ClassGroup::elem((1, 2, 3));
  }

  #[test]
  /*fn test_elem_from() {
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
  }*/
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
    /* let mut to_reduce = construct_raw_elem_from_strings(
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

    let mut reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );*/
    let mut to_reduce = construct_raw_elem_from_strings(
      "59162244921619725812008939143220718157267937427074598447911241410131470159247784852210767449\
      675610037288729551814191198624164179866076352187405442496568188988272422133088755036699145362\
      385840772236403043664778415471196678638241785773530531198720497580622741709880533724904220122\
      358854068046553219863419609777498761804625479650772123754523807001976654588225908928022367436\
      8",
      "18760351095004839755193532164856605650590306627169248964100884295652838905828158941233738613\
      175821849253748329102319504958410190952820220503570113920576542676928659211807590199941027958\
      195895385446372444261885022800653454209101497963588809819572703579484085278913354621371362285\
      341138299691587953249270188429393417132110841259813122945626515477865766896056280729710478647\
      13",
      "14872270891432803054791175727694631095755964943358394411314110783404577714102170379700365256\
      599679049493824862742803590079461712691146098397470840896560034332315858221821103076776907123\
      277315116632337385101204055232891361405428635972040596205450316747012080794838691280547894128\
      246741601088755087359234554141346980837292342320288111397175220296098629890108459305643419353\
      36"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "26888935961824081232597112540509824504614070059776273347136888921115497522070287009841688662\
      983066376019079593372296556420848446780369918809384119124783870290778875424468497961559643807\
      918398860928578027038014112641529893817109240852544158309292025321122680747989987560029531021\
      808743313150630063377037854944",
      "14529985196481999393995154363327100184407232892559561136140792409262328867440167480822808496\
      853924547751298342980606034124112579835255733824790020119078588372593288210628255956605240171\
      744703418426092073347584357826862813733154338737148962212641444735717023402201569115323580814\
      54099903972209626147819759991",
      "28467266502267127591420289007165819749231433586093061478772560429058231137856046130384492811\
      816456933286039468940950129263300933723839212086399375780796041634531383342902918719073416087\
      614456845205980227091403964285870107268917183244016635907926846271829374679124848388403486656\
      1564478239095738726823372184204"
    );
    assert_eq!(to_reduce, reduced_ground_truth);

    let already_reduced = reduced_ground_truth.clone();

    /*let mut already_reduced = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );*/
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
    let mut d = new_mpz();
    &g.discriminant(&mut d);
    println!("{:?}", d);
    println!("{:?}", &ClassGroup::rep().inner);
    //  let mut d_str = CString::default().as_ptr();
    //  mpz_get_str(d_str, 10, &ClassGroup::rep().inner);
    //println!("{}", d_str.to_str().unwrap());
    assert!(unsafe { mpz_cmp(&d, &ClassGroup::rep().inner) == 0 });
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

    let x = ClassGroup::op(&a, &b);

    assert_eq!(x, ground_truth);
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
