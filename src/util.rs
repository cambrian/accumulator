use crate::group::Group;
use gmp_mpfr_sys::gmp::{
  mpz_add, mpz_cmp, mpz_cmp_si, mpz_cmp_ui, mpz_fdiv_q, mpz_fdiv_q_ui, mpz_fdiv_qr, mpz_gcd,
  mpz_gcdext, mpz_get_str, mpz_init, mpz_mod, mpz_mul, mpz_mul_ui, mpz_neg, mpz_set, mpz_set_str,
  mpz_set_ui, mpz_sub, mpz_t,
};
use rug::Integer;
use std::cmp::Ordering;
use std::ffi::CString;
use std::hash::{Hash, Hasher};
use std::mem::uninitialized;
use std::num::ParseIntError;
use std::slice;
use std::str::FromStr;

/// Poor man's type-level programming.
/// This trait allows us to reflect "type-level" (i.e. static) information at runtime.
pub trait TypeRep: 'static {
  type Rep: 'static;
  fn rep() -> &'static Self::Rep;
}

#[derive(Debug)]
pub enum UtilityError {
  NoSolutionToLinearCongruence,
}

pub fn int<T>(val: T) -> Integer
where
  Integer: From<T>,
{
  Integer::from(val)
}

#[derive(Debug)]
#[cfg_attr(repr_transparent, repr(transparent))]
pub struct Mpz {
  inner: mpz_t,
}

impl Default for Mpz {
  fn default() -> Self {
    let inner = unsafe {
      let mut ret = uninitialized();
      mpz_init(&mut ret);
      ret
    };
    Self { inner }
  }
}

impl Clone for Mpz {
  fn clone(&self) -> Self {
    let mut ret = Mpz::default();
    ret.set(&self);
    ret
  }
}

impl PartialEq for Mpz {
  fn eq(&self, other: &Mpz) -> bool {
    self.cmp(&other) == 0
  }
}

impl Eq for Mpz {}

impl PartialOrd for Mpz {
  fn partial_cmp(&self, other: &Mpz) -> Option<Ordering> {
    match self.cmp(&other) {
      x if x < 0 => Some(Ordering::Less),
      0 => Some(Ordering::Equal),
      _ => Some(Ordering::Greater),
    }
  }
}

impl Ord for Mpz {
  fn cmp(&self, other: &Mpz) -> Ordering {
    match self.cmp(&other) {
      x if x < 0 => Ordering::Less,
      0 => Ordering::Equal,
      _ => Ordering::Greater,
    }
  }
}

impl Hash for Mpz {
  fn hash<H: Hasher>(&self, state: &mut H) {
    let size = self.inner.size;
    size.hash(state);
    if size != 0 {
      let limbs = size.checked_abs().expect("overflow") as usize;
      let slice = unsafe { slice::from_raw_parts(self.inner.d, limbs) };
      slice.hash(state);
    }
  }
}

impl FromStr for Mpz {
  type Err = std::ffi::NulError;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    let mut ret = Mpz::default();
    let c_str = CString::new(s)?;
    ret.set_cstr(&c_str);
    Ok(ret)
  }
}

// TODO: Make functions inline?
impl Mpz {
  pub fn add(&mut self, x: &Mpz, y: &Mpz) {
    unsafe {
      mpz_add(&mut self.inner, &x.inner, &y.inner);
    }
  }
  pub fn add_mut(&mut self, x: &Mpz) {
    unsafe { mpz_add(&mut self.inner, &self.inner, &x.inner) }
  }
  pub fn cmp(&self, other: &Mpz) -> i32 {
    unsafe { mpz_cmp(&self.inner, &other.inner) }
  }
  pub fn cmp_si(&self, val: i64) -> i32 {
    unsafe { mpz_cmp_si(&self.inner, val) }
  }
  pub fn floor_div(&mut self, x: &Mpz, y: &Mpz) {
    unsafe {
      mpz_fdiv_q(&mut self.inner, &x.inner, &y.inner);
    }
  }
  pub fn floor_div_mut(&mut self, x: &Mpz) {
    unsafe { mpz_fdiv_q(&mut self.inner, &self.inner, &x.inner) }
  }
  pub fn floor_div_rem(&mut self, r: &mut Mpz, x: &Mpz, y: &Mpz) {
    unsafe { mpz_fdiv_qr(&mut self.inner, &mut r.inner, &x.inner, &y.inner) }
  }
  pub fn floor_div_ui(&mut self, x: &Mpz, val: u64) {
    unsafe {
      mpz_fdiv_q_ui(&mut self.inner, &x.inner, val);
    }
  }
  pub fn floor_div_ui_mut(&mut self, val: u64) {
    unsafe {
      mpz_fdiv_q_ui(&mut self.inner, &self.inner, val);
    }
  }
  pub fn gcd(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_gcd(&mut self.inner, &x.inner, &y.inner) }
  }
  pub fn gcd_mut(&mut self, x: &Mpz) {
    unsafe { mpz_gcd(&mut self.inner, &self.inner, &x.inner) }
  }
  pub fn gcd_cofactors(&mut self, d: &mut Mpz, e: &mut Mpz, a: &Mpz, m: &Mpz) {
    unsafe {
      mpz_gcdext(
        &mut self.inner,
        &mut d.inner,
        &mut e.inner,
        &a.inner,
        &m.inner,
      )
    }
  }
  pub fn modulo(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_mod(&mut self.inner, &x.inner, &y.inner) }
  }
  pub fn modulo_mut(&mut self, x: &Mpz) {
    unsafe { mpz_mod(&mut self.inner, &self.inner, &x.inner) }
  }
  pub fn mul(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_mul(&mut self.inner, &x.inner, &y.inner) }
  }
  pub fn mul_mut(&mut self, x: &Mpz) {
    unsafe { mpz_mul(&mut self.inner, &self.inner, &x.inner) }
  }
  pub fn mul_ui(&mut self, x: &Mpz, val: u64) {
    unsafe { mpz_mul_ui(&mut self.inner, &x.inner, val) }
  }
  pub fn mul_ui_mut(&mut self, val: u64) {
    unsafe { mpz_mul_ui(&mut self.inner, &self.inner, val) }
  }
  pub fn neg(&mut self, x: &Mpz) {
    unsafe { mpz_neg(&mut self.inner, &x.inner) }
  }
  pub fn neg_mut(&mut self) {
    unsafe { mpz_neg(&mut self.inner, &self.inner) }
  }
  pub fn set(&mut self, x: &Mpz) {
    unsafe { mpz_set(&mut self.inner, &x.inner) }
  }
  pub fn set_cstr(&mut self, cs: &CString) {
    unsafe {
      mpz_set_str(&mut self.inner, cs.as_ptr(), 10);
    }
  }
  pub fn set_ui(&mut self, val: u64) {
    unsafe { mpz_set_ui(&mut self.inner, val) }
  }
  pub fn sub(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_sub(&mut self.inner, &x.inner, &y.inner) }
  }
  pub fn sub_mut(&mut self, x: &Mpz) {
    unsafe { mpz_sub(&mut self.inner, &self.inner, &x.inner) }
  }
  pub fn square_mut(&mut self) {
    unsafe { mpz_mul(&mut self.inner, &self.inner, &self.inner) }
  }
}

unsafe impl Send for Mpz {}
unsafe impl Sync for Mpz {}

/// Computes the `(xy)`th root of `g` given the `x`th and `y`th roots of `g` and `(x, y)` coprime.
/// Consider moving this to accumulator?
pub fn shamir_trick<G: Group>(
  xth_root: &G::Elem,
  yth_root: &G::Elem,
  x: &Integer,
  y: &Integer,
) -> Option<G::Elem> {
  if G::exp(xth_root, x) != G::exp(yth_root, y) {
    return None;
  }

  let (gcd, a, b) = <(Integer, Integer, Integer)>::from(x.gcd_cofactors_ref(&y));

  if gcd != int(1) {
    return None;
  }

  Some(G::op(&G::exp(xth_root, &b), &G::exp(yth_root, &a)))
}

// Solve a linear congruence of form `ax = b mod m` for the set of solutions x,
// characterized by integers mu and v such that x = mu + vn where n is any integer.
pub fn solve_linear_congruence(
  a: &Integer,
  b: &Integer,
  m: &Integer,
) -> Result<(Integer, Integer), UtilityError> {
  // g = gcd(a, m) => da + em = g
  let (g, d, _) = a.clone().gcd_cofactors(m.clone(), Integer::new());

  // q = floor_div(b, g)
  // r = b % g
  let (q, r) = b.clone().div_rem_floor(g.clone());
  if r != Integer::from(0) {
    return Err(UtilityError::NoSolutionToLinearCongruence);
  }

  // mu = (q * d) % m
  // v = m / g
  let mu = (q * d) % m;
  let (v, _) = m.clone().div_rem_floor(g);
  Ok((mu, v))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{Group, Rsa2048, UnknownOrderGroup};
  use crate::util::int;

  #[test]
  fn test_linear_congruence_solver() {
    assert_eq!(
      (Integer::from(-2), Integer::from(4)),
      solve_linear_congruence(&Integer::from(3), &Integer::from(2), &Integer::from(4)).unwrap()
    );

    assert_eq!(
      (Integer::from(-2), Integer::from(4)),
      solve_linear_congruence(&Integer::from(3), &Integer::from(2), &Integer::from(4)).unwrap()
    );

    assert_eq!(
      (Integer::from(1), Integer::from(2)),
      solve_linear_congruence(&Integer::from(5), &Integer::from(1), &Integer::from(2)).unwrap()
    );

    assert_eq!(
      (Integer::from(-3), Integer::from(5)),
      solve_linear_congruence(&Integer::from(2), &Integer::from(4), &Integer::from(5)).unwrap()
    );

    assert_eq!(
      (Integer::from(2491), Integer::from(529)),
      solve_linear_congruence(
        &Integer::from(230),
        &Integer::from(1081),
        &Integer::from(12167)
      )
      .unwrap()
    );
  }

  #[test]
  fn test_linear_congruence_solver_no_solution() {
    // Let g = gcd(a, m). If b is not divisible by g, there are no solutions. If b is divisible by
    // g, there are g solutions.
    let result =
      solve_linear_congruence(&Integer::from(33), &Integer::from(7), &Integer::from(143));
    assert!(result.is_err());

    let result =
      solve_linear_congruence(&Integer::from(13), &Integer::from(14), &Integer::from(39));
    assert!(result.is_err());
  }

  #[test]
  fn test_shamir_trick() {
    let (x, y, z) = (&int(13), &int(17), &int(19));
    let xth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(y * z));
    let yth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(x * z));
    let xyth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), z);
    assert!(shamir_trick::<Rsa2048>(&xth_root, &yth_root, x, y) == Some(xyth_root));
  }

  #[test]
  fn test_shamir_trick_failure() {
    let (x, y, z) = (&int(7), &int(14), &int(19)); // Inputs not coprime.
    let xth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(y * z));
    let yth_root = Rsa2048::exp(&Rsa2048::unknown_order_elem(), &int(x * z));
    assert!(shamir_trick::<Rsa2048>(&xth_root, &yth_root, x, y) == None);
  }
}
