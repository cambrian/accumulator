use gmp_mpfr_sys::gmp::{
  mpz_add, mpz_cmp, mpz_cmp_si, mpz_cmpabs, mpz_fdiv_q, mpz_fdiv_q_ui, mpz_fdiv_qr, mpz_gcd,
  mpz_gcdext, mpz_getlimbn, mpz_init, mpz_mod, mpz_mul, mpz_mul_si, mpz_mul_ui, mpz_neg, mpz_set,
  mpz_set_str, mpz_set_ui, mpz_sgn, mpz_size, mpz_sub, mpz_swap, mpz_t,
};
use std::cmp::Ordering;
use std::ffi::CString;
use std::hash::{Hash, Hasher};
use std::mem::uninitialized;
use std::slice;
use std::str::FromStr;

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

impl From<u64> for Mpz {
  fn from(x: u64) -> Self {
    let mut ret = Mpz::default();
    unsafe { mpz_set_ui(&mut ret.inner, x) };
    ret
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

unsafe impl Send for Mpz {}
unsafe impl Sync for Mpz {}

impl Mpz {
  #[inline]
  pub fn add(&mut self, x: &Mpz, y: &Mpz) {
    unsafe {
      mpz_add(&mut self.inner, &x.inner, &y.inner);
    }
  }
  #[inline]
  pub fn add_mut(&mut self, x: &Mpz) {
    unsafe { mpz_add(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
  pub fn cmp(&self, other: &Mpz) -> i32 {
    unsafe { mpz_cmp(&self.inner, &other.inner) }
  }
  pub fn cmp_abs(&self, other: &Mpz) -> i32 {
    unsafe { mpz_cmpabs(&self.inner, &other.inner) }
  }
  #[inline]
  pub fn cmp_si(&self, val: i64) -> i32 {
    unsafe { mpz_cmp_si(&self.inner, val) }
  }
  #[inline]
  pub fn floor_div(&mut self, x: &Mpz, y: &Mpz) {
    unsafe {
      mpz_fdiv_q(&mut self.inner, &x.inner, &y.inner);
    }
  }
  #[inline]
  pub fn floor_div_mut(&mut self, x: &Mpz) {
    unsafe { mpz_fdiv_q(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
  pub fn floor_div_rem(&mut self, r: &mut Mpz, x: &Mpz, y: &Mpz) {
    unsafe { mpz_fdiv_qr(&mut self.inner, &mut r.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn floor_div_ui(&mut self, x: &Mpz, val: u64) {
    unsafe {
      mpz_fdiv_q_ui(&mut self.inner, &x.inner, val);
    }
  }
  #[inline]
  pub fn floor_div_ui_mut(&mut self, val: u64) {
    unsafe {
      mpz_fdiv_q_ui(&mut self.inner, &self.inner, val);
    }
  }
  #[inline]
  pub fn gcd(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_gcd(&mut self.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn gcd_mut(&mut self, x: &Mpz) {
    unsafe { mpz_gcd(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
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
  #[inline]
  pub fn limbn(&self, size: i64) -> u64 {
    unsafe { mpz_getlimbn(&self.inner, size) as u64 }
  }
  #[inline]
  pub fn modulo(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_mod(&mut self.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn modulo_mut(&mut self, x: &Mpz) {
    unsafe { mpz_mod(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
  pub fn mul(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_mul(&mut self.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn mul_mut(&mut self, x: &Mpz) {
    unsafe { mpz_mul(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
  pub fn mul_si(&mut self, x: &Mpz, val: i64) {
    unsafe { mpz_mul_si(&mut self.inner, &x.inner, val) }
  }
  #[inline]
  pub fn mul_ui(&mut self, x: &Mpz, val: u64) {
    unsafe { mpz_mul_ui(&mut self.inner, &x.inner, val) }
  }
  #[inline]
  pub fn mul_ui_mut(&mut self, val: u64) {
    unsafe { mpz_mul_ui(&mut self.inner, &self.inner, val) }
  }
  #[inline]
  pub fn neg(&mut self, x: &Mpz) {
    unsafe { mpz_neg(&mut self.inner, &x.inner) }
  }
  #[inline]
  pub fn neg_mut(&mut self) {
    unsafe { mpz_neg(&mut self.inner, &self.inner) }
  }
  #[inline]
  pub fn is_neg(&self) -> bool {
    unsafe { mpz_sgn(&self.inner) == -1 }
  }
  #[inline]
  pub fn set(&mut self, x: &Mpz) {
    unsafe { mpz_set(&mut self.inner, &x.inner) }
  }
  #[inline]
  pub fn set_cstr(&mut self, cs: &CString) {
    unsafe {
      mpz_set_str(&mut self.inner, cs.as_ptr(), 10);
    }
  }
  #[inline]
  pub fn set_ui(&mut self, val: u64) {
    unsafe { mpz_set_ui(&mut self.inner, val) }
  }
  #[inline]
  pub fn size(&self) -> u64 {
    unsafe { mpz_size(&self.inner) as u64 }
  }
  #[inline]
  pub fn sub(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_sub(&mut self.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn sub_mut(&mut self, x: &Mpz) {
    unsafe { mpz_sub(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
  pub fn square_mut(&mut self) {
    unsafe { mpz_mul(&mut self.inner, &self.inner, &self.inner) }
  }
  #[inline]
  pub fn swap(&mut self, other: &mut Mpz) {
    unsafe { mpz_swap(&mut self.inner, &mut other.inner) }
  }
}
