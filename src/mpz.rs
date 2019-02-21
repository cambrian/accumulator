use gmp_mpfr_sys::gmp::{
  mpz_abs, mpz_add, mpz_add_ui, mpz_cdiv_q, mpz_cdiv_r, mpz_cmp, mpz_cmp_si, mpz_cmpabs,
  mpz_divexact, mpz_fdiv_q, mpz_fdiv_q_ui, mpz_fdiv_qr, mpz_fdiv_r, mpz_fits_slong_p, mpz_gcd,
  mpz_gcdext, mpz_get_si, mpz_init, mpz_mod, mpz_mul, mpz_mul_ui, mpz_neg, mpz_odd_p, mpz_root,
  mpz_set, mpz_set_si, mpz_set_str, mpz_set_ui, mpz_sgn, mpz_sub, mpz_submul, mpz_t,
};

use std::cmp::Ordering;
use std::ffi::CString;
use std::hash::{Hash, Hasher};
use std::mem::uninitialized;
use std::ops::Neg;
use std::slice;
use std::str::FromStr;

pub type mp_limb_t = ::std::os::raw::c_ulong;

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct __mpz_struct {
  pub _mp_alloc: ::std::os::raw::c_int,
  pub _mp_size: ::std::os::raw::c_int,
  pub _mp_d: *mut mp_limb_t,
}

#[derive(Debug)]
#[cfg_attr(repr_transparent, repr(transparent))]
pub struct Mpz {
  pub inner: mpz_t,
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

impl From<__mpz_struct> for Mpz {
  fn from(x: __mpz_struct) -> Self {
    let mut ret = Mpz::default();
    ret.inner.alloc = x._mp_alloc;
    ret.inner.size = x._mp_size;
    ret.inner.d = x._mp_d;
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

impl Mpz {
  #[inline]
  pub fn abs(&mut self, x: &Mpz) {
    unsafe { mpz_abs(&mut self.inner, &x.inner) }
  }
  #[inline]
  pub fn add(&mut self, x: &Mpz, y: &Mpz) {
    unsafe {
      mpz_add(&mut self.inner, &x.inner, &y.inner);
    }
  }
  #[inline]
  pub fn add_mut_ui(&mut self, x: u64) {
    unsafe { mpz_add_ui(&mut self.inner, &self.inner, x) }
  }
  #[inline]
  pub fn add_mut(&mut self, x: &Mpz) {
    unsafe { mpz_add(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
  pub fn cmp(&self, other: &Mpz) -> i32 {
    unsafe { mpz_cmp(&self.inner, &other.inner) }
  }
  #[inline]
  pub fn cmp_abs(&self, other: &Mpz) -> i32 {
    unsafe { mpz_cmpabs(&self.inner, &other.inner) }
  }
  #[inline]
  pub fn cmp_si(&self, val: i64) -> i32 {
    unsafe { mpz_cmp_si(&self.inner, val) }
  }
  #[inline]
  pub fn ceil_div(&mut self, x: &Mpz, y: &Mpz) {
    unsafe {
      mpz_cdiv_q(&mut self.inner, &x.inner, &y.inner);
    }
  }
  #[inline]
  pub fn ceil_div_rem(&mut self, x: &Mpz, y: &Mpz) {
    unsafe {
      mpz_cdiv_r(&mut self.inner, &x.inner, &y.inner);
    }
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
  pub fn floor_div_rem(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_fdiv_r(&mut self.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn floor_div_qrem(&mut self, r: &mut Mpz, x: &Mpz, y: &Mpz) {
    unsafe { mpz_fdiv_qr(&mut self.inner, &mut r.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn floor_div_ui(&mut self, x: &Mpz, val: u64) {
    unsafe {
      mpz_fdiv_q_ui(&mut self.inner, &x.inner, val);
    }
  }
  #[inline]
  pub fn get_si(&mut self) -> i64 {
    unsafe { mpz_get_si(&mut self.inner) }
  }
  #[inline]
  pub fn div_exact(&mut self, n: &Mpz, d: &Mpz) {
    unsafe { mpz_divexact(&mut self.inner, &n.inner, &d.inner) }
  }
  #[inline]
  pub fn div_exact_mut(&mut self, d: &Mpz) {
    unsafe {
      mpz_divexact(&mut self.inner, &self.inner, &d.inner);
    }
  }
  #[inline]
  pub fn fits_slong_p(&mut self) -> i32 {
    unsafe { mpz_fits_slong_p(&mut self.inner) }
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
  pub fn odd(&mut self) -> i32 {
    unsafe { mpz_odd_p(&mut self.inner) }
  }
  #[inline]
  pub fn root_mut(&mut self, x: u64) -> i32 {
    unsafe { mpz_root(&mut self.inner, &self.inner, x) }
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
  pub fn set_si(&mut self, val: i64) {
    unsafe { mpz_set_si(&mut self.inner, val) }
  }
  #[inline]
  pub fn sgn(&self) -> i32 {
    unsafe { mpz_sgn(&self.inner) }
  }
  #[inline]
  pub fn sub(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_sub(&mut self.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn sub_mul(&mut self, x: &Mpz, y: &Mpz) {
    unsafe { mpz_submul(&mut self.inner, &x.inner, &y.inner) }
  }
  #[inline]
  pub fn sub_mut(&mut self, x: &Mpz) {
    unsafe { mpz_sub(&mut self.inner, &self.inner, &x.inner) }
  }
  #[inline]
  pub fn square_mut(&mut self) {
    unsafe { mpz_mul(&mut self.inner, &self.inner, &self.inner) }
  }
}

unsafe impl Send for Mpz {}
unsafe impl Sync for Mpz {}
