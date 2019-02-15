//! TODO: reduce I256/I512 duplication with a macro
use gmp_mpfr_sys::gmp;
use gmp_mpfr_sys::gmp::mpz_t;
use std::cmp::{min, Ord, Ordering, PartialOrd};
use std::convert::From;
use std::fmt::Debug;
use std::mem::transmute;
use std::ops;

#[derive(PartialEq, Eq, Hash, Debug, Clone, Copy)]
pub struct I512 {
  size: i64,       // Number of limbs in use. Negative size represents a negative number.
  limbs: [u64; 8], // GMP limbs are lower-endian
}

impl I512 {
  fn ptr(&self) -> *const u64 {
    &self.limbs as *const u64
  }
  fn ptr_mut(&mut self) -> *mut u64 {
    &self.limbs as *const u64 as *mut u64
  }
  fn as_mpz(&self) -> mpz_t {
    mpz_t {
      size: self.size as i32,
      d: self.ptr() as *mut u64,
      alloc: 8,
    }
  }
  pub fn zero() -> Self {
    I512 {
      size: 0,
      limbs: [0; 8],
    }
  }

  pub fn one() -> Self {
    let mut limbs = [0; 8];
    limbs[0] = 1;
    I512 { size: 1, limbs }
  }

  pub fn minus_one() -> Self {
    let mut limbs = [0; 8];
    limbs[0] = 1;
    I512 { size: -1, limbs }
  }

  fn normalize_size(&mut self) {
    self.size = 0;
    for i in (0..8).rev() {
      if self.limbs[i] != 0 {
        self.size = (i + 1) as i64;
        break;
      }
    }
  }

  /// Returns the lower half of this I512 as a I256
  pub fn low_i256(self) -> I256 {
    let mut x = unsafe { transmute::<I512, (I256, [u64; 4])>(self) }.0;
    x.normalize_size();
    x
  }

  /// Returns (result, remainder)
  pub fn div_rem(self, x: Self) -> (Self, Self) {
    let (mut y, mut rem) = (I512::zero(), I512::zero());
    unsafe {
      gmp::mpn_tdiv_qr(
        y.ptr_mut(),
        rem.ptr_mut(),
        0,
        self.ptr(),
        self.size,
        x.ptr(),
        x.size,
      )
    };
    y.normalize_size();
    rem.normalize_size();
    // dbg!((self, x, y, rem));
    (y, rem)
  }

  /// mutates self to the remainder, returning result;
  pub fn div_rem_mut(&mut self, x: Self) -> Self {
    let mut y = I512::zero();
    unsafe {
      gmp::mpn_tdiv_qr(
        y.ptr_mut(),
        self.ptr_mut(),
        0,
        self.ptr(),
        self.size,
        x.ptr(),
        x.size,
      )
    };
    self.normalize_size();
    y.normalize_size();
    y
  }
}

/// Lower-endian u64s
impl From<[u64; 8]> for I512 {
  fn from(limbs: [u64; 8]) -> Self {
    let mut x = I512 { size: 0, limbs };
    x.normalize_size();
    x
  }
}

impl From<u64> for I512 {
  fn from(x: u64) -> Self {
    let limbs = [x, 0, 0, 0, 0, 0, 0, 0];
    Self::from(limbs)
  }
}

// TODO: more efficient implementation via memcpy?
impl From<I256> for I512 {
  fn from(x: I256) -> Self {
    let mut limbs = [0; 8];
    limbs[..4].copy_from_slice(&x.limbs);
    I512 {
      size: x.size,
      limbs,
    }
  }
}

impl ops::ShlAssign<u32> for I512 {
  fn shl_assign(&mut self, mut x: u32) {
    loop {
      let sz = min(gmp::LIMB_BITS as u32, x);
      if sz == 0 {
        break;
      };
      x -= sz;
      unsafe { gmp::mpn_lshift(self.ptr_mut(), self.ptr(), 8, sz) };
    }
    self.normalize_size();
  }
}

impl ops::Shl<u32> for I512 {
  type Output = I512;
  fn shl(self, x: u32) -> I512 {
    let mut y = self;
    y <<= x;
    y
  }
}

impl ops::ShrAssign<u32> for I512 {
  fn shr_assign(&mut self, mut x: u32) {
    loop {
      let sz = min(gmp::LIMB_BITS as u32, x);
      if sz == 0 {
        break;
      };
      x -= sz;
      unsafe { gmp::mpn_rshift(self.ptr_mut(), self.ptr(), 8, sz) };
    }
    self.normalize_size();
  }
}

impl ops::Shr<u32> for I512 {
  type Output = I512;
  fn shr(self, x: u32) -> I512 {
    let mut y = self;
    y >>= x;
    y
  }
}

impl ops::AddAssign for I512 {
  /// panics if result overflows.
  fn add_assign(&mut self, x: Self) {
    let carry = unsafe { gmp::mpn_add_n(self.ptr_mut(), self.ptr(), x.ptr(), 8) };
    assert!(carry == 0);
    self.normalize_size();
  }
}

impl ops::Add for I512 {
  /// panics if result overflows.
  type Output = Self;
  fn add(self, x: Self) -> Self {
    let mut y = self;
    y += x;
    y
  }
}

impl ops::SubAssign for I512 {
  /// panics if result is negative.
  fn sub_assign(&mut self, x: Self) {
    let borrow = unsafe { gmp::mpn_sub_n(self.ptr_mut(), self.ptr(), x.ptr(), 8) };
    assert!(borrow == 0);
    self.normalize_size();
  }
}

impl ops::Sub for I512 {
  type Output = Self;
  /// panics if result is negative.
  fn sub(self, x: Self) -> Self {
    let mut y = self;
    y -= x;
    y
  }
}

impl ops::Sub<u64> for I512 {
  type Output = Self;
  /// panics if result is negative.
  fn sub(self, x: u64) -> Self {
    self - I512::from(x)
  }
}

impl ops::Mul for I512 {
  type Output = I512;
  fn mul(self, x: Self) -> I512 {
    let mut y = I512::zero();
    if self.size >= x.size {
      unsafe { gmp::mpn_mul(y.ptr_mut(), self.ptr(), self.size, x.ptr(), x.size) };
    } else {
      unsafe { gmp::mpn_mul(y.ptr_mut(), x.ptr(), x.size, self.ptr(), self.size) };
    }
    y.normalize_size();
    y
  }
}

impl ops::Mul for &I512 {
  type Output = I512;
  fn mul(self, x: Self) -> I512 {
    let mut y = I512::zero();
    if self.size >= x.size {
      unsafe { gmp::mpn_mul(y.ptr_mut(), self.ptr(), self.size, x.ptr(), x.size) };
    } else {
      unsafe { gmp::mpn_mul(y.ptr_mut(), x.ptr(), x.size, self.ptr(), self.size) };
    }
    y.normalize_size();
    y
  }
}

impl ops::Div for I512 {
  type Output = Self;
  fn div(self, x: Self) -> Self {
    self.div_rem(x).0
  }
}

impl ops::Rem for I512 {
  type Output = Self;
  fn rem(self, x: Self) -> Self {
    self.div_rem(x).1
  }
}

impl ops::RemAssign for I512 {
  fn rem_assign(&mut self, x: Self) {
    self.div_rem_mut(x);
  }
}

impl ops::Rem<I256> for I512 {
  type Output = I256;
  fn rem(self, x: I256) -> I256 {
    self.div_rem(i512(x)).1.low_i256()
  }
}

#[derive(PartialEq, Eq, Hash, Debug, Clone, Copy)]
pub struct I256 {
  size: i64,       // Number of limbs in use. Negative size represents a negative number.
  limbs: [u64; 4], // GMP limbs are lower-endian
}

#[allow(unused_mut)]
fn mut_ptr<T>(mut t: &T) -> *mut T {
  t as *const T as *mut T
}

impl I256 {
  fn ptr(&self) -> *const u64 {
    &self.limbs as *const u64
  }
  fn ptr_mut(&mut self) -> *mut u64 {
    &self.limbs as *const u64 as *mut u64
  }
  fn as_mpz(&self) -> mpz_t {
    mpz_t {
      size: self.size as i32,
      d: self.ptr() as *mut u64,
      alloc: 4,
    }
  }
  pub fn zero() -> Self {
    Self {
      size: 0,
      limbs: [0; 4],
    }
  }

  pub fn one() -> Self {
    let mut limbs = [0; 4];
    limbs[0] = 1;
    Self { size: 1, limbs }
  }

  pub fn minus_one() -> Self {
    let mut limbs = [0; 4];
    limbs[0] = 1;
    Self { size: -1, limbs }
  }

  fn normalize_size(&mut self) {
    self.size = 0;
    for i in (0..4).rev() {
      if self.limbs[i] != 0 {
        self.size = (i + 1) as i64;
        break;
      }
    }
  }

  pub fn is_odd(&self) -> bool {
    self.limbs[0] & 1 == 1
  }

  /// Returns (result, remainder)
  pub fn div_rem(self, x: Self) -> (Self, Self) {
    let (mut y, mut rem) = (Self::zero(), Self::zero());
    unsafe {
      gmp::mpn_tdiv_qr(
        y.ptr_mut(),
        rem.ptr_mut(),
        0,
        self.ptr(),
        self.size,
        x.ptr(),
        x.size,
      )
    };
    y.normalize_size();
    rem.normalize_size();
    (y, rem)
  }

  /// mutates self to the remainder, returning result;
  pub fn div_rem_mut(&mut self, x: Self) -> Self {
    let mut y = Self::zero();
    unsafe {
      gmp::mpn_tdiv_qr(
        y.ptr_mut(),
        self.ptr_mut(),
        0,
        self.ptr(),
        self.size,
        x.ptr(),
        x.size,
      )
    };
    self.normalize_size();
    y.normalize_size();
    y
  }

  /// returns (result of removing all fs, number of fs removed)
  pub fn remove_factor(self, f: Self) -> (Self, u64) {
    // for some reason this needs extra scratch space
    let mut out = I512::zero();
    let outmpz = out.as_mpz();
    let s = self.as_mpz();
    let f = f.as_mpz();
    let c = unsafe { gmp::mpz_remove(mut_ptr(&outmpz), mut_ptr(&s), mut_ptr(&f)) };
    out.size = i64::from(outmpz.size);
    (out.low_i256(), c)
  }

  pub fn mod_inv(self, m: Self) -> Self {
    let mut out = I512::zero();
    let outmpz = out.as_mpz();
    let s = self.as_mpz();
    let m = m.as_mpz();
    let exists = unsafe { gmp::mpz_invert(mut_ptr(&outmpz), mut_ptr(&s), mut_ptr(&m)) };
    assert!(exists != 0);
    out.size = i64::from(outmpz.size);
    out.low_i256()
  }

  pub fn pow_mod(self, e: Self, m: Self) -> Self {
    let mut out = I512::zero();
    let outmpz = out.as_mpz();
    let s = self.as_mpz();
    let e = e.as_mpz();
    let m = m.as_mpz();
    unsafe { gmp::mpz_powm(mut_ptr(&outmpz), mut_ptr(&s), mut_ptr(&e), mut_ptr(&m)) };
    out.size = i64::from(outmpz.size);
    out.low_i256()
  }
}

/// Lower-endian bytes
impl From<[u8; 32]> for I256 {
  fn from(bytes: [u8; 32]) -> Self {
    let chunks = unsafe { transmute::<[u8; 32], [[u8; 8]; 4]>(bytes) };
    let mut limbs = [0; 4];
    for i in 0..4 {
      limbs[i] = u64::from_le_bytes(chunks[i]);
    }
    Self::from(limbs)
  }
}

/// Lower-endian bytes
impl From<&[u8; 32]> for I256 {
  fn from(bytes: &[u8; 32]) -> Self {
    let chunks = unsafe { transmute::<[u8; 32], [[u8; 8]; 4]>(*bytes) };
    let mut limbs = [0; 4];
    for i in 0..4 {
      limbs[i] = u64::from_le_bytes(chunks[i]);
    }
    Self::from(limbs)
  }
}

/// Lower-endian u64s
impl From<[u64; 4]> for I256 {
  fn from(limbs: [u64; 4]) -> Self {
    let mut x = I256 { size: 0, limbs };
    x.normalize_size();
    x
  }
}

impl From<u64> for I256 {
  fn from(x: u64) -> Self {
    let limbs = [x, 0, 0, 0];
    Self::from(limbs)
  }
}

impl PartialOrd for I256 {
  fn partial_cmp(&self, x: &I256) -> Option<Ordering> {
    let x = unsafe { gmp::mpn_cmp(self.ptr(), x.ptr(), 4) };
    Some({
      if x < 0 {
        Ordering::Less
      } else if x == 0 {
        Ordering::Equal
      } else {
        Ordering::Greater
      }
    })
  }
}

impl Ord for I256 {
  fn cmp(&self, x: &I256) -> Ordering {
    let x = unsafe { gmp::mpn_cmp(self.ptr(), x.ptr(), 4) };
    if x < 0 {
      Ordering::Less
    } else if x == 0 {
      Ordering::Equal
    } else {
      Ordering::Greater
    }
  }
}

impl ops::ShlAssign<u32> for I256 {
  fn shl_assign(&mut self, mut x: u32) {
    loop {
      let sz = min(gmp::LIMB_BITS as u32, x);
      if sz == 0 {
        break;
      };
      x -= sz;
      unsafe { gmp::mpn_lshift(self.ptr_mut(), self.ptr(), 4, sz) };
    }
    self.normalize_size();
  }
}

impl ops::Shl<u32> for I256 {
  type Output = I256;
  fn shl(self, x: u32) -> I256 {
    let mut y = self;
    y <<= x;
    y
  }
}

impl ops::ShrAssign<u32> for I256 {
  fn shr_assign(&mut self, mut x: u32) {
    loop {
      let sz = min(gmp::LIMB_BITS as u32, x);
      if sz == 0 {
        break;
      };
      x -= sz;
      unsafe { gmp::mpn_rshift(self.ptr_mut(), self.ptr(), 4, sz) };
    }
    self.normalize_size();
  }
}

impl ops::Shr<u32> for I256 {
  type Output = I256;
  fn shr(self, x: u32) -> I256 {
    let mut y = self;
    y >>= x;
    y
  }
}

impl ops::AddAssign for I256 {
  fn add_assign(&mut self, x: Self) {
    let carry = unsafe { gmp::mpn_add_n(self.ptr_mut(), self.ptr(), x.ptr(), 4) };
    assert!(carry == 0);
    self.normalize_size();
  }
}

impl ops::Add for I256 {
  type Output = Self;
  fn add(self, x: Self) -> Self {
    let mut y = self;
    y += x;
    y
  }
}

impl ops::Add<u64> for I256 {
  type Output = Self;
  /// panics if result overflows.
  fn add(self, x: u64) -> Self {
    self + I256::from(x)
  }
}

impl ops::SubAssign for I256 {
  /// panics if result is negative.
  fn sub_assign(&mut self, x: Self) {
    let borrow = unsafe { gmp::mpn_sub_n(self.ptr_mut(), self.ptr(), x.ptr(), 4) };
    assert!(borrow == 0);
    self.normalize_size();
  }
}

impl ops::Sub for I256 {
  type Output = Self;
  /// panics if result is negative.
  fn sub(self, x: Self) -> Self {
    let mut y = self;
    y -= x;
    y
  }
}

impl ops::Sub for &I256 {
  type Output = I256;
  /// panics if result is negative.
  fn sub(self, x: Self) -> I256 {
    let mut y = I256::zero();
    let borrow = unsafe { gmp::mpn_sub_n(y.ptr_mut(), self.ptr(), x.ptr(), 4) };
    assert!(borrow == 0);
    y.normalize_size();
    y
  }
}

impl ops::Sub<u64> for I256 {
  type Output = Self;
  /// panics if result is negative.
  fn sub(self, x: u64) -> Self {
    self - I256::from(x)
  }
}

impl ops::Sub<u64> for &I256 {
  type Output = I256;
  /// panics if result is negative.
  fn sub(self, x: u64) -> I256 {
    self - &I256::from(x)
  }
}

impl ops::Mul for I256 {
  type Output = I512;
  fn mul(self, x: Self) -> I512 {
    let mut y = I512::zero();
    unsafe { gmp::mpn_mul_n(y.ptr_mut(), self.ptr(), x.ptr(), 4) };
    y.normalize_size();
    y
  }
}

impl ops::Mul for &I256 {
  type Output = I512;
  fn mul(self, x: Self) -> I512 {
    let mut y = I512::zero();
    unsafe { gmp::mpn_mul_n(y.ptr_mut(), self.ptr(), x.ptr(), 4) };
    y.normalize_size();
    y
  }
}

impl ops::Div for I256 {
  type Output = Self;
  fn div(self, x: Self) -> Self {
    self.div_rem(x).0
  }
}

impl ops::Rem for I256 {
  type Output = Self;
  fn rem(self, x: Self) -> Self {
    self.div_rem(x).1
  }
}

impl ops::RemAssign for I256 {
  fn rem_assign(&mut self, x: Self) {
    self.div_rem_mut(x);
  }
}

pub fn i256<T>(t: T) -> I256
where
  I256: From<T>,
{
  I256::from(t)
}

pub fn i512<T>(t: T) -> I512
where
  I512: From<T>,
{
  I512::from(t)
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct Reduced<T: PartialEq + Eq + Debug + Clone + Copy>(T);

/// We choose r = 2^256 for simplicity. Therefore this reducer only works for odd m where
/// 3 <= m < 2^256.
pub struct MontgomeryReducer {
  pub m: I256,                          // modulus
  r_inv: I256,                          // r_inv = r^(-1) (mod m)
  k: I512,                              // k = (r*r_inv - 1) / m
  pub one_reduced: Reduced<I256>,       // the value of 1, reduced
  pub minus_one_reduced: Reduced<I256>, // the value of -1, reduced
}

impl MontgomeryReducer {
  /// m must be odd and >= 3
  /// We choose r to be 2^256 for simplicity. For correctness, it need only be greater than m and
  /// coprime to m.
  /// This lets us simplify x & r_mask to x.low_I256().
  pub fn new(m: &I256) -> Self {
    assert!(m.is_odd() && *m >= i256(3));
    let r = i512(1) << 256;
    let r_inv = (r % *m).mod_inv(*m);
    let k = ((i512(r_inv) << 256) - 1) / i512(*m);
    let one_reduced = Reduced(r % *m);
    let minus_one_reduced = Reduced(((i512(*m) - 1) << 256) % *m);
    MontgomeryReducer {
      m: *m,
      r_inv,
      k,
      one_reduced,
      minus_one_reduced,
    }
  }

  pub fn reduce(&self, a: I256) -> Reduced<I256> {
    Reduced((i512(a) << 256) % self.m)
  }
  pub fn unreduce(&self, a: Reduced<I256>) -> I256 {
    (a.0 * self.r_inv) % self.m
  }

  pub fn mul_mod_(&self, a: Reduced<I256>, b: Reduced<I256>) -> Reduced<I256> {
    let c = a.0 * b.0;
    let temp = (i512(c.low_i256()) * self.k).low_i256();
    let reduced = ((c + temp * self.m) >> 256).low_i256();
    if reduced < self.m {
      Reduced(reduced)
    } else {
      Reduced(reduced - self.m)
    }
  }

  /// n must be non-negative
  /// a must be reduced
  pub fn exp_mod_(&self, a: Reduced<I256>, n: I256) -> Reduced<I256> {
    let mut a = a;
    let mut out = self.one_reduced;
    let mut n = n;
    loop {
      if n == I256::zero() {
        return out;
      }
      if n.is_odd() {
        out = self.mul_mod_(out, a);
      }
      a = self.mul_mod_(a, a);
      n >>= 1;
    }
  }

  /// n must be non-negative
  pub fn exp_mod(&self, a: I256, n: I256) -> I256 {
    self.unreduce(self.exp_mod_(self.reduce(a), n))
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_add() {
    assert!(i256(1) + i256(2) == i256(3));
  }

  #[test]
  fn test_add_big() {
    assert!(i256([0, 1, 0, 0]) + i256([0, 1, 0, 0]) == i256([0, 2, 0, 0]));
  }

  #[test]
  fn test_add_different_sizes() {
    assert!(i256([0, 1, 0, 0]) + i256([0, 1, 1, 1]) == i256([0, 2, 1, 1]));
  }

  #[test]
  fn test_mul() {
    assert!(i256(2) * i256(3) == i512(6));
  }

  #[test]
  fn test_mul_big() {
    assert!(i256([0, 1, 0, 0]) * i256([0, 1, 0, 0]) == i512([0, 0, 1, 0, 0, 0, 0, 0]));
  }

  #[test]
  fn test_mul_different_sizes() {
    assert!(i256([0, 2, 0, 0]) * i256([0, 1, 0, 1]) == i512([0, 0, 2, 0, 2, 0, 0, 0]));
  }
}
