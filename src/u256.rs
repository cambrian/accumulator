//! TODO: reduce U256/U512 duplication with a macro
use gmp_mpfr_sys::gmp;
use gmp_mpfr_sys::gmp::mpz_t;
use std::cmp::{min, Ord, Ordering, PartialOrd};
use std::convert::From;
use std::mem::transmute;
use std::ops;

macro_rules! u_types {
  ($($t:ident,$size:expr),+) => {
    $(
      #[derive(PartialEq, Eq, Hash, Debug, Clone, Copy)]
      pub struct $t {
        size: i64,
        limbs: [u64; $size],
      }

      impl $t {
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
            alloc: $size,
          }
        }
        pub fn zero() -> Self {
          Self {
            size: 0,
            limbs: [0; $size],
          }
        }
        pub fn one() -> Self {
          let mut limbs = [0; $size];
          limbs[0] = 1;
          Self { size: 1, limbs }
        }
        pub fn minus_one() -> Self {
          let mut limbs = [0; $size];
          limbs[0] = 1;
          Self { size: -1, limbs }
        }
        fn normalize_size(&mut self) {
          self.size = 0;
          for i in (0..$size).rev() {
            if self.limbs[i] != 0 {
              self.size = (i + 1) as i64;
              break;
            }
          }
        }
        /// Returns (result, remainder)
        pub fn div_rem(self, x: &Self) -> (Self, Self) {
          if x.size > self.size {
            return (Self::zero(), self);
          }
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
        pub fn div_rem_mut(&mut self, x: &Self) -> Self {
          if x.size > self.size {
            return Self::zero();
          }
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
      }

      impl From<[u64; $size]> for $t {
        fn from(limbs: [u64; $size]) -> Self {
          let mut x = Self { size: 0, limbs };
          x.normalize_size();
          x
        }
      }

      impl From<u64> for $t {
        fn from(x: u64) -> Self {
          let mut limbs = [0; $size];
          limbs[0] = x;
          Self::from(limbs)
        }
      }

      impl PartialOrd for $t {
        fn partial_cmp(&self, x: &Self) -> Option<Ordering> {
          let x = unsafe { gmp::mpn_cmp(self.ptr(), x.ptr(), $size) };
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

      impl Ord for $t {
        fn cmp(&self, x: &Self) -> Ordering {
          let x = unsafe { gmp::mpn_cmp(self.ptr(), x.ptr(), $size) };
          if x < 0 {
            Ordering::Less
          } else if x == 0 {
            Ordering::Equal
          } else {
            Ordering::Greater
          }
        }
      }

      impl ops::ShlAssign<u32> for $t {
        fn shl_assign(&mut self, mut x: u32) {
          while x != 0 {
            let sz = min(gmp::LIMB_BITS as u32, x);
            x -= sz;
            unsafe { gmp::mpn_lshift(self.ptr_mut(), self.ptr(), $size, sz) };
          }
          self.normalize_size();
        }
      }

      impl ops::Shl<u32> for $t {
        type Output = Self;
        fn shl(self, x: u32) -> Self {
          let mut y = self;
          y <<= x;
          y
        }
      }

      impl ops::ShrAssign<u32> for $t {
        fn shr_assign(&mut self, mut x: u32) {
          while x != 0 {
            let sz = min(gmp::LIMB_BITS as u32, x);
            x -= sz;
            unsafe { gmp::mpn_rshift(self.ptr_mut(), self.ptr(), $size, sz) };
          }
          self.normalize_size();
        }
      }

      impl ops::Shr<u32> for $t {
        type Output = Self;
        fn shr(self, x: u32) -> Self {
          let mut y = self;
          y >>= x;
          y
        }
      }

      impl ops::AddAssign for $t {
        /// panics if result overflows.
        fn add_assign(&mut self, x: Self) {
          let carry = unsafe { gmp::mpn_add_n(self.ptr_mut(), self.ptr(), x.ptr(), $size) };
          assert!(carry == 0);
          self.normalize_size();
        }
      }

      impl ops::Add for $t {
        /// panics if result overflows.
        type Output = Self;
        fn add(self, x: Self) -> Self {
          let mut y = self;
          y += x;
          y
        }
      }

      impl ops::SubAssign for $t {
        /// panics if result is negative.
        fn sub_assign(&mut self, x: Self) {
          let borrow = unsafe { gmp::mpn_sub_n(self.ptr_mut(), self.ptr(), x.ptr(), $size) };
          assert!(borrow == 0);
          self.normalize_size();
        }
      }

      impl ops::Sub for $t {
        type Output = Self;
        /// panics if result is negative.
        fn sub(self, x: Self) -> Self {
          let mut y = self;
          y -= x;
          y
        }
      }

      impl ops::Sub<u64> for $t {
        type Output = Self;
        /// panics if result is negative.
        fn sub(self, x: u64) -> Self {
          self - Self::from(x)
        }
      }

      impl ops::Div for $t {
        type Output = Self;
        fn div(self, x: Self) -> Self {
          self.div_rem(&x).0
        }
      }

      impl ops::RemAssign for $t {
        fn rem_assign(&mut self, x: Self) {
          self.div_rem_mut(&x);
        }
      }
    )+
  }
}

u_types!(U256, 4, U512, 8);

impl U512 {
  /// Returns the lower half of this U512 as a U256
  pub fn low_u256(self) -> U256 {
    let mut x = unsafe { transmute::<U512, (U256, [u64; 4])>(self) }.0;
    x.normalize_size();
    x
  }
}

// TODO: more efficient implementation via memcpy?
impl From<U256> for U512 {
  fn from(x: U256) -> Self {
    let mut limbs = [0; 8];
    limbs[..4].copy_from_slice(&x.limbs);
    U512 {
      size: x.size,
      limbs,
    }
  }
}

impl ops::Mul for U512 {
  type Output = U512;
  fn mul(self, x: Self) -> U512 {
    let mut y = U512::zero();
    if self.size >= x.size {
      unsafe { gmp::mpn_mul(y.ptr_mut(), self.ptr(), self.size, x.ptr(), x.size) };
    } else {
      unsafe { gmp::mpn_mul(y.ptr_mut(), x.ptr(), x.size, self.ptr(), self.size) };
    }
    y.normalize_size();
    y
  }
}

impl ops::Mul for &U512 {
  type Output = U512;
  fn mul(self, x: Self) -> U512 {
    let mut y = U512::zero();
    if self.size >= x.size {
      unsafe { gmp::mpn_mul(y.ptr_mut(), self.ptr(), self.size, x.ptr(), x.size) };
    } else {
      unsafe { gmp::mpn_mul(y.ptr_mut(), x.ptr(), x.size, self.ptr(), self.size) };
    }
    y.normalize_size();
    y
  }
}

impl ops::Rem for U512 {
  type Output = Self;
  fn rem(self, x: Self) -> Self {
    self.div_rem(&x).1
  }
}

// This gets its own implementation because it needs to be fast.
impl ops::Rem<&U256> for U512 {
  type Output = U256;
  fn rem(self, x: &U256) -> U256 {
    if x.size > self.size {
      return self.low_u256();
    }
    let (mut y, mut rem) = (U512::zero(), U256::zero());
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
    rem.normalize_size();
    rem
  }
}

impl ops::Rem<U256> for U512 {
  type Output = U256;
  fn rem(self, x: U256) -> U256 {
    #![allow(clippy::op_ref)]
    self % &x
  }
}

#[allow(unused_mut)]
fn mut_ptr<T>(mut t: &T) -> *mut T {
  t as *const T as *mut T
}

impl U256 {
  pub fn is_odd(&self) -> bool {
    self.limbs[0] & 1 == 1
  }

  /// returns (result of removing all fs, number of fs removed)
  pub fn remove_factor(self, f: Self) -> (Self, u64) {
    // for some reason this needs extra scratch space
    let mut out = U512::zero();
    let outmpz = out.as_mpz();
    let s = self.as_mpz();
    let f = f.as_mpz();
    let c = unsafe { gmp::mpz_remove(mut_ptr(&outmpz), mut_ptr(&s), mut_ptr(&f)) };
    out.size = i64::from(outmpz.size);
    (out.low_u256(), c)
  }

  pub fn mod_inv(self, m: Self) -> Self {
    let mut out = U512::zero();
    let outmpz = out.as_mpz();
    let s = self.as_mpz();
    let m = m.as_mpz();
    let exists = unsafe { gmp::mpz_invert(mut_ptr(&outmpz), mut_ptr(&s), mut_ptr(&m)) };
    assert!(exists != 0);
    out.size = i64::from(outmpz.size);
    out.low_u256()
  }

  pub fn pow_mod(self, e: Self, m: Self) -> Self {
    let mut out = U512::zero();
    let outmpz = out.as_mpz();
    let s = self.as_mpz();
    let e = e.as_mpz();
    let m = m.as_mpz();
    unsafe { gmp::mpz_powm(mut_ptr(&outmpz), mut_ptr(&s), mut_ptr(&e), mut_ptr(&m)) };
    out.size = i64::from(outmpz.size);
    out.low_u256()
  }

  pub fn is_perfect_square(&self) -> bool {
    let issqr = unsafe { gmp::mpn_perfect_square_p(self.ptr(), self.size) };
    issqr != 0
  }

  // TODO: this maybe doesn't belong here
  pub fn jacobi(a: i32, b: &Self) -> i32 {
    let mut a_data = 0;
    let a = i32_to_mpz(a, &mut a_data);
    let b = b.as_mpz();
    unsafe { gmp::mpz_jacobi(&a as *const mpz_t, &b as *const mpz_t) }
  }

  pub fn is_congruent(self, i: i32, m: &U256) -> bool {
    let mut data = 0;
    let x = i32_to_mpz(i, &mut data);
    let s = self.as_mpz();
    let m = m.as_mpz();
    let res = unsafe { gmp::mpz_congruent_p(mut_ptr(&s), mut_ptr(&x), mut_ptr(&m)) };
    res != 0
  }

  pub fn write_binary(&self, buf: &mut [u8; 257]) -> usize {
    unsafe { gmp::mpn_get_str(mut_ptr(&buf[0]), 2, self.ptr() as *mut u64, self.size) }
  }
}

/// Lower-endian bytes
impl From<[u8; 32]> for U256 {
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
impl From<&[u8; 32]> for U256 {
  fn from(bytes: &[u8; 32]) -> Self {
    let chunks = unsafe { transmute::<[u8; 32], [[u8; 8]; 4]>(*bytes) };
    let mut limbs = [0; 4];
    for i in 0..4 {
      limbs[i] = u64::from_le_bytes(chunks[i]);
    }
    Self::from(limbs)
  }
}

impl ops::Add<u64> for U256 {
  type Output = Self;
  /// panics if result overflows.
  fn add(self, x: u64) -> Self {
    self + U256::from(x)
  }
}

impl ops::Sub for &U256 {
  type Output = U256;
  /// panics if result is negative.
  fn sub(self, x: Self) -> U256 {
    let mut y = U256::zero();
    let borrow = unsafe { gmp::mpn_sub_n(y.ptr_mut(), self.ptr(), x.ptr(), 4) };
    assert!(borrow == 0);
    y.normalize_size();
    y
  }
}

impl ops::Sub<u64> for &U256 {
  type Output = U256;
  /// panics if result is negative.
  fn sub(self, x: u64) -> U256 {
    self - &U256::from(x)
  }
}

impl ops::Mul for U256 {
  type Output = U512;
  fn mul(self, x: Self) -> U512 {
    let mut y = U512::zero();
    unsafe { gmp::mpn_mul_n(y.ptr_mut(), self.ptr(), x.ptr(), 4) };
    y.normalize_size();
    y
  }
}

impl ops::Mul for &U256 {
  type Output = U512;
  fn mul(self, x: Self) -> U512 {
    let mut y = U512::zero();
    unsafe { gmp::mpn_mul_n(y.ptr_mut(), self.ptr(), x.ptr(), 4) };
    y.normalize_size();
    y
  }
}

impl ops::Rem<&Self> for U256 {
  type Output = Self;
  fn rem(self, x: &Self) -> Self {
    self.div_rem(x).1
  }
}

impl ops::Rem for U256 {
  type Output = Self;
  fn rem(self, x: Self) -> Self {
    #![allow(clippy::op_ref)]
    self % &x
  }
}

pub fn u256<T>(t: T) -> U256
where
  U256: From<T>,
{
  U256::from(t)
}

pub fn u512<T>(t: T) -> U512
where
  U512: From<T>,
{
  U512::from(t)
}

fn i32_to_mpz(i: i32, data: &mut u64) -> mpz_t {
  *data = i.abs() as u64;
  mpz_t {
    size: i.signum(),
    d: mut_ptr(&data),
    alloc: 1,
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_add() {
    assert!(u256(1) + u256(2) == u256(3));
  }

  #[test]
  fn test_add_big() {
    assert!(u256([0, 1, 0, 0]) + u256([0, 1, 0, 0]) == u256([0, 2, 0, 0]));
  }

  #[test]
  fn test_add_different_sizes() {
    assert!(u256([0, 1, 0, 0]) + u256([0, 1, 1, 1]) == u256([0, 2, 1, 1]));
  }

  #[test]
  fn test_mul() {
    assert!(u256(2) * u256(3) == u512(6));
  }

  #[test]
  fn test_mul_big() {
    assert!(u256([0, 1, 0, 0]) * u256([0, 1, 0, 0]) == u512([0, 0, 1, 0, 0, 0, 0, 0]));
  }

  #[test]
  fn test_mul_different_sizes() {
    assert!(u256([0, 2, 0, 0]) * u256([0, 1, 0, 1]) == u512([0, 0, 2, 0, 2, 0, 0, 0]));
  }
}
