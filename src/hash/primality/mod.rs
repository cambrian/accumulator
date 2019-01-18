//use num::bigint::BigInt;
use num::bigint::BigUint;
use num_traits::cast::ToPrimitive;

mod utils;

// Baillie-PSW probabilistic primality test:
// 1. Filter composites with small divisors.
// 2. Do Miller-Rabin base 2.
// 3. Filter squares.
// 4. Do Lucas.
#[allow(dead_code)]
pub fn is_prob_prime(n: &BigUint) -> bool {
  !has_small_prime_factor(n) && passes_miller_rabin_base_2(n) && !is_square(n) && passes_lucas(n)
}

#[allow(dead_code)]
fn has_small_prime_factor(n: &BigUint) -> bool {
  for &divisor in utils::SMALL_PRIMES.iter() {
    let divisor = &BigUint::from(divisor);
    if divisor > n {
      break;
    }
    if n % divisor == BigUint::from(0u64) {
      return true;
    }
  }
  false
}

#[allow(dead_code)]
fn passes_miller_rabin_base_2(n: &BigUint) -> bool {
  let one = BigUint::from(1u64);
  let two = BigUint::from(2u64);

  // write n-1 = 2^r * d
  let mut d = n - 1u64;
  let mut r = 0;
  while &d % 2u64 == BigUint::from(0u64) {
    d /= 2u64;
    r += 1;
  }
  // println!("{} = 2^{} * {}", n, r, d);
  let mut x = two.modpow(&d, n);
  if x == one || x == n - &one {
    return true;
  }
  for _ in 0..(r - 1) {
    x = x.modpow(&two, n);
    if x == one {
      return false;
    }
    if x == n - &one {
      return true;
    }
  }
  false
}

fn is_square(n: &BigUint) -> bool {
  // Step 1
  let zero = BigUint::from(0 as u8);
  let one = BigUint::from(1 as u8);
  if n & BigUint::from(2 as u8) != zero
    || n & BigUint::from(7 as u8) == BigUint::from(5 as u8)
    || n & BigUint::from(11 as u8) == BigUint::from(8 as u8)
  {
    return false;
  }
  // Maybe unneccessary
  if *n == zero {
    return true;
  }

  println!("Step 2");

  // Step 2
  let copy = n.clone();
  let copy = (copy.clone() & BigUint::from(4_294_967_295 as u64)) + (copy >> 32);
  let copy = (copy.clone() & BigUint::from(65535 as u16)) + (copy >> 16);
  let copy = (copy.clone() & BigUint::from(255 as u8))
    + ((copy.clone() >> 8) & BigUint::from(255 as u8))
    + (copy >> 16);
  // println!("{}", n.to_u64().unwrap());
  if utils::BAD_255[copy.to_u64().unwrap() as usize] {
    return false;
  }

  println!("Step 3");

  let mut x = n.clone();
  if x.clone() & BigUint::from(4_294_967_295 as u64) == zero {
    x >>= 32;
  }
  if x.clone() & BigUint::from(65535 as u64) == zero {
    x >>= 16;
  }
  if x.clone() & BigUint::from(255 as u16) == zero {
    x >>= 8;
  }
  if x.clone() & BigUint::from(15 as u8) == zero {
    x >>= 4;
  }
  if x.clone() & BigUint::from(3 as u8) == zero {
    x >>= 2;
  }
  if x.clone() & BigUint::from(7 as u8) != one {
    return false;
  }

  println!("Step 4");

  // let mut r: i64 = start[((n >> 3) & BigUint::from(1023 as u64)).to_u64().unwrap() as usize];
  // let mut t: BigInt;
  // let mut z: BigInt;
  // let zero_i = BigInt::from(0 as i8);
  // while {
  //   z = BigInt::from(x.clone()) - BigInt::from(r * r);
  //   if z == zero_i {
  //     return true;
  //   }
  //   t = z.clone() & -z.clone();
  //   r += ((z & t.clone()) >> 1).to_i64().unwrap();
  //   if r > (t.clone() >> 1).to_i64().unwrap() {
  //     r = t.to_i64().unwrap() - r;
  //   }
  //   t <= BigInt::from(1 << 33)
  // } {}
  // println!("All else fails");

  //0xC840C04048404040
  // let inbase16 = &[12, 8, 4, 0, 12, 0, 4, 0, 4, 8, 4, 0, 4, 0, 4, 0];
  // let good_mask = BigUint::from_radix_be(inbase16, 16).unwrap();
  // if good_mask << n >= zero {
  //   return false;
  // }
  true
}

#[allow(dead_code)]
fn passes_lucas(_n: &BigUint) -> bool {
  false
}

#[cfg(test)]
mod tests {
  use super::*;

  fn isqrt(n: usize) -> usize {
    n == 0 && return n;
    let mut s = (n as f64).sqrt() as usize;
    s = (s + n / s) >> 1;
    if s * s > n {
      s - 1
    } else {
      s
    }
  }

  fn perfect_sqrt(n: usize) -> isize {
    match n & 0xf {
      0 | 1 | 4 | 9 => {
        let t = isqrt(n);
        if t * t == n {
          t as isize
        } else {
          -1
        }
      }
      _ => -1,
    }
  }
  #[test]
  fn test_is_square() {
    let mut squares = vec![];
    for i in 0..200 {
      if perfect_sqrt(i) >= 0 {
        squares.push(true);
      } else {
        squares.push(false);
      }
    }
    for i in 0..squares.len() {
      println!("{}", i);
      let val = BigUint::from(i as u8);
      assert!(is_square(&val) == squares[i]);
    }
  }
  #[test]
  fn test_small_prime_factor() {
    let n_prime = BigUint::from(233u64);
    let n_composite = BigUint::from(50_621u64);
    let n_composite_large = BigUint::from(104_927u64);

    assert!(n_composite == BigUint::from(223u64) * BigUint::from(227u64));
    assert!(n_composite_large == BigUint::from(317u64) * BigUint::from(331u64));

    assert!(!has_small_prime_factor(&n_prime));
    assert!(has_small_prime_factor(&n_composite));
    assert!(!has_small_prime_factor(&n_composite_large));
  }

  #[test]
  fn test_miller_rabin() {
    assert!(passes_miller_rabin_base_2(&BigUint::from(13u64)));
    assert!(!passes_miller_rabin_base_2(&BigUint::from(65u64)));
  }
}
