use num::bigint::{BigInt, BigUint, Sign, ToBigInt, ToBigUint};
use num_traits::cast::ToPrimitive;

mod utils;

macro_rules! bu {
  ($x:expr) => {
    //BigUint::from($x as u64)
    ToBigUint::to_biguint(&($x as u64)).unwrap()
  };
}

macro_rules! bi {
  ($x:expr) => {
    //BigInt::from($x as i64)
    //($x as i64).to_bigint()
    ToBigInt::to_bigint(&($x as i64)).unwrap()
  };
}

const MAX_JACOBI_ITERS: u64 = 500;

// Baillie-PSW probabilistic primality test:
// 1. Filter composites with small divisors.
// 2. Do Miller-Rabin base 2.
// 3. Filter squares.
// 4. Do Lucas.
#[allow(dead_code)]
pub fn is_prob_prime(n: &BigUint) -> bool {
  for &p in utils::SMALL_PRIMES.iter() {
    if n == &bu!(p) {
      return true;
    }
  }

  if has_small_prime_factor(n) {
    return false;
  }
  //println!("no small prime factors...");
  if !passes_miller_rabin_base_2(n) {
    return false;
  }
  //println!("passes Miller-Rabin base 2...");
  if is_prob_square(n) {
    return false;
  }
  //println!("probably not a square...");
  let n = &BigInt::from_biguint(Sign::Plus, n.clone());
  match choose_d(n, MAX_JACOBI_ITERS) {
    Some(d) => passes_lucas(n, &d),
    None => false,
  }
}

#[allow(dead_code)]
fn has_small_prime_factor(n: &BigUint) -> bool {
  for &divisor in utils::SMALL_PRIMES.iter() {
    let divisor = &bu!(divisor);
    if divisor == n {
      break;
    }
    if n % divisor == bu!(0) {
      return true;
    }
  }
  false
}

#[allow(dead_code)]
fn passes_miller_rabin_base_2(n: &BigUint) -> bool {
  // write n-1 = 2^r * d
  let mut d = n - 1u64;
  let mut r = 0;
  while &d % 2u64 == bu!(0) {
    d /= 2u64;
    r += 1;
  }
  // println!("{} = 2^{} * {}", n, r, d);
  let mut x = bu!(2).modpow(&d, n);
  if x == bu!(1) || x == n - &bu!(1) {
    return true;
  }
  for _ in 0..(r - 1) {
    x = x.modpow(&bu!(2), n);
    if x == bu!(1) {
      return false;
    }
    if x == n - &bu!(1) {
      return true;
    }
  }
  false
}

fn is_prob_square(n: &BigUint) -> bool {
  // Step 1
  let zero = bu!(0);
  let one = bu!(1);
  if n & bu!(2) != zero || n & bu!(7) == bu!(5) || n & bu!(11) == bu!(8) {
    return false;
  }
  // Maybe unneccessary
  if *n == zero {
    return true;
  }

  // println!("Step 2");

  // Step 2
  let copy = n.clone();
  let copy = (copy.clone() & bu!(4_294_967_295)) + (copy >> 32);
  let copy = (copy.clone() & bu!(65535)) + (copy >> 16);
  let copy = (copy.clone() & bu!(255)) + ((copy.clone() >> 8) & bu!(255)) + (copy >> 16);
  // println!("{}", n.to_u64().unwrap());
  if utils::BAD_255[copy.to_u64().unwrap() as usize] {
    return false;
  }

  // println!("Step 3");

  let mut x = n.clone();
  if x.clone() & bu!(4_294_967_295) == zero {
    x >>= 32;
  }
  if x.clone() & bu!(65535) == zero {
    x >>= 16;
  }
  if x.clone() & bu!(255) == zero {
    x >>= 8;
  }
  if x.clone() & bu!(15) == zero {
    x >>= 4;
  }
  if x.clone() & bu!(3) == zero {
    x >>= 2;
  }
  if x.clone() & bu!(7) != one {
    return false;
  }

  // println!("Step 4");

  // let mut r: i64 = start[((n >> 3) & bu!(1023 as u64)).to_u64().unwrap() as usize];
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

// Finds and returns first D in [5, -7, 9, ..., 5 + 2*max_iter] for which Jacobi symbol (D/n) = -1,
// or None if no such D exists. In the case that n is square, there is no such D even with max_iter
// infinite. Hence if you are not precisely sure that n is nonsquare, you should pass a low value
// to max_iter to avoid wasting too much time. Note that the average number of iterations required
// for nonsquare n is 1.8.
fn choose_d(n: &BigInt, max_iter: u64) -> Option<BigInt> {
  let mut d = bi!(5);

  for _ in 0..max_iter {
    if jacobi_symbol(&d, n) == -1 {
      return Some(d);
    }
    if d > bi!(0) {
      d += bi!(2);
    } else {
      d -= bi!(2);
    }
    d *= bi!(-1);
  }
  None
}

// Jacobi symbol (a/n)
fn jacobi_symbol(a: &BigInt, n: &BigInt) -> i64 {
  // base cases
  if n == &bi!(1) {
    return 1;
  }
  if a == &bi!(0) {
    return 0;
  } else if a == &bi!(1) {
    return 1;
  } else if a == &bi!(2) {
    let n_mod_8 = n % 8;
    if n_mod_8 == bi!(3) || n_mod_8 == bi!(5) {
      return -1;
    } else if n_mod_8 == bi!(1) || n_mod_8 == bi!(7) {
      return 1;
    }
  }
  // recursive cases
  if *a < bi!(0) {
    // symbol is (-1)^((n-1)/2) (-a/n)
    let j = jacobi_symbol(&(a * &bi!(-1)), n);
    let exp_mod_2 = ((n - bi!(1)) / bi!(2)) % 2;
    if exp_mod_2 == bi!(0) {
      return j;
    } else {
      return -j;
    }
  }
  if a % 2 == bi!(0) {
    jacobi_symbol(&bi!(2), n) * jacobi_symbol(&(a / &bi!(2)), n)
  } else if a % n != *a {
    jacobi_symbol(&(a % n), n)
  } else if a % 4 == bi!(3) && n % 4 == bi!(3) {
    -jacobi_symbol(n, a)
  } else {
    jacobi_symbol(n, a)
  }
}

#[allow(dead_code)]
// strong Lucas probable prime test (NOT the more common Lucas primality test which requires
// factorization of n-1).
fn passes_lucas(n: &BigInt, d: &BigInt) -> bool {
  let p = bi!(1);
  let q = (bi!(1) - d) / bi!(4);
  let delta = n + &bi!(1);

  // println!("Lucas test: (n, d, p, q) = ({}, {}, {}, {})", n, d, p, q);

  let (u_delta, _v_delta) = compute_u_and_v_k(&delta, n, &bi!(1), &p, &p, &q, d);
  // u_delta % n != 0 proves n composite
  u_delta == bi!(0)
  // EXTEND TO STRONG TEST
}

// Computes the Lucas sequences {u_i(p, q)} and {v_i(p, q)} up to a specified index k in log(k)
// time by recursively calculating only the (2i)th and (2i+1)th elements in an order determined by
// the binary expansion of k. In the Lucas probabilistic primality test we specify that d = p^2 - 4q
// and are generally concerned with the case that k = delta = n+1.
// Cf. https://en.wikipedia.org/wiki/Lucas_pseudoprime
fn compute_u_and_v_k(
  k: &BigInt,
  n: &BigInt,
  u_1: &BigInt,
  v_1: &BigInt,
  p: &BigInt,
  q: &BigInt,
  d: &BigInt,
) -> (BigInt, BigInt) {
  let k_bits = binary_rep(k);
  let mut u_k = u_1.clone();
  let mut v_k = v_1.clone();
  let mut q_k = q.clone();

  let mod_n = |x: &BigInt| {
    if *x < bi!(0) {
      n - (-x % n)
    } else {
      x % n
    }
  };

  let half = |x: &BigInt| {
    if x % 2 == bi!(1) {
      mod_n(&((x + n) / 2))
    } else {
      mod_n(&(x / 2))
    }
  };

  // Write binary expansion of k as [x_1, ..., x_l], e.g. [1, 0, 1, 1] for 11. x_1 is always
  // 1. For i = 2, 3, ..., l, do the following: if x_i = 0 then update u_k and v_k to u_{2k} and
  // v_{2k}, respectively. Else if x_i = 1, update to u_{2k+1} and v_{2k+1}. At the end of the loop
  // we will have computed u_k and v_k, with k as given, in log(delta) time.
  let mut k = bi!(1);
  for &bit in k_bits[1..].iter() {
    // compute (u, v)_{2k} from (u, v)_k
    u_k = mod_n(&(u_k.clone() * v_k.clone()));
    v_k = mod_n(&(v_k.modpow(&bi!(2), n) - bi!(2) * &q_k));
    q_k = mod_n(&(&q_k * &q_k));
    k *= bi!(2);
    if bit == 1 {
      // compute (u, v)_{2k+1} from (u, v)_{2k}
      let pu_plus_v = mod_n(&(p * u_k.clone() + v_k.clone()));
      let du_plus_pv = mod_n(&(d * u_k.clone() + p * v_k.clone()));
      // if &pu_plus_v % 2 == bi!(0) {
      //   if &du_plus_pv % 2 == bi!(0) {
      //     u_k = half(&pu_plus_v);
      //     v_k = half(&du_plus_pv);
      //   } else {
      //     u_k = half(&pu_plus_v);
      //     v_k = half(&(du_plus_pv + n));
      //   }
      // } else if &du_plus_pv % 2 == bi!(0) {
      //   u_k = half(&(pu_plus_v + n));
      //   v_k = half(&du_plus_pv);
      // } else {
      //   u_k = half(&(pu_plus_v + n));
      //   v_k = half(&(du_plus_pv + n));
      // }
      k += bi!(1);
      // u_k = (u_k % n + n) % n;
      // v_k = (v_k % n + n) % n;
      u_k = half(&pu_plus_v);
      v_k = half(&du_plus_pv);
      q_k = mod_n(&(q_k.clone() * q.clone()));
      assert!(u_k >= bi!(0) && v_k >= bi!(0));
    }
    // println!("(u, v)_{} = ({}, {})", k, u_k, v_k);
  }
  (u_k, v_k)
}

fn binary_rep(n: &BigInt) -> Vec<u8> {
  let mut bits_rev = Vec::new();
  let mut n_mod_n = n.clone();
  while n_mod_n > bi!(0) {
    bits_rev.push((n_mod_n.clone() % 2 == bi!(1)) as u8);
    n_mod_n /= 2;
  }
  bits_rev.reverse();
  bits_rev
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
      let val = bu!(i);
      assert!(is_prob_square(&val) == squares[i]);
    }
  }
  #[test]
  fn test_small_prime_factor() {
    let n_prime = bu!(233u64);
    let n_composite = bu!(50_621u64);
    let n_composite_large = bu!(104_927u64);

    assert!(n_composite == bu!(223u64) * bu!(227u64));
    assert!(n_composite_large == bu!(317u64) * bu!(331u64));

    assert!(!has_small_prime_factor(&n_prime));
    assert!(has_small_prime_factor(&n_composite));
    assert!(!has_small_prime_factor(&n_composite_large));
  }

  #[test]
  fn test_miller_rabin() {
    assert!(passes_miller_rabin_base_2(&bu!(13u64)));
    assert!(!passes_miller_rabin_base_2(&bu!(65u64)));
    for &p in utils::LARGE_PRIMES.iter() {
      assert!(passes_miller_rabin_base_2(&bu!(p)));
      assert!(!passes_miller_rabin_base_2(&(bu!(p) * bu!(7))));
    }
  }

  #[test]
  fn test_jacobi() {
    assert_eq!(jacobi_symbol(&bi!(0), &bi!(1)), 1);
    assert_eq!(jacobi_symbol(&bi!(15), &bi!(17)), 1);
    assert_eq!(jacobi_symbol(&bi!(14), &bi!(17)), -1);
    assert_eq!(jacobi_symbol(&bi!(30), &bi!(59)), -1);
    assert_eq!(jacobi_symbol(&bi!(27), &bi!(57)), 0);
  }

  #[test]
  fn test_binary_rep() {
    assert_eq!(binary_rep(&bi!(1)), [1]);
    assert_eq!(binary_rep(&bi!(44)), [1, 0, 1, 1, 0, 0]);
  }

  #[test]
  fn test_lucas() {
    // fair inputs are odd, nonsquare, and have no small prime factors
    for &n in utils::SMALL_PRIMES.iter() {
      match choose_d(&bi!(n), MAX_JACOBI_ITERS) {
        Some(d) => assert!(passes_lucas(&bi!(n), &d)),
        None => assert!(false),
      }
    }

    for &n in utils::LARGE_PRIMES.iter() {
      match choose_d(&bi!(n), MAX_JACOBI_ITERS) {
        Some(d) => assert!(passes_lucas(&bi!(n), &d)),
        None => assert!(false),
      }
    }
    // for &p in utils::SMALL_PRIMES.iter() {
    //   assert!(is_prob_prime(&bu!(p)));
    // }
    // for &p in utils::LARGE_PRIMES.iter() {
    //   assert!(is_prob_prime(&bu!(p)));
    // }
  }

  #[test]
  fn test_is_prob_prime() {
    // assert!(is_prob_prime(&bu!(2)));
    // assert!(is_prob_prime(&bu!(5)));
    // assert!(is_prob_prime(&bu!(7)));
    // assert!(is_prob_prime(&bu!(241)));
    // assert!(is_prob_prime(&bu!(7919)));
    // assert!(is_prob_prime(&bu!(48131)));
    // assert!(is_prob_prime(&bu!(75913)));
    // assert!(is_prob_prime(&bu!(76463)));
    // assert!(is_prob_prime(&bu!(115_547)));
    // for &p in utils::LARGE_PRIMES.iter() {
    //   for &q in utils::LARGE_PRIMES.iter() {
    //     assert!(!is_prob_prime(&(bu!(p) * bu!(q))));
    //   }
    // }
    let mut tps = 0f64;
    for &p in utils::MED_PRIMES.iter() {
      let d = choose_d(&bi!(p), 500);
      if is_prob_prime(&bu!(p)) {
        tps += 1f64;
        // println!("{}: success", d.unwrap().to_str_radix(10));
      } else {
        println!("{}, {}", d.unwrap().to_str_radix(10), p);
      }
    }
    println!("Small prime true positive rate: {}", tps / 456f64);
    // tested on ~10k 32-bit composites with no failures
  }

  #[test]
  fn jkdsaflkaflds() {
    assert!(is_prob_prime(&bu!(48131)));
    assert!(is_prob_prime(&bu!(106957)));
  }
}
