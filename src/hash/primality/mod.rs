// TODO: Add reading links.
use crate::util::bi;
use num::bigint::{BigInt, BigUint, ToBigInt};
mod constants;

const MAX_JACOBI_ITERS: u64 = 500;

// Baillie-PSW probabilistic primality test:
// 1. Filter composites with small divisors.
// 2. Do Miller-Rabin base 2.
// 3. Filter squares.
// 4. Do Lucas.
// Many BPSW implementations filter out square values of x after the base-2 Miller-Rabin test. We
// do not do this for several reasons: (1) See comments above choose_d. (2) num::bigint does not
// have a native sqrt function using e.g. Newton's method. (3) Since this implementation is geared
// toward very large (256-bit and up) values of n, for which squares are far sparser than primes,
// the expected marginal utility of catching squares before running out MAX_JACOBI_ITERS is
// extremely low.
pub fn is_prob_prime(n: &BigUint) -> bool {
  let n = n.to_bigint().expect("Could not convert BigUint to BigInt!");
  for &p in constants::SMALL_PRIMES.iter() {
    if &n % p == bi(0) {
      return n == bi(p);
    }
  }
  if !passes_miller_rabin_base_2(&n) {
    return false;
  }
  match choose_d(&n, MAX_JACOBI_ITERS) {
    Some(d) => passes_lucas(&n, &d),
    None => false,
  }
}

fn passes_miller_rabin_base_2(n: &BigInt) -> bool {
  // Write n-1 = 2^r * d.
  let mut d = n - 1;
  let mut r = 0;
  while &d % 2 == bi(0) {
    d /= 2;
    r += 1;
  }
  let mut x = bi(2).modpow(&d, n);
  if x == bi(1) || x == n - &bi(1) {
    return true;
  }
  for _ in 0..(r - 1) {
    x = x.modpow(&bi(2), n);
    if x == bi(1) {
      return false;
    }
    if x == n - &bi(1) {
      return true;
    }
  }
  false
}

// Finds and returns first D in [5, -7, 9, ..., 5 + 2 * max_iter] for which Jacobi symbol (D/n) =
// -1, or None if no such D exists. In the case that n is square, there is no such D even with
// max_iter infinite. Hence if you are not precisely sure that n is nonsquare, you should pass a
// low value to max_iter to avoid wasting too much time. Note that the average number of iterations
// required for nonsquare n is 1.8, and empirically we find it is extremely rare that |d| > 13.
fn choose_d(n: &BigInt, max_iter: u64) -> Option<BigInt> {
  let mut d = bi(5);

  for _ in 0..max_iter {
    if jacobi_symbol(&d, n) == -1 {
      return Some(d);
    }
    if d > bi(0) {
      d += bi(2);
    } else {
      d -= bi(2);
    }
    d *= bi(-1);
  }

  None
}

/// Jacobi symbol (a / n).
fn jacobi_symbol(a: &BigInt, n: &BigInt) -> i64 {
  // Base cases.
  if n == &bi(1) {
    return 1;
  }
  if a == &bi(0) {
    return 0;
  } else if a == &bi(1) {
    return 1;
  } else if a == &bi(2) {
    let n_mod_8 = n % 8;
    if n_mod_8 == bi(3) || n_mod_8 == bi(5) {
      return -1;
    } else if n_mod_8 == bi(1) || n_mod_8 == bi(7) {
      return 1;
    }
  }

  // Recursive cases.
  if *a < bi(0) {
    // Symbol is (-1)^((n-1)/2) (-a/n).
    let j = jacobi_symbol(&(a * &bi(-1)), n);
    let exp_mod_2 = ((n - bi(1)) / bi(2)) % 2;
    if exp_mod_2 == bi(0) {
      return j;
    } else {
      return -j;
    }
  }
  if a % 2 == bi(0) {
    jacobi_symbol(&bi(2), n) * jacobi_symbol(&(a / &bi(2)), n)
  } else if a % n != *a {
    jacobi_symbol(&(a % n), n)
  } else if a % 4 == bi(3) && n % 4 == bi(3) {
    -jacobi_symbol(n, a)
  } else {
    jacobi_symbol(n, a)
  }
}

// Strong Lucas probable prime test (NOT the more common Lucas primality test which requires
// factorization of n-1).
fn passes_lucas(n: &BigInt, d: &BigInt) -> bool {
  let p = bi(1);
  let q = (bi(1) - d) / bi(4);
  let delta = n + &bi(1);

  // TODO: Extend to stronger test.
  let (u_delta, _v_delta) = compute_u_and_v_k(&delta, n, &bi(1), &p, &p, &q, d);
  u_delta == bi(0) // u_delta % n != 0 proves n composite.
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
  let k_bits = to_binary(k);
  let mut u_k = u_1.clone();
  let mut v_k = v_1.clone();
  let mut q_k = q.clone();

  let mod_n = |x: &BigInt| {
    if *x < bi(0) {
      n - (-x % n)
    } else {
      x % n
    }
  };

  let half = |x: &BigInt| {
    if x % 2 == bi(1) {
      mod_n(&((x + n) / 2))
    } else {
      mod_n(&(x / 2))
    }
  };

  // Write binary expansion of k as [x_1, ..., x_l], e.g. [1, 0, 1, 1] for 11. x_1 is always
  // 1. For i = 2, 3, ..., l, do the following: if x_i = 0 then update u_k and v_k to u_{2k} and
  // v_{2k}, respectively. Else if x_i = 1, update to u_{2k+1} and v_{2k+1}. At the end of the loop
  // we will have computed u_k and v_k, with k as given, in log(delta) time.
  let mut k = bi(1);
  for bit in k_bits[1..].chars() {
    // Compute (u, v)_{2k} from (u, v)_k.
    u_k = mod_n(&(u_k.clone() * v_k.clone()));
    v_k = mod_n(&(v_k.modpow(&bi(2), n) - bi(2) * &q_k));
    q_k = mod_n(&(q_k.clone() * q_k.clone()));
    k *= bi(2);
    if bit == '1' {
      // Compute (u, v)_{2k+1} from (u, v)_{2k}.
      // TODO: Why is mod_n necessary here?
      let pu_plus_v = p * u_k.clone() + v_k.clone();
      let du_plus_pv = mod_n(&(d * u_k.clone() + p * v_k.clone()));
      u_k = half(&pu_plus_v);
      v_k = half(&du_plus_pv);
      q_k = mod_n(&(q_k.clone() * q.clone()));
      k += bi(1);
    }
  }

  (u_k, v_k)
}

fn to_binary(n: &BigInt) -> String {
  format!("{:b}", n)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::util::bu;

  #[test]
  fn test_miller_rabin() {
    assert!(passes_miller_rabin_base_2(&bi(13)));
    assert!(!passes_miller_rabin_base_2(&bi(65)));
    for &p in constants::LARGE_PRIMES.iter() {
      assert!(passes_miller_rabin_base_2(&bi(p)));
      assert!(!passes_miller_rabin_base_2(&(bi(p) * bi(7))));
    }
  }

  #[test]
  fn test_jacobi() {
    assert_eq!(jacobi_symbol(&bi(0), &bi(1)), 1);
    assert_eq!(jacobi_symbol(&bi(15), &bi(17)), 1);
    assert_eq!(jacobi_symbol(&bi(14), &bi(17)), -1);
    assert_eq!(jacobi_symbol(&bi(30), &bi(59)), -1);
    assert_eq!(jacobi_symbol(&bi(27), &bi(57)), 0);
  }

  #[test]
  fn test_to_binary() {
    assert_eq!(to_binary(&bi(1)), "1");
    assert_eq!(to_binary(&bi(44)), "101100");
  }

  #[test]
  fn test_is_prob_prime() {
    // Sanity checks.
    assert!(is_prob_prime(&bu(2u64)));
    assert!(is_prob_prime(&bu(5u64)));
    assert!(is_prob_prime(&bu(7u64)));
    assert!(is_prob_prime(&bu(241u64)));
    assert!(is_prob_prime(&bu(7919u64)));
    assert!(is_prob_prime(&bu(48131u64)));
    assert!(is_prob_prime(&bu(75913u64)));
    assert!(is_prob_prime(&bu(76463u64)));
    assert!(is_prob_prime(&bu(115_547u64)));

    // Medium primes.
    for &p in constants::MED_PRIMES.iter() {
      assert!(is_prob_prime(&bu(p as u64)));
    }

    // Large primes.
    for &p in constants::LARGE_PRIMES.iter() {
      assert!(is_prob_prime(&bu(p as u64)));
    }

    // Large, difficult-to-factor composites.
    for &p in constants::LARGE_PRIMES.iter() {
      for &q in constants::LARGE_PRIMES.iter() {
        assert!(!is_prob_prime(&(bu(p as u64) * bu(q as u64))));
      }
    }
  }
}
