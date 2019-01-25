// TODO: Add reading links.
use crate::util::{int, modpow_inplace};
use rug::Integer;
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
pub fn is_prob_prime(n: &Integer) -> bool {
  for &p in constants::SMALL_PRIMES.iter() {
    if n.is_congruent_u(0, p) {
      return *n == p;
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

pub fn passes_miller_rabin_base_2(n: &Integer) -> bool {
  // Write n-1 = 2^r * d.
  let mut d: Integer = n.clone() - 1;
  let mut r = 0;
  while d.is_congruent_u(0, 2) {
    d /= 2;
    r += 1;
  }
  // compute 2^d (mod n)
  let mut x = int(2);
  modpow_inplace(&mut x, &d, n);
  if x == 1 || x == n.clone() - 1 {
    return true;
  }
  for _ in 1..r {
    modpow_inplace(&mut x, &int(2), n);
    if x == 1 {
      return false;
    }
    if x == n.clone() - 1 {
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
fn choose_d(n: &Integer, max_iter: u64) -> Option<Integer> {
  let mut d = int(5);

  for _ in 0..max_iter {
    if d.jacobi(&n) == -1 {
      return Some(d);
    }
    if d < 0 {
      d = -d + 2;
    } else {
      d = -d - 2;
    }
  }
  None
}

// Strong Lucas probable prime test (NOT the more common Lucas primality test which requires
// factorization of n-1).
fn passes_lucas(n: &Integer, d: &Integer) -> bool {
  let p = int(1);
  let q = (1 - d.clone()) / 4;
  let delta = n.clone() + 1;

  // TODO: Extend to stronger test.
  let (u_delta, _v_delta) = compute_u_and_v_k(&delta, n, &int(1), &p, &p, &q, d);
  // u_delta % n != 0 proves n composite.
  u_delta == 0
}

// Computes the Lucas sequences {u_i(p, q)} and {v_i(p, q)} up to a specified index k in log(k)
// time by recursively calculating only the (2i)th and (2i+1)th elements in an order determined by
// the binary expansion of k. In the Lucas probabilistic primality test we specify that d = p^2 - 4q
// and are generally concerned with the case that k = delta = n+1.
// Cf. https://en.wikipedia.org/wiki/Lucas_pseudoprime
fn compute_u_and_v_k(
  k: &Integer,
  n: &Integer,
  u_1: &Integer,
  v_1: &Integer,
  p: &Integer,
  q: &Integer,
  d: &Integer,
) -> (Integer, Integer) {
  let k_bits = to_binary(k);
  let mut u_k = u_1.clone();
  let mut v_k = v_1.clone();
  let mut q_k = q.clone();

  let mod_n = |x: &Integer| {
    if *x < 0 {
      n - (-x.clone() % n)
    } else {
      x.clone() % n
    }
  };

  let half = |x: &Integer| {
    if x.is_congruent_u(1, 2) {
      mod_n(&((x.clone() + n) / 2))
    } else {
      mod_n(&(x.clone() / 2))
    }
  };

  // Write binary expansion of k as [x_1, ..., x_l], e.g. [1, 0, 1, 1] for 11. x_1 is always
  // 1. For i = 2, 3, ..., l, do the following: if x_i = 0 then update u_k and v_k to u_{2k} and
  // v_{2k}, respectively. Else if x_i = 1, update to u_{2k+1} and v_{2k+1}. At the end of the loop
  // we will have computed u_k and v_k, with k as given, in log(delta) time.
  let mut k = int(1);
  for bit in k_bits[1..].chars() {
    // Compute (u, v)_{2k} from (u, v)_k.
    // u_2k = u_k * v_k % n
    // v_2k = v_k^2 - 2*q^k
    u_k = mod_n(&(u_k.clone() * v_k.clone()));
    modpow_inplace(&mut v_k, &int(2), n);
    v_k = mod_n(&(v_k - 2 * q_k.clone()));
    q_k = mod_n(&(q_k.clone() * q_k.clone()));
    k *= 2;
    if bit == '1' {
      // Compute (u, v)_{2k+1} from (u, v)_{2k}.
      // TODO: Why is mod_n necessary here?
      let pu_plus_v = p * u_k.clone() + v_k.clone();
      let du_plus_pv = mod_n(&(d * u_k.clone() + p * v_k.clone()));
      u_k = half(&pu_plus_v);
      v_k = half(&du_plus_pv);
      q_k = mod_n(&(q_k.clone() * q.clone()));
      k += 1;
    }
  }
  (u_k, v_k)
}

fn to_binary(n: &Integer) -> String {
  format!("{:b}", n)
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_miller_rabin() {
    assert!(passes_miller_rabin_base_2(&int(13)));
    assert!(!passes_miller_rabin_base_2(&int(65)));
    for &p in constants::LARGE_PRIMES.iter() {
      assert!(passes_miller_rabin_base_2(&int(p)));
      assert!(!passes_miller_rabin_base_2(&(int(p) * int(106957))));
    }
  }

  #[test]
  fn test_jacobi() {
    assert_eq!(int(0).jacobi(&int(1)), 1);
    assert_eq!(int(15).jacobi(&int(17)), 1);
    assert_eq!(int(14).jacobi(&int(17)), -1);
    assert_eq!(int(30).jacobi(&int(59)), -1);
    assert_eq!(int(27).jacobi(&int(57)), 0);
  }

  #[test]
  fn test_to_binary() {
    assert_eq!(to_binary(&int(1)), "1");
    assert_eq!(to_binary(&int(44)), "101100");
  }

  #[test]
  fn test_is_prob_prime() {
    // Sanity checks.
    assert!(is_prob_prime(&int(2)));
    assert!(is_prob_prime(&int(5)));
    assert!(is_prob_prime(&int(7)));
    assert!(is_prob_prime(&int(241)));
    assert!(is_prob_prime(&int(7919)));
    assert!(is_prob_prime(&int(48131)));
    assert!(is_prob_prime(&int(75913)));
    assert!(is_prob_prime(&int(76463)));
    assert!(is_prob_prime(&int(115_547)));

    // Medium primes.
    for &p in constants::MED_PRIMES.iter() {
      assert!(is_prob_prime(&int(p)));
    }

    // Large primes.
    for &p in constants::LARGE_PRIMES.iter() {
      assert!(is_prob_prime(&int(p)));
    }

    // Large, difficult-to-factor composites.
    for &p in constants::LARGE_PRIMES.iter() {
      for &q in constants::LARGE_PRIMES.iter() {
        assert!(!is_prob_prime(&(int(p) * int(q))));
      }
    }
  }
}
