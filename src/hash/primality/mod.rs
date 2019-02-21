use crate::util::int;
use rug::integer::Order;
use rug::{Assign, Integer};

mod constants;
use constants::{D_VALUES, SMALL_PRIMES};

/// Implements the Baillie-PSW probabilistic primality test, which is known to be deterministic over
/// all integers up to 64 bits (u64). Offers more bang for your buck than Miller-Rabin (i.e.
/// iterated Fermat tests of random base) at wide n since Fermat and Lucas pseudoprimes have been
/// shown to be anticorrelated. Steps of BPSW are as follows:
/// 1. Accept small primes and reject multiples of them.
/// 2. Do a single iteration of Miller-Rabin (base-2 Fermat test).
/// 4. Do a strong probabilistic Lucas test (squares filtered during test initialization).
pub fn is_prob_prime(n: &Integer) -> bool {
  for &p in SMALL_PRIMES.iter() {
    if n.is_divisible_u(p) {
      return *n == p;
    }
  }
  passes_miller_rabin_base_2(&n) && passes_lucas(&n)
}

pub fn passes_miller_rabin_base_2(n: &Integer) -> bool {
  let (d, r) = int(n - 1).remove_factor(&int(2));
  let mut x = int(2);
  x.pow_mod_mut(&d, n).unwrap();
  if x == 1 || x == int(n - 1) {
    return true;
  }
  for _ in 1..r {
    x *= x.clone();
    x %= n;
    if x == 1 {
      return false;
    }
    if x == int(n - 1) {
      return true;
    }
  }
  false
}

/// Strong Lucas probable prime test (NOT the more common Lucas primality test which requires
/// factorization of `n-1`). Selects parameters `d`, `p`, `q` according to Selfridge's method.
/// Cf. [Lucas pseudoprime](https://en.wikipedia.org/wiki/Lucas_pseudoprime) on Wikipedia
/// If `n` passes, it is either prime or a "strong" Lucas pseudoprime. (The precise meaning of
/// "strong" is not fixed in the literature.) Procedure can be further strengthened by implementing
/// more tests in section 6 of [Baillie & Wagstaff 1980], but for now this is TODO.
/// Filters perfect squares as part of `choose_d`.
fn passes_lucas(n: &Integer) -> bool {
  let d_res = choose_d(&n);
  if d_res.is_err() {
    return false;
  }
  let d = d_res.unwrap();
  let q = (int(1) - &d) / 4;
  let delta = int(n + 1);

  let (u_delta, v_delta, q_delta_over_2) =
    compute_lucas_sequences(&delta, n, &int(1), &int(1), &int(1), &q, &d);
  // u_delta % n != 0 proves n composite.
  u_delta == 0
    // Additional check which is not strictly part of Lucas test but nonetheless filters some
    // composite n for free. See section "Checking additional congruence conditions" on Wikipedia.
    && v_delta.is_congruent(&int(2 * &q), &n)
    // Congruence check which holds for prime n by Euler's criterion.
    && q_delta_over_2.is_congruent(&int(&q * q.jacobi(&n)), &n)
}

#[derive(Debug)]
struct IsPerfectSquare();

/// Finds and returns first `D` in `[5, -7, 9, ..., 5 + 2 * max_iter]` for which Jacobi symbol
/// `(D/n) = -1`, or None if no such `D` exists. In the case that `n` is square, there is no such
/// `D` even with `max_iter` infinite. Hence if you are not precisely sure that `n` is nonsquare,
/// you should pass a low value to `max_iter` to avoid wasting too much time. Note that the average
/// number of iterations required for nonsquare `n` is 1.8, and empirically we find it is extremely
/// rare that `|d| > 13`.
///
/// We experimented with postponing the `is_perfect_square` check until after some number of
/// iterations but ultimately found no performance gain. It is likely that most perfect squares
/// are caught by the Miller-Rabin test.
fn choose_d(n: &Integer) -> Result<Integer, IsPerfectSquare> {
  if n.is_perfect_square() {
    return Err(IsPerfectSquare());
  }
  for &d in D_VALUES.iter() {
    let d_ = int(d);
    if d_.jacobi(&n) == -1 {
      return Ok(d_);
    }
  }
  panic!("n is not square but we still couldn't find a d value!")
}

/// Computes the Lucas sequences `{u_i(p, q)}` and `{v_i(p, q)}` up to a specified index `k_target`
/// in O(log(`k_target`)) time by recursively calculating only the `(2i)`th and `(2i+1)`th elements
/// in an order determined by the binary expansion of `k`. Also returns `q^{k/2} (mod n)`, which is
/// used in a stage of the strong Lucas test. In the Lucas case we specify that `d = p^2 - 4q` and
/// set `k_target = delta = n - (d/n) = n + 1`.
fn compute_lucas_sequences(
  k_target: &Integer,
  n: &Integer,
  u_1: &Integer,
  v_1: &Integer,
  p: &Integer,
  q: &Integer,
  d: &Integer,
) -> (Integer, Integer, Integer) {
  let mut u_k = u_1.clone();
  let mut v_k = v_1.clone();
  let mut q_k = q.clone();
  let mut q_k_over_2 = q.clone();
  let mut u_old = Integer::new(); // Ugly performance hack.

  // Finds t in Z_n with 2t = x (mod n).
  // Assumes x in 0..n
  let half = |x: Integer| {
    if x.is_odd() {
      (x + n) / 2
    } else {
      x / 2
    }
  };

  // Write binary expansion of k as [x_1, ..., x_l], e.g. [1, 0, 1, 1] for 11. x_1 is always
  // 1. For i = 2, 3, ..., l, do the following: if x_i = 0 then update u_k and v_k to u_{2k} and
  // v_{2k}, respectively. Else if x_i = 1, update to u_{2k+1} and v_{2k+1}. At the end of the loop
  // we will have computed u_k and v_k, with k as given, in log(delta) time.
  let k_target_bits = k_target.to_digits::<bool>(Order::MsfBe);
  for &bit in k_target_bits.iter().skip(1) {
    // Compute (u, v)_{2k} from (u, v)_k according to the following:
    // u_2k = u_k * v_k (mod n)
    // v_2k = v_k^2 - 2*q^k (mod n)
    u_k = (u_k * &v_k) % n;
    // We use *= for squaring to avoid the performance penalty of unboxing a MulIncomplete.
    v_k *= v_k.clone();
    v_k = (v_k - 2 * &q_k) % n;
    // Continuously maintain q_k = q^k (mod n) and q_k_over_2 = q^{k/2} (mod n).
    q_k_over_2.assign(&q_k);
    q_k *= q_k.clone();
    q_k %= n;
    if bit {
      // Compute (u, v)_{2k+1} from (u, v)_{2k} according to the following:
      // u_{2k+1} = 1/2 * (p*u_{2k} + v_{2k}) (mod n)
      // v_{2k+1} = 1/2 * (d*u_{2k} + p*v_{2k}) (mod n)
      u_old.assign(&u_k);
      u_k = half((p * u_k + &v_k) % n);
      v_k = half((d * &u_old + p * v_k) % n);
      q_k = (q_k * q) % n;
    }
  }
  (u_k, v_k, q_k_over_2)
}

#[cfg(test)]
mod tests {
  use self::constants::*;
  use super::*;

  #[test]
  fn test_miller_rabin() {
    assert!(passes_miller_rabin_base_2(&int(13)));
    assert!(!passes_miller_rabin_base_2(&int(65)));
    for &p in LARGE_PRIMES.iter() {
      assert!(passes_miller_rabin_base_2(&int(p)));
      assert!(!passes_miller_rabin_base_2(&(int(p) * int(106_957))));
    }
    for &n in STRONG_BASE_2_PSEUDOPRIMES.iter() {
      assert!(passes_miller_rabin_base_2(&int(n)));
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
  fn test_lucas() {
    // Should fail on p = 2.
    for &sp in SMALL_PRIMES[1..].iter() {
      assert!(passes_lucas(&int(sp)));
      // Note: The factor cannot be in `SMALL_PRIMES` or this test will fail because `choose_d`
      // fails on square numbers.
      assert!(!passes_lucas(&(int(sp) * 2047)));
    }
    for &mp in MED_PRIMES.iter() {
      assert!(passes_lucas(&int(mp)));
      assert!(!passes_lucas(&(int(mp) * 5)));
    }
    for &lp in LARGE_PRIMES.iter() {
      assert!(passes_lucas(&int(lp)));
      assert!(!passes_lucas(&(int(lp) * 7)));
    }
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
    assert!(is_prob_prime(&int(76463)));
    assert!(is_prob_prime(&int(115_547)));

    // Medium primes.
    for &p in MED_PRIMES.iter() {
      assert!(is_prob_prime(&int(p)));
    }

    // Large primes.
    for &p in LARGE_PRIMES.iter() {
      assert!(is_prob_prime(&int(p)));
    }

    // Large, difficult-to-factor composites.
    for &p in LARGE_PRIMES.iter() {
      for &q in LARGE_PRIMES.iter() {
        assert!(!is_prob_prime(&(int(p) * int(q))));
      }
    }
  }
}
