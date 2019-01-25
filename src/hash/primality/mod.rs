use crate::util::int;
use rug::Integer;
mod constants;

pub const MAX_JACOBI_ITERS: u32 = 500;

/// Implements the Baillie-PSW probabilistic primality test, which is known to be deterministic over
/// all integers up to 64 bits (u64). Offers more bang for your buck than Miller-Rabin (i.e.
/// iterated Fermat tests of random base) at wide n since Fermat and Lucas pseudoprimes have been
/// shown to be anticorrelated. Steps of BPSW are as follows:
/// 1. Accept small primes and reject multiples of them.
/// 2. Do a single iteration of Miller-Rabin (base-2 Fermat test).
/// 3. Filter squares.
/// 4. Do a strong probabilistic Lucas test.
pub fn is_prob_prime(n: &Integer) -> bool {
  for &p in constants::SMALL_PRIMES.iter() {
    if n.is_congruent_u(0, p) {
      return *n == p;
    }
  }
  passes_miller_rabin_base_2(&n) && !n.is_perfect_square() && passes_lucas(&n)
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
  x.pow_mod_mut(&d, n).unwrap();
  if x == 1 || x == n.clone() - 1 {
    return true;
  }
  for _ in 1..r {
    x.pow_mod_mut(&int(2), n).unwrap();
    if x == 1 {
      return false;
    }
    if x == n.clone() - 1 {
      return true;
    }
  }
  false
}

/// Strong Lucas probable prime test (NOT the more common Lucas primality test which requires
/// factorization of n-1). Selects parameters d, p, q according to Selfridge's method.
/// Cf. https://en.wikipedia.org/wiki/Lucas_pseudoprime
fn passes_lucas(n: &Integer) -> bool {
  let d = choose_d(&n);
  let p = int(1);
  let q = (1 - d.clone()) / 4;
  let delta = n.clone() + 1;

  // TODO: Extend to stronger test.
  let (u_delta, v_delta, q_delta_over_2) =
    compute_lucas_sequences(&delta, n, &int(1), &p, &p, &q, &d);
  // u_delta % n != 0 proves n composite.
  if u_delta != 0 {
    return false;
  }
  // Additional check which is not strictly part of Lucas test but nonetheless filters some
  // composite n for free. See section "Checking additional congruence conditions" on Wikipedia.
  if !v_delta.is_congruent(&(2 * q.clone()), &n) {
    return false;
  }
  // Congruence check which holds for prime n by Euler's criterion.
  if !q_delta_over_2.is_congruent(&(q.clone() * q.jacobi(&n)), &n) {
    return false;
  }
  // n is either prime or a "strong" Lucas pseudoprime. (The precise meaning of "strong" is not
  // fixed in the literature.) Procedure can be further strengthened by implementing more tests in
  // section 6 of [Baillie & Wagstaff 1980], but for now this is TODO.
  true
}

/// Finds and returns first D in [5, -7, 9, ..., 5 + 2 * max_iter] for which Jacobi symbol (D/n) =
/// -1, or None if no such D exists. In the case that n is square, there is no such D even with
/// max_iter infinite. Hence if you are not precisely sure that n is nonsquare, you should pass a
/// low value to max_iter to avoid wasting too much time. Note that the average number of iterations
/// required for nonsquare n is 1.8, and empirically we find it is extremely rare that |d| > 13.
fn choose_d(n: &Integer) -> Integer {
  let mut d = int(5);
  for _ in 0..MAX_JACOBI_ITERS {
    if d.jacobi(&n) == -1 {
      return d;
    }
    if d < 0 {
      d = -d + 2;
    } else {
      d = -d - 2;
    }
  }
  panic!("Could not find d with (d/n) = -1! Perhaps n is square?");
}

/// Computes the Lucas sequences {u_i(p, q)} and {v_i(p, q)} up to a specified index k_target in
/// O(log(k_target)) time by recursively calculating only the (2i)th and (2i+1)th elements in an
/// order determined by the binary expansion of k. Also returns q^{k/2} (mod n), which is used in
/// a stage of the strong Lucas test. In the Lucas case we specify that d = p^2 - 4q and set
/// k_target = delta = n - (d/n) = n + 1.
fn compute_lucas_sequences(
  k_target: &Integer,
  n: &Integer,
  u_1: &Integer,
  v_1: &Integer,
  p: &Integer,
  q: &Integer,
  d: &Integer,
) -> (Integer, Integer, Integer) {
  let k_target_bits = to_binary(k_target);
  let mut u_k = u_1.clone();
  let mut v_k = v_1.clone();
  let mut q_k = q.clone();
  let mut q_k_over_2 = q.clone();

  let mod_n = |x: &Integer| {
    let (_quotient, rem) = x.clone().div_rem_euc(n.clone());
    rem
  };

  // Finds t in Z_n with 2t = x (mod n)
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
  for bit in k_target_bits[1..].chars() {
    // Compute (u, v)_{2k} from (u, v)_k according to the following:
    // u_2k = u_k * v_k (mod n)
    // v_2k = v_k^2 - 2*q^k (mod n)
    u_k = mod_n(&(u_k.clone() * v_k.clone()));
    v_k.pow_mod_mut(&int(2), n).unwrap();
    v_k = mod_n(&(v_k - 2 * q_k.clone()));
    // Continuously maintain q_k = q^k (mod n) and q_k_over_2 = q^{k/2} (mod n).
    q_k_over_2 = q_k.clone();
    q_k = mod_n(&(q_k.clone() * q_k.clone()));
    k *= 2;
    if bit == '1' {
      // Compute (u, v)_{2k+1} from (u, v)_{2k} according to the following:
      // u_{2k+1} = 1/2 * (p*u_{2k} + v_{2k}) (mod n)
      // v_{2k+1} = 1/2 * (d*u_{2k} + p*v_{2k}) (mod n)
      // TODO: Why is mod_n necessary here?
      let pu_plus_v = p * u_k.clone() + v_k.clone();
      let du_plus_pv = mod_n(&(d * u_k.clone() + p * v_k.clone()));
      u_k = half(&pu_plus_v);
      v_k = half(&du_plus_pv);
      q_k = mod_n(&(q_k.clone() * q.clone()));
      k += 1;
    }
  }
  (u_k, v_k, q_k_over_2)
}

fn to_binary(n: &Integer) -> String {
  format!("{:b}", n)
}

#[cfg(test)]
mod tests {
  use super::*;
  use constants::*;

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
  fn test_to_binary() {
    assert_eq!(to_binary(&int(1)), "1");
    assert_eq!(to_binary(&int(44)), "101100");
  }

  #[test]
  fn test_lucas() {
    // should fail on p = 2
    for &sp in SMALL_PRIMES[1..].iter() {
      assert!(passes_lucas(&int(sp)));
      assert!(!passes_lucas(&(int(sp) * 241)));
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
