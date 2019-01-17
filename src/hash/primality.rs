use num::bigint::BigUint;

const SMALL_PRIMES: [u64; 50] = [
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
  101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
  197, 199, 211, 223, 227, 229,
];

// Baillie-PSW probabilistic primality test:
// 1. filter composites with small divisors
// 2. do Miller-Rabin base 2
// 3. filter squares
// 4. do Lucas
pub fn is_prob_prime(n: &BigUint) -> bool {
  !has_small_prime_factor(n) && passes_miller_rabin_base_2(n) && !is_square(n) && passes_lucas(n)
}

fn has_small_prime_factor(n: &BigUint) -> bool {
  for &divisor in SMALL_PRIMES.iter() {
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

// WIP
fn passes_miller_rabin_base_2(n: &BigUint) -> bool {
  //
  false
}

fn is_square(n: &BigUint) -> bool {
  true
}

fn passes_lucas(n: &BigUint) -> bool {
  false
}

#[test]
fn test_trial_div() {
  let n_prime = BigUint::from(233u64);
  let n_composite = BigUint::from(50621u64);
  assert!(n_composite == BigUint::from(223u64) * BigUint::from(227u64));
  assert!(!has_small_prime_factor(&n_prime));
  assert!(has_small_prime_factor(&n_composite));
}
