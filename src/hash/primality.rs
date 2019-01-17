use num::bigint::BigUint;

const SMALL_PRIMES: [u64; 50] = [
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
  101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
  197, 199, 211, 223, 227, 229,
];

// Baillie-PSW probabilistic primality test:
// 1. Filter composites with small divisors.
// 2. Do Miller-Rabin base 2.
// 3. Filter squares.
// 4. Do Lucas.
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

// Threefold integer-only square test for BigUint, after https://stackoverflow.com/a/424936:
// 1. Rule out n for which hex representation does not end in 0, 1, 4 or 9
// 2. ...
fn is_square(n: &BigUint) -> bool {
  true
}

fn passes_lucas(n: &BigUint) -> bool {
  false
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
