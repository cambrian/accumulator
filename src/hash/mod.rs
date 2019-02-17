use crate::util::int;
use rug::Integer;
use std::hash::{Hash, Hasher};

mod blake2b;
pub use blake2b::Blake2b;
pub mod primality;

/// Like `std::hash::Hasher`, but general over output type.
pub trait GeneralHasher: Hasher {
  type Output;
  /// Similar to `Hasher::finish`, but consumes `self`.
  fn finalize(self) -> Self::Output;
}

// Note: We explicitly pass in the hasher constructor so we don't have to specify its type via
// generics. Rust has poor support for type applications, so if we wanted to pass `H` at the
// type-level, we'd need to fully specify `T` as well, which is a pain in the ass.
//
// Instead of writing:
// `hash::<Blake2b, (&G::Elem, &BigUint, &G::Elem)>(&(base, exp, result))`
//
// This lets us write:
// `hash(&Blake2b::default, &(base, exp, result))`
pub fn hash<H: GeneralHasher, T: Hash + ?Sized>(new_hasher: &Fn() -> H, t: &T) -> H::Output {
  let mut h = new_hasher();
  t.hash(&mut h);
  h.finalize()
}

/// Calls `hash` with Blake2b hasher.
pub fn blake2b<T: Hash + ?Sized>(t: &T) -> Integer {
  hash(&Blake2b::default, t)
}

/// Hashes `t` with an incrementing counter until a prime is found.
#[allow(clippy::stutter)]
pub fn hash_to_prime_<H: GeneralHasher, T: Hash + ?Sized>(new_hasher: &Fn() -> H, t: &T) -> Integer
where
  Integer: From<H::Output>,
{
  let mut counter = 0_u64;
  loop {
    let mut candidate_prime = int(hash(new_hasher, &(t, counter)));
    // Make the candidate prime odd. This gives ~7% performance gain on a 2018 Macbook Pro.
    candidate_prime.set_bit(0, true);
    if primality::is_prob_prime(&candidate_prime) {
      return candidate_prime;
    }
    counter += 1;
  }
}

/// Calls `hash_to_prime_` with Blake2b hasher.
#[allow(clippy::stutter)]
pub fn hash_to_prime<T: Hash + ?Sized>(t: &T) -> Integer {
  hash_to_prime_(&Blake2b::default, t)
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_blake2() {
    let data = b"martian cyborg gerbil attack";
    hash(&Blake2b::default, data);
  }

  #[test]
  fn test_hash_to_prime() {
    let b_1 = "boom i got ur boyfriend";
    let b_2 = "boom i got ur boyfriene";
    assert_ne!(b_1, b_2);
    let h_1 = hash_to_prime(b_1);
    let h_2 = hash_to_prime(b_2);
    assert_ne!(h_1, h_2);
    assert!(primality::is_prob_prime(&h_1));
    assert!(primality::is_prob_prime(&h_2));
  }
}
