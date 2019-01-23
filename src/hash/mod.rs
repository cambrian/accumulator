use blake2_rfc::blake2b::{blake2b, Blake2b as Blake2b_, Blake2bResult};
use num::bigint::BigUint;
use sha2::{Digest, Sha256};
use std::hash::{Hash, Hasher};

mod primality;

// 32 bytes = 256 bits.
const HASH_LENGTH_IN_BYTES: usize = 32;

/// Just like Hasher, but general over output type.
pub trait GeneralHasher: Hasher + Default {
  type Output;
  /// Similar to Hasher::finish, but consumes self.
  fn finalize(self) -> Self::Output;
}

pub struct Blake2b(pub Blake2b_);

impl Default for Blake2b {
  fn default() -> Self {
    Blake2b(Blake2b_::new(HASH_LENGTH_IN_BYTES))
  }
}

impl Hasher for Blake2b {
  /// We could return a truncated hash but it's easier just to not use this fn for now.
  fn finish(&self) -> u64 {
    unimplemented!()
  }
  fn write(&mut self, bytes: &[u8]) {
    Blake2b_::update(&mut self.0, bytes)
  }
}

impl GeneralHasher for Blake2b {
  type Output = BigUint;
  fn finalize(self) -> Self::Output {
    BigUint::from_bytes_be(self.0.finalize().as_bytes())
  }
}

// Optional key can be used as a nonce for data in hash function.
pub fn blake2(data: &[u8], key: Option<&[u8]>) -> BigUint {
  let key: &[u8] = match key {
    Some(bytes) => bytes,
    None => &[],
  };
  let res: Blake2bResult = blake2b(HASH_LENGTH_IN_BYTES, key, data);
  let hash_bytes = res.as_bytes();
  BigUint::from_bytes_be(hash_bytes)
}

pub fn sha256(data: &[u8], key: Option<&[u8]>) -> BigUint {
  let mut hasher = Sha256::new();
  if let Some(bytes) = key {
    hasher.input(bytes)
  };
  hasher.input(data);
  let hash_bytes = hasher.result();
  BigUint::from_bytes_be(&hash_bytes[..])
}

// Note: we explicitly pass in the hasher constructor so we don't have to specify its type via
// generics. Rust has poor support for type applications, so if we wanted to pass H at the
// type-level, we'd need to fully specify T as well, which is a pain in the ass.
//
// Instead of writing:
// hash_to_prime::<Blake2b, (&G::Elem, &BigUint, &G::Elem)>(&(base, exp, result))
//
// This lets us write:
// hash_to_prime(Blake2b::default, &(base, exp, result))
pub fn hash<H: GeneralHasher, T: Hash + ?Sized>(new_hasher: &Fn() -> H, t: &T) -> BigUint
where
  BigUint: From<H::Output>,
{
  let mut h = new_hasher();
  t.hash(&mut h);
  BigUint::from(h.finalize())
}

pub fn hash_to_prime<H: GeneralHasher, T: Hash + ?Sized>(new_hasher: &Fn() -> H, t: &T) -> BigUint
where
  BigUint: From<H::Output>,
{
  let mut counter = 0u64;
  loop {
    let mut h = new_hasher();
    counter.hash(&mut h);
    t.hash(&mut h);
    let candidate_prime = BigUint::from(h.finalize());
    if primality::is_prob_prime(&candidate_prime) {
      return candidate_prime;
    }
    counter += 1;
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_blake2() {
    let data = b"martian cyborg gerbil attack";
    assert_ne!(blake2(data, None), blake2(data, Some(&[1])));
  }

  #[test]
  fn test_sha256() {
    let data = b"frick the world shyt aint real i bend tha spoon with my mind";
    assert_ne!(sha256(data, None), sha256(data, Some(&[1])));
    assert_ne!(sha256(data, Some(&[1])), sha256(data, Some(&[2])));
  }

  #[test]
  fn test_hash_to_prime() {
    let b_1 = "boom i got ur boyfriend";
    let b_2 = "boom i got ur boyfriene";
    assert_ne!(b_1, b_2);
    let h_1 = hash_to_prime(Blake2b::default(), b_1);
    let h_2 = hash_to_prime(Blake2b::default(), b_2);
    assert_ne!(h_1, h_2);
    assert!(primality::is_prob_prime(&h_1));
    assert!(primality::is_prob_prime(&h_2));
  }
}
