use crate::util::bi;
use blake2_rfc::blake2b::{blake2b, Blake2bResult};
use num::bigint::BigUint;
use sha2::{Digest, Sha256};

mod primality;

// 32 bytes = 256 bits.
const HASH_LENGTH_IN_BYTES: usize = 32;

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

#[allow(dead_code)]
pub fn sha256(data: &[u8], key: Option<&[u8]>) -> BigUint {
  let mut hasher = Sha256::new();
  if let Some(bytes) = key {
    hasher.input(bytes)
  };
  hasher.input(data);
  let hash_bytes = hasher.result();
  BigUint::from_bytes_be(&hash_bytes[..])
}

type HashFn = Fn(&[u8], Option<&[u8]>) -> BigUint;

#[allow(dead_code)]
pub fn h_prime(h: &HashFn, data: &[u8]) -> BigUint {
  let mut counter = bi(0);
  loop {
    let hash_val = h(data, Some(&counter.to_bytes_be().1));
    if primality::is_prob_prime(&hash_val) {
      return hash_val;
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
  fn test_h_prime() {
    let b_1: &[u8] = b"boom i got ur boyfriend";
    let b_2: &[u8] = b"boom i got ur boyfriene";
    assert_ne!(b_1, b_2);
    assert_ne!(h_prime(&blake2, b_1), h_prime(&blake2, b_2));
    assert_ne!(h_prime(&sha256, b_1), h_prime(&sha256, b_2));
    assert_ne!(h_prime(&blake2, b_1), h_prime(&sha256, b_1));
  }
}
