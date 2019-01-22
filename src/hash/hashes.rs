use super::primality::is_prob_prime;
use blake2_rfc::blake2b::{blake2b, Blake2bResult};
use num::bigint::{BigUint};
use super::super::util::{bi};
use sha2::{Digest, Sha256};

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

pub fn h_prime(h: &HashFn, data: &[u8]) -> BigUint {
  let mut counter = bi(0);
  loop {
    let hash_val = h(data, Some(&counter.to_bytes_be().1));
    let hash_val_signed = bi(hash_val);
    if is_prob_prime(&hash_val_signed) {
      return hash_val_signed
        .to_biguint()
        .expect("positive BigInt expected");
    }
    counter += 1;
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_blake2() {
    let data = b"test";
    assert_ne!(blake2(data, None), blake2(data, Some(&[1])));
  }

  #[test]
  fn test_sha256() {
    let data = b"hello world";
    assert_ne!(sha256(data, None), sha256(data, Some(&[1])));
    assert_ne!(sha256(data, Some(&[1])), sha256(data, Some(&[2])));
  }

  #[test]
  fn test_h_prime() {
    unimplemented!();
  }
}
