use blake2_rfc::blake2b::{blake2b, Blake2bResult};
use num::bigint::BigUint;
// use sha2::Sha256;

// 32 bytes = 256 bits.
const HASH_LENGTH_IN_BYTES: usize = 32;

// optional key can be used as a nonce for data in hash function
pub fn blake2(data: &[u8], key: Option<&[u8]>) -> BigUint {
  let key: &[u8] = match key {
    Some(bytes) => bytes,
    None => &[],
  };
  let res: Blake2bResult = blake2b(HASH_LENGTH_IN_BYTES, key, data);
  let res_bytes: &[u8] = res.as_bytes();
  BigUint::from_bytes_be(res_bytes)
}

#[allow(dead_code)]
pub fn sha256(_data: &[u8], _key: Option<&[u8]>) -> BigUint {
  unimplemented!();
  // let mut hasher = Sha256::new();
  // hasher.input(data);
  // hasher.result();
}

pub fn blake2_prime(_data: &[u8]) -> BigUint {
  unimplemented!()
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_blake2() {
    let data = b"test";
    let hash = blake2(data, None);
    let hash2 = blake2(data, Some(&[1]));
    assert_ne!(hash, hash2);
  }

  #[test]
  fn test_sha256() {
    let _data = b"hello world";
    // let hash = sha256(data, None);
    // assert_eq!(hash, b"94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9");
  }
}
