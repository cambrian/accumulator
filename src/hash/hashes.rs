use super::primality::is_prob_prime;
use blake2_rfc::blake2b::{blake2b, Blake2bResult};
use num::bigint::BigUint;
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
  let res_bytes: &[u8] = res.as_bytes();
  BigUint::from_bytes_be(res_bytes)
}

#[allow(dead_code)]
pub fn sha256(data: &[u8], _key: Option<&[u8]>) -> BigUint {
  let mut hasher = Sha256::new();
  hasher.input(data);
  let res_bytes: &[u8] = &hasher.result()[..];
  BigUint::from_bytes_be(res_bytes)
}

#[allow(dead_code)]
pub fn blake2_prime(data: &[u8]) -> BigUint {
  let mut counter = BigUint::from(0u64);
  loop {
    let h = blake2(data, Some(&counter.to_bytes_be()));
    if is_prob_prime(&h) {
      return h;
    }
    counter += 1u64;
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
    let data2 = b"wut";
    assert_ne!(sha256(data, None), sha256(data2, None));
    // assert_eq!(
    //   sha256(data, None),
    //   b"94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
    // );
  }

  // WIP: benchmarking blake2, sha256, and eventually *_prime
  // #[bench]
  // fn bench_blake2(b: &mut Bencher) {
  //   unimplemented!()
  // }
}
