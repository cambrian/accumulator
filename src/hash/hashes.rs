use blake2_rfc::blake2b::{blake2b, Blake2bResult};
use num::bigint::BigUint;

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

pub fn blake2_prime(_data: &[u8]) -> BigUint {
  unimplemented!()
}

#[test]
fn test_blake2_key() {
  let data = b"test";
  let hash = blake2(data, None);
  let hash2 = blake2(data, Some(&[1]));
  assert_ne!(hash, hash2);
}
