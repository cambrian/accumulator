use super::U256;
use blake2_rfc::blake2b::{blake2b, Blake2bResult};


// 32 bytes = 256 bits.
const HASH_LENGTH_IN_BYTES: usize = 32;

// optional key can be used as a nonce for data in hash function
pub fn blake2(data: &[u8], key: Option<&[u8]>) -> U256 {
  let key: &[u8] = match key {
    Some(bytes) => bytes,
    None => &[],
  };
  let res: Blake2bResult = blake2b(HASH_LENGTH_IN_BYTES, key, data);
  let res_bytes: &[u8] = res.as_bytes();
  U256::from(res_bytes)
}
