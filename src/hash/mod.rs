mod hashes;
mod primality;

// define U256 type as 4 words of length 64
construct_uint! {
  pub struct U256(4);
}
