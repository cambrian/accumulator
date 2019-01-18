use num::BigInt;
pub mod poe;
pub mod poke2;

/// T is the type of the group element
pub struct PoE<T> {
  q: T,
}

pub struct PoKE2<T> {
  z: T,
  q: T,
  r: BigInt,
}
