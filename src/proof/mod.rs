use num::BigUint;
pub mod poe;
pub mod poke2;

pub struct PoE<G> {
  q: G,
}

pub struct PoKE2<G> {
  z: G,
  q: G,
  r: BigUint,
}
