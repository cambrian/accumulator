// TODO
use num::BigUint;

struct PoKE2i<G> {
  z: G,
  q: G,
  r: BigUint,
}
// pub fn compute_poke2<G: Group>(G base, BigInt exp, G result) -> PoKE2
