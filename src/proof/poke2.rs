// TODO
use super::PoKE2;
use num::BigUint;

pub fn compute_poke2<G>(base: &G, exp: &BigUint, result: &G) -> PoKE2<G> {
  unimplemented!()
  // PoKE2 {
  //   z: base,
  //   q: result,
  //   r: exp,
  // }
}

pub fn verify_poke2<G>(base: &G, result: &G, proof: &PoKE2<G>) -> bool {
  false
}
