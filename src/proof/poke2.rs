use crate::group::{Group, InvertibleGroup};
use crate::hash::{hash, hash_to_prime, Blake2b};
use crate::util;
use crate::util::bi;
use num::{BigInt, BigUint};
use num_integer::Integer;

#[allow(non_snake_case)]
#[derive(PartialEq, Eq)]
pub struct PoKE2<G: Group> {
  z: G::Elem,
  Q: G::Elem,
  r: BigUint,
}

impl<G: InvertibleGroup> PoKE2<G> {
  /// See page 16 of B&B.
  pub fn prove(base: &G::Elem, exp: &BigInt, result: &G::Elem) -> PoKE2<G> {
    let g = G::base_elem();
    let z = G::exp_signed(&g, exp);
    let l = hash_to_prime(&Blake2b::default, &(base, result, &z));
    let alpha = hash(&Blake2b::default, &(base, result, &z, &l));
    let q = exp.div_floor(&bi(l.clone()));
    let r = util::mod_euc_big(exp, &l);
    #[allow(non_snake_case)]
    let Q = G::exp_signed(&G::op(&base, &G::exp(&g, &alpha)), &q);
    PoKE2 { z, Q, r }
  }
}

impl<G: Group> PoKE2<G> {
  /// See page 16 of B&B.
  pub fn verify(base: &G::Elem, result: &G::Elem, proof: &PoKE2<G>) -> bool {
    #[allow(non_snake_case)]
    let PoKE2 { z, Q, r } = proof;
    let g = G::base_elem();
    let l = hash_to_prime(&Blake2b::default, &(base, result, &z));
    let alpha = hash(&Blake2b::default, &(base, result, &z, &l));
    let lhs = G::op(
      &G::exp(Q, &l),
      &G::exp(&G::op(&base, &G::exp(&g, &alpha)), &r),
    );
    let rhs = G::op(result, &G::exp(&z, &alpha));
    lhs == rhs
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::RSA2048;
  use crate::util::bu;

  #[test]
  fn test_poke2() {
    // 2^20 = 1048576
    let base = RSA2048::base_elem();
    let exp = bi(20);
    let result = RSA2048::elem_of(1_048_576);
    let proof = PoKE2::<RSA2048>::prove(&base, &exp, &result);
    assert!(PoKE2::verify(&base, &result, &proof));
    // Must compare entire structs since elements z, Q, and r are private.
    assert!(
      proof
        == PoKE2 {
          z: RSA2048::elem_of(1_048_576),
          Q: RSA2048::elem_of(1),
          r: bu(20u8)
        }
    );

    // 2^35 = 34359738368
    let exp_2 = bi(35);
    let result_2 = RSA2048::elem_of(34_359_738_368);
    let proof_2 = PoKE2::<RSA2048>::prove(&base, &exp_2, &result_2);
    assert!(PoKE2::verify(&base, &result_2, &proof_2));
    // Cannot verify wrong base/exp/result triple with wrong pair.
    assert!(!PoKE2::verify(&base, &result_2, &proof));
    assert!(
      proof_2
        == PoKE2 {
          z: RSA2048::elem_of(34_359_738_368),
          Q: RSA2048::elem_of(1),
          r: bu(35u8)
        }
    );
  }

  #[test]
  fn test_poke2_negative() {
    let base = RSA2048::elem_of(2);
    let exp = bi(-5);
    let result = RSA2048::exp_signed(&base, &exp);
    let proof = PoKE2::<RSA2048>::prove(&base, &exp, &result);
    assert!(PoKE2::verify(&base, &result, &proof));
  }
}
