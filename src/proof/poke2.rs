use crate::group::UnknownOrderGroup;
use crate::hash::{hash, hash_to_prime, Blake2b};
use rug::Integer;

#[allow(non_snake_case)]
#[derive(PartialEq, Eq)]
pub struct PoKE2<G: UnknownOrderGroup> {
  z: G::Elem,
  Q: G::Elem,
  r: Integer,
}

impl<G: UnknownOrderGroup> PoKE2<G> {
  /// See page 16 of B&B.
  pub fn prove(base: &G::Elem, exp: &Integer, result: &G::Elem) -> PoKE2<G> {
    let g = G::unknown_order_elem();
    let z = G::exp(&g, exp);
    let l = hash_to_prime(&Blake2b::default, &(base, result, &z));
    let alpha = hash(&Blake2b::default, &(base, result, &z, &l));
    let (q, r) = exp.clone().div_rem_euc(l);
    #[allow(non_snake_case)]
    let Q = G::exp(&G::op(&base, &G::exp(&g, &alpha)), &q);
    PoKE2 { z, Q, r }
  }

  /// See page 16 of B&B.
  #[allow(non_snake_case)]
  pub fn verify(base: &G::Elem, result: &G::Elem, PoKE2 { z, Q, r }: &PoKE2<G>) -> bool {
    let g = G::unknown_order_elem();
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
  use crate::group::{Group, GroupElemFrom, RSA2048};
  use crate::util::int;

  #[test]
  fn test_poke2() {
    // 2^20 = 1048576
    let base = RSA2048::unknown_order_elem();
    let exp = int(20);
    let result = RSA2048::elem(1_048_576);
    let proof = PoKE2::<RSA2048>::prove(&base, &exp, &result);
    assert!(PoKE2::verify(&base, &result, &proof));
    // Must compare entire structs since elements z, Q, and r are private.
    assert!(
      proof
        == PoKE2 {
          z: RSA2048::elem(1_048_576),
          Q: RSA2048::elem(1),
          r: int(20)
        }
    );

    // 2^35 = 34359738368
    let exp_2 = int(35);
    let result_2 = RSA2048::elem(34_359_738_368u64);
    let proof_2 = PoKE2::<RSA2048>::prove(&base, &exp_2, &result_2);
    assert!(PoKE2::verify(&base, &result_2, &proof_2));
    // Cannot verify wrong base/exp/result triple with wrong pair.
    assert!(!PoKE2::verify(&base, &result_2, &proof));
    assert!(
      proof_2
        == PoKE2 {
          z: RSA2048::elem(34_359_738_368u64),
          Q: RSA2048::elem(1),
          r: int(35)
        }
    );
  }

  #[test]
  fn test_poke2_negative() {
    let base = RSA2048::elem(2);
    let exp = int(-5);
    let result = RSA2048::exp(&base, &exp);
    let proof = PoKE2::<RSA2048>::prove(&base, &exp, &result);
    assert!(PoKE2::verify(&base, &result, &proof));
  }
}
