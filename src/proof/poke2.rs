//! Non-Interactive Proofs of Knowledge of Exponent (NI-PoKE2). See BBF (pages 10 and 42) for
//! details.
use crate::group::UnknownOrderGroup;
use crate::hash::{blake2b, hash_to_prime};
use rug::Integer;

#[allow(non_snake_case)]
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
/// Struct for NI-PoKE2.
pub struct Poke2<G: UnknownOrderGroup> {
  z: G::Elem,
  Q: G::Elem,
  r: Integer,
}

impl<G: UnknownOrderGroup> Poke2<G> {
  /// Computes a proof that you know `exp` s.t. `base ^ exp = result`.
  pub fn prove(base: &G::Elem, exp: &Integer, result: &G::Elem) -> Self {
    let g = G::unknown_order_elem();
    let z = G::exp(&g, exp).unwrap();
    let l = hash_to_prime(&(base, result, &z));
    let alpha = blake2b(&(base, result, &z, &l));
    let (q, r) = <(Integer, Integer)>::from(exp.div_rem_euc_ref(&l));
    #[allow(non_snake_case)]
    let Q = G::exp(&G::op(&base, &G::exp(&g, &alpha).unwrap()), &q).unwrap();
    Self { z, Q, r }
  }

  /// Verifies that the prover knows `exp` s.t. `base ^ exp = result`.
  #[allow(non_snake_case)]
  pub fn verify(base: &G::Elem, result: &G::Elem, Self { z, Q, r }: &Self) -> bool {
    let g = G::unknown_order_elem();
    let l = hash_to_prime(&(base, result, &z));
    let alpha = blake2b(&(base, result, &z, &l));
    let lhs = G::op(
      &G::exp(Q, &l).unwrap(),
      &G::exp(&G::op(&base, &G::exp(&g, &alpha).unwrap()), &r).unwrap(),
    );
    let rhs = G::op(result, &G::exp(&z, &alpha).unwrap());
    lhs == rhs
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{ElemFrom, Group, Rsa2048};
  use crate::util::int;

  #[test]
  fn test_poke2() {
    use std::str::FromStr;
    // sage: w = power_mod(2,exp, modulus)
    let base = Rsa2048::unknown_order_elem();
     let exp = Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110",
    )
    .unwrap();
    let w = Integer::from_str("15237009150211370041572066643854992199159670014401836849321696862635102033487835342310727017245109132166684919786539411147576425083300413858833269356670380323733544946009726244587299888075528737163608201739522141432863879185104979614488213225007619266202959930396741246840028757785072423669876995919918707162762105031693124069429835211177047936412676083097631109112467835488434055566930455343640875193245804869807246696358272733220445826908935579926381184476706321520364895733176015236667338933737155347587968575990509888262873494415904958502766481314251287061092434837169635961698728491245532926158449261934834101518").unwrap();
    let result = Rsa2048::elem(w);
    let proof = Poke2::<Rsa2048>::prove(&base, &exp, &result);
    assert!(Poke2::verify(&base, &result, &proof));


    // 2^35 = 34359738368
    let exp_2 = int(35);
    let result_2 = Rsa2048::elem(34_359_738_368u64);
    let proof_2 = Poke2::<Rsa2048>::prove(&base, &exp_2, &result_2);
    assert!(Poke2::verify(&base, &result_2, &proof_2));
    // Cannot verify wrong base/exp/result triple with wrong pair.
    assert!(!Poke2::verify(&base, &result_2, &proof));

  }

  #[test]
  fn test_poke2_negative() {
    let base = Rsa2048::elem(2);
    let exp = int(-5);
    let result = Rsa2048::exp(&base, &exp).unwrap();
    let proof = Poke2::<Rsa2048>::prove(&base, &exp, &result);
    assert!(Poke2::verify(&base, &result, &proof));
  }
}
