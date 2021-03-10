//! Non-Interactive Proofs of Exponentiation (NI-PoE). See BBF (pages 8 and 42) for details.
use crate::group::Group;
use crate::hash::hash_to_prime;
use crate::util::int;
use rug::Integer;

#[allow(non_snake_case)]
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
/// Struct for NI-PoE.
pub struct Poe<G: Group> {
  Q: G::Elem,
}

impl<G: Group> Poe<G> {
  /// Computes a proof that `base ^ exp` was performed to derive `result`.
  pub fn prove(base: &G::Elem, exp: &Integer, result: &G::Elem) -> Self {
    let l = hash_to_prime(&(base, exp, result));
    let q = exp / l;
    Self {
      Q: G::exp(&base, &q).unwrap(),
    }
  }

  /// Verifies that `base ^ exp = result` using the given proof to avoid computation.
  pub fn verify(base: &G::Elem, exp: &Integer, result: &G::Elem, proof: &Self) -> bool {
    let l = hash_to_prime(&(base, exp, result));
    let r = int(exp % &l);
    // w = Q^l * u^r
    let w = G::op(&G::exp(&proof.Q, &l).unwrap(), &G::exp(&base, &r).unwrap());
    w == *result
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{ElemFrom, Rsa2048, UnknownOrderGroup};
  use crate::util::int;

  #[test]
  fn test_poe_small_exp() {
    use std::str::FromStr;
    // sage: w = power_mod(2,exp, modulus)
    let base = Rsa2048::unknown_order_elem();
    let exp = Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110",
    )
    .unwrap();
    let w = Integer::from_str("15237009150211370041572066643854992199159670014401836849321696862635102033487835342310727017245109132166684919786539411147576425083300413858833269356670380323733544946009726244587299888075528737163608201739522141432863879185104979614488213225007619266202959930396741246840028757785072423669876995919918707162762105031693124069429835211177047936412676083097631109112467835488434055566930455343640875193245804869807246696358272733220445826908935579926381184476706321520364895733176015236667338933737155347587968575990509888262873494415904958502766481314251287061092434837169635961698728491245532926158449261934834101518").unwrap();
    let result = Rsa2048::elem(w);
    let proof = Poe::<Rsa2048>::prove(&base, &exp, &result);
    assert!(Poe::verify(&base, &exp, &result, &proof));


    // 2^35 = 34359738368
    let exp_2 = int(35);
    let result_2 = Rsa2048::elem(34_359_738_368u64);
    let proof_2 = Poe::<Rsa2048>::prove(&base, &exp_2, &result_2);
    assert!(Poe::verify(&base, &exp_2, &result_2, &proof_2));

  }
}
