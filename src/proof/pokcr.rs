use crate::group::{multi_exp, Group};
use rug::Integer;

#[allow(non_snake_case)]
#[derive(PartialEq, Eq)]
pub struct PoKCR<G: Group> {
  w: G::Elem,
}

impl<G: Group> PoKCR<G> {
  pub fn prove(witnesses: &[G::Elem]) -> PoKCR<G> {
    PoKCR {
      w: witnesses.iter().fold(G::id(), |a, b| G::op(&a, b)),
    }
  }

  pub fn verify(alphas: &[G::Elem], x: &[Integer], proof: &PoKCR<G>) -> bool {
    let y = multi_exp::<G>(alphas, x);
    let lhs = G::exp(&proof.w, &x.iter().product());
    lhs == y
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::{ElemFrom, RSA2048};
  use crate::util::int;

  #[test]
  fn test_pokcr() {
    let witnesses = [RSA2048::elem(2), RSA2048::elem(3)];
    let x = [int(2), int(2)];
    let alphas = [RSA2048::elem(4), RSA2048::elem(9)];
    let proof = PoKCR::<RSA2048>::prove(&witnesses);
    assert!(proof.w == RSA2048::elem(6));
    assert!(PoKCR::verify(&alphas, &x, &proof));
  }
}
