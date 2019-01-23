use crate::group::{multi_exp, InvertibleGroup};
use crate::util::product;
use num::BigInt;

#[allow(non_snake_case)]
#[derive(PartialEq, Eq)]
pub struct PoKCR<G: InvertibleGroup> {
  w: G::Elem,
}

impl<G: InvertibleGroup> PoKCR<G> {
  pub fn prove(witnesses: &[&G::Elem]) -> PoKCR<G> {
    PoKCR {
      w: witnesses.iter().fold(G::id(), |a, b| G::op(&a, b)),
    }
  }

  pub fn verify(alphas: &[&G::Elem], x: &[&BigInt], proof: &PoKCR<G>) -> bool {
    let x_star = product(x);
    let y = multi_exp::<G>(alphas, x);
    let lhs = G::exp_signed(&proof.w, &x_star);
    lhs == y
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::group::DummyRSA;
  use crate::util::bi;

  #[test]
  fn test_pokcr() {
    let witnesses = [&DummyRSA::elem_of(2), &DummyRSA::elem_of(3)];
    let x = [&bi(2), &bi(2)];
    let alphas = [&DummyRSA::elem_of(4), &DummyRSA::elem_of(9)];
    let proof = PoKCR::<DummyRSA>::prove(&witnesses);
    assert!(proof.w == DummyRSA::elem_of(6));
    assert!(PoKCR::verify(&alphas, &x, &proof));
  }
}
