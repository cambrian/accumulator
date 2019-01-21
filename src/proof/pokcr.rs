// TODO
use super::super::group::multi_exp;
use super::super::group::Group;
use super::super::util::ConvertBytes;
use num::BigInt;
use num_traits::cast::ToPrimitive;
use num_traits::pow;

pub struct PoKCR {
  w: BigInt,
}

pub fn prove_pokcr<G: Group>(_alphas: &[G::Elem], _x: &[BigInt], witness: &[BigInt]) -> PoKCR {
  let w = witness.iter().fold(BigInt::from(1 as u8), |a, b| a * b);
  PoKCR { w }
}

fn exp(base: BigInt, exp: &BigInt) -> BigInt {
  pow(base, exp.to_usize().unwrap())
}

pub fn verify_pokcr<G: Group + ConvertBytes>(
  alphas: &[G::Elem],
  x: &[BigInt],
  _witness: &[BigInt],
  proof: PoKCR,
) -> bool {
  let x_star = x.iter().fold(BigInt::from(1 as u8), |a, b| a * b);
  let y = G::to_le_bytes(&multi_exp::<G>(alphas.len(), alphas, x));
  let lhs = exp(proof.w, &x_star).to_bytes_le().1;
  // Workaround for differently sized vecs from to_bytes functions
  for i in 0..lhs.len() {
    if lhs[i] != y[i] {
      return false;
    }
  }
  true
}

#[cfg(test)]
mod tests {
  use super::super::super::group::dummy::DummyRSA;
  use super::*;

  #[test]
  fn test_pokcr() {
    let alpha_1 = DummyRSA::elem_of(4);
    let alpha_2 = DummyRSA::elem_of(9);
    let x_1 = BigInt::from(2 as u8);
    let x_2 = BigInt::from(2 as u8);
    let w_1 = BigInt::from(2 as u8);
    let w_2 = BigInt::from(3 as u8);
    let alphas = [alpha_1, alpha_2];
    let x = [x_1, x_2];
    let witness = [w_1, w_2];
    let proof = prove_pokcr::<DummyRSA>(&alphas, &x, &witness);
    assert!(proof.w == BigInt::from(6 as u8));
    assert!(verify_pokcr::<DummyRSA>(&alphas, &x, &witness, proof));
  }
}
