use super::super::util;
use super::{Group, InvertibleGroup};
use core::clone::Clone;
use core::str::FromStr;
use num::BigUint;
use num_traits::identities::One;
use ring::arithmetic::montgomery::{Unencoded, R};
use ring::rsa::bigint::{elem_exp_consttime, elem_mul, Elem as RingElem, Modulus, PrivateExponent};
use serde::ser::{Serialize, Serializer};
use std::result::Result;
use untrusted::Input;

/// Type parameter for ring's modulus. Kind of misleading since it doesn't encode any info itself,
/// but we're forced into this style.
#[derive(PartialEq, Eq, Debug)]
enum M {}

struct RSA2048 {
  m: Modulus<M>,
}

const ELEM_BYTES: usize = 256;

/// RSA-2048 modulus, taken from https://en.wikipedia.org/wiki/RSA_numbers#RSA-2048
const RSA2048_MODULUS_DECIMAL: &str = "25195908475657893494027183240048398571429282126204\
                                       03202777713783604366202070759555626401852588078440\
                                       69182906412495150821892985591491761845028084891200\
                                       72844992687392807287776735971418347270261896375014\
                                       97182469116507761337985909570009733045974880842840\
                                       17974291006424586918171951187461215151726546322822\
                                       16869987549182422433637259085141865462043576798423\
                                       38718477444792073993423658482382428119816381501067\
                                       48104516603773060562016196762561338441436038339044\
                                       14952634432190114657544454178424020924616515723350\
                                       77870774981712577246796292638635637328991215483143\
                                       81678998850404453640235273819513786365643912120103\
                                       97122822120720357";

lazy_static! {
  pub static ref RSA2048_MODULUS: BigUint = BigUint::from_str(RSA2048_MODULUS_DECIMAL).unwrap();
  static ref RSA2048_: RSA2048 = {
    RSA2048 {
      m: Modulus::<M>::from_be_bytes_with_bit_length(Input::from(
        RSA2048_MODULUS.to_bytes_be().as_slice(),
      ))
      .unwrap()
      .0,
    }
  };
}

/// We restrict our group to singly-montgomery-encoded values, to have proper closure.
/// As a result, it's a pain to do some of the operations (see Group impl). We may want to
/// make more substantial changes to our Ring fork to alleviate this.
#[derive(PartialEq, Eq, Clone, Debug)]
struct RSA2048Elem(RingElem<M, R>);

impl Serialize for RSA2048Elem {
  /// TODO: the copying involved in fill_be_bytes is not strictly necessary
  /// In fact the copying involded in serialize isn't strictly necessary either...
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: Serializer,
  {
    let RSA2048Elem(e) = self.clone();
    let mut bytes = [0; ELEM_BYTES];
    e.into_unencoded(&RSA2048::get().m)
      .fill_be_bytes(&mut bytes);
    serializer.serialize_bytes(&bytes)
  }
}

/// TODO: implement {x, -x} cosets
impl Group for RSA2048 {
  type Elem = RSA2048Elem;
  fn get() -> &'static Self {
    &RSA2048_
  }
  fn id_(&self) -> RSA2048Elem {
    RSA2048Elem(self.m.oneR_elem())
  }
  // We multiply an unencoded 2 by a doubly-encoded 1 to get a singly-encoded 2. There is no other
  // way (using only public functions from ring) to get singly-encoded values.
  fn base_elem_(&self) -> RSA2048Elem {
    let unencoded_2 =
      RingElem::from_be_bytes_padded(Input::from(&[2 as u8]), &Self::get().m).unwrap();
    encode(unencoded_2)
  }
  fn op_(&self, RSA2048Elem(a): &RSA2048Elem, RSA2048Elem(b): &RSA2048Elem) -> RSA2048Elem {
    RSA2048Elem(elem_mul(&a, b.clone(), &self.m))
  }
  /// Constant-time exponentiation, via montgomery-multiplication
  /// Can we avoid needing to re-encode the result?
  fn exp(RSA2048Elem(a): &RSA2048Elem, n: &BigUint) -> RSA2048Elem {
    let exponent = PrivateExponent::from_be_bytes_padded(
      Input::from(n.to_bytes_be().as_slice()),
      &Self::get().m,
    )
    .unwrap();
    let unencoded = elem_exp_consttime(a.clone(), &exponent, &Self::get().m).unwrap();
    encode(unencoded)
  }
}

impl InvertibleGroup for RSA2048 {
  fn inv_(&self, x: &RSA2048Elem) -> RSA2048Elem {
    let x_big = biguint_from_elem(x);
    let (a, _, gcd) = util::bezout(&x_big, &RSA2048_MODULUS);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    elem_from_biguint(&util::mod_euc_big::<BigUint>(&a, &RSA2048_MODULUS))
  }
}

fn biguint_from_elem(RSA2048Elem(a): &RSA2048Elem) -> BigUint {
  let mut bytes = [0; ELEM_BYTES];
  a.clone()
    .into_unencoded(&RSA2048::get().m)
    .fill_be_bytes(&mut bytes);
  BigUint::from_bytes_be(&bytes)
}

fn elem_from_biguint(a: &BigUint) -> RSA2048Elem {
  let unencoded =
    RingElem::from_be_bytes_padded(Input::from(a.to_bytes_be().as_slice()), &RSA2048::get().m)
      .unwrap();
  encode(unencoded)
}

/// Performs Montgomery encoding via multiplication with a doubly-encoded 1.
fn encode(a: RingElem<M, Unencoded>) -> RSA2048Elem {
  RSA2048Elem(elem_mul(
    RSA2048::get().m.oneRR().as_ref(),
    a,
    &RSA2048::get().m,
  ))
}

#[cfg(test)]
mod tests {
  use super::*;

  fn elem_of<U>(n: U) -> RSA2048Elem
  where
    BigUint: From<U>,
  {
    elem_from_biguint(&BigUint::from(n))
  }

  #[test]
  fn test_init() {
    let _x = &RSA2048::get().m;
  }

  #[test]
  fn test_op() {
    let a = RSA2048::op(&elem_of(2u32), &elem_of(3u32));
    assert!(a == elem_of(6u32));
    let b = RSA2048::op(
      &elem_of(RSA2048_MODULUS.clone() - BigUint::from(2u32)),
      &elem_of(RSA2048_MODULUS.clone() - BigUint::from(3u32)),
    );
    assert!(b == elem_of(BigUint::from(6u32)));
  }

  /// Tests that -x and x are treated as the same element.
  #[test]
  fn test_cosets() {
    unimplemented!();
    // assert!(elem_of(2) == elem_of(2_794_712_452_613_795));
    // let r = RSA2048::op(&elem_of(931_570_817_537_932), &elem_of(2));
    // assert!(r == elem_of(931_570_817_537_933));
  }

  #[test]
  fn test_exp() {
    let a = RSA2048::exp(&elem_of(2u32), &BigUint::from(3u32));
    assert!(a == elem_of(8u32));
    let b = RSA2048::exp(&elem_of(2u32), &BigUint::from(4096u32));
    assert!(b == elem_of(BigUint::from_str("21720738995539542858936915878186921869751915989840152165899303861582487240810878492659751749672737203717627738047648700009977053044057502917091973287111671693426065546612150833232954361536709981055037121764270784874720971933716065574032615073613728454497477072129686538873330572773963696018637078230885896090312654536801520372853122471254294946328305929844982319416384204134056551840145916685870951507887895129356414704422748714217113880489703934147612551938082501753055296801829703017260731439871110215618988509545129088484396848644805730347466581515692959313583208325725034506693916571047785061884094866050395109710").unwrap()));
  }

  #[test]
  fn test_inv() {
    unimplemented!();
    // let r = RSA2048::inv(&elem_of(2));
    // assert!(r == elem_of(1_397_356_226_306_899));
    // let r = RSA2048::inv(&elem_of(32_416_188_490));
    // assert!(r == elem_of(173_039_603_491_119));
  }
}
