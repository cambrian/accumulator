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

#[derive(PartialEq, Eq)]
enum RSA2048 {}

const ELEM_BYTES: usize = 256;

///
const RSA_2048_MODULUS_DECIMAL: &str = "25195908475657893494027183240048398571429282126204
  03202777713783604366202070759555626401852588078440
  69182906412495150821892985591491761845028084891200
  72844992687392807287776735971418347270261896375014
  97182469116507761337985909570009733045974880842840
  17974291006424586918171951187461215151726546322822
  16869987549182422433637259085141865462043576798423
  38718477444792073993423658482382428119816381501067
  48104516603773060562016196762561338441436038339044
  14952634432190114657544454178424020924616515723350
  77870774981712577246796292638635637328991215483143
  81678998850404453640235273819513786365643912120103
  97122822120720357";

lazy_static! {
  static ref RSA_2048_MODULUS: BigUint = BigUint::from_str(RSA_2048_MODULUS_DECIMAL).unwrap();
  static ref RSA_2048: Modulus<RSA2048> = {
    Modulus::<RSA2048>::from_be_bytes_with_bit_length(Input::from(
      RSA_2048_MODULUS.to_bytes_be().as_slice(),
    ))
    .unwrap()
    .0
  };
}

/// We restrict our group to singly-montgomery-encoded values, to have proper closure.
/// As a result, it's a pain to do some of the operations (see Group impl). We may want to
/// make more substantial changes to our Ring fork to alleviate this.
#[derive(PartialEq, Eq, Clone)]
struct RSA2048Elem(RingElem<RSA2048, R>);

impl Serialize for RSA2048Elem {
  /// TODO: the copying involved in fill_be_bytes is not strictly necessary
  /// In fact the copying involded in serialize isn't strictly necessary either...
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: Serializer,
  {
    let RSA2048Elem(e) = self.clone();
    let mut bytes = [0; ELEM_BYTES];
    e.into_unencoded(Modulus::<RSA2048>::get())
      .fill_be_bytes(&mut bytes);
    serializer.serialize_bytes(&bytes)
  }
}

impl Group for Modulus<RSA2048> {
  type Elem = RSA2048Elem;
  fn get() -> &'static Self {
    &RSA_2048
  }
  fn id_(&self) -> RSA2048Elem {
    RSA2048Elem(self.oneR_elem())
  }
  // We multiply an unencoded 2 by a doubly-encoded 1 to get a singly-encoded 2. There is no other
  // way (using only public functions from ring) to get singly-encoded values.
  fn base_elem_(&self) -> RSA2048Elem {
    let unencoded_2 = RingElem::from_be_bytes_padded(Input::from(&[2 as u8]), Self::get()).unwrap();
    RSA2048Elem(encode(unencoded_2))
  }
  fn op_(&self, RSA2048Elem(a): &RSA2048Elem, RSA2048Elem(b): &RSA2048Elem) -> RSA2048Elem {
    RSA2048Elem(elem_mul(&a, b.clone(), &self))
  }
  /// Constant-time exponentiation, via montgomery-multiplication
  /// Can we avoid needing to re-encode the result?
  fn exp(RSA2048Elem(a): &RSA2048Elem, n: &BigUint) -> RSA2048Elem {
    let exponent =
      PrivateExponent::from_be_bytes_padded(Input::from(n.to_bytes_be().as_slice()), Self::get())
        .unwrap();
    let unencoded = elem_exp_consttime(a.clone(), &exponent, Self::get()).unwrap();
    RSA2048Elem(encode(unencoded))
  }
}

impl InvertibleGroup for Modulus<RSA2048> {
  fn inv_(&self, x: &RSA2048Elem) -> RSA2048Elem {
    let x_big = biguint_from_elem(x);
    let (a, _, gcd) = util::bezout(&x_big, &RSA_2048_MODULUS);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    elem_from_biguint(&util::mod_euc_big::<BigUint>(&a, &RSA_2048_MODULUS))
  }
}

fn biguint_from_elem(RSA2048Elem(a): &RSA2048Elem) -> BigUint {
  let mut bytes = [0; ELEM_BYTES];
  a.clone()
    .into_unencoded(Modulus::<RSA2048>::get())
    .fill_be_bytes(&mut bytes);
  BigUint::from_bytes_be(&bytes)
}

fn elem_from_biguint(a: &BigUint) -> RSA2048Elem {
  let unencoded = RingElem::from_be_bytes_padded(
    Input::from(a.to_bytes_be().as_slice()),
    Modulus::<RSA2048>::get(),
  )
  .unwrap();
  RSA2048Elem(encode(unencoded))
}

/// Performs Montgomery encoding via multiplication with a doubly-encoded 1.
fn encode(a: RingElem<RSA2048, Unencoded>) -> RingElem<RSA2048, R> {
  elem_mul(
    Modulus::<RSA2048>::get().oneRR().as_ref(),
    a,
    Modulus::<RSA2048>::get(),
  )
}
