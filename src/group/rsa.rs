use super::{Group, InvertibleGroup};
use num::BigInt;
use num::BigUint;
use ring::rsa::bigint::{Modulus, Elem as RingElem, elem_mul, PrivateExponent, elem_exp_consttime};
use ring::arithmetic::montgomery::{Unencoded, R};
use serde::ser::{Serialize, Serializer};
use untrusted::Input;
use std::result::Result;

const P: u64 = 226_022_213;
const Q: u64 = 12_364_769;

#[derive (PartialEq, Eq)]
enum RSA2048 {}

lazy_static! {
  static ref RSA_2048: Modulus<RSA2048> = {
    // TODO
    unimplemented!()
  };
}

/// We restrict our group to singly-montgomery-encoded values, to have proper closure.
/// As a result, it's a pain to do some of the operations (see Group impl). We may want to
/// make more substantial changes to our Ring fork to alleviate this.
#[derive (PartialEq, Eq, Clone)]
struct RSA2048Elem (RingElem<RSA2048, R>);

impl Serialize for RSA2048Elem {
  /// TODO: the copying involved in fill_be_bytes is not strictly necessary
  /// In fact the copying involded in serialize isn't strictly necessary either...
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: Serializer,
  {
    let RSA2048Elem(e) = self.clone();
    let mut bytes = [0; 256]; // 256 bytes, or 2048 bits
    e.into_unencoded(Modulus::<RSA2048>::get()).fill_be_bytes(&mut bytes);
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
  fn exp(RSA2048Elem(a): &RSA2048Elem, n: &BigUint) -> RSA2048Elem {
    let exponent = PrivateExponent::from_be_bytes_padded(Input::from(n.to_bytes_be().as_slice()), Self::get()).unwrap();
    let unencoded = elem_exp_consttime(a.clone(), &exponent, Self::get()).unwrap();
    RSA2048Elem(encode(unencoded))
  }
}

fn encode(a: RingElem<RSA2048, Unencoded>) -> RingElem<RSA2048, R> {
  elem_mul(Modulus::<RSA2048>::get().oneRR().as_ref(), a, Modulus::<RSA2048>::get())
}