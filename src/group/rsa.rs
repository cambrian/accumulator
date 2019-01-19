use super::{Group, InvertibleGroup};
use num::BigInt;
use num::BigUint;
use ring::rsa::bigint::{Modulus, Elem as RingElem, elem_mul};
use ring::arithmetic::montgomery::R;
use serde::ser::{Serialize, Serializer};

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
  fn base_elem_(&self) -> RSA2048Elem {
    unimplemented!()
    //self.oneRR().decode_once(&self)
  }
  fn op_(&self, RSA2048Elem(a): &RSA2048Elem, RSA2048Elem(b): &RSA2048Elem) -> RSA2048Elem {
    RSA2048Elem(elem_mul(&a, b.clone(), &self))
  }
}
