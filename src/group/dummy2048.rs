//! Dummy RSA group for 64-bit numbers.
//! Use this group for testing while we figure out ring integration.

use super::{ElemFromUnsigned, Group, UnknownOrderGroup};
use crate::util;
use crate::util::{bi, bu, Singleton};
use num::BigUint;
use num_traits::identities::One;
use num_traits::Unsigned;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub enum DummyRSA2048 {}

/// RSA-2048 modulus, taken from https://en.wikipedia.org/wiki/RSA_numbers#RSA-2048.
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
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct DummyRSA2048Elem {
  val: BigUint,
}

impl Singleton for DummyRSA2048 {
  type Rep = BigUint;
  fn rep() -> &'static Self::Rep {
    &RSA2048_MODULUS
  }
}

pub fn elem_from_biguint(a: &BigUint) -> DummyRSA2048Elem {
  let n = DummyRSA2048::rep();
  let val = a % n;
  if val > n / bu(2u8) {
    DummyRSA2048Elem {
      val: util::mod_euc_big(&-bi(val), n),
    }
  } else {
    DummyRSA2048Elem { val }
  }
}

impl Group for DummyRSA2048 {
  type Elem = DummyRSA2048Elem;
  fn op_(
    modulus: &BigUint,
    a_elem: &DummyRSA2048Elem,
    b_elem: &DummyRSA2048Elem,
  ) -> DummyRSA2048Elem {
    // Note: This is a pretty naive implementation of op.
    let (a, b) = (&a_elem.val, &b_elem.val);
    let op_result = (a * b) % modulus;
    DummyRSA2048::elem_of(op_result)
  }
  fn id_(_: &BigUint) -> DummyRSA2048Elem {
    DummyRSA2048::elem_of(bu(1u8))
  }
  fn inv_(modulus: &BigUint, x: &DummyRSA2048Elem) -> DummyRSA2048Elem {
    let x_big = bu(x.val.clone());
    let (a, _, gcd) = util::bezout(&x_big, modulus);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    DummyRSA2048::elem_of(util::mod_euc_big(&a, modulus))
  }
}

impl ElemFromUnsigned for DummyRSA2048 {
  fn elem_of<U: Unsigned>(n: U) -> DummyRSA2048Elem
  where
    BigUint: From<U>,
  {
    elem_from_biguint(&bu(n))
  }
}

impl UnknownOrderGroup for DummyRSA2048 {
  fn unknown_order_elem_(_: &BigUint) -> DummyRSA2048Elem {
    DummyRSA2048::elem_of(bu(2u8))
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::accumulator;
  use crate::hash::{hash_to_prime, Blake2b};
  use rand::Rng;

  #[test]
  fn test_init() {
    let _x = &DummyRSA2048::rep();
  }

  #[test]
  fn test_op() {
    let a = DummyRSA2048::op(&DummyRSA2048::elem_of(2u32), &DummyRSA2048::elem_of(3u32));
    assert!(a == DummyRSA2048::elem_of(6u32));
    let b = DummyRSA2048::op(
      &DummyRSA2048::elem_of(RSA2048_MODULUS.clone() - bu(2u32)),
      &DummyRSA2048::elem_of(RSA2048_MODULUS.clone() - bu(3u32)),
    );
    assert!(b == DummyRSA2048::elem_of(bu(6u32)));
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
    let a = DummyRSA2048::exp(&DummyRSA2048::elem_of(2u32), &bu(3u32));
    assert!(a == DummyRSA2048::elem_of(8u32));
    let b = DummyRSA2048::exp(&DummyRSA2048::elem_of(2u32), &bu(4096u32));
    assert!(b == DummyRSA2048::elem_of(BigUint::from_str("217207389955395428589369158781869218697519159898401521658993038615824872408108784926597517\
        496727372037176277380476487000099770530440575029170919732871116716934260655466121508332329\
        543615367099810550371217642707848747209719337160655740326150736137284544974770721296865388\
        733305727739636960186370782308858960903126545368015203728531224712542949463283059298449823\
        194163842041340565518401459166858709515078878951293564147044227487142171138804897039341476\
        125519380825017530552968018297030172607314398711102156189885095451290884843968486448057303\
        47466581515692959313583208325725034506693916571047785061884094866050395109710").unwrap()));
    let c = DummyRSA2048::exp(&DummyRSA2048::elem_of(2u32), &RSA2048_MODULUS);
    dbg!(c);
    let d = DummyRSA2048::exp(
      &DummyRSA2048::elem_of(2u32),
      &(RSA2048_MODULUS.clone() * bu(2u32)),
    );
    dbg!(d);
  }

  #[test]
  fn test_inv() {
    let x = DummyRSA2048::elem_of(2u32);
    let inv = DummyRSA2048::inv(&x);
    assert!(DummyRSA2048::op(&x, &inv) == DummyRSA2048::id());
  }

  #[test]
  fn test_add_failure() {
    let mut elems = Vec::new();
    // Works fine for 8 elements, 9 elements causes a panic on .unwrap, indicating ring is likely
    // unable to process this exponent (product of 9 elements).
    for _ in 0..9 {
      let random_bytes = rand::thread_rng().gen::<[u8; 32]>();
      let prime = hash_to_prime(&Blake2b::default, &random_bytes);
      elems.push(prime);
    }
    let acc = accumulator::setup::<DummyRSA2048>();
    accumulator::add::<DummyRSA2048>(acc, &elems[..]);
  }
}
