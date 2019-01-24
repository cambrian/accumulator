//! Integration of Brian Smith's ring library into our group interface.
//! There are a lot of rough edges to the interface boundary. One thing to note in particular:
//!
//! We restrict our group to singly-encoded values (via Montgomery encoding), whereas ring
//! will use unencoded, inverse, or doubly-encoded values where they are more convenient. As a
//! result, we perform extra encodings/decodings that ring doesn't. The performance impact of
//! this should be minor, as it makes only the functions exp, id, and base_elem slower by a
//! constant factor. op is unaffected.
use super::{Group, UnknownOrderGroup};
use crate::util::{bezout, mod_euc_big, Singleton};
use num::BigUint;
use num_traits::identities::One;
use ring::arithmetic::montgomery::{Unencoded, R};
use ring::rsa::bigint::{elem_exp_consttime, elem_mul, Elem, Modulus, PrivateExponent};
use std::clone::Clone;
use std::str::FromStr;
use untrusted::Input;

/// Type parameter for ring's modulus. Kind of misleading since it doesn't encode any info itself,
/// but we're forced into this style.
#[derive(PartialEq, Eq, Debug, Hash)]
pub enum M {}

#[derive(Debug, PartialEq, Eq)]
pub enum RSA2048 {}

const ELEM_BYTES: usize = 256;

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
  static ref RSA2048_: Modulus<M> = {
    Modulus::<M>::from_be_bytes_with_bit_length(Input::from(
      RSA2048_MODULUS.to_bytes_be().as_slice(),
    ))
    .unwrap()
    .0
  };
}

/// We restrict our group to singly-montgomery-encoded values, to have proper closure.
/// As a result, it's a pain to do some of the operations (see Group impl). We may want to make
/// more substantial changes to our Ring fork to allgceviate this.
#[derive(PartialEq, Eq, Clone, Debug, Hash)]
pub struct RSA2048Elem(Elem<M, R>);

impl Singleton for RSA2048 {
  type Rep = Modulus<M>;
  fn rep() -> &'static Self::Rep {
    &RSA2048_
  }
}

/// TODO: Implement {x, -x} cosets.
impl Group for RSA2048 {
  type Elem = RSA2048Elem;

  fn id_(m: &Modulus<M>) -> RSA2048Elem {
    RSA2048Elem(m.oneR_elem())
  }

  fn op_(
    m: &Modulus<M>,
    RSA2048Elem(a): &RSA2048Elem,
    RSA2048Elem(b): &RSA2048Elem,
  ) -> RSA2048Elem {
    RSA2048Elem(elem_mul(&a, b.clone(), m))
  }

  /// Constant-time exponentiation via montgomery-multiplication.
  fn exp_(m: &Modulus<M>, RSA2048Elem(a): &RSA2048Elem, n: &BigUint) -> RSA2048Elem {
    let exponent =
      PrivateExponent::from_be_bytes_padded(Input::from(n.to_bytes_be().as_slice()), m).unwrap();
    let unencoded = elem_exp_consttime(a.clone(), &exponent, m).unwrap();
    encode(unencoded)
  }

  fn inv_(_m: &Modulus<M>, x: &RSA2048Elem) -> RSA2048Elem {
    let x_big = biguint_from_elem(x);
    let (a, _, gcd) = bezout(&x_big, &RSA2048_MODULUS);
    assert!(gcd.is_one()); // TODO: Handle this impossibly rare failure?
    elem_from_biguint(&mod_euc_big::<BigUint>(&a, &RSA2048_MODULUS))
  }
}

impl UnknownOrderGroup for RSA2048 {
  fn unknown_order_elem_(m: &Modulus<M>) -> RSA2048Elem {
    let unencoded_2 = Elem::from_be_bytes_padded(Input::from(&[2 as u8]), m).unwrap();
    encode(unencoded_2)
  }
}

fn biguint_from_elem(RSA2048Elem(a): &RSA2048Elem) -> BigUint {
  let mut bytes = [0; ELEM_BYTES];
  a.clone()
    .into_unencoded(&RSA2048::rep())
    .fill_be_bytes(&mut bytes);
  BigUint::from_bytes_be(&bytes)
}

fn elem_from_biguint(a: &BigUint) -> RSA2048Elem {
  let unencoded =
    Elem::from_be_bytes_padded(Input::from(a.to_bytes_be().as_slice()), &RSA2048::rep()).unwrap();
  encode(unencoded)
}

impl RSA2048 {
  pub fn elem_of(val_unbounded: u64) -> RSA2048Elem {
    let n = BigUint::from(val_unbounded);
    elem_from_biguint(&n)
  }
}

/// Performs Montgomery encoding via multiplication with a doubly-encoded 1. This is often necessary
/// because ring will give results in whatever encoding type is most convenient. For unencoded
/// values we have to encode the result explicitly.
fn encode(a: Elem<M, Unencoded>) -> RSA2048Elem {
  RSA2048Elem(elem_mul(
    RSA2048::rep().oneRR().as_ref(),
    a,
    &RSA2048::rep(),
  ))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::accumulator::{add, setup};
  use crate::hash::{hash_to_prime, Blake2b};
  use crate::util::bu;
  use num_traits::Unsigned;
  use rand::Rng;

  fn elem_of<U: Unsigned>(n: U) -> RSA2048Elem
  where
    BigUint: From<U>,
  {
    elem_from_biguint(&bu(n))
  }

  #[test]
  fn test_init() {
    let _x = &RSA2048::rep();
  }

  #[test]
  fn test_op() {
    let a = RSA2048::op(&elem_of(2u32), &elem_of(3u32));
    assert!(a == elem_of(6u32));
    let b = RSA2048::op(
      &elem_of(RSA2048_MODULUS.clone() - bu(2u32)),
      &elem_of(RSA2048_MODULUS.clone() - bu(3u32)),
    );
    assert!(b == elem_of(bu(6u32)));
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
    let a = RSA2048::exp(&elem_of(2u32), &bu(3u32));
    assert!(a == elem_of(8u32));
    let b = RSA2048::exp(&elem_of(2u32), &bu(4096u32));
    assert!(b == elem_of(BigUint::from_str("217207389955395428589369158781869218697519159898401521658993038615824872408108784926597517\
        496727372037176277380476487000099770530440575029170919732871116716934260655466121508332329\
        543615367099810550371217642707848747209719337160655740326150736137284544974770721296865388\
        733305727739636960186370782308858960903126545368015203728531224712542949463283059298449823\
        194163842041340565518401459166858709515078878951293564147044227487142171138804897039341476\
        125519380825017530552968018297030172607314398711102156189885095451290884843968486448057303\
        47466581515692959313583208325725034506693916571047785061884094866050395109710").unwrap()));
    let c = RSA2048::exp(&elem_of(2u32), &RSA2048_MODULUS);
    dbg!(c);
    let d = RSA2048::exp(&elem_of(2u32), &(RSA2048_MODULUS.clone() * bu(2u32)));
    dbg!(d);
  }

  #[test]
  fn test_inv() {
    let x = elem_of(2u32);
    let inv = RSA2048::inv(&x);
    assert!(RSA2048::op(&x, &inv) == RSA2048::id());
  }
}
