//! RSA group using rug's GMP integers.
use super::{GroupElemFrom, Group, UnknownOrderGroup};
use crate::util::{int, Singleton};
use rug::Integer;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub enum RSA2048 {}

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
  pub static ref RSA2048_MODULUS: Integer = Integer::from_str(RSA2048_MODULUS_DECIMAL).unwrap();
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct RSA2048Elem(Integer);

impl Singleton for RSA2048 {
  type Rep = Integer;
  fn rep() -> &'static Self::Rep {
    &RSA2048_MODULUS
  }
}

impl Group for RSA2048 {
  type Elem = RSA2048Elem;
  fn op_(
    modulus: &Integer,
    a: &RSA2048Elem,
    b: &RSA2048Elem,
  ) -> RSA2048Elem {
    RSA2048Elem((a.0.clone() * b.0.clone()) % modulus)
  }
  fn id_(_: &Integer) -> RSA2048Elem {
    RSA2048Elem(int(1))
  }
  fn inv_(modulus: &Integer, x: &RSA2048Elem) -> RSA2048Elem {
    RSA2048Elem(x.0.clone().invert(modulus).unwrap())
  }
  fn exp_(modulus: &Integer, x: &RSA2048Elem, n: &Integer) -> RSA2048Elem {
    RSA2048Elem(
        x.0.clone()
          .pow_mod(n, modulus)
          .unwrap()
    )
  }
}

impl<T> GroupElemFrom<T> for RSA2048
where
  Integer: From<T>,
{
  fn elem(t: T) -> RSA2048Elem {
    RSA2048Elem(Integer::from(t))
  }
}

impl UnknownOrderGroup for RSA2048 {
  fn unknown_order_elem_(_: &Integer) -> RSA2048Elem {
    RSA2048Elem(int(2))
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_init() {
    let _x = &RSA2048::rep();
  }

  #[test]
  fn test_op() {
    let a = RSA2048::op(&RSA2048::elem(2), &RSA2048::elem(3));
    assert!(a == RSA2048::elem(6));
    let b = RSA2048::op(
      &RSA2048::elem(RSA2048_MODULUS.clone() - int(2)),
      &RSA2048::elem(RSA2048_MODULUS.clone() - int(3)),
    );
    assert!(b == RSA2048::elem(6));
  }

  /// Tests that -x and x are treated as the same element.
  #[test]
  fn test_cosets() {
    unimplemented!();
    // assert!(elem(2) == elem(2_794_712_452_613_795));
    // let r = RSA2048::op(&elem(931_570_817_537_932), &elem(2));
    // assert!(r == elem(931_570_817_537_933));
  }

  #[test]
  fn test_exp() {
    let a = RSA2048::exp(&RSA2048::elem(2), &int(3));
    assert!(a == RSA2048::elem(8));
    let b = RSA2048::exp(&RSA2048::elem(2), &int(4096u16));
    assert!(b == RSA2048::elem(Integer::parse("217207389955395428589369158781869218697519159898401521658993038615824872408108784926597517\
        496727372037176277380476487000099770530440575029170919732871116716934260655466121508332329\
        543615367099810550371217642707848747209719337160655740326150736137284544974770721296865388\
        733305727739636960186370782308858960903126545368015203728531224712542949463283059298449823\
        194163842041340565518401459166858709515078878951293564147044227487142171138804897039341476\
        125519380825017530552968018297030172607314398711102156189885095451290884843968486448057303\
        47466581515692959313583208325725034506693916571047785061884094866050395109710").unwrap()));
    let c = RSA2048::exp(&RSA2048::elem(2), &RSA2048_MODULUS);
    dbg!(c);
    let d = RSA2048::exp(
      &RSA2048::elem(2),
      &(RSA2048_MODULUS.clone() * int(2)),
    );
    dbg!(d);
  }

  #[test]
  fn test_inv() {
    let x = RSA2048::elem(2);
    let inv = RSA2048::inv(&x);
    assert!(RSA2048::op(&x, &inv) == RSA2048::id());
  }
}
