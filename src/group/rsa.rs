//! RSA (2048) group using rug's GMP integers.
use super::{Group, GroupElemFrom, UnknownOrderGroup};
use crate::util::{int, Singleton};
use rug::Integer;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub enum RSA2048 {}

/// RSA-2048 modulus, taken from https://en.wikipedia.org/wiki/RSA_numbers#RSA-2048.
const RSA2048_MODULUS_DECIMAL: &str = "25195908475657893494027183240048398571429282126204032027777\
                                       13783604366202070759555626401852588078440691829064124951508\
                                       21892985591491761845028084891200728449926873928072877767359\
                                       71418347270261896375014971824691165077613379859095700097330\
                                       45974880842840179742910064245869181719511874612151517265463\
                                       22822168699875491824224336372590851418654620435767984233871\
                                       84774447920739934236584823824281198163815010674810451660377\
                                       30605620161967625613384414360383390441495263443219011465754\
                                       44541784240209246165157233507787077498171257724679629263863\
                                       56373289912154831438167899885040445364023527381951378636564\
                                       391212010397122822120720357";

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

/// Review: whoever changed all of the RSA2048Elem(x) constructors to RSA2048::Elem change them
/// back. In these cases we know that the int value satisfies 0 <= x < modulus so we don't need to
/// pay the price of checking.
impl Group for RSA2048 {
  type Elem = RSA2048Elem;
  fn op_(modulus: &Integer, a: &RSA2048Elem, b: &RSA2048Elem) -> RSA2048Elem {
    RSA2048::elem((a.0.clone() * b.0.clone()) % modulus)
  }
  fn id_(_: &Integer) -> RSA2048Elem {
    RSA2048::elem(1)
  }
  fn inv_(modulus: &Integer, x: &RSA2048Elem) -> RSA2048Elem {
    RSA2048::elem(x.0.invert_ref(modulus).unwrap())
  }
  fn exp_(modulus: &Integer, x: &RSA2048Elem, n: &Integer) -> RSA2048Elem {
    RSA2048::elem(x.0.pow_mod_ref(n, modulus).unwrap())
    // This is the code we'd use if we need side-channel attack resilience. It's 40% slower so not
    // using it til we decide we have to.
    // if *n < int(0) {
    //   RSA2048::elem(Integer::from(x.0.invert_ref(modulus).unwrap()).secure_pow_mod(&-n.clone(), modulus))
    // }
    // else {
    //   RSA2048::elem(x.0.secure_pow_mod_ref(n, modulus))
    // }
  }
}

// We would normally do this inline, but the Integer: From<T> constraint on GroupElemFrom creates
// problems with doing so.
fn half(x: &Integer) -> Integer {
  int(x / 2)
}

impl<T> GroupElemFrom<T> for RSA2048
where
  Integer: From<T>,
{
  fn elem(t: T) -> RSA2048Elem {
    let modulus = Self::rep();
    let val = int(t) % modulus;
    if val > half(modulus) {
      RSA2048Elem((-val).div_rem_euc(modulus.clone()).1)
    } else {
      RSA2048Elem(val)
    }
  }
}

impl UnknownOrderGroup for RSA2048 {
  fn unknown_order_elem_(_: &Integer) -> RSA2048Elem {
    RSA2048::elem(2)
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
    assert!(RSA2048::elem(3) == RSA2048::elem(RSA2048_MODULUS.clone() - 3));
    // TODO: Add a trickier coset test involving `op`.
  }

  #[test]
  fn test_exp() {
    let a = RSA2048::exp(&RSA2048::elem(2), &int(3));
    assert!(a == RSA2048::elem(8));
    let b = RSA2048::exp(&RSA2048::elem(2), &int(4096));
    assert!(
      b == RSA2048::elem(
        Integer::parse(
          "217207389955395428589369158781869218697519159898401521658993038615824872408108784926597\
           517496727372037176277380476487000099770530440575029170919732871116716934260655466121508\
           332329543615367099810550371217642707848747209719337160655740326150736137284544974770721\
           296865388733305727739636960186370782308858960903126545368015203728531224712542949463283\
           059298449823194163842041340565518401459166858709515078878951293564147044227487142171138\
           804897039341476125519380825017530552968018297030172607314398711102156189885095451290884\
           843968486448057303474665815156929593135832083257250345066939165710477850618840948660503\
           95109710"
        )
        .unwrap()
      )
    );
    let c = RSA2048::exp(&RSA2048::elem(2), &RSA2048_MODULUS);
    dbg!(c);
    let d = RSA2048::exp(&RSA2048::elem(2), &(RSA2048_MODULUS.clone() * int(2)));
    dbg!(d);
  }

  #[test]
  fn test_inv() {
    let x = RSA2048::elem(2);
    let inv = RSA2048::inv(&x);
    assert!(RSA2048::op(&x, &inv) == RSA2048::id());
  }
}
