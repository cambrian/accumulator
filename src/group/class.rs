//! Fixed-discriminant implementation of a form class group, with optimizations.
//!
//! Using a class group instead of an RSA group for accumulators or vector commitments eliminates
//! the need for a trusted setup, albeit at the expense of slower operations.
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::util;
use crate::util::{int, TypeRep};
use rug::{Assign, Integer};
use std::hash::{Hash, Hasher};
use std::str::FromStr;

#[allow(clippy::module_name_repetitions)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
/// Class group implementation, with future optimizations available via the `--features` flag.
pub enum ClassGroup {}

// 2048-bit prime, negated, congruent to `3 mod 4`. Generated using OpenSSL.
// According to "A Survey of IQ Cryptography" (Buchmann & Hamdy) Table 1, IQ-MPQS for computing
// discrete logarithms in class groups with a 2048-bit discriminant is comparable in complexity to
// GNFS for factoring a 4096-bit integer.
const DISCRIMINANT2048_DECIMAL: &str =
  "-30616069034807523947093657516320815215492876376165067902716988657802400037331914448218251590830\
  1102189519215849430413184776658192481976276720778009261808832630304841711366872161223643645001916\
  6969493423497224870506311710491233557329479816457723381368788734079933165653042145718668727765268\
  0575673207678516369650123480826989387975548598309959486361425021860161020248607833276306314923730\
  9854570972702350567411779734372573754840570138310317754359137013512655926325773048926718050691092\
  9453371727344087286361426404588335160385998280988603297435639020911295652025967761702701701471162\
  3966286152805654229445219531956098223";

lazy_static! {
  pub static ref CLASS_GROUP_DISCRIMINANT: Integer =
    Integer::from_str(DISCRIMINANT2048_DECIMAL).unwrap();
}

#[allow(clippy::module_name_repetitions)]
#[derive(Clone, Debug, Eq)]
/// A class group element, which wraps three `Rug` integers. You should never need to construct a
/// class group element yourself.
pub struct ClassElem {
  a: Integer,
  b: Integer,
  c: Integer,
}

// `ClassElem` and `ClassGroup` ops based on Chia's fantastic doc explaining applied class groups:
// https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf.
impl ClassGroup {
  /// This method is only public for benchmarking. You should not need to use it.
  pub fn normalize(a: Integer, b: Integer, c: Integer) -> (Integer, Integer, Integer) {
    if Self::is_normal(&a, &b, &c) {
      return (a, b, c);
    }
    // r = floor_div((a - b), 2a)
    // (a, b, c) = (a, b + 2ra, ar^2 + br + c)
    let (r, _) = int(&a - &b).div_rem_floor(int(2 * &a));
    let new_b = &b + 2 * int(&r * &a);
    let new_c = c + b * &r + &a * r.square();
    (a, new_b, new_c)
  }

  /// This method is only public for benchmarking. You should not need to use it.
  // Note: Does not return a `ClassElem` because the output is not guaranteed to be
  // a valid `ClassElem` for all inputs.
  pub fn reduce(mut a: Integer, mut b: Integer, mut c: Integer) -> (Integer, Integer, Integer) {
    while !Self::is_reduced(&a, &b, &c) {
      // s = floor_div(c + b, 2c)
      let (s, _) = int(&c + &b).div_rem_floor(int(2 * &c));

      // (a, b, c) = (c, −b + 2sc, cs^2 − bs + a)
      let old_a = a.clone();
      let old_b = b.clone();
      a = c.clone();
      b = -b + 2 * int(&s * &c);
      c = -int(&old_b * &s) + old_a + c * s.square();
    }
    Self::normalize(a, b, c)
  }

  #[allow(non_snake_case)]
  /// This method is only public for benchmarking. You should not need to use it.
  pub fn square(x: &ClassElem) -> ClassElem {
    // Solve `bk = c mod a` for `k`, represented by `mu`, `v` and any integer `n` s.t.
    // `k = mu + v * n`.
    let (mu, _) = util::solve_linear_congruence(&x.b, &x.c, &x.a).unwrap();

    // A = a^2
    // B = b - 2a * mu
    // tmp = (b * mu) / a
    // C = mu^2 - tmp
    let a = int(x.a.square_ref());
    let b = &x.b - int(2 * &x.a) * &mu;
    let (tmp, _) = <(Integer, Integer)>::from(int((&x.b * &mu) - &x.c).div_rem_floor_ref(&x.a));
    let c = mu.square() - tmp;

    Self::elem((a, b, c))
  }

  fn discriminant(a: &Integer, b: &Integer, c: &Integer) -> Integer {
    int(b.square_ref()) - int(4) * a * c
  }

  fn validate(a: &Integer, b: &Integer, c: &Integer) -> bool {
    Self::discriminant(a, b, c) == *Self::rep()
  }

  fn is_reduced(a: &Integer, b: &Integer, c: &Integer) -> bool {
    Self::is_normal(a, b, c) && (a <= c && !(a == c && *b < int(0)))
  }

  fn is_normal(a: &Integer, b: &Integer, _c: &Integer) -> bool {
    -int(a) < int(b) && b <= a
  }
}

impl TypeRep for ClassGroup {
  type Rep = Integer;
  fn rep() -> &'static Self::Rep {
    &CLASS_GROUP_DISCRIMINANT
  }
}

impl Group for ClassGroup {
  type Elem = ClassElem;

  #[allow(non_snake_case)]
  fn op_(_: &Integer, x: &ClassElem, y: &ClassElem) -> ClassElem {
    // g = (b1 + b2) / 2
    // h = (b2 - b1) / 2
    // w = gcd(a1, a2, g)
    let (g, _) = (int(&x.b) + &y.b).div_rem_floor(int(2));
    let (h, _) = (&y.b - int(&x.b)).div_rem_floor(int(2));
    let w = int(x.a.gcd_ref(&y.a)).gcd(&g);

    // j = w
    // s = a1 / w
    // t = a2 / w
    // u = g / ww
    // r = 0
    let j = int(&w);
    let (s, _) = <(Integer, Integer)>::from(x.a.div_rem_floor_ref(&w));
    let (t, _) = <(Integer, Integer)>::from(y.a.div_rem_floor_ref(&w));
    let (u, _) = g.div_rem_floor(w);

    // a = tu
    // b = hu + sc
    // m = st
    // Solve linear congruence `(tu)k = hu + sc mod st` or `ak = b mod m` for solutions `k`.
    let a = int(&t * &u);
    let b = int(&h * &u) + (&s * &x.c);
    let mut m = int(&s * &t);
    let (mu, v) = util::solve_linear_congruence(&a, &b, &m).unwrap();

    // a = tv
    // b = h - t * mu
    // m = s
    // Solve linear congruence `(tv)k = h - t * mu mod s` or `ak = b mod m` for solutions `k`.
    let a = int(&t * &v);
    let b = &h - int(&t * &mu);
    m.assign(&s);
    let (lambda, _) = util::solve_linear_congruence(&a, &b, &m).unwrap();

    // k = mu + v * lambda
    // l = (k * t - h) / s
    // m = (tuk - hu - cs) / st
    let k = &mu + int(&v * &lambda);
    let (l, _) = <(Integer, Integer)>::from((int(&k * &t) - &h).div_rem_floor_ref(&s));
    let (m, _) = (int(&t * &u) * &k - &h * &u - &x.c * &s).div_rem_floor(int(&s * &t));

    // A = st
    // B = ju - kt + ls
    // C = kl - jm
    let a = int(&s * &t);
    let b = int(&j * &u) - (int(&k * &t) + int(&l * &s));
    let c = int(&k * &l) - int(&j * &m);
    Self::elem((a, b, c))
  }

  // Constructs the reduced element directly instead of using `Self::Elem()`.
  fn id_(d: &Integer) -> ClassElem {
    let a = int(1);
    let b = int(1);

    // c = (b * b - d) / 4a
    let (c, _) = int(1 - d).div_rem_floor(int(4));
    ClassElem { a, b, c }
  }

  // Constructs the inverse directly instead of using `Self::Elem()`.
  fn inv_(_: &Integer, x: &ClassElem) -> ClassElem {
    ClassElem {
      a: int(&x.a),
      b: int(-(&x.b)),
      c: int(&x.c),
    }
  }

  fn exp_(_: &Integer, a: &ClassElem, n: &Integer) -> ClassElem {
    let (mut val, mut a, mut n) = {
      if *n < int(0) {
        (Self::id(), Self::inv(a), int(-n))
      } else {
        (Self::id(), a.clone(), n.clone())
      }
    };
    loop {
      if n == int(0) {
        return val;
      }
      if n.is_odd() {
        val = Self::op(&val, &a);
      }
      a = Self::square(&a);
      n >>= 1;
    }
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Integer) -> ClassElem {
    // a = 2
    // b = 1
    // c = (b * b - d) / 4a
    let a = int(2);
    let b = int(1);
    let c = int(1 - d) / int(8);
    ClassElem { a, b, c }
  }
}

impl Hash for ClassElem {
  // Assumes `ClassElem` is reduced and normalized, which will be the case unless a struct is
  // instantiated manually in this module.
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.a.hash(state);
    self.b.hash(state);
    self.c.hash(state);
  }
}

impl PartialEq for ClassElem {
  fn eq(&self, other: &Self) -> bool {
    self.a == other.a && self.b == other.b && self.c == other.c
  }
}

/// Panics if `(a, b, c)` cannot be reduced to a valid class element.
impl<A, B, C> ElemFrom<(A, B, C)> for ClassGroup
where
  Integer: From<A>,
  Integer: From<B>,
  Integer: From<C>,
{
  fn elem(abc: (A, B, C)) -> ClassElem {
    let (a, b, c) = Self::reduce(int(abc.0), int(abc.1), int(abc.2));

    // Ideally, this should return an error and the return type of `ElemFrom` should be
    // `Result<Self::Elem, Self:err>`, but this would require a lot of ugly `unwrap`s in the
    // accumulator library. Besides, users should not need to create new class group elements, so
    // an invalid `ElemFrom` here should signal a severe internal error.
    assert!(Self::validate(&a, &b, &c));

    ClassElem { a, b, c }
  }
}

// Caveat: Tests that use "ground truth" use outputs from Chia's sample implementation in python:
// https://github.com/Chia-Network/vdf-competition/blob/master/inkfish/classgroup.py.
#[cfg(test)]
mod tests {
  use super::*;
  use std::collections::hash_map::DefaultHasher;

  // Makes a class elem tuple but does not reduce.
  fn construct_raw_elem_from_strings(a: &str, b: &str, c: &str) -> ClassElem {
    ClassElem {
      a: Integer::from_str(a).unwrap(),
      b: Integer::from_str(b).unwrap(),
      c: Integer::from_str(c).unwrap(),
    }
  }

  #[should_panic]
  #[test]
  fn test_bad_elem() {
    let _ = ClassGroup::elem((1, 2, 3));
  }

  #[test]
  fn test_elem_from() {
    let a1 = Integer::from_str("16").unwrap();
    let b1 = Integer::from_str("105").unwrap();
    let c1 = Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",
    )
    .unwrap();

    let a2 = Integer::from_str("16").unwrap();
    let b2 = Integer::from_str("9").unwrap();
    let c2 = Integer::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036",
    )
    .unwrap();

    let reduced_elem = ClassGroup::elem((a1, b1, c1));
    let also_reduced_elem = ClassGroup::elem((a2, b2, c2));
    assert_eq!(reduced_elem, also_reduced_elem);
  }

  #[test]
  fn test_equality() {
    let not_reduced = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    let diff_elem = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    assert!(not_reduced != reduced_ground_truth);
    assert!(not_reduced == not_reduced.clone());
    assert!(reduced_ground_truth == reduced_ground_truth.clone());
    assert!(not_reduced != diff_elem);
    assert!(reduced_ground_truth != diff_elem);

    let reduced = ClassGroup::elem((not_reduced.a, not_reduced.b, not_reduced.c));
    assert!(reduced == reduced_ground_truth);
  }

  #[test]
  fn test_hash() {
    let not_reduced = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    let diff_elem = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    let mut hasher_lh = DefaultHasher::new();
    let mut hasher_rh = DefaultHasher::new();
    not_reduced.hash(&mut hasher_lh);
    reduced_ground_truth.hash(&mut hasher_rh);

    assert!(hasher_lh.finish() != hasher_rh.finish());
    assert!(hasher_lh.finish() == hasher_lh.finish());
    assert!(hasher_rh.finish() == hasher_rh.finish());

    hasher_lh = DefaultHasher::new();
    hasher_rh = DefaultHasher::new();
    let reduced = ClassGroup::elem((not_reduced.a, not_reduced.b, not_reduced.c));
    reduced.hash(&mut hasher_lh);
    reduced_ground_truth.hash(&mut hasher_rh);
    assert!(hasher_lh.finish() == hasher_rh.finish());

    hasher_lh = DefaultHasher::new();
    hasher_rh = DefaultHasher::new();
    reduced.hash(&mut hasher_lh);
    diff_elem.hash(&mut hasher_rh);
    assert!(hasher_lh.finish() != hasher_rh.finish());
  }

  #[test]
  fn test_reduce_basic() {
    // Unreduced element.
    let to_reduce = construct_raw_elem_from_strings(
      "59162244921619725812008939143220718157267937427074598447911241410131470159247784852210767449\
      675610037288729551814191198624164179866076352187405442496568188988272422133088755036699145362\
      385840772236403043664778415471196678638241785773530531198720497580622741709880533724904220122\
      358854068046553219863419609777498761804625479650772123754523807001976654588225908928022367436\
      8",
      "18760351095004839755193532164856605650590306627169248964100884295652838905828158941233738613\
      175821849253748329102319504958410190952820220503570113920576542676928659211807590199941027958\
      195895385446372444261885022800653454209101497963588809819572703579484085278913354621371362285\
      341138299691587953249270188429393417132110841259813122945626515477865766896056280729710478647\
      13",
      "14872270891432803054791175727694631095755964943358394411314110783404577714102170379700365256\
      599679049493824862742803590079461712691146098397470840896560034332315858221821103076776907123\
      277315116632337385101204055232891361405428635972040596205450316747012080794838691280547894128\
      246741601088755087359234554141346980837292342320288111397175220296098629890108459305643419353\
      36"
    );

    let reduced_ground_truth = construct_raw_elem_from_strings(
      "26888935961824081232597112540509824504614070059776273347136888921115497522070287009841688662\
      983066376019079593372296556420848446780369918809384119124783870290778875424468497961559643807\
      918398860928578027038014112641529893817109240852544158309292025321122680747989987560029531021\
      808743313150630063377037854944",
      "14529985196481999393995154363327100184407232892559561136140792409262328867440167480822808496\
      853924547751298342980606034124112579835255733824790020119078588372593288210628255956605240171\
      744703418426092073347584357826862813733154338737148962212641444735717023402201569115323580814\
      54099903972209626147819759991",
      "28467266502267127591420289007165819749231433586093061478772560429058231137856046130384492811\
      816456933286039468940950129263300933723839212086399375780796041634531383342902918719073416087\
      614456845205980227091403964285870107268917183244016635907926846271829374679124848388403486656\
      1564478239095738726823372184204"
    );

    let (a, b, c) = ClassGroup::reduce(to_reduce.a, to_reduce.b, to_reduce.c);
    assert_eq!(ClassElem { a, b, c }, reduced_ground_truth.clone());

    let reduced_ground_truth_ = reduced_ground_truth.clone();
    let (a, b, c) = ClassGroup::reduce(
      reduced_ground_truth_.a,
      reduced_ground_truth_.b,
      reduced_ground_truth_.c,
    );
    assert_eq!(ClassElem { a, b, c }, reduced_ground_truth);
  }

  #[test]
  // REVIEW: This test should be restructured to not construct `ClassElem`s but it will do for now.
  fn test_normalize_basic() {
    let unnormalized = construct_raw_elem_from_strings(
      "16",
      "105",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",

    );

    let normalized_ground_truth = construct_raw_elem_from_strings(
      "16",
      "9",
      "4783760786688675616733383986925127377420761933775791859799529477781625005833111632534101811\
       06720472171123774764735020601213528425753087932376215639471576300984851315174010737751911943\
       19531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526\
       85419887133231339948938699768182757831793879217091871179468485931169743972659665650159413844\
       97394942286170683296647767144847422761580905834957146491938390841109871491186151583613524884\
       88402038894799695420483272708933239751363849397287571692736881031223140446926522431859701738\
       9945629057462766047140854869124473221137588347335081555186814036",
    );

    let (a, b, c) = ClassGroup::normalize(unnormalized.a, unnormalized.b, unnormalized.c);
    assert_eq!(normalized_ground_truth, ClassElem { a, b, c });
  }

  #[test]
  // REVIEW: This test should be rewritten, because it may be broken by `unknown_order_elem` not
  // working correctly.
  fn test_discriminant_basic() {
    let g = ClassGroup::unknown_order_elem();
    assert_eq!(
      ClassGroup::discriminant(&g.a, &g.b, &g.c),
      *ClassGroup::rep()
    );
  }

  #[test]
  // REVIEW: This test should be rewritten. See review for `test_discriminant_basic`.
  fn test_discriminant_across_ops() {
    let id = ClassGroup::id();
    let g1 = ClassGroup::unknown_order_elem();
    let g2 = ClassGroup::op(&g1, &g1);
    let g3 = ClassGroup::op(&id, &g2);
    let g3_inv = ClassGroup::inv(&g3);

    assert!(ClassGroup::validate(&id.a, &id.b, &id.c));
    assert!(ClassGroup::validate(&g1.a, &g1.b, &g1.c));
    assert!(ClassGroup::validate(&g2.a, &g2.b, &g2.c));
    assert!(ClassGroup::validate(&g3.a, &g3.b, &g3.c));
    assert!(ClassGroup::validate(&g3_inv.a, &g3_inv.b, &g3_inv.c));
  }

  #[test]
  fn test_op_single() {
    let a = construct_raw_elem_from_strings(
      "4",
      "1",
      "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    let b = construct_raw_elem_from_strings(
      "16",
      "41",
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814061"
    );

    let ground_truth = construct_raw_elem_from_strings(
      "64",
      "9",
      "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      50972369992385512081817723330993784096234932189292318422025780578511173163060796492543474864\
      07264365691511785213717281118305284397086833770388796703509"
    );

    assert_eq!(ClassGroup::op(&a, &b), ground_truth);
  }

  #[test]
  fn test_op_alternating() {
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();
    let mut g_star = ClassGroup::id();

    // g
    g = ClassGroup::op(&g_anchor, &g);

    // g^2, g^* = g^2
    g = ClassGroup::op(&g_anchor, &g);
    g_star = ClassGroup::op(&g, &g_star);

    // g^3
    g = ClassGroup::op(&g_anchor, &g);

    // g^4, g^* = g^2 * g^4 = g^6
    g = ClassGroup::op(&g_anchor, &g);
    g_star = ClassGroup::op(&g, &g_star);

    let ground_truth = construct_raw_elem_from_strings(
      "64",
      "9",
      "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      509723699923855120818177233309937840962349321892923184220257805785111731630607964925434748640\
      7264365691511785213717281118305284397086833770388796703509"
    );

    assert_eq!(ground_truth, g_star);
  }

  #[test]
  fn test_op_complex() {
    // 1. Take g^100, g^200, ..., g^1000.
    // 2. Compute g^* = g^100 * ... * g^1000.
    // 3. For each of g^100, g^200, ..., g^1000 compute the inverse of that element and assert that
    //    g^* * current_inverse = product of g^100, g^200, ..., g^1000 without the inversed-out
    //    element.
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    let mut gs = vec![];
    let mut gs_invs = vec![];

    let mut g_star = ClassGroup::id();
    for i in 1..=1000 {
      g = ClassGroup::op(&g_anchor, &g);
      assert!(ClassGroup::validate(&g.a, &g.b, &g.c));
      if i % 100 == 0 {
        gs.push(g.clone());
        gs_invs.push(ClassGroup::inv(&g));
        g_star = ClassGroup::op(&g, &g_star);
        assert!(ClassGroup::validate(&g_star.a, &g_star.b, &g_star.c));
      }
    }

    let elems_n_invs = gs.iter().zip(gs_invs.iter());
    for (g_elem, g_inv) in elems_n_invs {
      assert!(ClassGroup::validate(&g_elem.a, &g_elem.b, &g_elem.c));
      assert!(ClassGroup::validate(&g_inv.a, &g_inv.b, &g_inv.c));
      let mut curr_prod = ClassGroup::id();
      for elem in &gs {
        if elem != g_elem {
          curr_prod = ClassGroup::op(&curr_prod, &elem);
          assert!(ClassGroup::validate(
            &curr_prod.a,
            &curr_prod.b,
            &curr_prod.c
          ));
        }
      }
      assert_eq!(ClassGroup::id(), ClassGroup::op(&g_inv, &g_elem));
      assert_eq!(curr_prod, ClassGroup::op(&g_inv, &g_star));
    }
  }

  #[test]
  fn test_id_basic() {
    let g = ClassGroup::unknown_order_elem();
    let id = ClassGroup::id();
    assert_eq!(g, ClassGroup::op(&g, &id));
    assert_eq!(g, ClassGroup::op(&id, &g));
    assert_eq!(id, ClassGroup::op(&id, &id));
  }

  #[test]
  fn test_id_repeated() {
    let mut id = ClassGroup::id();
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();
    for _ in 0..1000 {
      id = ClassGroup::op(&id, &id);
      assert_eq!(id, ClassGroup::id());
      g = ClassGroup::op(&g, &ClassGroup::id());
      assert_eq!(g, g_anchor);
    }
  }

  #[test]
  fn test_inv() {
    let id = ClassGroup::id();
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();

    for _ in 0..1000 {
      g = ClassGroup::op(&g, &g_anchor);
      let g_inv = ClassGroup::inv(&g);
      assert_eq!(id, ClassGroup::op(&g_inv, &g));
      assert_eq!(id, ClassGroup::op(&g, &g_inv));
      assert_eq!(g, ClassGroup::inv(&g_inv));
    }
  }

  #[test]
  fn test_exp_basic() {
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    for i in 1..=1000 {
      g = ClassGroup::op(&g, &g_anchor);
      assert_eq!(&g, &ClassGroup::exp(&g_anchor, &int(i)));
    }
  }

  #[test]
  fn test_square_basic() {
    let g = ClassGroup::unknown_order_elem();
    let mut g4 = ClassGroup::id();

    // g^4
    for _ in 0..4 {
      g4 = ClassGroup::op(&g, &g4);
    }

    // g^2
    let mut g2 = ClassGroup::op(&g, &g);

    // g^4
    g2 = ClassGroup::square(&g2);

    assert_eq!(&g2, &g4);
  }
}
