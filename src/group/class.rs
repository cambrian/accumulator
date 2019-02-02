//! Class Group implementation
use super::{Group, UnknownOrderGroup};
use crate::util;
use crate::util::{int, Singleton};
use rug::{Assign, Integer};
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
pub enum ClassGroup {}

// 2048-bit prime, negated, congruent to 3 mod 4..  Generated using OpenSSL.

// According to "A Survey of IQ Cryptography" (Buchmann & Hamdy) Table 1,
// IQ-MPQS for computing discrete logarithms in class groups with a 2048-bit discriminant is
// comparable in complexity to GNFS for factoring a 4096-bit integer.
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

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ClassElem {
  a: Integer,
  b: Integer,
  c: Integer,
}

// ClassElem and ClassGroup ops based on Chia's fantastic doc explaining applied class groups:
//  https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf.

impl ClassElem {
  fn normalize(&mut self) {
    if self._is_normalized() {
      return;
    }
    // r = floor_div((a - b), 2a)
    // (a, b, c) = (a, b + 2ra, ar^2 + br + c)
    let (r, _) = Integer::from(&self.a - &self.b).div_rem_floor(Integer::from(2 * &self.a));
    let new_b = &self.b + 2 * Integer::from(&r * &self.a);
    let new_c = &self.c + Integer::from(&self.b * &r) + &self.a * r.square();
    self.b = new_b;
    self.c = new_c;
  }

  fn reduce(&mut self) {
    self.normalize();
    while !self._is_reduced() {
      // s = floor_div(c + b, 2c)
      let (s, _) = Integer::from(&self.c + &self.b).div_rem_floor(Integer::from(2 * &self.c));
      let old_c = self.c.clone();

      // (a, b, c) = (c, −b + 2sc, cs^2 − bs + a)
      let new_a = old_c;
      let new_b = Integer::from(-&self.b) + 2 * Integer::from(&s * &new_a);
      let new_c = -Integer::from(&self.b * &s) + &self.a + &self.c * s.square();
      self.a = new_a;
      self.b = new_b;
      self.c = new_c;
    }
    self.normalize();
  }

  #[allow(non_snake_case)]
  fn square(&mut self) {
    // Solve `bk = c mod a` for k, represented by mu, v and any integer n s.t. k = mu + v * n
    //
    let (mu, _) = util::solve_linear_congruence(&self.b, &self.c, &self.a);

    // A = a^2
    // B = b - 2a * mu
    // tmp = (b * mu) / a
    // C = mu^2 - tmp
    let A = self.a.clone().square();
    let B = &self.b - Integer::from(2 * &self.a) * &mu;
    let (tmp, _) = Integer::from((&self.b * &mu) - &self.c).div_rem_floor(self.a.clone());
    let C = mu.square() - tmp;
    self.a = A;
    self.b = B;
    self.c = C;
    self.reduce();
  }

  fn _discriminant(&self) -> Integer {
    &self.b * &self.b - Integer::from(4) * &self.a * &self.c
  }

  fn _validate(&self) -> bool {
    &self._discriminant() == ClassGroup::rep()
  }

  fn _is_reduced(&self) -> bool {
    self._is_normalized() && (self.a <= self.c && !(self.a == self.c && self.b < 0))
  }

  fn _is_normalized(&self) -> bool {
    -Integer::from(&self.a) < self.b && self.b <= self.a
  }
}

impl Singleton for ClassGroup {
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
    let (g, _) = (Integer::from(&x.b) + &y.b).div_rem_floor(Integer::from(2));
    let (h, _) = (&y.b - Integer::from(&x.b)).div_rem_floor(Integer::from(2));
    let w = util::three_gcd(&x.a, &y.a, &g);

    // j = w
    // s = a1 / w
    // t = a2 / w
    // u = g / w
    // r = 0
    let j = Integer::from(&w);
    let (s, _) = x.a.clone().div_rem_floor(w.clone());
    let (t, _) = y.a.clone().div_rem_floor(w.clone());
    let (u, _) = g.div_rem_floor(w.clone());

    // a = tu
    // b = hu + sc
    // m = st
    // Solve linear congruence `(tu)k = hu + sc mod st` or `ak = b mod m` for solutions k.
    let a = Integer::from(&t * &u);
    let b = Integer::from(&h * &u) + (&s * &x.c);
    let mut m = Integer::from(&s * &t);
    let (mu, v) = util::solve_linear_congruence(&a, &b, &m);

    // a = tv
    // b = h - t * mu
    // m = s
    // Solve linear congruence `(tv)k = h - t * mu mod s` or `ak = b mod m` for solutions k
    let a = Integer::from(&t * &v);
    let b = &h - Integer::from(&t * &mu);
    m.assign(&s);
    let (lambda, _) = util::solve_linear_congruence(&a, &b, &m);

    // k = mu + v * lambda
    // l = (k * t - h) / s
    // m = (tuk - hu - cs) / st
    let k = &mu + Integer::from(&v * &lambda);
    let (l, _) = (Integer::from(&k * &t) - &h).div_rem_floor(s.clone());
    let (m, _) =
      (Integer::from(&t * &u) * &k - &h * &u - &x.c * &s).div_rem_floor(Integer::from(&s * &t));

    // A = st
    // B = ju - kt + ls
    // C = kl - jm
    let A = Integer::from(&s * &t);
    let B = Integer::from(&j * &u) - (Integer::from(&k * &t) + Integer::from(&l * &s));
    let C = Integer::from(&k * &l) - Integer::from(&j * &m);
    let mut ret = ClassElem { a: A, b: B, c: C };

    ret.reduce();
    ret
  }

  fn id_(d: &Integer) -> ClassElem {
    let a = Integer::from(1);
    let b = Integer::from(1);

    // c = (b * b - d) / 4a
    let (c, _) = Integer::from(1 - d).div_rem_floor(Integer::from(4));
    ClassElem { a, b, c }
  }

  fn inv_(_: &Integer, x: &ClassElem) -> ClassElem {
    ClassElem {
      a: Integer::from(&x.a),
      b: Integer::from(-(&x.b)),
      c: Integer::from(&x.c),
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
      a.square();
      n >>= 1;
    }
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Integer) -> ClassElem {
    // a = 2
    // b = 1
    // c = (b * b - d) / 4a
    let a = Integer::from(2);
    let b = Integer::from(1);
    let c = Integer::from(1 - d) / Integer::from(8);
    let mut ret = ClassElem { a, b, c };
    ret.reduce();
    ret
  }
}

// Caveat: tests that use "ground truth" use outputs from
//  Chia's sample implementation in python:
//    https://github.com/Chia-Network/vdf-competition/blob/master/inkfish/classgroup.py.
#[cfg(test)]
mod tests {
  use super::*;

  fn construct_elem_from_strings(a: &str, b: &str, c: &str) -> ClassElem {
    ClassElem {
      a: Integer::from_str(a).unwrap(),
      b: Integer::from_str(b).unwrap(),
      c: Integer::from_str(c).unwrap(),
    }
  }

  #[test]
  fn test_reduce_basic() {
    let mut to_reduce = construct_elem_from_strings(
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

    let reduced_ground_truth = construct_elem_from_strings(
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

    to_reduce.reduce();
    assert_eq!(to_reduce, reduced_ground_truth);

    let mut already_reduced = reduced_ground_truth.clone();
    already_reduced.reduce();
    assert_eq!(already_reduced, reduced_ground_truth);
  }

  #[test]
  fn test_normalize_basic() {
    let mut unnormalized = construct_elem_from_strings(
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

    let normalized_ground_truth = construct_elem_from_strings(
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

    unnormalized.normalize();
    assert_eq!(normalized_ground_truth, unnormalized);
  }

  #[test]
  fn test_discriminant_basic() {
    let g = ClassGroup::unknown_order_elem();
    assert_eq!(&g._discriminant(), ClassGroup::rep());
  }

  #[test]
  fn test_discriminant_across_ops() {
    let id = ClassGroup::id();
    let g1 = ClassGroup::unknown_order_elem();
    let g2 = ClassGroup::op(&g1, &g1);
    let g3 = ClassGroup::op(&id, &g2);
    let g3_inv = ClassGroup::inv(&g3);

    assert!(id._validate());
    assert!(g1._validate());
    assert!(g2._validate());
    assert!(g3._validate());
    assert!(g3_inv._validate());
  }

  #[test]
  fn test_op_single() {
    let a = construct_elem_from_strings(
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

    let b = construct_elem_from_strings(
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

    let ground_truth = construct_elem_from_strings(
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

    let ground_truth = construct_elem_from_strings(
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
    // 2. Compute g^* = g^100 * ... * g^1000
    // 3. For each of g^100, g^200, ..., g^1000
    //    compute the inverse of the element and
    //    assert that g^* * current_inverse = product of g^100, g^200, ..., g^1000
    //    without the inversed-out element.
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::id();

    let mut gs = vec![];
    let mut gs_invs = vec![];

    let mut g_star = ClassGroup::id();
    for i in 1..=1000 {
      g = ClassGroup::op(&g_anchor, &g);
      assert!(g._validate());
      if i % 100 == 0 {
        gs.push(g.clone());
        gs_invs.push(ClassGroup::inv(&g));
        g_star = ClassGroup::op(&g, &g_star);
        assert!(g_star._validate());
      }
    }

    let elems_n_invs = gs.iter().zip(gs_invs.iter());
    for (g_elem, g_inv) in elems_n_invs {
      assert!(g_elem._validate());
      assert!(g_inv._validate());
      let mut curr_prod = ClassGroup::id();
      for elem in &gs {
        if elem != g_elem {
          curr_prod = ClassGroup::op(&curr_prod, &elem);
          assert!(curr_prod._validate());
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
      assert_eq!(&g, &ClassGroup::exp(&g_anchor, &Integer::from(i)));
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
    g2.square();

    assert_eq!(&g2, &g4);
  }
}
