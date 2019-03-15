//! Class Group implementation, with optimizations. Class group ops are run within a
//! a dedicated Mpz context.
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::num::flint;
use crate::num::flint::fmpz;
use crate::num::mpz::{flinty_mpz, Mpz};
use crate::util::{int, TypeRep};
use rug::Integer;
use std::cell::RefCell;

mod class_ctx;
mod discriminant;
mod elem;
mod lin_congruence_ctx;

use class_ctx::ClassCtx;
use discriminant::CLASS_GROUP_DISCRIMINANT;
use elem::ClassElem;

thread_local! {
  // Thread-local context for class group operations.
  static CTX: RefCell<ClassCtx> = Default::default();
}

// Runs the given closure with the Class Context. The expression passed must be
// a closure that takes in an element of type &mut ClassElem. Furthermore, the lambda
// cannot contain subroutines which themselves call the `with_ctx` macro, or the
// compiler will not be happy.
macro_rules! with_ctx {
  ($logic:expr) => {
    CTX.with(|refcell| {
      let mut ctx_ = refcell.borrow_mut();
      $logic(&mut ctx_)
    })
  };
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ClassGroup {}

impl TypeRep for ClassGroup {
  type Rep = Mpz;
  fn rep() -> &'static Self::Rep {
    &CLASS_GROUP_DISCRIMINANT
  }
}

//  Class group operations based on Chia's fantastic doc explaining form class groups:
//  https://github.com/Chia-Network/vdf-competition/blob/master/classgroups.pdf,
//  hereafter refered to as "Binary Quadratic Forms".  Includes optimizations from Chia's VDF
//  competition and from Jacobson, Michael J., and Alfred J. Van Der Poorten.
//  "Computational aspects of NUCOMP."
impl Group for ClassGroup {
  type Elem = ClassElem;

  fn op_(_: &Mpz, x: &ClassElem, y: &ClassElem) -> ClassElem {
    let mut unreduced = with_ctx!(|ctx: &mut ClassCtx| {
      let (g, h, j, w, r, s, t, u, a, b, l, m, mut mu, mut v, mut lambda, mut sigma, k) =
        mut_tuple_elems!(ctx.op_ctx, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

      // Binary Quadratic Forms, 6.1.1
      g.add(&x.b, &y.b);
      g.fdiv_q_ui_mut(2);
      h.sub(&y.b, &x.b);
      h.fdiv_q_ui_mut(2);
      w.gcd(&x.a, &y.a);
      w.gcd_mut(&g);
      j.set(&w);
      r.set_ui(0);
      s.fdiv_q(&x.a, &w);
      t.fdiv_q(&y.a, &w);
      u.fdiv_q(&g, &w);
      a.mul(&t, &u);
      b.mul(&h, &u);
      m.mul(&s, &x.c);
      b.add_mut(&m);
      m.mul(&s, &t);
      ctx
        .lin_cong_ctx
        .solve_linear_congruence(&mut mu, &mut v, &a, &b, &m)
        .unwrap();

      a.mul(&t, &v);
      m.mul(&t, &mu);
      b.sub(&h, &m);
      m.set(&s);
      ctx
        .lin_cong_ctx
        .solve_linear_congruence(&mut lambda, &mut sigma, &a, &b, &m)
        .unwrap();

      a.mul(&v, &lambda);
      k.add(&mu, &a);
      l.mul(&k, &t);
      l.sub_mut(&h);
      l.fdiv_q_mut(&s);
      m.mul(&t, &u);
      m.mul_mut(&k);
      a.mul(&h, &u);
      m.sub_mut(&a);
      a.mul(&x.c, &s);
      m.sub_mut(&a);
      a.mul(&s, &t);
      m.fdiv_q_mut(&a);

      let mut ret = ClassElem::default();

      ret.a.mul(&s, &t);
      a.mul(&r, &u);
      ret.a.sub_mut(&a);

      ret.b.mul(&j, &u);
      a.mul(&m, &r);
      ret.b.add_mut(&a);
      a.mul(&k, &t);
      ret.b.sub_mut(&a);
      a.mul(&l, &s);
      ret.b.sub_mut(&a);

      ret.c.mul(&k, &l);
      a.mul(&j, &m);
      ret.c.sub_mut(&a);
      ret
    });

    Self::reduce_mut(&mut unreduced);
    unreduced
  }

  fn id_(d: &Mpz) -> ClassElem {
    with_ctx!(|ctx: &mut ClassCtx| {
      let (a,) = mut_tuple_elems!(ctx.op_ctx, 0);

      // Binary Quadratic Forms, Definition 5.4
      // The identity is the Principal Form of Discriminant d.
      let mut ret = ClassElem::default();
      ret.a.set_ui(1);
      ret.b.set_ui(1);
      a.sub(&ret.b, &d);
      ret.c.fdiv_q_ui(&a, 4);
      ret
    })
  }

  fn inv_(_: &Mpz, x: &ClassElem) -> ClassElem {
    let mut ret = ClassElem::default();
    ret.a.set(&x.a);
    ret.b.neg(&x.b);
    ret.c.set(&x.c);
    ret
  }

  fn exp_(_: &Mpz, a: &ClassElem, n: &Integer) -> ClassElem {
    let (mut val, mut a, mut n) = {
      if *n < int(0) {
        (Self::id(), Self::inv(&a), int(-n))
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

      ClassGroup::square(&mut a);
      n >>= 1;
    }
  }
}

impl UnknownOrderGroup for ClassGroup {
  fn unknown_order_elem_(d: &Mpz) -> ClassElem {
    // Binary Quadratic Forms, Definition 5.4
    let mut ret = ClassElem::default();
    ret.a.set_ui(2);
    ret.b.set_ui(1);
    ret.c.set_ui(1);
    ret.c.sub_mut(&d);
    ret.c.fdiv_q_ui_mut(8);

    let (a, b, c) = Self::reduce(ret.a, ret.b, ret.c);
    ClassElem { a, b, c }
  }
}

impl<A, B, C> ElemFrom<(A, B, C)> for ClassGroup
where
  Mpz: From<A>,
  Mpz: From<B>,
  Mpz: From<C>,
{
  fn elem(abc: (A, B, C)) -> ClassElem {
    let (a, b, c) = ClassGroup::reduce(Mpz::from(abc.0), Mpz::from(abc.1), Mpz::from(abc.2));

    // Ideally, this should return an error and the
    // return type of ElemFrom should be Result<Self::Elem, Self:err>,
    // but this would require a lot of ugly "unwraps" in the accumulator
    // library. Besides, users should not need to create new class group
    // elements, so an invalid ElemFrom here should signal a severe internal error.
    assert!(ClassGroup::validate(&a, &b, &c));

    ClassElem { a, b, c }
  }
}

impl ClassGroup {
  // Normalize, reduce, and square are public for benchmarking.
  pub fn normalize(mut a: Mpz, mut b: Mpz, mut c: Mpz) -> (Mpz, Mpz, Mpz) {
    let already_normal = with_ctx!(|ctx: &mut ClassCtx| {
      let (scratch,) = mut_tuple_elems!(ctx.op_ctx, 0);
      Self::elem_is_normal(scratch, &a, &b, &c)
    });

    if already_normal {
      (a, b, c)
    } else {
      ClassGroup::normalize_(&mut a, &mut b, &mut c);
      (a, b, c)
    }
  }

  pub fn reduce(a: Mpz, b: Mpz, c: Mpz) -> (Mpz, Mpz, Mpz) {
    let (a, b, c) = Self::normalize(a, b, c);
    let mut elem = ClassElem { a, b, c };
    Self::reduce_(&mut elem);
    Self::normalize(elem.a, elem.b, elem.c)
  }

  pub fn square(x: &mut ClassElem) {
    if cfg!(feature = "nudulp") {
      // NUDULP optimization, maybe using FLINT bindings.
      Self::square_nudulp(x);
    } else {
      Self::square_regular(x);
    }
  }

  // Private functions.

  fn validate(a: &Mpz, b: &Mpz, c: &Mpz) -> bool {
    ClassGroup::discriminant(a, b, c) == *ClassGroup::rep()
  }

  fn discriminant(a: &Mpz, b: &Mpz, c: &Mpz) -> Mpz {
    with_ctx!(|ctx: &mut ClassCtx| {
      let (scratch,) = mut_tuple_elems!(ctx.op_ctx, 0);

      let mut d = Mpz::default();
      d.mul(&b, &b);
      scratch.mul(&a, &c);
      scratch.mul_ui_mut(4);
      d.sub_mut(&scratch);
      d
    })
  }

  fn square_regular(x: &mut ClassElem) {
    with_ctx!(|ctx: &mut ClassCtx| {
      let (mut mu, mut v, m, a, old_a) = mut_tuple_elems!(ctx.op_ctx, 0, 1, 2, 3, 4);

      // Binary Quadratic Forms, 6.3.1
      ctx
        .lin_cong_ctx
        .solve_linear_congruence(&mut mu, &mut v, &x.b, &x.c, &x.a)
        .unwrap();

      m.mul(&x.b, &mu);
      m.sub_mut(&x.c);
      m.fdiv_q_mut(&x.a);

      old_a.set(&x.a);
      x.a.square_mut();

      a.mul(&mu, &old_a);
      a.mul_ui_mut(2);
      x.b.sub_mut(&a);

      x.c.mul(&mu, &mu);
      x.c.sub_mut(&m);
    });

    ClassGroup::reduce_mut(x);
  }

  #[allow(non_snake_case)]
  fn square_nudulp(x: &mut ClassElem) {
    // Jacobson, Michael J., and Alfred J. Van Der Poorten. "Computational aspects of NUCOMP."
    // Algorithm 2 (Alg 2).

    with_ctx!(|ctx: &mut ClassCtx| {
      let (
        G_sq_op,
        scratch,
        y_sq_op,
        By_sq_op,
        Dy_sq_op,
        bx_sq_op,
        by_sq_op,
        dx_sq_op,
        q_sq_op,
        t_sq_op,
        ax_sq_op,
        ay_sq_op,
        Q1_sq_op,
        x_sq_op,
        z_sq_op,
        dy_sq_op,
      ) = mut_tuple_elems!(ctx.op_ctx, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);

      let L_sq_op = &mut ctx.L;

      // Step 1 in Alg 2.
      G_sq_op.gcdext(scratch, y_sq_op, &x.a, &x.b);
      By_sq_op.divexact(&x.a, &G_sq_op);
      Dy_sq_op.divexact(&x.b, &G_sq_op);

      // Step 2 in Alg 2.
      bx_sq_op.mul(&y_sq_op, &x.c);
      bx_sq_op.modulo_mut(&By_sq_op);
      by_sq_op.set(&By_sq_op);

      if by_sq_op.cmpabs(&L_sq_op) <= 0 {
        // Step 4 in Alg 2.
        dx_sq_op.mul(&bx_sq_op, &Dy_sq_op);
        dx_sq_op.sub_mut(&x.c);
        dx_sq_op.divexact_mut(&By_sq_op);
        x.a.mul(&by_sq_op, &by_sq_op);
        x.c.mul(&bx_sq_op, &bx_sq_op);
        t_sq_op.add(&bx_sq_op, &by_sq_op);
        t_sq_op.square_mut();

        x.b.sub_mut(&t_sq_op);
        x.b.add_mut(&x.a);
        x.b.add_mut(&x.c);
        t_sq_op.mul(&G_sq_op, &dx_sq_op);
        x.c.sub_mut(&t_sq_op);
        return;
      }

      // Most of Step 3 in Alg 2.
      if cfg!(feature = "flint") {
        // Subroutine as handled by top entry to the Chia VDF competition "bulaiden."
        let mut fy: fmpz = 0;
        let mut fx: fmpz = 0;
        let mut fby: fmpz = 0;
        let mut fbx: fmpz = 0;
        let mut fL: fmpz = 0;

        // Convert our Mpz's into a form that Flint can understand as an Mpz.
        let mut y_square_clone = flinty_mpz::from(y_sq_op.clone());
        let mut x_square_clone = flinty_mpz::from(x_sq_op.clone());
        let mut by_square_clone = flinty_mpz::from(by_sq_op.clone());
        let mut bx_square_clone = flinty_mpz::from(bx_sq_op.clone());
        let mut L_square_clone = flinty_mpz::from(L_sq_op.clone());

        // Convert the Flint-style Mpz's into true fmpz types.
        flint::set_mpz(&mut fy, &mut y_square_clone);
        flint::set_mpz(&mut fx, &mut x_square_clone);
        flint::set_mpz(&mut fby, &mut by_square_clone);
        flint::set_mpz(&mut fbx, &mut bx_square_clone);
        flint::set_mpz(&mut fL, &mut L_square_clone);

        // Flint Lehmer partial extended GCD.
        flint::xgcd_partial(&mut fy, &mut fx, &mut fby, &mut fbx, &mut fL);

        // Unwrap the fmpz results back to our Mpz type.
        flint::get_mpz(&mut y_square_clone, &mut fy);
        flint::get_mpz(&mut x_square_clone, &mut fx);
        flint::get_mpz(&mut by_square_clone, &mut fby);
        flint::get_mpz(&mut bx_square_clone, &mut fbx);

        *y_sq_op = Mpz::from(y_square_clone);
        *x_sq_op = Mpz::from(x_square_clone);
        *by_sq_op = Mpz::from(by_square_clone);
        *bx_sq_op = Mpz::from(bx_square_clone);

        x_sq_op.neg_mut();
        if x_sq_op.sgn() > 0 {
          y_sq_op.neg_mut();
        } else {
          by_sq_op.neg_mut();
        }
      } else {
        // Subroutine as presented in "Computational aspects of NUCOMP.", Algorithm 2, most of step 3.
        x_sq_op.set_ui(1);
        y_sq_op.set_ui(0);
        z_sq_op.set_ui(0);
        while by_sq_op.cmpabs(&L_sq_op) > 0 && bx_sq_op.sgn() == 0 {
          q_sq_op.fdiv_q(&by_sq_op, &bx_sq_op);
          t_sq_op.fdiv_r(&by_sq_op, &bx_sq_op);

          by_sq_op.set(&bx_sq_op);
          bx_sq_op.set(&t_sq_op);
          y_sq_op.submul(&q_sq_op, &x_sq_op);
          t_sq_op.set(&y_sq_op);
          y_sq_op.set(&x_sq_op);
          x_sq_op.set(&t_sq_op);
          z_sq_op.add_ui_mut(1);
        }

        if z_sq_op.odd() != 0 {
          by_sq_op.neg_mut();
          y_sq_op.neg_mut();
        }
      }

      ax_sq_op.mul(&G_sq_op, &x_sq_op);
      ay_sq_op.mul(&G_sq_op, &y_sq_op);

      // Step 5 in Alg 2.
      t_sq_op.mul(&Dy_sq_op, &bx_sq_op);
      t_sq_op.submul(&x.c, &x_sq_op);
      dx_sq_op.divexact(&t_sq_op, &By_sq_op);
      Q1_sq_op.mul(&y_sq_op, &dx_sq_op);
      dy_sq_op.add(&Q1_sq_op, &Dy_sq_op);
      x.b.add(&dy_sq_op, &Q1_sq_op);
      x.b.mul_mut(&G_sq_op);
      dy_sq_op.divexact_mut(&x_sq_op);
      x.a.mul(&by_sq_op, &by_sq_op);
      x.c.mul(&bx_sq_op, &bx_sq_op);
      t_sq_op.add(&bx_sq_op, &by_sq_op);
      x.b.submul(&t_sq_op, &t_sq_op);
      x.b.add_mut(&x.a);
      x.b.add_mut(&x.c);
      x.a.submul(&ay_sq_op, &dy_sq_op);
      x.c.submul(&ax_sq_op, &dx_sq_op);
    });

    Self::reduce_mut(x);
  }

  fn elem_is_reduced(scratch: &mut Mpz, a: &Mpz, b: &Mpz, c: &Mpz) -> bool {
    Self::elem_is_normal(scratch, a, b, c) && (a <= c && !(a == c && b.cmp_si(0) < 0))
  }

  fn elem_is_normal(scratch: &mut Mpz, a: &Mpz, b: &Mpz, _c: &Mpz) -> bool {
    scratch.neg(&a);
    *scratch < *b && b <= a
  }

  fn reduce_mut(x: &mut ClassElem) {
    Self::normalize_mut(x);
    Self::reduce_(x);
    Self::normalize_mut(x);
  }

  fn reduce_(elem: &mut ClassElem) {
    with_ctx!(|ctx: &mut ClassCtx| {
      let (x, s, old_a, old_b, scratch) = mut_tuple_elems!(ctx.op_ctx, 0, 1, 2, 3, 4);

      // Binary Quadratic Forms, 5.2.1
      while !Self::elem_is_reduced(scratch, &elem.a, &elem.b, &elem.c) {
        s.add(&elem.c, &elem.b);
        x.mul_ui(&elem.c, 2);
        s.fdiv_q_mut(&x);
        old_a.set(&elem.a);
        old_b.set(&elem.b);
        elem.a.set(&elem.c);
        elem.b.neg_mut();
        x.mul(&s, &elem.c);
        x.mul_ui_mut(2);
        elem.b.add_mut(&x);

        elem.c.mul_mut(&s);
        elem.c.mul_mut(&s);
        x.mul(&old_b, &s);
        elem.c.sub_mut(&x);
        elem.c.add_mut(&old_a);
      }
    })
  }

  fn normalize_mut(x: &mut ClassElem) {
    let already_normal = with_ctx!(|ctx: &mut ClassCtx| {
      let (scratch,) = mut_tuple_elems!(ctx.op_ctx, 0);
      if Self::elem_is_normal(scratch, &x.a, &x.b, &x.c) {
        return true;
      }
      false
    });

    if !already_normal {
      ClassGroup::normalize_(&mut x.a, &mut x.b, &mut x.c);
    }
  }

  fn normalize_(a: &mut Mpz, b: &mut Mpz, c: &mut Mpz) {
    with_ctx!(|ctx: &mut ClassCtx| {
      let (r, denom, old_b, ra) = mut_tuple_elems!(ctx.op_ctx, 0, 1, 2, 3);

      // Binary Quadratic Forms, 5.1.1
      r.sub(&a, &b);
      denom.mul_ui(&a, 2);
      r.fdiv_q_mut(&denom);

      old_b.set(&b);

      ra.mul(&r, &a);
      b.add_mut(&ra);
      b.add_mut(&ra);

      ra.mul_mut(&r);
      c.add_mut(&ra);

      ra.set(&r);
      ra.mul_mut(&old_b);
      c.add_mut(&ra);
    })
  }
}

//  Caveat: tests that use "ground truth" use outputs from
//  Chia's sample implementation in python:
//    https://github.com/Chia-Network/vdf-competition/blob/master/inkfish/classgroup.py.
#[cfg(test)]
mod tests {
  use super::*;
  use crate::util::int;
  use std::collections::hash_map::DefaultHasher;
  use std::hash::{Hash, Hasher};
  use std::str::FromStr;
  use std::fs::File;
  use std::io::{BufRead, BufReader, Result};

  #[should_panic]
  #[test]
  fn test_bad_elem() {
    let _ = ClassGroup::elem((1, 2, 3));
  }

  #[test]
  fn test_elem_from() {
    let a1 = Mpz::from_str("16").unwrap();
    let b1 = Mpz::from_str("105").unwrap();
    let c1 = Mpz::from_str(
      "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207",
    )
    .unwrap();

    let a2 = Mpz::from_str("16").unwrap();
    let b2 = Mpz::from_str("9").unwrap();
    let c2 = Mpz::from_str(
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

    let not_reduced = ClassGroup::reduce(not_reduced.a, not_reduced.b, not_reduced.c);
    assert!(
      not_reduced
        == (
          reduced_ground_truth.a,
          reduced_ground_truth.b,
          reduced_ground_truth.c
        )
    );
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
    let reduced = ClassGroup::reduce(not_reduced.a, not_reduced.b, not_reduced.c);
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

    let already_reduced = reduced_ground_truth.clone();
    assert_eq!(already_reduced, reduced_ground_truth);

    assert_ne!(to_reduce, reduced_ground_truth);
    let reduced = ClassGroup::reduce(to_reduce.a, to_reduce.b, to_reduce.c);
    assert_eq!(
      reduced,
      (
        reduced_ground_truth.a,
        reduced_ground_truth.b,
        reduced_ground_truth.c
      )
    );
  }

  #[test]
  fn test_normalize_basic() {
    let unnorm_a = Mpz::from_str("16").unwrap();
    let unnorm_b = Mpz::from_str("105").unwrap();
    let unnorm_c = Mpz::from_str(
      "4783760786688675616733383986925127377420761933775791859799529477781625005833111632534101811\
       0672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194\
       3195315494838983347421441386016611204764255243332731221321519278338873239699989557133287835\
       2685419887133231339948938699768182757831793879217091871179468485931169743972659665650159413\
       8449739494228617068329664776714484742276158090583495714649193839084110987149118615158361352\
       4884884020388947996954204832727089332397513638493972875716927368810312231404469265224318597\
       017389945629057462766047140854869124473221137588347335081555186814207",
    )
    .unwrap();

    let norm_a = Mpz::from_str("16").unwrap();
    let norm_b = Mpz::from_str("9").unwrap();
    let norm_c = Mpz::from_str(
      "4783760786688675616733383986925127377420761933775791859799529477781625005833111632534101811\
       06720472171123774764735020601213528425753087932376215639471576300984851315174010737751911943\
       19531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526\
       85419887133231339948938699768182757831793879217091871179468485931169743972659665650159413844\
       97394942286170683296647767144847422761580905834957146491938390841109871491186151583613524884\
       88402038894799695420483272708933239751363849397287571692736881031223140446926522431859701738\
       9945629057462766047140854869124473221137588347335081555186814036",
    )
    .unwrap();

    let elem_tuple = ClassGroup::normalize(unnorm_a, unnorm_b, unnorm_c);
    assert_eq!((norm_a, norm_b, norm_c), elem_tuple);
  }

  #[test]
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
        assert!(ClassGroup::validate(&g.a, &g.b, &g.c));
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

    // Test with 1000 consecutive elements from the generator.
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();
    for _ in 0..1000 {
      g = ClassGroup::op(&g, &id);
      assert_eq!(g, g_anchor);
    }

    // Test with 1000 random elements.
    for elem in read_group_elems() {
      assert_eq!(ClassGroup::op(&elem, &id), elem);
    }
  }

  #[test]
  fn test_inv_repeated() {
    let id = ClassGroup::id();
    let g_anchor = ClassGroup::unknown_order_elem();
    let mut g = ClassGroup::unknown_order_elem();

    // Test with 1000 consecutive elements from generator.
    for _ in 0..1000 {
      g = ClassGroup::op(&g, &g_anchor);
      let g_inv = ClassGroup::inv(&g);
      assert_eq!(id, ClassGroup::op(&g_inv, &g));
      assert_eq!(id, ClassGroup::op(&g, &g_inv));
      assert_eq!(g, ClassGroup::inv(&g_inv));
    }

    // Test with 1000 random elements
    for elem in read_group_elems() {
      assert_eq!(ClassGroup::op(&elem, &ClassGroup::inv(&elem)), id);
      assert_eq!(ClassGroup::op(&ClassGroup::inv(&elem), &elem), id);
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
    let mut g2 = g.clone();
    // g^4
    ClassGroup::square(&mut g2);
    ClassGroup::square(&mut g2);

    assert_eq!(&g2, &g4);
  }

  #[test]
  fn test_square_repeated() {
    let mut g = ClassGroup::unknown_order_elem();
    let g_ = g.clone();

    for i in 0..12 {
      ClassGroup::square(&mut g);
      let mut base = ClassGroup::id();

      for _ in 0..(2i32.pow(i + 1)) {
        base = ClassGroup::op(&g_, &base);
      }

      assert_eq!(g, base);
    }
  }

  #[test]
  fn test_op_with_gt() {

    // Multiply one way...
    test_with_gt(
      "src/group/class/test/op.txt",
      |s: String| {
        let operands: Vec<&str> = s.split(":").collect();
        let rh = construct_raw_elem_from_string(operands[0]);
        let lh = construct_raw_elem_from_string(operands[1]);
        ClassGroup::op(&rh, &lh)
      }
    );

    // ...and the other way
    test_with_gt(
      "src/group/class/test/op.txt",
      |s: String| {
        let operands: Vec<&str> = s.split(":").collect();
        let rh = construct_raw_elem_from_string(operands[0]);
        let lh = construct_raw_elem_from_string(operands[1]);
        ClassGroup::op(&lh, &rh)
      }
    );
  }

  #[test]
  fn test_normalize_with_gt() {
    test_with_gt(
      "src/group/class/test/normalize.txt",
      |s: String| {
        let mut to_normalize = construct_raw_elem_from_string(&s);
        let (norm_a, norm_b, norm_c) = ClassGroup::normalize(
          to_normalize.a, to_normalize.b, to_normalize.c
        );
        ClassElem {a: norm_a, b: norm_b, c: norm_c}
      }
    );
  }

  #[test]
  fn test_reduce_with_gt() {
    test_with_gt(
      "src/group/class/test/reduce.txt",
      |s: String| {
        let mut to_reduce = construct_raw_elem_from_string(&s);
        let (reduced_a, reduced_b, reduced_c) = ClassGroup::reduce(
          to_reduce.a, to_reduce.b, to_reduce.c
        );
        ClassElem {a: reduced_a, b: reduced_b, c: reduced_c}
      }
    );
  }

  #[test]
  fn test_square_with_gt() {
    test_with_gt(
      "src/group/class/test/square.txt",
      |s: String| {
        let mut to_square = construct_raw_elem_from_string(&s);
        ClassGroup::square(&mut to_square);
        to_square
      }
    );
  }

  #[test]
  fn test_exp_with_gt() {
    test_with_gt(
      "src/group/class/test/exp.txt",
      |s: String| {
        let elem_and_exp: Vec<&str> = s.split(":").collect();
        let elem = construct_raw_elem_from_string(elem_and_exp[0]);
        let exp = Integer::from_str(elem_and_exp[1]).unwrap();
        ClassGroup::exp(&elem, &exp)
      }
    );
  }

  /* Helper functions */

  // Makes a class elem tuple but does not reduce.
  fn construct_raw_elem_from_strings(a: &str, b: &str, c: &str) -> ClassElem {
    ClassElem {
      a: Mpz::from_str(a).unwrap(),
      b: Mpz::from_str(b).unwrap(),
      c: Mpz::from_str(c).unwrap(),
    }
  }

  // Makes a class elem tuple but does not reduce, from a string of format "(a,b,c)"
  fn construct_raw_elem_from_string(elem: &str) -> ClassElem {
    let len = elem.len();
    let parts: Vec<&str> = elem[1..len - 1].split(',').collect();
    construct_raw_elem_from_strings(parts[0], parts[1], parts[2])
  }

  // Given a file path with per-line ground-truth test cases, returns a vector of triples (A, B, C)
  // where A is a colon-separated list of input strings corresponding to class elements of form
  // "(a,b,c)" or single integers (it is up to a user of this function to know which it is),
  // B is a class element corresponding to the ground-truth output, and C is true if this
  // is a positive (i.e. true) example.
  fn parse_gt_file(file_path: &str) -> Vec<(String, ClassElem, bool)> {
    let mut out: Vec<(String, ClassElem, bool)> = Vec::new();
    let file = File::open(file_path).unwrap();
    for line in BufReader::new(file).lines() {
      match line {
        Ok(inner_line) => {
          let main_parts: Vec<&str> = inner_line.split("|").collect();
          let inputs_str = &main_parts[0][1..main_parts[0].len() - 1];
          let output_str = &main_parts[1];
          let pos_or_neg_str = &main_parts[2];
          let true_example = pos_or_neg_str == &"+";
          let out_elem = construct_raw_elem_from_string(output_str);
          out.push((inputs_str.to_string(), out_elem, true_example));
        },
        Err(_) => ()
      }
    }
    out
  }

  // Reads static list of 1000 random group elements into a vector of ClassElems.
  fn read_group_elems() -> Vec<ClassElem> {
    let mut group_elems = Vec::new();
    let file = File::open("src/group/class/test/group_elems.txt").unwrap();
    for line in BufReader::new(file).lines() {
      match line {
        Ok(inner_line) => {
          group_elems.push(construct_raw_elem_from_string(&inner_line));
        },
        Err(_) => ()
      }
    }
    group_elems
  }

  // Generic test function for comparing ground-truth inputs to ground-truth outputs.
  // Expects to be called with the file path containing the test cases and a procedure that maps
  // the input string of the test cases to the final ClassElem to be tested.
  fn test_with_gt(file_path: &str, procedure: fn(String) -> ClassElem) {
    let parsed = parse_gt_file(file_path);
    for (inp_s, gt_out, pos_or_neg) in parsed {
      let to_test = procedure(inp_s);
      if pos_or_neg {
        assert_eq!(to_test, gt_out);
      } else {
        assert_ne!(to_test, gt_out);
      }
    }
  }
}
