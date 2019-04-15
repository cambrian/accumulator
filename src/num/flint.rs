//! Bindings to Flint 2.5.2. Mostly generated with rust-bindgen, modified for readability.
use crate::num::mpz::flinty_mpz;

#[allow(non_camel_case_types)]
pub type fmpz = ::std::os::raw::c_long;

extern "C" {
    #[link_name = "\u{1}_fmpz_get_mpz"]
    fn fmpz_get_mpz(x: *mut flinty_mpz, f: *mut fmpz);
}
extern "C" {
    #[link_name = "\u{1}_fmpz_set_mpz"]
    fn fmpz_set_mpz(f: *mut fmpz, x: *mut flinty_mpz);
}
extern "C" {
    #[link_name = "\u{1}_fmpz_xgcd_partial"]
    fn fmpz_xgcd_partial(
        co2: *mut fmpz,
        co1: *mut fmpz,
        r2: *mut fmpz,
        r1: *mut fmpz,
        L: *mut fmpz,
    );
}

#[inline]
pub fn get_mpz(x: &mut flinty_mpz, f: &mut fmpz) {
    unsafe { fmpz_get_mpz(x, f) }
}

#[inline]
pub fn set_mpz(f: &mut fmpz, x: &mut flinty_mpz) {
    unsafe { fmpz_set_mpz(f, x) }
}

#[inline]
pub fn xgcd_partial(co2: &mut fmpz, co1: &mut fmpz, r2: &mut fmpz, r1: &mut fmpz, l: &mut fmpz) {
    unsafe { fmpz_xgcd_partial(co2, co1, r2, r1, l) }
}

#[test]
fn bindgen_test_layout_flint_mpz_struct() {
    assert_eq!(
        ::std::mem::size_of::<flinty_mpz>(),
        16usize,
        concat!("Size of: ", stringify!(flint_mpz_struct))
    );
    assert_eq!(
        ::std::mem::align_of::<flinty_mpz>(),
        8usize,
        concat!("Alignment of ", stringify!(flint_mpz_struct))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<flinty_mpz>())).mp_alloc as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(flint_mpz_struct),
            "::",
            stringify!(_mp_alloc)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<flinty_mpz>())).mp_size as *const _ as usize },
        4usize,
        concat!(
            "Offset of field: ",
            stringify!(flint_mpz_struct),
            "::",
            stringify!(_mp_size)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<flinty_mpz>())).mp_d as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(flint_mpz_struct),
            "::",
            stringify!(_mp_d)
        )
    );
}
