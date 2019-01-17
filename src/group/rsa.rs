// use bigint::uint::U256;
// use alga::general::MultiplicativeGroup;
// use num::One;
// use std::ops::Mul;
use ring::rsa::bigint;

// TODO: use ring here

// struct RSA<T>( T); // represents RSA group elements of type T. Hard-coded modulus for now.

// impl Mul for RSA<u128> {
//   fn mul(self, rhs: RSA<u128>) -> RSA<u128> { RSA::<u128>(self * rhs) }
// }

// impl<T> One for RSA<T> where T: One, RSA<T>: Mul {
//   fn one() -> RSA<T> { RSA::<T>(T::one()) }
// }

// impl MultiplicativeGroup for RSA<u128> {
// }
