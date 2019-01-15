use bigint::uint::U128;
use alga::general::MultiplicativeGroup;

struct RSA<T>(T); // represents RSA group over elements of type T. Hard-coded modulus for now.

impl MultiplicativeGroup for RSA<U128> {}