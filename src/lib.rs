#![allow(clippy::unknown_clippy_lints)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::empty_enum)]

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate arrayref;

mod accumulator;
pub use accumulator::*;
mod vector_commitment;
pub use vector_commitment::*;

pub mod group;
pub mod hash;
pub mod proof;
pub mod uint;
pub mod util;
