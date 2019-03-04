#![allow(clippy::unknown_clippy_lints)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::empty_enum)]

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate arrayref;

pub mod util;

mod accumulator;
pub use crate::accumulator::*;
pub mod group;
pub mod hash;
pub mod num;
pub mod proof;

mod vector_commitment;
pub use vector_commitment::*;
mod uint;
pub use uint::*;
