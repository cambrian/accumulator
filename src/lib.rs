#![allow(clippy::unknown_clippy_lints)]
#![allow(clippy::many_single_char_names)]

#[macro_use]
extern crate lazy_static;

mod accumulator;
pub use crate::accumulator::*;
pub mod group;
pub mod hash;
pub mod proof;
pub mod simulation;
pub mod util;
mod vector_commitment;
pub use vector_commitment::*;
