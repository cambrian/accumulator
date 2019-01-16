// apparently this is needed before mod statements to properly import macros
#[macro_use]
extern crate uint;

pub mod accumulator;
mod group;
mod hash;
pub mod proof;
