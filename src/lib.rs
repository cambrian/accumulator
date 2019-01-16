// apparently this is needed before mod statements to properly import macros
#[macro_use]
extern crate uint;
#[macro_use]
extern crate serde_derive;

pub mod accumulator;
mod group;
mod hash;
pub mod proof;
