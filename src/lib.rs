#![allow(clippy::many_single_char_names)]

#[macro_use]
extern crate serde_derive;
extern crate lazy_static;

pub mod accumulator;
mod group;
mod hash;
pub mod proof;
mod util;
