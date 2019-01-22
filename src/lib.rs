#![allow(clippy::many_single_char_names)]

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate serde_derive;

pub mod accumulator;
/// Must be public to access from /benches to benchmark properly
pub mod group;
pub mod hash;
pub mod proof;
mod util;
pub mod vector;
