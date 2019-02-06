/// User executable for simulation.
extern crate accumulator;
use accumulator::simulation::user::User;
use std::collections::HashSet;

pub fn main() {
  let _user = User::setup(HashSet::new());
  println!("Hello!")
}
