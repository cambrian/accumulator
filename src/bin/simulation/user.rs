/// User executable for simulation.
extern crate crypto;
use crypto::simulation::user::User;
use std::collections::HashSet;

pub fn main() {
  let _user = User::setup(HashSet::new());
  println!("Hello!")
}
