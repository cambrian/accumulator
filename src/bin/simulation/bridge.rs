/// Bridge executable for simulation.
extern crate crypto;
use crypto::accumulator::Accumulator;
use crypto::group::Rsa2048;
use crypto::simulation::bridge::Bridge;

pub fn main() {
  let _bridge = Bridge::<Rsa2048>::setup(Accumulator::new(), 0);
  println!("Hello!")
}
