/// Bridge executable for simulation.
extern crate accumulator;
use accumulator::group::Rsa2048;
use accumulator::simulation::bridge::Bridge;
use accumulator::Accumulator;

pub fn main() {
  let _bridge = Bridge::<Rsa2048>::setup(Accumulator::new(), 0);
  println!("Hello!")
}
