/// Miner executable for simulation.
extern crate accumulator;
use accumulator::group::Rsa2048;
use accumulator::simulation::miner::Miner;
use accumulator::Accumulator;

pub fn main() {
  let _miner = Miner::<Rsa2048>::setup(true, Accumulator::new(), 0);
  println!("Hello!")
}
