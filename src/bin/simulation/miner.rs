/// Miner executable for simulation.
extern crate crypto;
use crypto::accumulator::Accumulator;
use crypto::group::Rsa2048;
use crypto::simulation::miner::Miner;

pub fn main() {
  let _miner = Miner::<Rsa2048>::setup(true, Accumulator::new(), 0);
  println!("Hello!")
}
