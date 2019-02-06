use super::state::Utxo;
use crate::accumulator::Accumulator;
use crate::group::UnknownOrderGroup;
use crate::util::int;
use rug::Integer;

#[allow(dead_code)]
struct Bridge<G: UnknownOrderGroup> {
  utxo_set: Vec<Utxo>,
  utxo_set_product: Integer,
  utxo_set_witness: Accumulator<G>,
}

impl<G: UnknownOrderGroup> Bridge<G> {
  #[allow(dead_code)]
  pub fn setup(acc: Accumulator<G>) -> Self {
    Bridge {
      utxo_set: Vec::new(),
      utxo_set_product: int(1),
      utxo_set_witness: acc,
    }
  }
}
