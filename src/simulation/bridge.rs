use super::state::{Block, Utxo};
use super::util;
use crate::accumulator::Accumulator;
use crate::group::UnknownOrderGroup;
use crate::hash::hash_to_prime;
use crate::util::int;
use rug::Integer;
use std::clone::Clone;

#[allow(dead_code)]
#[derive(Clone)]
pub struct Bridge<G: UnknownOrderGroup> {
  utxo_set_product: Integer,
  utxo_set_witness: Accumulator<G>,
  block_height: u64,
}

#[allow(dead_code)]
impl<G: UnknownOrderGroup> Bridge<G> {
  pub fn setup(acc: Accumulator<G>, block_height: u64) -> Self {
    Bridge {
      utxo_set_product: int(1),
      utxo_set_witness: acc,
      block_height,
    }
  }

  pub fn update(&mut self, block: Block<G>) {
    // Preserves idempotency if multiple miners are leaders.
    if block.height != self.block_height + 1 {
      return;
    }

    let (elems_added, elems_deleted) = util::elems_from_transactions(&block.transactions);
    let elems_added_product: Integer = elems_added.iter().product();
    let elems_deleted_product: Integer = elems_deleted.iter().map(|(u, _wit)| u).product();

    self.utxo_set_product *= elems_added_product;
    self.utxo_set_product /= elems_deleted_product;

    // TODO: Avoid clone.
    self.utxo_set_witness = self
      .utxo_set_witness
      .clone()
      .delete(&elems_deleted)
      .unwrap()
      .0;
    self.utxo_set_witness = self.utxo_set_witness.clone().add(&elems_added).0;
    self.block_height = block.height;
  }

  fn create_aggregate_membership_witness(self, utxos: Vec<Utxo>) -> Accumulator<G> {
    let subproduct: Integer = utxos.iter().map(|u| hash_to_prime(u)).product();
    let update_exponent = self.utxo_set_product / subproduct;
    Accumulator(G::exp(&self.utxo_set_witness.0, &update_exponent))
  }

  /// Slow O(N^2) algorithm for creating individual membership witnesses for several UTXOs.
  /// TODO: Implement O(N log N) RootFactor algorithm from BBF V3 p. 18.
  pub fn create_membership_witnesses(self, utxos: Vec<Utxo>) -> Vec<Accumulator<G>> {
    utxos
      .iter()
      .map(|u| Self::create_aggregate_membership_witness(self.clone(), vec![u.clone()]))
      .collect()
  }
}
