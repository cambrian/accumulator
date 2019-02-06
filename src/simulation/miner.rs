use super::state::{Block, Transaction};
use super::util::elems_from_transactions;
use crate::accumulator::Accumulator;
use crate::group::UnknownOrderGroup;
use rug::Integer;

struct Miner<G: UnknownOrderGroup> {
  acc: Accumulator<G>,
  block_height: u64,
  pending_transactions: Vec<Transaction<G>>,
}

#[allow(dead_code)]
impl<G: UnknownOrderGroup> Miner<G> {
  pub fn setup(acc: Accumulator<G>, block_height: u64) -> Self {
    Miner {
      acc,
      block_height,
      pending_transactions: Vec::new(),
    }
  }

  pub fn add_transaction(&mut self, transaction: Transaction<G>) {
    self.pending_transactions.push(transaction);
  }

  pub fn forge_block(&self) -> Block<G> {
    let (elems_added, elems_deleted) = elems_from_transactions(&self.pending_transactions);
    let (witness_deleted, proof_deleted) = self.acc.clone().delete(&elems_deleted).unwrap();
    let (new_acc, proof_added) = witness_deleted.clone().add(&elems_added);
    Block {
      height: self.block_height + 1,
      transactions: self.pending_transactions.clone(),
      new_acc,
      proof_added,
      proof_deleted,
    }
  }

  pub fn validate_block(&mut self, block: Block<G>) {
    let (elems_added, elem_witnesses_deleted) = elems_from_transactions(&block.transactions);
    let elems_deleted: Vec<Integer> = elem_witnesses_deleted
      .iter()
      .map(|(u, _wit)| u.clone())
      .collect();
    assert!(self
      .acc
      .verify_membership(&elems_deleted, &block.proof_deleted));
    assert!(block
      .new_acc
      .verify_membership(&elems_added, &block.proof_added));
    assert!(block.proof_deleted.witness == block.proof_added.witness);
    self.acc = block.new_acc;
    self.block_height = block.height;
    self.pending_transactions.clear();
  }
}
