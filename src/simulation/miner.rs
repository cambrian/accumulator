#![allow(dead_code)]

use super::state::{Block, Transaction};
use crate::accumulator::Accumulator;
use crate::group::UnknownOrderGroup;
use crate::hash::hash_to_prime;
use rug::Integer;

struct Miner<G: UnknownOrderGroup> {
  acc: Accumulator<G>,
  block_height: u64,
  pending_transactions: Vec<Transaction<G>>,
}

impl<G: UnknownOrderGroup> Miner<G> {
  fn setup(acc: Accumulator<G>, block_height: u64) -> Self {
    Miner {
      acc,
      block_height,
      pending_transactions: Vec::new(),
    }
  }

  fn add_transaction(&mut self, transaction: Transaction<G>) {
    self.pending_transactions.push(transaction);
  }

  fn elems_from_transactions(&self) -> (Vec<Integer>, Vec<(Integer, Accumulator<G>)>) {
    let mut elems_added = Vec::new();
    let mut elems_deleted = Vec::new();

    for tx in &self.pending_transactions {
      elems_added.extend(tx.utxos_added.iter().map(|u| hash_to_prime(u)));
      elems_deleted.extend(
        tx.utxos_deleted
          .iter()
          .map(|(u, wit)| (hash_to_prime(u), wit.clone())),
      );
    }

    (elems_added, elems_deleted)
  }

  fn forge_block(&self) -> Block<G> {
    let (elems_added, elems_deleted) = self.elems_from_transactions();
    let (witness_deleted, proof_deleted) = self.acc.clone().delete(&elems_deleted).unwrap();
    let (new_acc, proof_added) = witness_deleted.clone().add(&elems_added);
    Block {
      height: self.block_height + 1,
      _transactions: self.pending_transactions.clone(),
      new_acc,
      proof_added,
      proof_deleted,
    }
  }

  fn validate_block(&mut self, block: Block<G>) {
    let (elems_added, elem_witnesses_deleted) = self.elems_from_transactions();
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
