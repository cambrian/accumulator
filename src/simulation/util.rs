use super::state::Transaction;
use crate::accumulator::Accumulator;
use crate::group::UnknownOrderGroup;
use crate::hash::hash_to_prime;
use rug::Integer;

pub fn elems_from_transactions<G: UnknownOrderGroup>(
  transactions: &[Transaction<G>],
) -> (Vec<Integer>, Vec<(Integer, Accumulator<G>)>) {
  let mut elems_added = Vec::new();
  let mut elems_deleted = Vec::new();

  for tx in transactions {
    elems_added.extend(tx.utxos_added.iter().map(|u| hash_to_prime(u)));
    elems_deleted.extend(
      tx.utxos_deleted
        .iter()
        .map(|(u, wit)| (hash_to_prime(u), wit.clone())),
    );
  }

  (elems_added, elems_deleted)
}
