use super::state::Utxo;
use std::collections::HashSet;

pub struct User(HashSet<Utxo>);

#[allow(dead_code)]
impl User {
  pub fn setup(utxo_set: HashSet<Utxo>) -> Self {
    User(utxo_set)
  }

  // TODO: Maybe support more inputs than one.
  // Expects executable to call `update` to remove this UTXO when it is confirmed.
  pub fn get_input_for_transaction(&mut self) -> Utxo {
    self.0.iter().next().unwrap().clone()
  }

  pub fn update(&mut self, deleted_inputs: &[Utxo], added_outputs: &[Utxo]) {
    for del in deleted_inputs {
      self.0.remove(&del);
    }
    for add in added_outputs {
      self.0.insert(add.clone());
    }
  }
}
