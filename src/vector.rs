use super::accumulator;
use super::group::Group;
use super::proof::poe::PoE;
use bitvec::BitVec;
use num::BigUint;

pub enum UpdateError {
  IndexOutOfBounds,
  AccError,
}

pub enum UpdateResult<G: Group> {
  NoChange(G::Elem),
  Update(G::Elem, PoE<G>),
}

// struct T<G: Group> {
//   acc: G::Elem,
// }

// impl T {
//   pub fn setup<G: Group>() -> G::Elem {}
// }

pub fn setup<G: Group>() -> G::Elem {
  G::base_elem()
}

pub fn commit<G: Group>(m: BitVec) -> G::Elem {
  G::base_elem()
}

pub fn update<G: Group>(b: bool, i: BigUint) -> Result<UpdateResult<G>, UpdateError> {
  Err(UpdateError::AccError)
}
