use super::group::Generator;
use alga::general::AbstractGroup;
use alga::general::Operator;

pub fn setup<O, G: AbstractGroup<O> + Generator<O>>() -> G
where
  O: Operator,
{
  G::generator()
}
