use super::U256;
const MAX_TRIAL_DIVISION_DIVISOR: u64 = 100;

#[test]
fn make_uint() {
  let p = U256::from_dec_str(
    "38873241744847760218045702002058062581688990428170398542849190507947196700873",
  )
  .expect("lol");
  println!("{}", p);
}

pub fn is_prob_prime(n: U256) -> bool {
  if !passes_trial_division(n, MAX_TRIAL_DIVISION_DIVISOR) {
    return false;
  }
  true
}

// WIP
fn passes_trial_division(n: U256, max_divisor: u64) -> bool {
  for divisor in 1..=max_divisor {
    let divisor: U256 = divisor.into();
    if divisor > n {
      break;
    }
  }
  true
}

#[test]
fn test_trial_div() {
  let n_prime: U256 = 10000019.into();
  // 3209 * 3659
  let n_composite: U256 = 11741731.into();
  assert!(passes_trial_division(n_prime, 5000));
  assert!(!passes_trial_division(n_composite, 5000));
  assert!(passes_trial_division(n_composite, 3000));
}
