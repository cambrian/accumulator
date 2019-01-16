// TODO

// x, u, w: u^x = w
// pub fn compute_poe<G: Group> (G base, U256 exp, G result) -> PoE<G> {
//   let l = Hprime (exp, base, result);
//   let q = (exp / l);
//   PoE<G>(base^q)
// }

// pub fn verify_poe<G: Group> (G base, U256 exp, G result, PoE<G> proof) -> Bool{

// }

// impl<T> f
// Prove (x, u, w)
// l <- Hprime (x, u, w)
// q <- Int divide x / l
// r <- x mod l
// return PoE <- u^q

// Verify (x, u, w, Q)
// l <- Hprime (x, u, w)
// r <- x mod l
// Check Q^lu^r = w
// 1 if yes, else 0
