// Syntactic sugar for cloning elements into closures.
#[macro_export]
macro_rules! enclose {
  ( ($( $x:ident ),*) $y:expr ) => {
    {
      $(let $x = $x.clone();)*
      $y
    }
  };
}