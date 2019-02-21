fn main() {
  if cfg!(feature = "flint") {
    println!("cargo:rustc-link-lib=flint");
  }
}