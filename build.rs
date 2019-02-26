use std::env;
use std::fs::{self, File};
use std::ffi::{OsStr, OsString};
use std::io::{Result as IoResult};
use std::path::{Path, PathBuf};
use std::process::Command;

#[cfg(unix)]
use std::os::unix::fs as unix_fs;
#[cfg(linux)]
use std::os::linux::fs as linux_fs;

// Adapted from gmp-mpfr-sys/build.rs.
const FLINT_DIR: &'static str = "ext/flint-2.5.2";

struct Environment {
    out_dir: PathBuf,
    lib_dir: PathBuf,
    include_dir: PathBuf,
    build_dir: PathBuf
}

fn build_flint(env: &Environment, lib: &Path, header: &Path) {
  let src_dir = env.build_dir.join("flint-src");
  create_dir_or_panic(&src_dir);
  println!("$ cd {:?}", src_dir);

  // For now, we expect GMP and MPFR to be installed where they normally are.
  let mut conf = OsString::from("./configure --with-gmp=/usr --prefix=");
  conf.push(env.out_dir.clone().into_os_string());
  configure(&src_dir, &conf);
  make_and_install(&src_dir);
  panic!("done");
}

fn need_compile(
    flint_ah: &(PathBuf, PathBuf),
) -> bool {
    !(flint_ah.0.is_file() && flint_ah.1.is_file())
}

fn main() {
  let src_dir = PathBuf::from(cargo_env("CARGO_MANIFEST_DIR"));
  let out_dir = PathBuf::from(cargo_env("OUT_DIR"));

  let host = cargo_env("HOST");
  let target = cargo_env("TARGET");
  assert_eq!(host, target, "cross compilation is not supported");

  let env = Environment {
      out_dir: out_dir.clone(),
      lib_dir: out_dir.join("lib"),
      include_dir: out_dir.join("include"),
      build_dir: out_dir.join("build")
  };

  // Make sure we have target directories
  create_dir_or_panic(&env.lib_dir);
  create_dir_or_panic(&env.include_dir);

  // TODO: is the include dir fmpz.h or flint.h?
  let flint_ah = (env.lib_dir.join("libflint.a"), env.include_dir.join("fmpz.h"));

  // If flint needs compilation, compile it.
  if cfg!(feature = "flint") && need_compile(&flint_ah) {
    remove_dir_or_panic(&env.build_dir);
    create_dir_or_panic(&env.build_dir);
    link_dir(&src_dir.join(FLINT_DIR), &env.build_dir.join("flint-src"));
    let (ref a, ref h) = flint_ah;
    build_flint(&env, a, h);
  }

  if cfg!(feature = "flint") {
    println!("cargo:rustc-link-lib=flint");
  }
}

fn copy_file(src: &Path, dst: &Path) -> IoResult<u64> {
    println!("$ cp {:?} {:?}", src, dst);
    fs::copy(src, dst)
}

fn copy_file_or_panic(src: &Path, dst: &Path) {
    copy_file(src, dst).unwrap_or_else(|_| {
        panic!("Unable to copy {:?} -> {:?}", src, dst);
    });
}

fn configure(build_dir: &Path, conf_line: &OsStr) {
    let mut conf = Command::new("sh");
    conf.current_dir(build_dir).arg("-c").arg(conf_line);
    execute(conf);
}

fn execute(mut command: Command) {
    println!("$ {:?}", command);
    let status = command
        .status()
        .unwrap_or_else(|_| panic!("Unable to execute: {:?}", command));
    if !status.success() {
        if let Some(code) = status.code() {
            panic!("Program failed with code {}: {:?}", code, command);
        } else {
            panic!("Program failed: {:?}", command);
        }
    }
}

fn make_and_install(build_dir: &Path) {
    let mut make = Command::new("make");
    make.current_dir(build_dir);
    execute(make);

    let mut make_install = Command::new("make install");
    make_install.current_dir(build_dir);
    execute(make_install);
}

fn create_dir(dir: &Path) -> IoResult<()> {
    println!("$ mkdir -p {:?}", dir);
    fs::create_dir_all(dir)
}

fn create_dir_or_panic(dir: &Path) {
    create_dir(dir).unwrap_or_else(|_| panic!("Unable to create directory: {:?}", dir));
}

fn cargo_env(name: &str) -> OsString {
    env::var_os(name)
        .unwrap_or_else(|| panic!("Environment variable not found: {}, please use cargo.", name))
}

fn remove_dir(dir: &Path) -> IoResult<()> {
    if !dir.exists() {
        return Ok(());
    }
    assert!(dir.is_dir(), "Not a directory: {:?}", dir);
    println!("$ rm -r {:?}", dir);
    fs::remove_dir_all(dir)
}

fn remove_dir_or_panic(dir: &Path) {
    remove_dir(dir).unwrap_or_else(|_| panic!("Unable to remove directory: {:?}", dir));
}

#[cfg(unix)]
fn link_dir(src: &Path, dst: &Path) {
    println!("$ ln -s {:?} {:?}", src, dst);
    unix_fs::symlink(src, dst).unwrap_or_else(|_| {
        panic!("Unable to symlink {:?} -> {:?}", src, dst);
    });
}

#[cfg(linux)]
fn link_dir(src: &Path, dst: &Path) {
    println!("$ ln -s {:?} {:?}", src, dst);
    linux_fs::symlink(src, dst).unwrap_or_else(|_| {
        panic!("Unable to symlink {:?} -> {:?}", src, dst);
    });
}

// NOTES
/*
  - We're gonna need different installation instructions if you want to run Flint
  - Options: include Flint source code in this package or not.
  - I like this approach:
    - Have the source code in the package as an external git submodule. Then we can build and
      install FLINT from source if it's not already there on the system.
    - We're gonna require people to have GMP 5.1.1 or later installed, as well as MPFR and
      an implementation of pthread.
    -
*/