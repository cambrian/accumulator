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

const FLINT_DIR: &'static str = "ext/flint-2.5.2";

struct BuildEnvironment {
    out_dir: PathBuf,
    lib_dir: PathBuf,
    include_dir: PathBuf,
    build_dir: PathBuf,
    headers_dir: PathBuf,
    archive_file: PathBuf,
}

// Adapted from the build script from gmp-mpfr-sys.  We currently only care
// about building Flint.
fn main() {
    if !cfg!(feature = "flint") {
        return;
    }

    let host = cargo_env_or_panic("HOST");
    let target = cargo_env_or_panic("TARGET");
    assert_eq!(host, target, "Cross compilation is not supported with this feature.");

    let flint_src_dir = PathBuf::from(cargo_env_or_panic("CARGO_MANIFEST_DIR")).join(FLINT_DIR);
    let out_dir = PathBuf::from(cargo_env_or_panic("OUT_DIR"));

    let flint_out_dir = out_dir.join("flint");
    let flint_archive_file = flint_out_dir.join("lib").join("libflint.a");
    let flint_headers_dir = flint_out_dir.join("include").join("flint");

    let flint_env = BuildEnvironment {
        out_dir: flint_out_dir.clone(),
        lib_dir: flint_out_dir.join("lib"),
        include_dir: flint_out_dir.join("include"),
        build_dir: flint_out_dir.join("build"),
        archive_file: flint_archive_file,
        headers_dir: flint_headers_dir,
    };

    // Create target directories for Flint.
    create_dir_or_panic(&flint_env.out_dir);
    create_dir_or_panic(&flint_env.lib_dir);
    create_dir_or_panic(&flint_env.include_dir);

    if need_compile(&flint_env) {
        remove_dir_or_panic(&flint_env.build_dir);
        create_dir_or_panic(&flint_env.build_dir);
        link_dir(&flint_src_dir, &flint_env.build_dir.join("flint-src"));
        build_flint(&flint_env);
    }

    write_cargo_cmds(&flint_env);
}

fn build_flint(flint_env: &BuildEnvironment) {
    let src_dir = flint_env.build_dir.join("flint-src");
    create_dir_or_panic(&src_dir);
    println!("$ cd {:?}", src_dir);

    // For now, we expect GMP and MPFR to be installed on the target system. On MACOS they
    // can be installed with brew.  On linux with apt.
    let mut conf = OsString::from("./configure --disable-shared --with-gmp=/usr --prefix=");

    conf.push(flint_env.out_dir.clone().into_os_string());
    configure(&src_dir, &conf);
    make_and_install(&src_dir);
}

fn need_compile(
    env: &BuildEnvironment
) -> bool {
    !(env.archive_file.is_file() && env.headers_dir.join("flint.h").is_file())
}

fn write_cargo_cmds(
    env: &BuildEnvironment,
) {
    let out_str = env.out_dir.to_str().unwrap_or_else(|| {
        panic!(
            "Path contains unsupported characters, can only make {}",
            env.out_dir.display()
        );
    });
    let lib_str = env.lib_dir.to_str().unwrap_or_else(|| {
        panic!(
            "Path contains unsupported characters, can only make {}",
            env.lib_dir.display()
        )
    });
    let include_str = env.include_dir.to_str().unwrap_or_else(|| {
        panic!(
            "Path contains unsupported characters, can only make {}",
            env.include_dir.display()
        )
    });
    println!("cargo:out_dir={}", out_str);
    println!("cargo:lib_dir={}", lib_str);
    println!("cargo:include_dir={}", include_str);
    println!("cargo:rustc-link-search=native={}", lib_str);
    println!("cargo:rustc-link-lib=static=flint");
}

fn make_and_install(build_dir: &Path) {
    let mut make = Command::new("make");
    make.current_dir(build_dir);
    exec(make);

    let mut make_install = Command::new("make");
    make_install.current_dir(build_dir).arg("install");
    exec(make_install);
}

fn configure(build_dir: &Path, conf_line: &OsStr) {
    let mut conf = Command::new("sh");
    conf.current_dir(build_dir).arg("-c").arg(conf_line);
    exec(conf);
}

fn exec(mut command: Command) {
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

// Filesys utility functions

fn copy_file(src: &Path, dst: &Path) -> IoResult<u64> {
    println!("$ cp {:?} {:?}", src, dst);
    fs::copy(src, dst)
}

fn copy_file_or_panic(src: &Path, dst: &Path) {
    copy_file(src, dst).unwrap_or_else(|_| {
        panic!("Unable to copy {:?} -> {:?}", src, dst);
    });
}

fn cargo_env_or_panic(name: &str) -> OsString {
    env::var_os(name)
        .unwrap_or_else(|| panic!("Environment variable not found: {}, please use cargo.", name))
}

fn create_dir(dir: &Path) -> IoResult<()> {
    println!("$ mkdir -p {:?}", dir);
    fs::create_dir_all(dir)
}

fn create_dir_or_panic(dir: &Path) {
    create_dir(dir).unwrap_or_else(|_| panic!("Unable to create directory: {:?}", dir));
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
