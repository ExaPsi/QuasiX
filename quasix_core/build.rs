use std::env;
use std::path::PathBuf;

fn main() {
    // Only configure libcint if the feature is enabled
    if cfg!(feature = "real_libcint") {
        configure_libcint();
    }
}

fn configure_libcint() {
    // Try to find libcint using multiple methods

    // Method 1: Check CINT_DIR environment variable
    if let Ok(cint_dir) = env::var("CINT_DIR") {
        // Try both "build" and "lib" subdirectories
        let build_dir = PathBuf::from(&cint_dir).join("build");
        let lib_dir = PathBuf::from(&cint_dir).join("lib");

        // Check build directory first (for in-source builds)
        if build_dir.join("libcint.so").exists() {
            println!("cargo:rustc-link-search=native={}", build_dir.display());
            println!("cargo:rustc-link-lib=cint");
            println!("cargo:rustc-link-arg=-Wl,-rpath,{}", build_dir.display());
            println!("cargo:warning=Found libcint in {}", build_dir.display());
            return;
        }

        // Fall back to lib directory (for installed builds)
        if lib_dir.join("libcint.so").exists() {
            println!("cargo:rustc-link-search=native={}", lib_dir.display());
            println!("cargo:rustc-link-lib=cint");
            println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_dir.display());
            println!("cargo:warning=Found libcint in {}", lib_dir.display());
            return;
        }
    }

    // Method 2: Check common installation paths
    let home = env::var("HOME").unwrap_or_else(|_| String::from("/home/vyv"));
    let common_paths = vec![
        format!("{}/.local", home),
        String::from("/usr/local"),
        String::from("/usr"),
        String::from("/opt/libcint"),
    ];

    for path in common_paths {
        let lib_path = PathBuf::from(&path).join("lib");
        let lib_file = lib_path.join("libcint.so");

        if lib_file.exists() {
            println!("cargo:rustc-link-search=native={}", lib_path.display());
            println!("cargo:rustc-link-lib=cint");
            println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib_path.display());
            println!("cargo:warning=Found libcint at {}", lib_path.display());
            return;
        }
    }

    // If we can't find libcint, warn but continue
    // This allows tests without real_libcint to pass
    println!("cargo:warning=libcint not found in standard locations");
    println!("cargo:warning=Set CINT_DIR to point to libcint installation");
    println!("cargo:warning=Or install libcint to $HOME/.local");
}
