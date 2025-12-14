//! Utility functions for the QuasiX CLI

use tracing_subscriber::{EnvFilter, FmtSubscriber};

/// Print the QuasiX banner
#[allow(dead_code)]
pub fn print_banner(title: &str) {
    println!();
    println!("╔═══════════════════════════════════════════════════════════╗");
    println!("║                       QuasiX                              ║");
    println!("║       High-Performance GW/BSE Implementation              ║");
    println!("╠═══════════════════════════════════════════════════════════╣");
    println!("║  {}  ║", center_text(title, 57));
    println!("╚═══════════════════════════════════════════════════════════╝");
    println!();
}

#[allow(dead_code)]
fn center_text(text: &str, width: usize) -> String {
    if text.len() >= width {
        text[..width].to_string()
    } else {
        let _padding = (width - text.len()) / 2;
        format!("{:^width$}", text, width = width)
    }
}

/// Print a success message
#[allow(dead_code)]
pub fn print_success(message: &str) {
    println!("✅ {}", message);
}

/// Print a warning message
#[allow(dead_code)]
pub fn print_warning(message: &str) {
    println!("⚠️  {}", message);
}

/// Print an error message
#[allow(dead_code)]
pub fn print_error(message: &str) {
    eprintln!("❌ {}", message);
}

/// Setup the logger based on verbosity level
#[allow(dead_code)]
pub fn setup_logger(verbosity: u8) {
    let level = match verbosity {
        0 => "error",
        1 => "warn",
        2 => "info",
        3 => "debug",
        _ => "trace",
    };

    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(level));

    let subscriber = FmtSubscriber::builder()
        .with_env_filter(filter)
        .with_target(false)
        .with_thread_ids(false)
        .with_thread_names(false)
        .finish();

    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");
}
