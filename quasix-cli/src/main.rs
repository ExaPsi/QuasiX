// Clippy configuration - enforce important lints, allow some pedantic ones
#![warn(clippy::all)]
#![warn(clippy::correctness)]
#![warn(clippy::suspicious)]
#![warn(clippy::complexity)]
#![warn(clippy::perf)]
#![warn(clippy::style)]
// Allow some pedantic lints that are too strict
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::uninlined_format_args)]
#![allow(clippy::unnecessary_wraps)]
#![allow(clippy::upper_case_acronyms)] // BSE, GW are standard acronyms in quantum chemistry

use anyhow::Result;
use clap::{Parser, Subcommand};
use quasix_core::logging;
use tracing::{debug, error, info, warn};

mod commands;
mod gw;
mod utils;

#[derive(Parser)]
#[command(
    name = "quasix",
    version,
    about = "QuasiX: High-performance GW/BSE calculations for molecules and periodic systems",
    long_about = "QuasiX is a high-performance implementation of the GW approximation with \
                  Bethe-Salpeter Equation (BSE) for accurate quasiparticle and optical \
                  excitation calculations.",
    author = "QuasiX Development Team"
)]
struct Cli {
    /// Set log level (TRACE, DEBUG, INFO, WARN, ERROR)
    #[arg(long, env = "QUASIX_LOG", global = true)]
    log_level: Option<String>,

    /// Set log format (json, pretty, compact)
    #[arg(long, env = "QUASIX_LOG_FORMAT", global = true, default_value = "json")]
    log_format: Option<String>,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Perform GW calculations
    #[command(name = "gw")]
    GW {
        /// Input HDF5 file containing molecular data from PySCF
        #[arg(
            value_name = "INPUT",
            help = "Path to input HDF5 file with DFT/HF results"
        )]
        input: String,

        /// Output HDF5 file path (default: <input>_gw.h5)
        #[arg(short = 'o', long, value_name = "OUTPUT")]
        output: Option<String>,

        /// Maximum number of evGW iterations
        #[arg(
            long,
            value_name = "N",
            help = "Maximum iterations for self-consistency"
        )]
        max_iterations: Option<usize>,

        /// Energy convergence threshold in Hartree
        #[arg(
            long,
            value_name = "TOL",
            help = "Energy convergence threshold (default: 1e-4 Ha)"
        )]
        energy_tolerance: Option<f64>,

        /// Damping factor for self-consistency (0 < α < 1)
        #[arg(long, value_name = "ALPHA", help = "Damping factor (default: 0.5)")]
        damping: Option<f64>,

        /// Number of frequency points for contour deformation
        #[arg(long, value_name = "N", help = "Number of frequency points")]
        n_frequency: Option<usize>,

        /// Use DIIS acceleration
        #[arg(long, help = "Enable DIIS acceleration for convergence")]
        diis: bool,

        /// Verbosity level (0=silent, 1=summary, 2=detailed, 3=debug)
        #[arg(short = 'v', long, value_name = "LEVEL", default_value = "1")]
        verbose: u32,
    },

    /// Compare CD and AC methods for validation
    #[command(name = "compare")]
    Compare(commands::CompareCommand),

    /// Perform BSE calculations (placeholder for future implementation)
    #[command(name = "bse")]
    BSE {
        /// Input file path
        #[arg(value_name = "INPUT")]
        input: Option<String>,
    },

    /// Show version information
    #[command(name = "version")]
    Version,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Set up logging environment before initializing
    if let Some(ref level) = cli.log_level {
        std::env::set_var("QUASIX_LOG", level);
    }
    if let Some(ref format) = cli.log_format {
        std::env::set_var("QUASIX_LOG_FORMAT", format);
    }

    // Initialize logging
    logging::init_logger()?;

    info!("QuasiX CLI started");
    debug!("Command line arguments processed");

    match cli.command {
        Some(Commands::GW {
            input,
            output,
            max_iterations,
            energy_tolerance,
            damping,
            n_frequency,
            diis,
            verbose,
        }) => {
            info!("Starting GW calculation");
            debug!("Input file: {}", input);

            // Build GW command configuration
            let gw_config = gw::GWCommandConfig {
                input_file: input,
                output_file: output,
                max_iterations,
                energy_tolerance,
                damping_factor: damping,
                n_frequency_points: n_frequency,
                use_diis: diis,
                print_level: verbose,
            };

            // Run GW calculation
            match gw::run_gw_calculation(gw_config) {
                Ok(()) => {
                    info!("GW calculation completed successfully");
                }
                Err(e) => {
                    error!("GW calculation failed: {}", e);
                    eprintln!("Error: {:#}", e);
                    std::process::exit(1);
                }
            }
        }
        Some(Commands::Compare(cmd)) => {
            info!("Starting CD vs AC comparison");
            match cmd.execute() {
                Ok(()) => {
                    info!("Comparison completed successfully");
                }
                Err(e) => {
                    error!("Comparison failed: {}", e);
                    eprintln!("Error: {:#}", e);
                    std::process::exit(1);
                }
            }
        }
        Some(Commands::BSE { input }) => {
            warn!("BSE calculations not yet implemented. Coming in Sprint 6.");
            if let Some(path) = input {
                debug!("Input file specified: {}", path);
            }
            eprintln!("BSE calculations not yet implemented. Coming in Sprint 6.");
            std::process::exit(1);
        }
        Some(Commands::Version) => {
            info!("Version information requested");
            println!("QuasiX version {}", env!("CARGO_PKG_VERSION"));
            println!("Core library version: {}", quasix_core::VERSION);
            debug!("Version display complete");
        }
        None => {
            // When no subcommand is provided, show brief help
            info!("No subcommand provided, showing help");
            println!("QuasiX v{}", env!("CARGO_PKG_VERSION"));
            println!("High-performance GW/BSE calculations");
            println!();
            println!("Run 'quasix --help' for usage information.");
        }
    }

    Ok(())
}
