//! GW calculation command implementation
//!
//! This module provides the implementation for GW calculations
//! via the command-line interface.
//!
//! Note: CLI is currently a placeholder pending API stabilization.
//! Use the Python bindings for GW calculations.

use anyhow::Result;

/// GW calculation parameters from command line
#[allow(dead_code)] // Fields used when CLI is fully implemented
#[derive(Debug, Clone)]
pub struct GWCommandConfig {
    /// Input HDF5 file path
    pub input_file: String,
    /// Output HDF5 file path (defaults to `input_gw.h5`)
    pub output_file: Option<String>,
    /// Maximum iterations for evGW
    pub max_iterations: Option<usize>,
    /// Energy convergence threshold in Ha
    pub energy_tolerance: Option<f64>,
    /// Damping factor for self-consistency
    pub damping_factor: Option<f64>,
    /// Number of frequency points for contour deformation
    pub n_frequency_points: Option<usize>,
    /// Use DIIS acceleration
    pub use_diis: bool,
    /// Print level (0=silent, 1=summary, 2=detailed, 3=debug)
    pub print_level: u32,
}

/// Run GW calculation from command-line parameters
///
/// # Note
/// This CLI interface is currently a placeholder. The quasix_core API is under
/// active development as part of the M1 milestone. For production GW calculations,
/// please use the Python bindings (`quasix.evgw`) which provide a stable interface.
///
/// # Errors
/// Currently returns an error indicating CLI is not yet available.
pub fn run_gw_calculation(_config: GWCommandConfig) -> Result<()> {
    eprintln!("╔════════════════════════════════════════════════════════════════╗");
    eprintln!("║                    QuasiX GW CLI Notice                        ║");
    eprintln!("╠════════════════════════════════════════════════════════════════╣");
    eprintln!("║ The CLI interface is currently under development.              ║");
    eprintln!("║                                                                ║");
    eprintln!("║ For GW calculations, please use the Python interface:          ║");
    eprintln!("║                                                                ║");
    eprintln!("║   from quasix import evgw                                      ║");
    eprintln!("║   result = evgw.run_evgw(mol, mf, ...)                         ║");
    eprintln!("║                                                                ║");
    eprintln!("║ See: https://github.com/ExaPsi/QuasiX for documentation.       ║");
    eprintln!("╚════════════════════════════════════════════════════════════════╝");

    Err(anyhow::anyhow!(
        "CLI GW calculation not yet implemented. Use Python bindings instead."
    ))
}
