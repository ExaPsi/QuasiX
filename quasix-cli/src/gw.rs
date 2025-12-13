//! GW calculation command implementation
//!
//! This module provides the implementation for GW calculations
//! via the command-line interface.

#![warn(clippy::all, clippy::pedantic, clippy::perf)]
#![warn(missing_docs)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::module_name_repetitions)]

use anyhow::{Context, Result};
use ndarray::{Array1, Array2, Array3};
use quasix_core::{
    common::QuasixError,
    io::{
        hdf5_io,
        schema::{CalculationType, GWResults, QuasixData, SelfEnergyData},
    },
    qp::evgw::{EvGWConfig, EvGWDriver, EvGWResult},
};
use std::path::Path;
use tracing::{debug, error, info, warn};

/// GW calculation parameters from command line
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
pub fn run_gw_calculation(config: GWCommandConfig) -> Result<()> {
    info!("Starting GW calculation");
    debug!("Configuration: {:?}", config);

    // Validate input file exists
    let input_path = Path::new(&config.input_file);
    if !input_path.exists() {
        error!("Input file not found: {}", config.input_file);
        return Err(anyhow::anyhow!(
            "Input file not found: {}",
            config.input_file
        ));
    }

    // Determine output file path
    let output_path = config.output_file.unwrap_or_else(|| {
        let stem = input_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
        format!("{stem}_gw.h5")
    });
    info!("Output will be written to: {}", output_path);

    // Read input data from HDF5
    info!("Reading input data from: {}", config.input_file);
    let input_data = hdf5_io::read_hdf5(input_path)
        .with_context(|| format!("Failed to read input file: {}", config.input_file))?;

    // Validate input data has required fields
    validate_input_data(&input_data)?;

    // Extract molecular orbital data
    let mo_data = extract_mo_data(&input_data)?;
    info!(
        "Loaded molecular system: {} atoms, {} MOs, {} auxiliary basis functions",
        input_data.molecule.natoms, mo_data.n_mo, mo_data.n_aux
    );

    // Configure evGW calculation
    let mut evgw_config = EvGWConfig::default();

    // Apply command-line overrides
    if let Some(max_iter) = config.max_iterations {
        evgw_config.max_iterations = max_iter;
    }
    if let Some(tol) = config.energy_tolerance {
        evgw_config.energy_tolerance = tol;
    }
    if let Some(damping) = config.damping_factor {
        evgw_config.damping_factor = damping;
    }
    if let Some(n_freq) = config.n_frequency_points {
        evgw_config.cd_config.n_imag_points = n_freq;
    }
    evgw_config.use_diis = config.use_diis;
    evgw_config.print_level = config.print_level;

    info!(
        "evGW configuration: max_iter={}, tol={:.2e}, damping={:.2}",
        evgw_config.max_iterations, evgw_config.energy_tolerance, evgw_config.damping_factor
    );

    // Create evGW driver
    let mut driver = EvGWDriver::new(mo_data.n_mo, mo_data.n_aux, mo_data.n_occ, evgw_config);

    // Run evGW calculation
    info!("Starting evGW calculation...");
    let result = match driver.run_evgw(
        &mo_data.orbital_energies,
        &mo_data.occupations,
        &mo_data.mo_coefficients,
        &mo_data.df_integrals,
        &mo_data.coulomb_metric,
        &mo_data.vxc_diagonal,
    ) {
        Ok(res) => res,
        Err(QuasixError::ConvergenceFailed(msg)) => {
            warn!("evGW convergence failed: {}", msg);
            return Err(anyhow::anyhow!("evGW convergence failed: {}", msg));
        }
        Err(e) => {
            error!("evGW calculation error: {}", e);
            return Err(anyhow::anyhow!("evGW calculation failed: {}", e));
        }
    };

    // Report results
    if result.converged {
        info!(
            "evGW converged in {} iterations (wall time: {:.2} s)",
            result.n_iterations, result.wall_time
        );
    } else {
        warn!(
            "evGW did not converge within {} iterations",
            result.n_iterations
        );
    }

    // Print summary of results
    print_results_summary(&result, &mo_data);

    // Prepare output data
    let output_data = prepare_output_data(input_data, &result, &mo_data);

    // Write output to HDF5
    info!("Writing results to: {}", output_path);
    hdf5_io::write_hdf5(&output_data, Path::new(&output_path))
        .with_context(|| format!("Failed to write output file: {output_path}"))?;

    info!("GW calculation completed successfully");
    Ok(())
}

/// Molecular orbital data extracted from input
struct MOData {
    /// Number of molecular orbitals
    n_mo: usize,
    /// Number of auxiliary basis functions
    n_aux: usize,
    /// Number of occupied orbitals
    n_occ: usize,
    /// Orbital energies from DFT/HF
    orbital_energies: Array1<f64>,
    /// Orbital occupations
    occupations: Array1<f64>,
    /// MO coefficients [`n_ao`, `n_mo`]
    mo_coefficients: Array2<f64>,
    /// Density fitting integrals (ia|P)
    df_integrals: Array3<f64>,
    /// Coulomb metric matrix `V_PQ`
    coulomb_metric: Array2<f64>,
    /// Exchange-correlation potential diagonal
    vxc_diagonal: Array1<f64>,
}

/// Validate that input data has all required fields for GW
fn validate_input_data(data: &QuasixData) -> Result<()> {
    // Check for required reference data
    if data.results.reference.orbital_energies.is_empty() {
        return Err(anyhow::anyhow!(
            "Input file missing orbital data. Please provide DFT/HF results."
        ));
    }

    // Check basis set information
    if data.basis.n_ao == 0 || data.basis.n_aux == 0 {
        return Err(anyhow::anyhow!(
            "Invalid basis set information. Both AO and auxiliary basis required."
        ));
    }

    Ok(())
}

/// Extract molecular orbital data from input
fn extract_mo_data(data: &QuasixData) -> Result<MOData> {
    let reference = &data.results.reference;

    // Convert vectors to ndarray
    let n_mo = reference.orbital_energies.len();
    let n_aux = data.basis.n_aux;
    let n_basis_ao = data.basis.n_ao;

    // Count occupied orbitals
    let n_occ = reference
        .occupations
        .iter()
        .filter(|&&occ| occ > 0.5)
        .count();

    // Create arrays from data
    let orbital_energies = Array1::from_vec(reference.orbital_energies.clone());
    let occupations = Array1::from_vec(reference.occupations.clone());

    // Reshape MO coefficients
    let mo_coefficients = if let Some(ref mo_coeff) = reference.mo_coefficients {
        // Flatten nested Vec<Vec<f64>> to Vec<f64>
        let flat_coeffs: Vec<f64> = mo_coeff.iter().flat_map(Clone::clone).collect();
        Array2::from_shape_vec((n_basis_ao, n_mo), flat_coeffs)
            .context("Failed to reshape MO coefficients")?
    } else {
        // Create placeholder MO coefficients if not provided
        warn!("MO coefficients not found in input, using identity-like placeholder");
        Array2::eye(n_basis_ao.min(n_mo))
    };

    // Create placeholder DF integrals array
    // In a real implementation, these would be read from the HDF5 file
    // Shape should be (n_mo, n_mo, n_aux) for three-center integrals in MO basis
    let df_integrals = Array3::zeros((n_mo, n_mo, n_aux));

    // Create placeholder Coulomb metric
    // In a real implementation, this would be computed from auxiliary basis
    let coulomb_metric = Array2::eye(n_aux);

    // Create placeholder Vxc diagonal
    // In a real implementation, this would come from the DFT calculation
    let vxc_diagonal = Array1::zeros(n_mo);

    warn!("Using placeholder data for DF integrals, Coulomb metric, and Vxc. Real implementation would read these from the input file.");

    Ok(MOData {
        n_mo,
        n_aux,
        n_occ,
        orbital_energies,
        occupations,
        mo_coefficients,
        df_integrals,
        coulomb_metric,
        vxc_diagonal,
    })
}

/// Print summary of GW results
fn print_results_summary(result: &EvGWResult, mo_data: &MOData) {
    println!("\n========== evGW Results ==========");
    println!("Converged: {}", result.converged);
    println!("Iterations: {}", result.n_iterations);
    println!("Wall time: {:.2} s", result.wall_time);

    // Find HOMO and LUMO indices
    let homo_idx = mo_data.n_occ.saturating_sub(1);
    let lumo_idx = mo_data.n_occ;

    if homo_idx < mo_data.n_mo && lumo_idx < mo_data.n_mo {
        let homo_dft = mo_data.orbital_energies[homo_idx];
        let lumo_dft = mo_data.orbital_energies[lumo_idx];
        let homo_gw = result.qp_energies[homo_idx];
        let lumo_gw = result.qp_energies[lumo_idx];

        println!("\nFrontier orbitals:");
        println!(
            "  HOMO: DFT = {:.4} Ha, GW = {:.4} Ha, shift = {:.4} Ha ({:.3} eV)",
            homo_dft,
            homo_gw,
            homo_gw - homo_dft,
            (homo_gw - homo_dft) * 27.211
        );
        println!(
            "  LUMO: DFT = {:.4} Ha, GW = {:.4} Ha, shift = {:.4} Ha ({:.3} eV)",
            lumo_dft,
            lumo_gw,
            lumo_gw - lumo_dft,
            (lumo_gw - lumo_dft) * 27.211
        );
        println!(
            "  Gap:  DFT = {:.3} eV, GW = {:.3} eV, change = {:.3} eV",
            (lumo_dft - homo_dft) * 27.211,
            (lumo_gw - homo_gw) * 27.211,
            ((lumo_gw - homo_gw) - (lumo_dft - homo_dft)) * 27.211
        );
    }

    // Diagnostics
    println!("\nDiagnostics:");
    println!(
        "  Z-factor range: [{:.3}, {:.3}]",
        result.diagnostics.z_range.0, result.diagnostics.z_range.1
    );
    println!("  Average Z-factor: {:.3}", result.diagnostics.avg_z_factor);
    println!(
        "  Physical Z-factors: {}/{}",
        result.diagnostics.physical_z_count, mo_data.n_mo
    );
    println!(
        "  Max QP shift: {:.3} eV",
        result.diagnostics.max_qp_shift * 27.211
    );
    println!(
        "  Average QP shift: {:.3} eV",
        result.diagnostics.avg_qp_shift * 27.211
    );

    if result.diagnostics.oscillating_states > 0 {
        warn!(
            "  Warning: {} oscillating states detected",
            result.diagnostics.oscillating_states
        );
    }

    println!("==================================\n");
}

/// Prepare output data structure with GW results
fn prepare_output_data(
    mut input_data: QuasixData,
    result: &EvGWResult,
    mo_data: &MOData,
) -> QuasixData {
    // Update metadata
    input_data.metadata.timestamp = chrono::Utc::now();
    input_data.metadata.quasix_version = quasix_core::VERSION.to_string();

    // Update calculation type
    input_data.parameters.calculation_type = CalculationType::EvGW;

    // Store GW results
    input_data.results.gw = Some(GWResults {
        qp_energies: result.qp_energies.to_vec(),
        z_factors: result.z_factors.to_vec(),
        self_energy: Some(SelfEnergyData {
            sigma_x: result.sigma_x.diag().to_vec(),
            sigma_c_real: result.sigma_c.diag().iter().map(|c| c.re).collect(),
            sigma_c_imag: result.sigma_c.diag().iter().map(|c| c.im).collect(),
            matrix_elements: None,
        }),
        spectral_functions: None,
        converged: result.converged,
        iterations: result.n_iterations,
    });

    // Add convergence history to custom metadata
    input_data.metadata.custom.insert(
        "convergence_history".to_string(),
        serde_json::json!({
            "max_energy_changes": result.convergence_history.max_energy_changes,
            "rms_energy_changes": result.convergence_history.rms_energy_changes,
            "gaps": result.convergence_history.gaps,
            "z_ranges": result.convergence_history.z_ranges,
        }),
    );

    // Update HOMO/LUMO indices in reference results
    input_data.results.reference.homo = mo_data.n_occ.saturating_sub(1);
    input_data.results.reference.lumo = mo_data.n_occ;

    input_data
}
