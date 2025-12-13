//! Example of using the evGW integration with GW100 benchmark
//!
//! This example demonstrates how to use the real evGW calculation
//! implementation for GW100 benchmarking instead of mock results.

use quasix_core::validation::gw100::{
    run_evgw_calculation, BenchmarkExecutor, BenchmarkTask, BenchmarkResult,
};
use quasix_core::common::Result;
use ndarray::{Array1, Array2, Array3};

/// Example function that would load molecular data from PySCF or HDF5
fn load_molecular_data(
    molecule: &str,
    basis_set: &str,
) -> Result<(
    Array1<f64>,  // mo_energies
    Array1<f64>,  // mo_occ
    Array3<f64>,  // ia_occ_virt
    Array3<f64>,  // ij_occ_occ
    Array2<f64>,  // chol_v
    Array1<f64>,  // vxc_dft
    Vec<f64>,     // reference_energies
)> {
    // In a real implementation, this would:
    // 1. Load precomputed DFT results from HDF5 files
    // 2. Or interface with PySCF to compute them on-the-fly
    // 3. Load GW100 reference values from database

    // For this example, return mock data
    let nbasis = 20;
    let nocc = 10;
    let nvirt = 10;
    let naux = 100;

    let mo_energies = Array1::linspace(-20.0, 10.0, nbasis);
    let mut mo_occ = Array1::zeros(nbasis);
    mo_occ.slice_mut(ndarray::s![..nocc]).fill(2.0);

    let ia_occ_virt = Array3::from_elem((nocc, nvirt, naux), 0.1);
    let ij_occ_occ = Array3::from_elem((nocc, nocc, naux), 0.1);
    let chol_v = Array2::eye(naux) * 0.5;
    let vxc_dft = Array1::from_elem(nbasis, -0.3);

    // Mock reference energies from GW100 database
    let reference_energies = mo_energies.iter().map(|&e| e * 27.211386245988).collect();

    Ok((mo_energies, mo_occ, ia_occ_virt, ij_occ_occ, chol_v, vxc_dft, reference_energies))
}

/// Create an executor function that uses real evGW calculations
fn create_real_evgw_executor() -> impl Fn(&BenchmarkTask) -> Result<BenchmarkResult> {
    move |task: &BenchmarkTask| {
        // Load molecular data for this task
        let (mo_energies, mo_occ, ia_occ_virt, ij_occ_occ, chol_v, vxc_dft, reference_energies) =
            load_molecular_data(&task.molecule, &task.basis_set)?;

        // Run the actual evGW calculation
        run_evgw_calculation(
            task,
            &mo_energies,
            &mo_occ,
            &ia_occ_virt,
            &ij_occ_occ,
            &chol_v,
            &vxc_dft,
            reference_energies,
        )
    }
}

fn main() -> Result<()> {
    // Initialize logging
    env_logger::init();

    // Define benchmark tasks
    let tasks = vec![
        BenchmarkTask {
            molecule: "H2O".to_string(),
            basis_set: "def2-svp".to_string(),
            aux_basis: "def2-svp-jkfit".to_string(),
            conv_tol: 1e-4,
            max_iterations: 20,
            freq_integration: "cd".to_string(),
            n_frequencies: 32,
        },
        BenchmarkTask {
            molecule: "NH3".to_string(),
            basis_set: "def2-svp".to_string(),
            aux_basis: "def2-svp-jkfit".to_string(),
            conv_tol: 1e-4,
            max_iterations: 20,
            freq_integration: "cd".to_string(),
            n_frequencies: 32,
        },
        BenchmarkTask {
            molecule: "CO".to_string(),
            basis_set: "def2-tzvp".to_string(),
            aux_basis: "def2-tzvp-jkfit".to_string(),
            conv_tol: 1e-5,
            max_iterations: 30,
            freq_integration: "ac".to_string(),
            n_frequencies: 48,
        },
    ];

    // Create executor with 4 threads
    let executor = BenchmarkExecutor::new(4);

    // Use the real evGW executor instead of mock
    let executor_fn = create_real_evgw_executor();

    // Run benchmarks in parallel
    println!("Running evGW calculations for {} molecules...", tasks.len());
    let results = executor.run_parallel(tasks, executor_fn);

    // Print results
    println!("\nResults:");
    println!("{:-<80}", "");
    for result in &results {
        println!("Molecule: {}", result.molecule);
        println!("  Basis: {}", result.basis_set);
        println!("  Converged: {}", result.converged);
        if let Some(iterations) = result.convergence_iterations {
            println!("  Iterations: {}", iterations);
        }
        println!("  IP error: {:.3} eV", result.ip_error);
        println!("  EA error: {:.3} eV", result.ea_error);
        println!("  MAD: {:.3} eV", result.mad());
        println!("  Wall time: {:.2} s", result.wall_time);

        // Check Z-factors are physical
        let z_valid = result.z_factors.iter().all(|&z| z > 0.0 && z <= 1.0);
        println!("  Z-factors valid: {}", z_valid);
        println!("{:-<80}", "");
    }

    // Compute overall statistics
    use quasix_core::validation::gw100::compute_validation_stats;
    let stats = compute_validation_stats(&results)?;
    println!("\nOverall Statistics:");
    println!("  MAD: {:.3} eV", stats.mad);
    println!("  RMSE: {:.3} eV", stats.rmse);
    println!("  R²: {:.4}", stats.correlation);
    println!("  Max error: {:.3} eV", stats.max_error);

    // Check if validation passes GW100 thresholds
    let mad_threshold = 0.3; // ~0.3 eV is typical for good evGW
    let r2_threshold = 0.95;

    if stats.passes_validation(mad_threshold, r2_threshold) {
        println!("\n✓ Validation PASSED GW100 benchmarks!");
    } else {
        println!("\n✗ Validation did not meet GW100 thresholds");
        println!("  Expected: MAD < {:.2} eV, R² > {:.2}", mad_threshold, r2_threshold);
    }

    Ok(())
}