//! Debug test for W screening calculation using H2/STO-3G reference data
//!
//! This test loads PySCF reference data and runs W calculation with extensive
//! diagnostics to identify the source of the 3.7× over-normalization bug.

use ndarray::Array2;
use num_complex::Complex64;
use std::env;

#[test]
fn test_screening_h2_sto3g_debug() {
    // Set debug environment variable
    env::set_var("QUASIX_DEBUG_SCREENING", "1");

    // Load H2/STO-3G reference data from HDF5
    let h5_path = "/home/vyv/Working/QuasiX/pyscf_s3_3_W_references/H2_sto-3g_W_reference.h5";

    println!("\n{}", "=".repeat(80));
    println!("H2/STO-3G W Screening Debug Test");
    println!("{}", "=".repeat(80));
    println!("Loading reference data from: {}", h5_path);

    // Load reference data using hdf5 crate
    let file = match hdf5::File::open(h5_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("ERROR: Cannot open HDF5 file: {}", e);
            eprintln!("Skipping test (reference data not available)");
            return;
        }
    };

    // Load matrices - they're stored as real, convert to complex where needed
    let p0_ref_real: Array2<f64> = match load_real_matrix(&file, "P0_matrix") {
        Ok(m) => m,
        Err(e) => {
            eprintln!("ERROR loading P0: {}", e);
            return;
        }
    };
    let p0_ref = p0_ref_real.mapv(|x| Complex64::new(x, 0.0));

    let v_sqrt_ref: Array2<f64> = match load_real_matrix(&file, "v_sqrt") {
        Ok(m) => m,
        Err(e) => {
            eprintln!("ERROR loading v_sqrt: {}", e);
            return;
        }
    };

    let w_ref_real: Array2<f64> = match load_real_matrix(&file, "W_matrix") {
        Ok(m) => m,
        Err(e) => {
            eprintln!("ERROR loading W reference: {}", e);
            return;
        }
    };
    let w_ref = w_ref_real.mapv(|x| Complex64::new(x, 0.0));

    let naux = p0_ref.nrows();
    println!("Loaded reference data: naux = {}", naux);

    // Print reference values
    println!("\nReference (PySCF) values:");
    let w_ref_trace: f64 = (0..naux).map(|i| w_ref[[i, i]].re).sum();
    println!("  W trace = {:.4} Ha", w_ref_trace);
    println!("  W[0,0]  = {:.4} Ha", w_ref[[0, 0]].re);
    println!("  W[1,1]  = {:.4} Ha", w_ref[[1, 1]].re);

    // Create DielectricSolver
    use quasix_core::dielectric::screening::{DielectricSolver, SolverType, SolverBackend};

    let solver = DielectricSolver::with_backend(
        naux,
        SolverType::Direct,
        SolverBackend::LU,
    );

    println!("\nRunning QuasiX W calculation with DEBUG output:");
    println!("{}", "-".repeat(80));

    // Compute W using QuasiX implementation
    let w_quasix = match solver.compute_screened_interaction(&p0_ref, &v_sqrt_ref) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("ERROR computing W: {}", e);
            panic!("W calculation failed");
        }
    };

    println!("{}", "-".repeat(80));

    // Compare results
    let w_quasix_trace: f64 = (0..naux).map(|i| w_quasix[[i, i]].re).sum();
    println!("\nQuasiX results:");
    println!("  W trace = {:.4} Ha", w_quasix_trace);
    println!("  W[0,0]  = {:.4} Ha", w_quasix[[0, 0]].re);
    println!("  W[1,1]  = {:.4} Ha", w_quasix[[1, 1]].re);

    println!("\nComparison:");
    let trace_diff = (w_quasix_trace - w_ref_trace).abs();
    let trace_ratio = w_quasix_trace / w_ref_trace;
    println!("  Trace difference: {:.4} Ha", trace_diff);
    println!("  Trace ratio (QuasiX/PySCF): {:.4}x", trace_ratio);

    // Element-wise comparison
    let mut max_abs_diff = 0.0;
    let mut max_rel_diff = 0.0;
    for i in 0..naux {
        for j in 0..naux {
            let diff = (w_quasix[[i, j]] - w_ref[[i, j]]).norm();
            let rel = if w_ref[[i, j]].norm() > 1e-10 {
                diff / w_ref[[i, j]].norm()
            } else {
                diff
            };

            if diff > max_abs_diff {
                max_abs_diff = diff;
            }
            if rel > max_rel_diff {
                max_rel_diff = rel;
            }
        }
    }

    println!("  Max absolute difference: {:.2e}", max_abs_diff);
    println!("  Max relative difference: {:.2e}", max_rel_diff);

    // Check tolerance
    if max_rel_diff > 0.01 {
        println!("\n❌ FAIL: W differs from PySCF by {:.1}%", max_rel_diff * 100.0);
        println!("Expected: W trace ≈ {:.2} Ha", w_ref_trace);
        println!("Got:      W trace = {:.2} Ha ({:.2}× larger)", w_quasix_trace, trace_ratio);
    } else {
        println!("\n✓ PASS: W matches PySCF within 1%");
    }

    println!("{}", "=".repeat(80));

    // Clean up environment
    env::remove_var("QUASIX_DEBUG_SCREENING");
}

/// Load a complex matrix from HDF5 file
fn load_complex_matrix(file: &hdf5::File, name: &str) -> Result<Array2<Complex64>, String> {
    let dataset = file
        .dataset(name)
        .map_err(|e| format!("Dataset '{}' not found: {}", name, e))?;

    // Read as f64 array (real and imaginary interleaved)
    let data: Vec<f64> = dataset
        .read_raw()
        .map_err(|e| format!("Failed to read dataset: {}", e))?;

    let shape = dataset.shape();
    if shape.len() != 2 {
        return Err(format!("Expected 2D array, got {}D", shape.len()));
    }

    let nrows = shape[0];
    let ncols = shape[1];

    // Complex data is stored as [real_00, imag_00, real_01, imag_01, ...]
    let mut matrix = Array2::<Complex64>::zeros((nrows, ncols));
    let mut idx = 0;
    for i in 0..nrows {
        for j in 0..ncols {
            let re = data[idx];
            let im = data[idx + 1];
            matrix[[i, j]] = Complex64::new(re, im);
            idx += 2;
        }
    }

    Ok(matrix)
}

/// Load a real matrix from HDF5 file
fn load_real_matrix(file: &hdf5::File, name: &str) -> Result<Array2<f64>, String> {
    let dataset = file
        .dataset(name)
        .map_err(|e| format!("Dataset '{}' not found: {}", name, e))?;

    let data: Vec<f64> = dataset
        .read_raw()
        .map_err(|e| format!("Failed to read dataset: {}", e))?;

    let shape = dataset.shape();
    if shape.len() != 2 {
        return Err(format!("Expected 2D array, got {}D", shape.len()));
    }

    let nrows = shape[0];
    let ncols = shape[1];

    let matrix = Array2::from_shape_vec((nrows, ncols), data)
        .map_err(|e| format!("Failed to reshape: {}", e))?;

    Ok(matrix)
}
