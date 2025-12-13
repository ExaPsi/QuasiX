// Demonstration of SIMD-optimized W(ω) interpolation performance
//
// This example compares the performance of standard vs SIMD-optimized
// W interpolation for typical GW calculations.

use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;
use quasix_core::selfenergy::{
    WScreenedStorage, WScreenedStorageSimd, InterpolationMethod,
    EnergyDependentSelfEnergy, EnergyDependentSelfEnergySimd,
    ContourDeformationConfig, DFTensorBlocked,
    benchmark::{compare_implementations, print_comparison},
};
use std::time::Instant;

/// Generate a test W matrix with typical properties
fn generate_hermitian_w(n_aux: usize, decay_factor: f64) -> Array2<Complex64> {
    let mut w = Array2::<Complex64>::zeros((n_aux, n_aux));

    for i in 0..n_aux {
        for j in i..n_aux {
            let distance = ((i as f64 - j as f64).abs() + 1.0);
            let value = Complex64::new(
                (-decay_factor * distance).exp(),
                0.1 * (-decay_factor * distance * 2.0).exp(),
            );
            w[[i, j]] = value;
            if i != j {
                w[[j, i]] = value.conj();
            }
        }
    }

    w
}

/// Generate typical GW100 test system
fn generate_test_system(n_mo: usize, n_occ: usize, n_aux: usize, n_freq: usize)
    -> (Vec<Array2<Complex64>>, Array1<f64>, Array1<f64>, Array1<f64>, Array3<f64>)
{
    // Generate W matrices on imaginary frequency grid
    let mut w_matrices = Vec::new();
    for freq_idx in 0..n_freq {
        let decay = 0.1 + 0.05 * freq_idx as f64;
        w_matrices.push(generate_hermitian_w(n_aux, decay));
    }

    // Generate frequency grid (typical for GW)
    let xi_grid = Array1::linspace(0.1, 100.0, n_freq);

    // Generate MO energies (occupied: negative, virtual: positive)
    let mut mo_energy = Array1::<f64>::zeros(n_mo);
    for i in 0..n_occ {
        mo_energy[i] = -20.0 + i as f64 * 2.0; // Occupied
    }
    for i in n_occ..n_mo {
        mo_energy[i] = 1.0 + (i - n_occ) as f64 * 3.0; // Virtual
    }

    // Generate MO occupations
    let mut mo_occ = Array1::<f64>::zeros(n_mo);
    for i in 0..n_occ {
        mo_occ[i] = 2.0; // Doubly occupied
    }

    // Generate DF tensor (simplified)
    let mut df_tensor = Array3::<f64>::zeros((n_mo, n_mo, n_aux));
    for i in 0..n_mo {
        for j in 0..n_mo {
            for p in 0..n_aux {
                let value = (-(((i as f64 - j as f64).powi(2) +
                              (p as f64 / 10.0).powi(2)) / 100.0)).exp();
                df_tensor[[i, j, p]] = value * 0.1;
            }
        }
    }

    (w_matrices, xi_grid, mo_energy, mo_occ, df_tensor)
}

fn main() {
    println!("===========================================");
    println!("QuasiX W(ω) Interpolation Performance Demo");
    println!("===========================================\n");

    // Test different system sizes
    let test_cases = vec![
        ("Small molecule (H2O-like)", 20, 10, 50, 20),
        ("Medium molecule (Benzene-like)", 60, 30, 200, 30),
        ("Large molecule (Fullerene-like)", 120, 60, 400, 40),
    ];

    for (name, n_mo, n_occ, n_aux, n_freq) in test_cases {
        println!("\n### {} ###", name);
        println!("MOs: {}, Occupied: {}, Aux basis: {}, Freq points: {}",
                 n_mo, n_occ, n_aux, n_freq);
        println!("-------------------------------------------");

        // Generate test system
        let (w_matrices, xi_grid, mo_energy, mo_occ, df_tensor) =
            generate_test_system(n_mo, n_occ, n_aux, n_freq);

        // Test interpolation at multiple frequencies
        let test_frequencies = vec![
            Complex64::new(0.0, 1.0),   // Low imaginary frequency
            Complex64::new(0.0, 10.0),  // Medium imaginary frequency
            Complex64::new(0.0, 50.0),  // High imaginary frequency
        ];

        // Compare standard vs SIMD for interpolation
        println!("\n1. W Interpolation Performance:");

        let start = Instant::now();
        let w_storage_std = WScreenedStorage::new(
            w_matrices.clone(),
            xi_grid.clone(),
            InterpolationMethod::Cubic,
        ).unwrap();
        let setup_time_std = start.elapsed();

        let start = Instant::now();
        for &omega in &test_frequencies {
            let _ = w_storage_std.evaluate(omega);
        }
        let interp_time_std = start.elapsed();

        let start = Instant::now();
        let w_storage_simd = WScreenedStorageSimd::from_matrices(
            w_matrices.clone(),
            xi_grid.clone(),
            None,
        ).unwrap();
        let setup_time_simd = start.elapsed();

        let start = Instant::now();
        for &omega in &test_frequencies {
            let _ = w_storage_simd.cubic_interpolation_simd(omega);
        }
        let interp_time_simd = start.elapsed();

        println!("  Standard implementation:");
        println!("    Setup: {:.3} ms", setup_time_std.as_secs_f64() * 1000.0);
        println!("    Interpolation (3 points): {:.3} ms", interp_time_std.as_secs_f64() * 1000.0);
        println!("  SIMD implementation:");
        println!("    Setup: {:.3} ms", setup_time_simd.as_secs_f64() * 1000.0);
        println!("    Interpolation (3 points): {:.3} ms", interp_time_simd.as_secs_f64() * 1000.0);
        println!("  Speedup: {:.2}x", interp_time_std.as_secs_f64() / interp_time_simd.as_secs_f64());

        // Test self-energy evaluation
        println!("\n2. Self-Energy Evaluation Performance:");

        let config = ContourDeformationConfig {
            n_imag_points: 20,
            xi_max: 100.0,
            eta: 1e-3,
            pole_threshold: 0.1,
            convergence_tol: 1e-6,
            use_simd: false,
            n_threads: None,
            compute_spectral: false,
        };

        // Standard evaluator
        let w_storage_std = WScreenedStorage::new(
            w_matrices.clone(),
            xi_grid.clone(),
            InterpolationMethod::Cubic,
        ).unwrap();

        let evaluator_std = EnergyDependentSelfEnergy::new(
            w_storage_std,
            mo_energy.clone(),
            mo_occ.clone(),
            df_tensor.clone(),
            config.clone(),
        ).unwrap();

        // SIMD evaluator
        let w_storage_simd = WScreenedStorageSimd::from_matrices(
            w_matrices,
            xi_grid,
            None,
        ).unwrap();

        let evaluator_simd = EnergyDependentSelfEnergySimd::new(
            w_storage_simd,
            mo_energy,
            mo_occ,
            df_tensor,
            config.eta,
        ).unwrap();

        // Test evaluation for HOMO orbital
        let homo_idx = n_occ - 1;
        let test_energy = 0.0;

        let start = Instant::now();
        let _sigma_std = evaluator_std.evaluate_orbital(homo_idx, test_energy).unwrap();
        let eval_time_std = start.elapsed();

        let start = Instant::now();
        let _sigma_simd = evaluator_simd.evaluate_orbital_simd(homo_idx, test_energy).unwrap();
        let eval_time_simd = start.elapsed();

        println!("  Standard: {:.3} ms", eval_time_std.as_secs_f64() * 1000.0);
        println!("  SIMD: {:.3} ms", eval_time_simd.as_secs_f64() * 1000.0);
        println!("  Speedup: {:.2}x", eval_time_std.as_secs_f64() / eval_time_simd.as_secs_f64());

        // Test parallel evaluation
        println!("\n3. Parallel Orbital Evaluation:");
        let orbital_indices: Vec<usize> = (0..n_occ.min(8)).collect();

        let start = Instant::now();
        let _results_seq: Vec<_> = orbital_indices
            .iter()
            .map(|&idx| evaluator_simd.evaluate_orbital_simd(idx, test_energy))
            .collect();
        let seq_time = start.elapsed();

        let start = Instant::now();
        let _results_par = evaluator_simd.evaluate_orbitals_parallel(&orbital_indices, test_energy).unwrap();
        let par_time = start.elapsed();

        println!("  Sequential ({} orbitals): {:.3} ms", orbital_indices.len(), seq_time.as_secs_f64() * 1000.0);
        println!("  Parallel ({} orbitals): {:.3} ms", orbital_indices.len(), par_time.as_secs_f64() * 1000.0);
        println!("  Speedup: {:.2}x", seq_time.as_secs_f64() / par_time.as_secs_f64());
    }

    println!("\n===========================================");
    println!("Summary: SIMD optimizations provide significant");
    println!("performance improvements for W interpolation and");
    println!("self-energy evaluation, especially for larger systems.");
    println!("===========================================");
}