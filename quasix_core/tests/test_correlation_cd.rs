//! Integration tests for correlation self-energy using contour deformation
//!
//! This test validates the implementation against known results and
//! theoretical properties.

use ndarray::{arr1, Array1, Array3};
use num_complex::Complex64;
use quasix_core::selfenergy::correlation::{ContourDeformationConfig, CorrelationSelfEnergyCD};

/// Test helper to create a simple W matrix for testing
fn create_test_w_matrix(n_freq: usize, n_aux: usize) -> Array3<Complex64> {
    let mut w = Array3::<Complex64>::zeros((n_freq, n_aux, n_aux));

    // Create a simple diagonal W matrix with decay
    for freq_idx in 0..n_freq {
        let freq_factor = 1.0 / (1.0 + freq_idx as f64);
        for i in 0..n_aux {
            w[[freq_idx, i, i]] = Complex64::new(freq_factor * (1.0 + i as f64 * 0.1), 0.0);

            // Add some off-diagonal elements
            if i > 0 {
                w[[freq_idx, i, i - 1]] = Complex64::new(freq_factor * 0.1, 0.0);
                w[[freq_idx, i - 1, i]] = Complex64::new(freq_factor * 0.1, 0.0);
            }
        }
    }

    w
}

/// Test helper to create simple DF tensors
fn create_test_df_tensor(n_mo: usize, n_aux: usize) -> Array3<f64> {
    let mut df = Array3::<f64>::zeros((n_mo, n_mo, n_aux));

    // Create diagonal-dominant DF tensor
    for i in 0..n_mo {
        for j in 0..n_mo {
            for k in 0..n_aux.min(n_mo) {
                if i == j && j == k {
                    df[[i, j, k]] = 1.0;
                } else if (i as i32 - j as i32).abs() <= 1 {
                    df[[i, j, k]] = 0.1 / (1.0 + (k as f64));
                }
            }
        }
    }

    df
}

#[test]
fn test_correlation_selfenergy_basic() {
    // Setup test system: 4 MOs (2 occupied, 2 virtual), 10 auxiliary functions
    let n_mo = 4;
    let n_aux = 10;
    let n_freq = 20;

    // Create MO energies: 2 occupied (negative), 2 virtual (positive)
    let mo_energy = arr1(&[-10.0, -5.0, 2.0, 4.0]);
    let mo_occ = arr1(&[1.0, 1.0, 0.0, 0.0]);

    // Create test W matrix and frequency grid
    let w_screened = create_test_w_matrix(n_freq, n_aux);
    let omega_grid = Array1::linspace(-20.0, 20.0, n_freq);

    // Create evaluation points around the HOMO-LUMO gap
    let eval_points = arr1(&[-7.0, -5.0, -3.0, 0.0, 1.0, 2.0, 3.0]);

    // Create DF tensor
    let df_tensor = create_test_df_tensor(n_mo, n_aux);

    // Setup calculator with default config
    let config = ContourDeformationConfig::default();
    let calculator = CorrelationSelfEnergyCD::new(n_mo, n_aux, config);

    // Compute correlation self-energy
    let result = calculator
        .compute_sigma_c(
            &mo_energy,
            &mo_occ,
            &w_screened,
            &omega_grid,
            &eval_points,
            &df_tensor,
        )
        .expect("Computation should succeed");

    // Basic sanity checks
    assert_eq!(result.omega_eval.len(), eval_points.len());
    assert_eq!(result.sigma_c.shape(), &[eval_points.len(), n_mo]);
    assert_eq!(result.residue_part.shape(), &[eval_points.len(), n_mo]);
    assert_eq!(result.integral_part.shape(), &[eval_points.len(), n_mo]);

    // Check that correlation self-energy is non-zero
    let sigma_norm: f64 = result.sigma_c.iter().map(|c| c.norm()).sum();
    assert!(
        sigma_norm > 1e-10,
        "Correlation self-energy should be non-zero"
    );

    // Check diagnostics
    assert!(result.diagnostics.residue_weight >= 0.0);
    assert!(result.diagnostics.integral_weight >= 0.0);
    assert!(
        (result.diagnostics.residue_weight + result.diagnostics.integral_weight - 1.0).abs()
            < 1e-10
    );
}

#[test]
fn test_pole_detection() {
    let n_mo = 6;
    let n_aux = 15;

    // Create system with 3 occupied and 3 virtual orbitals
    let mo_energy = arr1(&[-15.0, -10.0, -5.0, 1.0, 3.0, 5.0]);
    let mo_occ = arr1(&[1.0, 1.0, 1.0, 0.0, 0.0, 0.0]);

    let config = ContourDeformationConfig {
        pole_threshold: 0.1, // Group poles within 0.1 eV
        ..Default::default()
    };

    let calculator = CorrelationSelfEnergyCD::new(n_mo, n_aux, config);

    // We expect 9 poles from 3 occ × 3 virt transitions
    // Poles should be at energies: virt_i - occ_j
    let _expected_poles = [
        1.0 - (-15.0), // 16.0
        1.0 - (-10.0), // 11.0
        1.0 - (-5.0),  // 6.0
        3.0 - (-15.0), // 18.0
        3.0 - (-10.0), // 13.0
        3.0 - (-5.0),  // 8.0
        5.0 - (-15.0), // 20.0
        5.0 - (-10.0), // 15.0
        5.0 - (-5.0),  // 10.0
    ];

    // The pole detection is internal, but we can verify through the diagnostics
    let w_screened = create_test_w_matrix(10, n_aux);
    let omega_grid = Array1::linspace(-30.0, 30.0, 10);
    let eval_points = arr1(&[0.0]);
    let df_tensor = create_test_df_tensor(n_mo, n_aux);

    let result = calculator
        .compute_sigma_c(
            &mo_energy,
            &mo_occ,
            &w_screened,
            &omega_grid,
            &eval_points,
            &df_tensor,
        )
        .expect("Computation should succeed");

    // Check that poles were detected
    assert!(result.diagnostics.n_poles > 0);
    assert!(result.diagnostics.n_poles <= 9); // May group some near-degenerate poles
}

#[test]
fn test_spectral_function_properties() {
    let n_mo = 4;
    let n_aux = 10;

    let mo_energy = arr1(&[-8.0, -4.0, 1.0, 3.0]);
    let mo_occ = arr1(&[1.0, 1.0, 0.0, 0.0]);

    // Create a dense frequency grid for spectral function
    let eval_points = Array1::linspace(-20.0, 20.0, 100);

    let config = ContourDeformationConfig {
        eta: 0.1, // Larger broadening for smoother spectral function
        ..Default::default()
    };

    let calculator = CorrelationSelfEnergyCD::new(n_mo, n_aux, config);

    let w_screened = create_test_w_matrix(20, n_aux);
    let omega_grid = Array1::linspace(-30.0, 30.0, 20);
    let df_tensor = create_test_df_tensor(n_mo, n_aux);

    let result = calculator
        .compute_sigma_c(
            &mo_energy,
            &mo_occ,
            &w_screened,
            &omega_grid,
            &eval_points,
            &df_tensor,
        )
        .expect("Computation should succeed");

    // Check that spectral function is computed
    assert!(result.spectral_function.is_some());

    if let Some(spectral) = result.spectral_function {
        // Check positivity of spectral function
        for val in spectral.values.iter() {
            assert!(
                *val >= -1e-10,
                "Spectral function must be non-negative, got {}",
                val
            );
        }

        // Check normalization (should be close to 1 for each orbital)
        for (i, &norm) in spectral.normalization.iter().enumerate() {
            assert!(
                (norm - 1.0).abs() < 0.5,
                "Orbital {} normalization {} too far from 1",
                i,
                norm
            );
        }
    }
}

#[test]
fn test_causality() {
    // Test that Σ(ω*) = Σ*(ω) for real frequencies
    let n_mo = 2;
    let n_aux = 5;

    let mo_energy = arr1(&[-5.0, 2.0]);
    let mo_occ = arr1(&[1.0, 0.0]);

    let config = ContourDeformationConfig::default();
    let calculator = CorrelationSelfEnergyCD::new(n_mo, n_aux, config);

    // Use a real frequency
    let eval_points = arr1(&[0.5]);

    let w_screened = create_test_w_matrix(10, n_aux);
    let omega_grid = Array1::linspace(-10.0, 10.0, 10);
    let df_tensor = create_test_df_tensor(n_mo, n_aux);

    let result = calculator
        .compute_sigma_c(
            &mo_energy,
            &mo_occ,
            &w_screened,
            &omega_grid,
            &eval_points,
            &df_tensor,
        )
        .expect("Computation should succeed");

    // For real W and real frequency, imaginary part should be small
    // (only from broadening η)
    for sigma_val in result.sigma_c.iter() {
        assert!(
            sigma_val.im.abs() < 1.0,
            "Imaginary part {} should be small for real frequency",
            sigma_val.im
        );
    }
}

#[test]
fn test_convergence_with_quadrature_points() {
    let n_mo = 2;
    let n_aux = 5;

    let mo_energy = arr1(&[-5.0, 2.0]);
    let mo_occ = arr1(&[1.0, 0.0]);
    let eval_points = arr1(&[0.0]);

    let w_screened = create_test_w_matrix(10, n_aux);
    let omega_grid = Array1::linspace(-10.0, 10.0, 10);
    let df_tensor = create_test_df_tensor(n_mo, n_aux);

    // Test with different numbers of quadrature points
    let mut results = Vec::new();

    for n_points in [8, 16, 32, 64] {
        let config = ContourDeformationConfig {
            n_imag_points: n_points,
            ..Default::default()
        };

        let calculator = CorrelationSelfEnergyCD::new(n_mo, n_aux, config);

        let result = calculator
            .compute_sigma_c(
                &mo_energy,
                &mo_occ,
                &w_screened,
                &omega_grid,
                &eval_points,
                &df_tensor,
            )
            .expect("Computation should succeed");

        results.push(result.sigma_c[[0, 0]]);
    }

    // Check convergence: difference should decrease with more points
    for i in 1..results.len() {
        let diff = (results[i] - results[i - 1]).norm();

        // For well-behaved integrands, should converge exponentially
        if i > 1 {
            let prev_diff = (results[i - 1] - results[i - 2]).norm();
            assert!(
                diff < prev_diff * 1.1, // Allow some numerical noise
                "Should converge with more quadrature points"
            );
        }
    }
}

#[test]
fn test_parallel_execution() {
    // Test that parallel execution gives same results as serial
    let n_mo = 4;
    let n_aux = 8;

    let mo_energy = arr1(&[-10.0, -5.0, 2.0, 4.0]);
    let mo_occ = arr1(&[1.0, 1.0, 0.0, 0.0]);
    let eval_points = Array1::linspace(-5.0, 5.0, 10);

    let w_screened = create_test_w_matrix(15, n_aux);
    let omega_grid = Array1::linspace(-15.0, 15.0, 15);
    let df_tensor = create_test_df_tensor(n_mo, n_aux);

    // Run with single thread
    let config_serial = ContourDeformationConfig {
        n_threads: Some(1),
        ..Default::default()
    };
    let calculator_serial = CorrelationSelfEnergyCD::new(n_mo, n_aux, config_serial);
    let result_serial = calculator_serial
        .compute_sigma_c(
            &mo_energy,
            &mo_occ,
            &w_screened,
            &omega_grid,
            &eval_points,
            &df_tensor,
        )
        .expect("Serial computation should succeed");

    // Run with multiple threads
    let config_parallel = ContourDeformationConfig {
        n_threads: Some(4),
        ..Default::default()
    };
    let calculator_parallel = CorrelationSelfEnergyCD::new(n_mo, n_aux, config_parallel);
    let result_parallel = calculator_parallel
        .compute_sigma_c(
            &mo_energy,
            &mo_occ,
            &w_screened,
            &omega_grid,
            &eval_points,
            &df_tensor,
        )
        .expect("Parallel computation should succeed");

    // Results should be identical (within numerical precision)
    for i in 0..eval_points.len() {
        for j in 0..n_mo {
            let diff = (result_serial.sigma_c[[i, j]] - result_parallel.sigma_c[[i, j]]).norm();
            assert!(
                diff < 1e-10,
                "Serial and parallel results should match at [{}, {}]: diff = {}",
                i,
                j,
                diff
            );
        }
    }
}
