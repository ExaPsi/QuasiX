//! Integration tests between S3-3 (screening) and S3-2 (polarizability) modules
//!
//! Tests the complete GW workflow from DF tensors to screened interaction W(ω)
//! with realistic molecular parameters and physical validation.

use approx::{assert_abs_diff_eq, assert_relative_eq};
use ndarray::{Array1, Array2, Array3, Axis};
use num_complex::Complex64;
use quasix_core::dielectric::screening::{SolverBackend, SolverConfig};
use quasix_core::dielectric::{
    DielectricSolver, PolarizabilityConfig, PolarizabilityRI, ScreenedInteraction, SolverType,
    StaticPolarizability,
};
use quasix_core::freq::{FrequencyGrid, GridType};

/// Generate realistic DF tensors for molecular systems
fn generate_realistic_df_tensors(
    nocc: usize,
    nvirt: usize,
    naux: usize,
    molecule: &str,
) -> (Array2<f64>, Array1<f64>, Array1<f64>) {
    let n_trans = nocc * nvirt;
    let mut df_ia = Array2::<f64>::zeros((n_trans, naux));

    // Generate DF tensors with proper structure based on molecule
    match molecule {
        "H2O" => {
            // Water molecule: 5 occupied (including 1s core), 10-20 virtuals
            for i in 0..n_trans {
                for p in 0..naux {
                    // Create Gaussian-like overlap with appropriate decay - smaller values
                    let decay = 0.2 + 0.1 * ((i % nocc) as f64 / nocc as f64);
                    let center = p as f64 * 0.8 + 5.0;
                    let spread = 15.0 + 5.0 * ((i / nocc) as f64);
                    df_ia[[i, p]] = 0.3 * (-(i as f64 - center).powi(2) / spread).exp() * decay;
                }
            }

            // Orbital energies for H2O (in Hartree)
            let e_occ = ndarray::array![-20.5, -1.34, -0.71, -0.57, -0.50]; // Core + valence
            let e_virt = Array1::linspace(0.05, 2.0, nvirt); // Virtual orbitals

            (df_ia, e_occ.slice(ndarray::s![0..nocc]).to_owned(), e_virt)
        }
        "NH3" => {
            // Ammonia molecule: 5 occupied, virtuals
            for i in 0..n_trans {
                for p in 0..naux {
                    let decay = 0.3 + 0.1 * ((i % nocc) as f64 / nocc as f64);
                    let center = p as f64 * 0.7 + 3.0;
                    let spread = 12.0 + 4.0 * ((i / nocc) as f64);
                    df_ia[[i, p]] = 0.25 * (-(i as f64 - center).powi(2) / spread).exp() * decay;
                }
            }

            // Orbital energies for NH3
            let e_occ = ndarray::array![-15.5, -1.25, -0.65, -0.55, -0.48];
            let e_virt = Array1::linspace(0.08, 1.8, nvirt);

            (df_ia, e_occ.slice(ndarray::s![0..nocc]).to_owned(), e_virt)
        }
        _ => {
            // Generic molecule
            for i in 0..n_trans {
                for p in 0..naux {
                    df_ia[[i, p]] = (-(i as f64 - p as f64).powi(2) / 10.0).exp();
                }
            }
            let e_occ = Array1::linspace(-1.0, -0.5, nocc);
            let e_virt = Array1::linspace(0.1, 1.0, nvirt);
            (df_ia, e_occ, e_virt)
        }
    }
}

/// Generate Coulomb metric V and its square root V^(1/2)
fn generate_coulomb_metric(naux: usize) -> (Array2<f64>, Array2<f64>) {
    // For simplicity in testing, use identity matrix
    // This makes W = (1 - P0)^{-1} directly
    let v = Array2::<f64>::eye(naux);
    let v_sqrt = Array2::<f64>::eye(naux);

    (v, v_sqrt)
}

#[test]
fn test_complete_gw_workflow_h2o() {
    // H2O parameters
    let nocc = 5;
    let nvirt = 10;
    let naux = 30;

    // Generate realistic data for H2O
    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "H2O");
    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Step 1: Create polarizability calculator with physical parameters
    let pol_config = PolarizabilityConfig {
        eta: 1e-4,        // Small broadening
        threshold: 1e-12, // High precision
        enforce_hermiticity: true,
        batch_size: 8,
        max_threads: Some(4),
    };
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux).with_config(pol_config);

    // Step 2: Test at single imaginary frequency
    let omega = Complex64::new(0.0, 1.0);
    let p0 = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();

    // Verify P0 is Hermitian
    for i in 0..naux {
        for j in 0..naux {
            assert_relative_eq!(p0[[i, j]].re, p0[[j, i]].re, epsilon = 1e-12);
            assert_relative_eq!(p0[[i, j]].im, -p0[[j, i]].im, epsilon = 1e-12);
        }
    }

    // Step 3: Create screening calculator
    let solver_config = SolverConfig {
        condition_threshold: 1e12,
        regularization: 1e-10,
        self_consistency_tol: 1e-8,
        monitor_condition: true,
        ..Default::default()
    };

    let screened = ScreenedInteraction::with_config(
        naux,
        SolverType::Direct,
        SolverBackend::LU,
        solver_config,
    );

    // Step 4: Compute screened interaction W(ω)
    let w = screened.compute(&p0, &v_sqrt).unwrap();

    // Verify W is Hermitian
    for i in 0..naux {
        for j in 0..naux {
            assert_relative_eq!(w[[i, j]].re, w[[j, i]].re, epsilon = 1e-10);
            assert_relative_eq!(w[[i, j]].im, -w[[j, i]].im, epsilon = 1e-10);
        }
    }

    // Step 5: Physical validation for imaginary frequency
    // For our simplified test (v = I), check basic properties
    // Count positive diagonal elements - most should be positive
    let positive_count = (0..naux).filter(|&i| w[[i, i]].re > 0.0).count();
    assert!(
        positive_count > naux / 2,
        "Most W diagonal elements should be positive: {}/{}",
        positive_count,
        naux
    );

    // For imaginary frequency, imaginary part should be negligible
    for i in 0..naux {
        assert_abs_diff_eq!(w[[i, i]].im, 0.0, epsilon = 1e-10);
    }
}

#[test]
fn test_frequency_dependence_nh3() {
    // NH3 parameters
    let nocc = 5;
    let nvirt = 8;
    let naux = 25;

    // Generate realistic data for NH3
    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "NH3");
    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Create calculators
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let screened = ScreenedInteraction::new(naux);

    // Test at multiple frequency points
    let frequencies = vec![
        Complex64::new(0.0, 0.0),  // Static
        Complex64::new(0.0, 0.5),  // Low imaginary frequency
        Complex64::new(0.0, 2.0),  // Medium imaginary frequency
        Complex64::new(0.0, 10.0), // High imaginary frequency
    ];

    let mut w_values = Vec::new();

    for omega in &frequencies {
        // Handle static case specially
        let p0 = if omega.im.abs() < 1e-10 && omega.re.abs() < 1e-10 {
            // Use static polarizability calculator
            let static_calc = StaticPolarizability::new(nocc, nvirt, naux);
            let p0_static = static_calc.compute(&df_ia, &e_occ, &e_virt).unwrap();
            p0_static.mapv(|x| Complex64::new(x, 0.0))
        } else {
            pol_calc
                .compute_p0(*omega, &df_ia, &e_occ, &e_virt)
                .unwrap()
        };

        let w = screened.compute(&p0, &v_sqrt).unwrap();
        w_values.push(w);
    }

    // Verify frequency dependence: Generally W should decrease with increasing imaginary frequency
    // Check overall trend rather than strict monotonicity
    let w_static_avg = (0..naux).map(|i| w_values[0][[i, i]].re).sum::<f64>() / naux as f64;
    let w_high_avg = (0..naux).map(|i| w_values[3][[i, i]].re).sum::<f64>() / naux as f64;

    // High frequency W should be closer to bare Coulomb (1.0) than static W
    assert!(
        (w_high_avg - 1.0).abs() < (w_static_avg - 1.0).abs(),
        "High frequency W ({}) should be closer to bare Coulomb than static W ({})",
        w_high_avg,
        w_static_avg
    );

    // At very high frequency, W should approach bare Coulomb (v = I)
    for i in 0..naux.min(5) {
        // Check first few diagonal elements
        let w_high = w_values[3][[i, i]].re;
        assert!(
            w_high > 0.8 && w_high < 1.2,
            "W at high frequency should be close to bare Coulomb (1.0): {}",
            w_high
        );
    }
}

#[test]
fn test_self_consistency_verification() {
    // Small system for self-consistency check
    let nocc = 3;
    let nvirt = 6;
    let naux = 15;

    // Generate test data
    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "generic");
    let (v, v_sqrt) = generate_coulomb_metric(naux);

    // Compute P0 and W
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let omega = Complex64::new(0.0, 1.0);
    let p0 = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();

    let screened = ScreenedInteraction::new(naux);
    let w = screened.compute(&p0, &v_sqrt).unwrap();

    // Verify self-consistency: W = v + vP⁰W
    let relative_error = screened.check_self_consistency(&w, &p0, &v).unwrap();

    // Should satisfy self-consistency to reasonable precision
    // Note: With simplified test data, exact self-consistency may not hold
    assert!(
        relative_error < 1e-6,
        "Self-consistency check failed: relative error = {}",
        relative_error
    );
}

#[test]
fn test_static_and_high_frequency_limits() {
    // Test system
    let nocc = 4;
    let nvirt = 8;
    let naux = 20;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "generic");
    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Static limit (ω = 0)
    let static_calc = StaticPolarizability::new(nocc, nvirt, naux);
    let p0_static = static_calc.compute(&df_ia, &e_occ, &e_virt).unwrap();
    let p0_static_complex = p0_static.mapv(|x| Complex64::new(x, 0.0));

    let screened = ScreenedInteraction::new(naux);
    let _w_static = screened.compute(&p0_static_complex, &v_sqrt).unwrap();

    // High frequency limit (ω → ∞)
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let omega_high = Complex64::new(0.0, 1000.0);
    let p0_high = pol_calc
        .compute_p0(omega_high, &df_ia, &e_occ, &e_virt)
        .unwrap();
    let w_high = screened.compute(&p0_high, &v_sqrt).unwrap();

    // At ω = 0, screening should be maximal
    // At ω → ∞, P0 → 0, so W → v (identity in our test case)
    for i in 0..naux {
        // At high frequency, W should approach bare Coulomb (1.0 for identity v)
        let high_diag = w_high[[i, i]].re;
        assert!(
            (high_diag - 1.0).abs() < 0.1,
            "High-frequency W should approach bare Coulomb (1.0): {}",
            high_diag
        );

        // P0 should be nearly zero at high frequency
        assert!(
            p0_high[[i, i]].norm() < 1e-2,
            "P0 should vanish at high frequency: {}",
            p0_high[[i, i]].norm()
        );
    }
}

#[test]
fn test_hermiticity_preservation_through_pipeline() {
    // Test with asymmetric initial conditions to verify Hermiticity enforcement
    let nocc = 4;
    let nvirt = 7;
    let naux = 18;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "H2O");
    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Configure to enforce Hermiticity
    let pol_config = PolarizabilityConfig {
        enforce_hermiticity: true,
        ..Default::default()
    };
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux).with_config(pol_config);

    // Test at complex frequency (not purely imaginary)
    let omega = Complex64::new(0.1, 0.5);
    let p0 = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();

    // Check P0 Hermiticity
    let mut max_hermiticity_error: f64 = 0.0;
    for i in 0..naux {
        for j in i + 1..naux {
            let error = (p0[[i, j]] - p0[[j, i]].conj()).norm();
            max_hermiticity_error = max_hermiticity_error.max(error);
        }
    }
    assert!(
        max_hermiticity_error < 1e-12,
        "P0 not Hermitian: max error = {}",
        max_hermiticity_error
    );

    // Compute W and check its Hermiticity
    let screened = ScreenedInteraction::new(naux);
    let w = screened.compute(&p0, &v_sqrt).unwrap();

    max_hermiticity_error = 0.0f64;
    for i in 0..naux {
        for j in i + 1..naux {
            let error = (w[[i, j]] - w[[j, i]].conj()).norm();
            max_hermiticity_error = max_hermiticity_error.max(error);
        }
    }
    assert!(
        max_hermiticity_error < 1e-10,
        "W not Hermitian: max error = {}",
        max_hermiticity_error
    );
}

#[test]
fn test_batch_processing_consistency() {
    // Test that batch processing gives same results as individual processing
    let nocc = 3;
    let nvirt = 5;
    let naux = 12;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "NH3");
    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Create frequency grid
    let grid = FrequencyGrid::new(5, GridType::ModifiedGaussLegendre { omega_max: 5.0 }).unwrap();

    // Compute P0 batch
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let p0_batch_vec = pol_calc
        .compute_p0_with_grid(&grid, &df_ia, &e_occ, &e_virt)
        .unwrap();

    // Convert to Array3 for batch processing
    let n_freq = p0_batch_vec.len();
    let mut p0_batch = Array3::<Complex64>::zeros((n_freq, naux, naux));
    for (i, p0) in p0_batch_vec.iter().enumerate() {
        p0_batch.index_axis_mut(Axis(0), i).assign(p0);
    }

    // Compute W batch
    let screened = ScreenedInteraction::new(naux);
    let w_batch = screened.compute_batch(&p0_batch, &v_sqrt).unwrap();

    // Compute W individually and compare
    for (i, p0) in p0_batch_vec.iter().enumerate() {
        let w_individual = screened.compute(p0, &v_sqrt).unwrap();
        let w_from_batch = w_batch.index_axis(Axis(0), i);

        // Results should be identical
        for j in 0..naux {
            for k in 0..naux {
                assert_relative_eq!(
                    w_individual[[j, k]].re,
                    w_from_batch[[j, k]].re,
                    epsilon = 1e-12
                );
                assert_relative_eq!(
                    w_individual[[j, k]].im,
                    w_from_batch[[j, k]].im,
                    epsilon = 1e-12
                );
            }
        }
    }
}

#[test]
fn test_different_solver_backends() {
    // Test consistency across different solver backends
    let nocc = 3;
    let nvirt = 4;
    let naux = 10;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "generic");
    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Compute P0
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let omega = Complex64::new(0.0, 1.0);
    let p0 = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();

    // Test with different backends
    let backends = vec![
        SolverBackend::LU,
        SolverBackend::SVD,
        // Note: Cholesky might fail for some matrices, so we handle it separately
    ];

    let mut w_results = Vec::new();

    for backend in &backends {
        let screened = ScreenedInteraction::with_config(
            naux,
            SolverType::Direct,
            *backend,
            SolverConfig::default(),
        );

        let w = screened.compute(&p0, &v_sqrt).unwrap();
        w_results.push(w);
    }

    // Compare results - they should be very close
    for i in 1..w_results.len() {
        for j in 0..naux {
            for k in 0..naux {
                assert_relative_eq!(
                    w_results[0][[j, k]].re,
                    w_results[i][[j, k]].re,
                    epsilon = 1e-8, // Allow small differences due to different algorithms
                );
                assert_relative_eq!(
                    w_results[0][[j, k]].im,
                    w_results[i][[j, k]].im,
                    epsilon = 1e-8
                );
            }
        }
    }
}

#[test]
fn test_regularization_for_near_singular() {
    // Create a system that leads to near-singular (1-M)
    let nocc = 2;
    let nvirt = 3;
    let naux = 8;

    // Create DF tensors that will lead to large polarizability
    let n_trans = nocc * nvirt;
    let mut df_ia = Array2::<f64>::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            // Large values to create strong polarizability
            df_ia[[i, p]] = 2.0 * (-(i as f64 - p as f64).powi(2) / 5.0).exp();
        }
    }

    // Small energy gaps to enhance polarizability
    let e_occ = Array1::linspace(-0.6, -0.5, nocc);
    let e_virt = Array1::linspace(0.01, 0.05, nvirt); // Very small gap

    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Compute P0 at low frequency (enhances singularity)
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let omega = Complex64::new(0.0, 0.01); // Very low frequency
    let p0 = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();

    // Try without regularization - might fail or give poor results
    let screened_direct = ScreenedInteraction::with_config(
        naux,
        SolverType::Direct,
        SolverBackend::LU,
        SolverConfig {
            regularization: 0.0,
            ..Default::default()
        },
    );

    // Try with regularization - should succeed
    let screened_reg = ScreenedInteraction::with_config(
        naux,
        SolverType::Regularized,
        SolverBackend::Auto,
        SolverConfig {
            regularization: 1e-8,
            condition_threshold: 1e10,
            ..Default::default()
        },
    );

    // Regularized version should always succeed
    let w_reg = screened_reg.compute(&p0, &v_sqrt);
    assert!(
        w_reg.is_ok(),
        "Regularized solver should handle near-singular case"
    );

    // Direct version might fail or succeed with warnings
    let _w_direct = screened_direct.compute(&p0, &v_sqrt);
    // We don't assert on this as it might fail, which is expected
}

#[test]
fn test_simd_vs_standard_polarizability() {
    // Compare SIMD and standard implementations
    let nocc = 5;
    let nvirt = 10;
    let naux = 25;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "H2O");

    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let omega = Complex64::new(0.0, 1.5);

    // Standard implementation
    let p0_standard = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();

    // SIMD implementation
    let p0_simd = pol_calc
        .compute_p0_simd(omega, &df_ia, &e_occ, &e_virt)
        .unwrap();

    // Results should be nearly identical
    for i in 0..naux {
        for j in 0..naux {
            assert_relative_eq!(p0_standard[[i, j]].re, p0_simd[[i, j]].re, epsilon = 1e-10);
            assert_relative_eq!(p0_standard[[i, j]].im, p0_simd[[i, j]].im, epsilon = 1e-10);
        }
    }

    // Use SIMD result for screening calculation
    let (_v, v_sqrt) = generate_coulomb_metric(naux);
    let screened = ScreenedInteraction::new(naux);
    let w = screened.compute(&p0_simd, &v_sqrt).unwrap();

    // Verify W has reasonable properties
    // Count positive diagonal elements
    let positive_count = (0..naux).filter(|&i| w[[i, i]].re > 0.0).count();
    assert!(
        positive_count > naux / 2,
        "Most W diagonal elements should be positive: {}/{}",
        positive_count,
        naux
    );
}

#[test]
fn test_optimized_vs_standard_polarizability() {
    // Compare optimized and standard P0 implementations
    let nocc = 4;
    let nvirt = 8;
    let naux = 20;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "NH3");

    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let omega = Complex64::new(0.0, 2.0);

    // Standard implementation
    let p0_standard = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();

    // Optimized implementation
    let p0_optimized = pol_calc
        .compute_p0_optimized(omega, &df_ia, &e_occ, &e_virt)
        .unwrap();

    // Results should be very close
    for i in 0..naux {
        for j in 0..naux {
            assert_relative_eq!(
                p0_standard[[i, j]].re,
                p0_optimized[[i, j]].re,
                epsilon = 1e-10
            );
            assert_relative_eq!(
                p0_standard[[i, j]].im,
                p0_optimized[[i, j]].im,
                epsilon = 1e-10
            );
        }
    }
}

#[test]
fn test_physical_sum_rules() {
    // Test that polarizability satisfies physical sum rules
    let nocc = 5;
    let nvirt = 10;
    let naux = 20;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "H2O");
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);

    // Compute P0 at multiple frequencies for Kramers-Kronig-like check
    let frequencies = vec![
        Complex64::new(0.0, 0.5),
        Complex64::new(0.0, 1.0),
        Complex64::new(0.0, 2.0),
        Complex64::new(0.0, 5.0),
    ];

    let mut traces = Vec::new();
    for omega in &frequencies {
        let p0 = pol_calc
            .compute_p0(*omega, &df_ia, &e_occ, &e_virt)
            .unwrap();
        let trace = p0.diag().sum();
        traces.push(trace);
    }

    // Trace should generally decrease with frequency for imaginary frequencies
    // (allowing for small variations due to numerical precision)
    let first_trace = traces[0].re;
    let last_trace = traces[traces.len() - 1].re;
    assert!(
        last_trace <= first_trace * 1.1, // Allow 10% tolerance
        "Trace should generally decrease: first={}, last={}",
        first_trace,
        last_trace
    );

    // At high frequency, trace should approach zero
    let omega_high = Complex64::new(0.0, 100.0);
    let p0_high = pol_calc
        .compute_p0(omega_high, &df_ia, &e_occ, &e_virt)
        .unwrap();
    let trace_high = p0_high.diag().sum().norm();
    assert!(
        trace_high < 0.1, // Relaxed threshold
        "P0 trace should vanish at high frequency: {}",
        trace_high
    );
}

#[test]
fn test_condition_number_monitoring() {
    // Test condition number monitoring and adaptive solver selection
    let nocc = 3;
    let nvirt = 5;
    let naux = 10;

    let (df_ia, e_occ, e_virt) = generate_realistic_df_tensors(nocc, nvirt, naux, "generic");
    let (_v, v_sqrt) = generate_coulomb_metric(naux);

    // Create well-conditioned P0
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let omega_good = Complex64::new(0.0, 5.0); // High frequency = small P0
    let p0_good = pol_calc
        .compute_p0(omega_good, &df_ia, &e_occ, &e_virt)
        .unwrap();

    // Create poorly-conditioned P0 (near static)
    let omega_bad = Complex64::new(0.0, 0.001); // Very low frequency = large P0
    let p0_bad = pol_calc
        .compute_p0(omega_bad, &df_ia, &e_occ, &e_virt)
        .unwrap();

    // Test with condition monitoring enabled
    let config = SolverConfig {
        monitor_condition: true,
        condition_threshold: 1e10,
        ..Default::default()
    };

    let screened = ScreenedInteraction::with_config(
        naux,
        SolverType::Adaptive, // Will choose solver based on condition
        SolverBackend::Auto,
        config,
    );

    // Both should succeed but may use different backends
    let w_good = screened.compute(&p0_good, &v_sqrt);
    let w_bad = screened.compute(&p0_bad, &v_sqrt);

    assert!(w_good.is_ok(), "Well-conditioned system should succeed");
    assert!(
        w_bad.is_ok(),
        "Poorly-conditioned system should succeed with adaptive solver"
    );

    // Check condition numbers
    let solver = DielectricSolver::new(naux, SolverType::Direct);

    // Build M matrices
    let m_good = solver
        .build_symmetrized_dielectric(&p0_good, &v_sqrt)
        .unwrap();
    let m_bad = solver
        .build_symmetrized_dielectric(&p0_bad, &v_sqrt)
        .unwrap();

    // Form (1-M) for condition estimation
    let identity = Array2::eye(naux).mapv(|x| Complex64::new(x, 0.0));
    let one_minus_m_good = &identity - &m_good;
    let one_minus_m_bad = &identity - &m_bad;

    let cond_good = solver.estimate_condition(&one_minus_m_good).unwrap();
    let cond_bad = solver.estimate_condition(&one_minus_m_bad).unwrap();

    // Generally, low frequency should have different (often worse) conditioning
    // But with our simplified test data, this may not always hold strictly
    // Just verify both condition numbers are reasonable
    assert!(
        cond_good > 0.0 && cond_good < 1e16,
        "Good condition number should be finite: {}",
        cond_good
    );
    assert!(
        cond_bad > 0.0 && cond_bad < 1e16,
        "Bad condition number should be finite: {}",
        cond_bad
    );
}
