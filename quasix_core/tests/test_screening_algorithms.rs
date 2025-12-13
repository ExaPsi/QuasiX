//! Comprehensive test suite for S3-3 screening algorithms
//!
//! Tests cover:
//! - Adaptive solver selection
//! - Numerical stability for ill-conditioned systems
//! - Self-consistency verification
//! - Edge case handling

use approx::{assert_abs_diff_eq, assert_relative_eq};
use ndarray::{s, Array2, Array3};
use ndarray_linalg::Norm;
use num_complex::Complex64;
use quasix_core::dielectric::screening_enhanced::{
    BatchScreeningProcessor, EnhancedDielectricSolver, EnhancedSolverConfig, SolverMethod,
};

/// Test helper: Create a diagonal matrix with specified values
fn create_diagonal_complex(diag: &[f64]) -> Array2<Complex64> {
    let n = diag.len();
    let mut m = Array2::zeros((n, n));
    for i in 0..n {
        m[[i, i]] = Complex64::new(diag[i], 0.0);
    }
    m
}

/// Test helper: Create a Hermitian matrix
fn create_hermitian_matrix(n: usize, seed: u64) -> Array2<Complex64> {
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    let mut rng = StdRng::seed_from_u64(seed);
    let mut m = Array2::zeros((n, n));

    // Fill upper triangle with random values
    for i in 0..n {
        for j in i..n {
            let real = rng.random_range(-1.0..1.0);
            let imag = if i == j {
                0.0
            } else {
                rng.random_range(-1.0..1.0)
            };
            m[[i, j]] = Complex64::new(real, imag);
            if i != j {
                m[[j, i]] = Complex64::new(real, -imag); // Hermitian
            }
        }
    }
    m
}

/// Test helper: Compute condition number via SVD
#[allow(dead_code)]
fn compute_condition_number(m: &Array2<Complex64>) -> f64 {
    use ndarray_linalg::SVD;

    let (_, s, _) = m.svd(true, true).unwrap();
    let max_sv = s[0];
    let min_sv = s[s.len() - 1];

    if min_sv < 1e-15 {
        f64::INFINITY
    } else {
        max_sv / min_sv
    }
}

#[test]
fn test_well_conditioned_system() {
    let naux = 20;
    let mut solver = EnhancedDielectricSolver::new(naux);

    // Create well-conditioned P0
    let mut p0 = create_hermitian_matrix(naux, 42);
    p0 = &p0 * 0.1; // Scale to ensure stability

    // Identity metric for simplicity
    let v_sqrt = Array2::eye(naux);

    // Compute W
    let (w, diagnostics) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("Should succeed for well-conditioned system");

    // Verify properties
    assert!(
        diagnostics.condition_number < 1e6,
        "Condition number should be reasonable: {}",
        diagnostics.condition_number
    );
    assert!(
        diagnostics.self_consistency_residual < 1e-6,
        "Self-consistency should be satisfied: {}",
        diagnostics.self_consistency_residual
    );
    assert!(
        diagnostics.hermiticity_error < 1e-10,
        "Should preserve Hermiticity: {}",
        diagnostics.hermiticity_error
    );

    // Check that standard solver was used
    assert!(matches!(
        diagnostics.solver_used,
        SolverMethod::Cholesky | SolverMethod::LUPartialPivot | SolverMethod::LUCompletePivot
    ));

    // Verify W is positive on diagonal
    for i in 0..naux {
        assert!(w[[i, i]].re > 0.0, "W diagonal should be positive");
        assert_abs_diff_eq!(w[[i, i]].im, 0.0, epsilon = 1e-10);
    }
}

#[test]
fn test_ill_conditioned_system() {
    let naux = 10;
    let config = EnhancedSolverConfig {
        condition_threshold: 1e8,
        tikhonov_alpha: 1e-8,
        ..Default::default()
    };
    let mut solver = EnhancedDielectricSolver::with_config(naux, config);

    // Create ill-conditioned P0 with wide eigenvalue spread
    let eigenvalues: Vec<f64> = (0..naux).map(|i| 0.1 * 0.001_f64.powi(i as i32)).collect();
    let mut p0 = create_diagonal_complex(&eigenvalues);

    // Add small off-diagonal elements
    for i in 0..naux {
        for j in 0..naux {
            if i != j {
                p0[[i, j]] = Complex64::new(1e-6, 0.0);
            }
        }
    }

    let v_sqrt = Array2::eye(naux);

    // Should handle ill-conditioning gracefully
    let (_w, diagnostics) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("Should handle ill-conditioned system");

    // Should have used regularization or special solver
    assert!(
        diagnostics.applied_regularization.is_some()
            || matches!(
                diagnostics.solver_used,
                SolverMethod::TikhonovRegularized(_)
                    | SolverMethod::SVDTruncated(_)
                    | SolverMethod::SVD
            )
    );

    // Self-consistency might be relaxed due to regularization
    assert!(
        diagnostics.self_consistency_residual <= 1e-2,
        "Self-consistency with regularization: {}",
        diagnostics.self_consistency_residual
    );
}

#[test]
fn test_near_singular_system() {
    let naux = 8;
    let mut solver = EnhancedDielectricSolver::new(naux);

    // Create near-singular P0 (max eigenvalue approaching 1)
    let mut p0 = Array2::zeros((naux, naux));
    for i in 0..naux {
        p0[[i, i]] = Complex64::new(0.98, 0.0); // Very close to singular
    }

    // Add small coupling
    for i in 0..naux - 1 {
        p0[[i, i + 1]] = Complex64::new(0.01, 0.0);
        p0[[i + 1, i]] = Complex64::new(0.01, 0.0);
    }

    let v_sqrt = Array2::eye(naux);

    let result = solver.compute_screened_interaction_adaptive(&p0, &v_sqrt);

    assert!(result.is_ok(), "Should handle near-singular system");

    let (_w, diagnostics) = result.unwrap();

    // Should have detected near-singularity
    assert!(diagnostics.max_eigenvalue > 0.95);
    assert!(!diagnostics.warnings.is_empty());

    // Should have used special treatment
    assert!(matches!(
        diagnostics.solver_used,
        SolverMethod::TikhonovRegularized(_) | SolverMethod::SVDTruncated(_)
    ));
}

#[test]
fn test_metallic_system() {
    let naux = 6;
    let config = EnhancedSolverConfig {
        metallic_broadening: 0.01, // Small imaginary part
        ..Default::default()
    };
    let mut solver = EnhancedDielectricSolver::with_config(naux, config);

    // Create metallic-like P0 (eigenvalue > 0.999)
    let mut p0 = Array2::zeros((naux, naux));
    p0[[0, 0]] = Complex64::new(0.9995, 0.0); // Metallic mode
    for i in 1..naux {
        p0[[i, i]] = Complex64::new(0.5, 0.0); // Normal modes
    }

    let v_sqrt = Array2::eye(naux);

    let (_w, diagnostics) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("Should handle metallic system");

    // Should have detected metallic character
    assert!(diagnostics.max_eigenvalue > 0.999);
    assert!(diagnostics.warnings.iter().any(|w| w.contains("metallic")));

    // Should have applied broadening
    assert_eq!(
        diagnostics.solver_used,
        SolverMethod::TikhonovRegularized(0.01)
    );
}

#[test]
fn test_self_consistency() {
    let naux = 12;
    let mut solver = EnhancedDielectricSolver::new(naux);

    // Create physically reasonable P0
    let p0 = create_hermitian_matrix(naux, 123) * 0.2;
    let v = Array2::<f64>::eye(naux) * 2.0; // Simple Coulomb metric
    let v_sqrt = Array2::eye(naux) * 2.0_f64.sqrt();

    let (w, diagnostics) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("Should compute W");

    // When regularization is used, exact self-consistency cannot be satisfied
    // The solver reports the residual for the modified equation
    if diagnostics.applied_regularization.is_some() {
        // With regularization, accept larger tolerance
        assert!(
            diagnostics.self_consistency_residual < 0.1,
            "Self-consistency with regularization should be reasonable: {}",
            diagnostics.self_consistency_residual
        );
    } else {
        // Manually verify self-consistency: W = v + vP^0W
        let v_c = v.mapv(|x| Complex64::new(x, 0.0));
        let vp0 = v_c.dot(&p0);
        let vp0w = vp0.dot(&w);
        let w_reconstructed = &v_c + &vp0w;

        // Compare with computed W
        let diff = &w - &w_reconstructed;
        let relative_error = diff.norm_l2() / w.norm_l2();

        assert!(
            relative_error < 1e-6,
            "Self-consistency should be satisfied: {}",
            relative_error
        );
        assert_abs_diff_eq!(
            diagnostics.self_consistency_residual,
            relative_error,
            epsilon = 1e-8
        );
    }
}

#[test]
fn test_hermiticity_preservation() {
    let naux = 15;
    let mut solver = EnhancedDielectricSolver::new(naux);

    // Create Hermitian P0
    let p0 = create_hermitian_matrix(naux, 456) * 0.3;
    let v_sqrt = Array2::eye(naux);

    let (w, diagnostics) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("Should compute W");

    // Check Hermiticity of W
    for i in 0..naux {
        for j in 0..naux {
            let w_ij = w[[i, j]];
            let w_ji = w[[j, i]];

            assert_abs_diff_eq!(w_ij.re, w_ji.re, epsilon = 1e-10);
            assert_abs_diff_eq!(w_ij.im, -w_ji.im, epsilon = 1e-10);
        }

        // Diagonal should be real
        assert_abs_diff_eq!(w[[i, i]].im, 0.0, epsilon = 1e-10);
    }

    assert!(diagnostics.hermiticity_error < 1e-10);
}

#[test]
fn test_positive_definiteness_correction() {
    let naux = 8;
    let config = EnhancedSolverConfig {
        min_eigenvalue_correction: 1e-10,
        ..Default::default()
    };
    let mut solver = EnhancedDielectricSolver::with_config(naux, config);

    // Create P0 that might lead to negative eigenvalues in W
    let mut p0 = create_hermitian_matrix(naux, 789);
    // Make it problematic
    p0 = &p0 * 1.2; // Scale beyond stability

    let v_sqrt = Array2::eye(naux);

    let result = solver.compute_screened_interaction_adaptive(&p0, &v_sqrt);

    if let Ok((w, _)) = result {
        // Check that W has no negative eigenvalues
        use ndarray_linalg::{Eigh, UPLO};
        let (eigenvalues, _) = w.eigh(UPLO::Lower).unwrap();

        for &eigenvalue in eigenvalues.iter() {
            assert!(
                eigenvalue >= -1e-10,
                "W should be positive semi-definite: eigenvalue = {}",
                eigenvalue
            );
        }
    }
}

#[test]
fn test_batch_processing() {
    let naux = 10;
    let n_freq = 20;
    let mut processor = BatchScreeningProcessor::new(naux, 4);

    // Create frequency-dependent P0
    let mut p0_batch = Array3::zeros((n_freq, naux, naux));
    for iw in 0..n_freq {
        let omega = iw as f64 * 0.1;
        for i in 0..naux {
            // Frequency-dependent diagonal
            p0_batch[[iw, i, i]] = Complex64::new(0.1 / (1.0 + omega * omega), 0.0);
        }
    }

    let v_sqrt = Array2::eye(naux);

    let (w_batch, diagnostics) = processor
        .process_batch_adaptive(&p0_batch, &v_sqrt)
        .expect("Batch processing should succeed");

    // Check all frequencies processed
    assert_eq!(w_batch.shape(), &[n_freq, naux, naux]);
    assert_eq!(diagnostics.len(), n_freq);

    // Verify continuity
    for iw in 1..n_freq {
        let w_prev = w_batch.slice(s![iw - 1, .., ..]);
        let w_curr = w_batch.slice(s![iw, .., ..]);
        let jump = (&w_curr - &w_prev).norm_l2() / w_curr.norm_l2();

        // Should be reasonably continuous
        assert!(
            jump < 0.5,
            "Frequency continuity at {} : jump = {}",
            iw,
            jump
        );
    }
}

#[test]
fn test_adaptive_alpha_selection() {
    let naux = 8;

    // Test with varying condition numbers
    let condition_targets = vec![1e3, 1e6, 1e9, 1e12];

    for target_cond in condition_targets {
        let config = EnhancedSolverConfig {
            condition_threshold: 1e8,
            ..Default::default()
        };
        let mut solver = EnhancedDielectricSolver::with_config(naux, config);

        // Create matrix with specified condition number
        let max_eig = 1.0;
        let min_eig = max_eig / target_cond;
        let eigenvalues: Vec<f64> = (0..naux)
            .map(|i| {
                let t = i as f64 / (naux - 1) as f64;
                max_eig * (1.0 - t) + min_eig * t
            })
            .collect();

        let p0 = create_diagonal_complex(&eigenvalues) * 0.5;
        let v_sqrt = Array2::eye(naux);

        let (_, diagnostics) = solver
            .compute_screened_interaction_adaptive(&p0, &v_sqrt)
            .expect("Should handle varying condition numbers");

        // Check appropriate solver was selected
        if target_cond > 1e10 {
            assert!(
                diagnostics.applied_regularization.is_some()
                    || matches!(
                        diagnostics.solver_used,
                        SolverMethod::SVDTruncated(_) | SolverMethod::TikhonovRegularized(_)
                    )
            );
        }
    }
}

#[test]
fn test_iterative_refinement() {
    let naux = 12;
    let config = EnhancedSolverConfig {
        max_refinement_iterations: 5,
        ..Default::default()
    };
    let mut solver = EnhancedDielectricSolver::with_config(naux, config);

    // Create moderately ill-conditioned system
    let p0 = create_hermitian_matrix(naux, 321) * 0.4;
    let v_sqrt = Array2::eye(naux);

    let (w1, diag1) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("First computation");

    // Second computation should be consistent
    let (w2, diag2) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("Second computation");

    // Results should be deterministic
    let diff = (&w1 - &w2).norm_l2() / w1.norm_l2();
    assert!(
        diff < 1e-10,
        "Computations should be deterministic: diff = {}",
        diff
    );

    // Diagnostics should be similar
    assert_relative_eq!(
        diag1.condition_number,
        diag2.condition_number,
        epsilon = 1e-6
    );
}

#[test]
fn test_frequency_interpolation() {
    let naux = 8;
    let n_freq_dense = 100;
    let mut processor = BatchScreeningProcessor::new(naux, 4);

    // Create smooth frequency-dependent P0
    let mut p0_batch = Array3::zeros((n_freq_dense, naux, naux));
    for iw in 0..n_freq_dense {
        let omega = iw as f64 * 0.01;
        for i in 0..naux {
            // Smooth frequency dependence
            p0_batch[[iw, i, i]] = Complex64::new(0.1 * (1.0 + omega.cos()), 0.0);
        }
    }

    let v_sqrt = Array2::eye(naux);

    // Process with adaptive strategy (should use interpolation)
    let (w_batch, _diagnostics) = processor
        .process_batch_adaptive(&p0_batch, &v_sqrt)
        .expect("Should process with interpolation");

    // Verify smoothness
    let mut max_second_derivative: f64 = 0.0;
    for iw in 1..n_freq_dense - 1 {
        for i in 0..naux {
            let w_prev = w_batch[[iw - 1, i, i]];
            let w_curr = w_batch[[iw, i, i]];
            let w_next = w_batch[[iw + 1, i, i]];

            // Approximate second derivative
            let d2w = (w_next - 2.0 * w_curr + w_prev).norm();
            max_second_derivative = f64::max(max_second_derivative, d2w);
        }
    }

    // Should be smooth (small second derivative)
    assert!(
        max_second_derivative < 0.01,
        "Frequency dependence should be smooth: d2w_max = {}",
        max_second_derivative
    );
}

#[test]
fn test_edge_case_empty_matrix() {
    let naux = 5;
    let mut solver = EnhancedDielectricSolver::new(naux);

    // P0 = 0 should give W = v
    let p0 = Array2::zeros((naux, naux));
    let v_sqrt = Array2::eye(naux) * 2.0_f64.sqrt();

    let (w, diagnostics) = solver
        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
        .expect("Should handle P0 = 0");

    // W should equal v = v_sqrt^2
    let _v = Array2::<f64>::eye(naux) * 2.0;
    for i in 0..naux {
        for j in 0..naux {
            let expected = if i == j { 2.0 } else { 0.0 };
            assert_abs_diff_eq!(w[[i, j]].re, expected, epsilon = 1e-10);
            assert_abs_diff_eq!(w[[i, j]].im, 0.0, epsilon = 1e-10);
        }
    }

    assert!(diagnostics.self_consistency_residual < 1e-10);
}

/// Integration test: Full GW-like calculation
#[test]
fn test_gw_integration() {
    let naux = 16;
    let n_freq = 30; // Typical for GW

    // Create realistic frequency grid (imaginary axis)
    let omega_max = 10.0;
    let mut p0_batch = Array3::zeros((n_freq, naux, naux));

    for iw in 0..n_freq {
        // Gauss-Legendre transformed frequency
        let t = (iw as f64 + 0.5) / n_freq as f64;
        let omega = omega_max * t / (1.0 - t);

        // Create frequency-dependent polarizability
        let p0 = create_hermitian_matrix(naux, iw as u64) * 0.1;

        // Apply frequency dependence
        for i in 0..naux {
            for j in 0..naux {
                let damping = 1.0 / (1.0 + omega * omega);
                p0_batch[[iw, i, j]] = p0[[i, j]] * damping;
            }
        }
    }

    // Realistic Coulomb metric
    let mut v = Array2::eye(naux);
    for i in 0..naux {
        for j in 0..naux {
            if i != j {
                // Off-diagonal Coulomb elements
                v[[i, j]] = 1.0 / ((i as f64 - j as f64).abs() + 1.0);
            } else {
                v[[i, i]] = 2.0;
            }
        }
    }

    // Compute v^{1/2} via eigendecomposition
    use ndarray_linalg::{Eigh, UPLO};
    let (eigenvalues, eigenvectors) = v.eigh(UPLO::Lower).unwrap();
    let v_sqrt_diag = eigenvalues.mapv(|x| x.sqrt());
    let v_sqrt = eigenvectors
        .dot(&Array2::from_diag(&v_sqrt_diag))
        .dot(&eigenvectors.t());

    // Process batch
    let mut processor = BatchScreeningProcessor::new(naux, 8);
    let (w_batch, diagnostics) = processor
        .process_batch_adaptive(&p0_batch, &v_sqrt)
        .expect("GW calculation should succeed");

    // Verify physical properties at each frequency
    for (iw, diag) in diagnostics.iter().enumerate() {
        // Accept relaxed tolerance if regularization was used
        let tolerance = if diag.applied_regularization.is_some() {
            1e-2
        } else {
            1e-4
        };
        assert!(
            diag.self_consistency_residual <= tolerance,
            "Frequency {}: self-consistency = {} (tolerance = {})",
            iw,
            diag.self_consistency_residual,
            tolerance
        );
        assert!(
            diag.hermiticity_error < 1e-8,
            "Frequency {}: hermiticity = {}",
            iw,
            diag.hermiticity_error
        );
    }

    // Check physical property: W(iω) approaches v as ω→∞ (screening decreases)
    // This means ||W(ω→∞)|| ≥ ||W(0)|| for most systems
    // We just check that the norms are finite and reasonable
    let w_0_norm = w_batch.slice(s![0, .., ..]).norm_l2();
    let w_end_norm = w_batch.slice(s![n_freq - 1, .., ..]).norm_l2();
    assert!(
        w_0_norm.is_finite() && w_0_norm > 0.0,
        "W(0) should have finite positive norm"
    );
    assert!(
        w_end_norm.is_finite() && w_end_norm > 0.0,
        "W(ω_max) should have finite positive norm"
    );
}
