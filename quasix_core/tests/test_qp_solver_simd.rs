//! Test SIMD optimizations in the QP solver

use approx::assert_relative_eq;
use ndarray::Array1;
use num_complex::Complex64 as Complex;
use quasix_core::common::Result;
use quasix_core::qp::solver::{QPEquationSolver, QPSolverConfig};

/// Test that SIMD derivative computation matches scalar version
#[test]
fn test_simd_derivative_accuracy() {
    let config = QPSolverConfig::default();
    let solver = QPEquationSolver::new(config);

    // Test function: f(x) = x³ + 2x² - 5x + 3
    // Derivative: f'(x) = 3x² + 4x - 5
    let f = |x: f64| -> Result<f64> { Ok(x.powi(3) + 2.0 * x.powi(2) - 5.0 * x + 3.0) };

    // Test at various points
    let test_points = vec![0.5, 1.0, 1.5, 2.0, 3.0, 5.0];

    for x in test_points {
        // Compute derivative using our solver
        let computed = solver.compute_derivative_simple(&f, x).unwrap();

        // Analytical derivative
        let expected = 3.0 * x.powi(2) + 4.0 * x - 5.0;

        // Check accuracy (should be within numerical error)
        assert_relative_eq!(computed, expected, epsilon = 1e-6);
    }
}

/// Test Richardson extrapolation with SIMD
#[test]
fn test_simd_richardson_accuracy() {
    let config = QPSolverConfig {
        use_richardson: true,
        derivative_delta: 1e-5,
        ..Default::default()
    };
    let solver = QPEquationSolver::new(config);

    // Test function: f(x) = sin(x) * exp(x/2)
    // Derivative: f'(x) = cos(x) * exp(x/2) + sin(x) * exp(x/2) / 2
    let f = |x: f64| -> Result<f64> { Ok(x.sin() * (x / 2.0).exp()) };

    let test_points = vec![0.1, 0.5, 1.0, 1.5, 2.0];

    for x in test_points {
        // Compute derivative using Richardson extrapolation
        let computed = solver.compute_derivative_richardson(&f, x).unwrap();

        // Analytical derivative
        let expected = x.cos() * (x / 2.0).exp() + x.sin() * (x / 2.0).exp() / 2.0;

        // Richardson should give better accuracy
        assert_relative_eq!(computed, expected, epsilon = 1e-8);
    }
}

/// Test Z-factor computation with SIMD
#[test]
fn test_simd_z_factor_computation() {
    let config = QPSolverConfig::default();
    let solver = QPEquationSolver::new(config);

    // Simple self-energy model for testing
    let sigma_func = |e: f64| -> Result<Complex> {
        // Σ(ω) = Σ_x + α·ω/(ω² + β²)
        let sigma_x = -0.5;
        let alpha = 0.1;
        let beta = 1.0;
        let sigma_c = alpha * e / (e * e + beta * beta);
        Ok(Complex::new(sigma_x + sigma_c, 0.0))
    };

    // Test at different QP energies
    let test_energies = vec![-0.5, -0.3, 0.0, 0.3, 0.5];

    for qp_energy in test_energies {
        let z_factor = solver.compute_z_factor(0, qp_energy, &sigma_func).unwrap();

        // Z-factor should be physical
        assert!(
            z_factor > 0.0 && z_factor < 1.0,
            "Unphysical Z-factor {} at E={}",
            z_factor,
            qp_energy
        );

        // Compute derivative manually for verification
        let h = 1e-5;
        let sigma_plus = sigma_func(qp_energy + h).unwrap();
        let sigma_minus = sigma_func(qp_energy - h).unwrap();
        let d_sigma_d_omega = (sigma_plus.re - sigma_minus.re) / (2.0 * h);
        let expected_z = 1.0 / (1.0 - d_sigma_d_omega);

        // Account for bounds enforcement
        let expected_z_bounded = expected_z.clamp(0.01, 0.99);

        assert_relative_eq!(z_factor, expected_z_bounded, epsilon = 1e-4);
    }
}

/// Test parallel orbital processing
#[test]
fn test_parallel_orbital_processing() {
    // Small test system
    let n_mo = 10;
    let mo_energies = Array1::from_vec((0..n_mo).map(|i| -1.0 + 0.1 * i as f64).collect());
    let mo_occ = Array1::from_vec(
        (0..n_mo)
            .map(|i| if i < n_mo / 2 { 2.0 } else { 0.0 })
            .collect(),
    );
    let vxc_diagonal = Array1::from_vec(vec![-0.4; n_mo]);

    // Simple self-energy function
    let sigma_func = |orbital_idx: usize, energy: f64| -> Result<Complex> {
        let sigma_x = -0.5 - 0.01 * orbital_idx as f64;
        let alpha = 0.1;
        let beta = 1.0;
        let sigma_c = alpha * energy / (energy * energy + beta * beta);
        Ok(Complex::new(sigma_x + sigma_c, 0.0))
    };

    // Test with different thread counts
    for n_threads in &[1, 2, 4] {
        let config = QPSolverConfig {
            n_threads: Some(*n_threads),
            max_newton_iterations: 20,
            ..Default::default()
        };
        let solver = QPEquationSolver::new(config);

        let solution = solver
            .solve_all_orbitals(&mo_energies, &mo_occ, sigma_func, &vxc_diagonal)
            .expect("Solver should succeed");

        // Check that all orbitals were processed
        assert_eq!(solution.qp_energies.len(), n_mo);
        assert_eq!(solution.z_factors.len(), n_mo);

        // Check physical constraints
        for i in 0..n_mo {
            assert!(
                solution.z_factors[i] > 0.0 && solution.z_factors[i] < 1.0,
                "Unphysical Z-factor for orbital {}",
                i
            );

            // QP energies should be reasonably close to MO energies
            assert!(
                (solution.qp_energies[i] - mo_energies[i]).abs() < 1.0,
                "QP energy too far from MO energy for orbital {}",
                i
            );
        }

        // Check convergence
        let converged_count = solution.convergence_flags.iter().filter(|&&c| c).count();
        assert!(
            converged_count > n_mo / 2,
            "Too few orbitals converged: {}/{}",
            converged_count,
            n_mo
        );
    }
}

/// Performance comparison test (not a benchmark, just correctness)
#[test]
fn test_simd_vs_scalar_consistency() {
    let config_scalar = QPSolverConfig {
        use_richardson: false, // Force simple derivative
        ..Default::default()
    };
    let config_simd = QPSolverConfig {
        use_richardson: true, // Use SIMD-optimized Richardson
        ..Default::default()
    };

    let solver_scalar = QPEquationSolver::new(config_scalar);
    let solver_simd = QPEquationSolver::new(config_simd);

    // Test function
    let f = |x: f64| -> Result<f64> { Ok(x.exp() * x.sin()) };

    let test_points = vec![0.5, 1.0, 1.5, 2.0];

    for x in test_points {
        let deriv_scalar = solver_scalar.compute_derivative_simple(&f, x).unwrap();
        let deriv_simd = solver_simd.compute_derivative_richardson(&f, x).unwrap();

        // SIMD should give similar or better accuracy
        // Analytical: f'(x) = exp(x) * (sin(x) + cos(x))
        let expected = x.exp() * (x.sin() + x.cos());

        let error_scalar = (deriv_scalar - expected).abs();
        let error_simd = (deriv_simd - expected).abs();

        // Both methods should be accurate, but Richardson may have different error characteristics
        // We just check that both are reasonably accurate
        assert!(
            error_simd < 1e-6,
            "SIMD error too large at x={}: error {:.2e}",
            x,
            error_simd
        );
        assert!(
            error_scalar < 1e-4,
            "Scalar error too large at x={}: error {:.2e}",
            x,
            error_scalar
        );
    }
}

#[test]
fn test_simd_workspace_alignment() {
    use quasix_core::qp::solver::SimdWorkspace;

    // Test that workspace is properly aligned
    let workspace = SimdWorkspace::new(100);

    // Check alignment (should be at least 8-byte aligned for f64)
    let ptr_omega = workspace.omega_points.as_ptr() as usize;
    let ptr_real = workspace.sigma_real.as_ptr() as usize;
    let ptr_imag = workspace.sigma_imag.as_ptr() as usize;

    assert_eq!(ptr_omega % 8, 0, "omega_points not aligned");
    assert_eq!(ptr_real % 8, 0, "sigma_real not aligned");
    assert_eq!(ptr_imag % 8, 0, "sigma_imag not aligned");

    // For AVX, we'd prefer 32-byte alignment
    if cfg!(target_feature = "avx") {
        assert_eq!(ptr_omega % 32, 0, "omega_points not AVX aligned");
        assert_eq!(ptr_real % 32, 0, "sigma_real not AVX aligned");
        assert_eq!(ptr_imag % 32, 0, "sigma_imag not AVX aligned");
    }
}
