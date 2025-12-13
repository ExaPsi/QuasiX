//! Demonstration of QuasiX condition monitoring capabilities
//!
//! This example shows how the new O(N²) condition monitoring prevents
//! silent numerical failures in GW calculations.

use ndarray::{Array2, Array3};
use num_complex::Complex64;
use quasix_core::dielectric::screening::{DielectricSolver, SolverConfig, SolverType, SolverBackend};
use quasix_core::linalg::conditioning::{
    estimate_condition_power_iteration, ConditionConfig, analyze_condition
};
use std::time::Instant;

fn main() {
    println!("QuasiX Condition Monitoring Demo");
    println!("=================================\n");

    // Test different matrix sizes
    let sizes = vec![50, 100, 200];

    for naux in sizes {
        println!("Testing matrix size: {}x{}", naux, naux);

        // Create test matrices with varying condition numbers
        let matrices = vec![
            ("Well-conditioned", create_well_conditioned_matrix(naux)),
            ("Moderately ill-conditioned", create_moderate_ill_matrix(naux)),
            ("Severely ill-conditioned", create_severe_ill_matrix(naux)),
        ];

        for (name, matrix) in matrices {
            println!("\n  {}:", name);

            // Time the condition estimation
            let start = Instant::now();
            let condition = estimate_condition_power_iteration(&matrix, 30, 1e-6);
            let elapsed = start.elapsed();

            println!("    Condition number: {:.2e}", condition);
            println!("    Estimation time: {:.3}ms", elapsed.as_secs_f64() * 1000.0);

            // Test with full analysis
            let config = ConditionConfig::default();
            let analysis = analyze_condition(&matrix, &config);

            if let Some(reg) = analysis.suggested_regularization {
                println!("    ⚠️  Regularization suggested: λ={:.2e}", reg);
            } else {
                println!("    ✓  Matrix is well-conditioned");
            }
        }
        println!();
    }

    // Demonstrate in actual GW context
    println!("GW Screening Matrix Test");
    println!("------------------------\n");

    demo_gw_screening();
}

fn create_well_conditioned_matrix(n: usize) -> Array2<Complex64> {
    let mut matrix = Array2::zeros((n, n));
    for i in 0..n {
        matrix[[i, i]] = Complex64::new(1.0 + 0.1 * (i as f64 / n as f64), 0.0);
        if i > 0 {
            matrix[[i, i-1]] = Complex64::new(0.05, 0.0);
            matrix[[i-1, i]] = Complex64::new(0.05, 0.0);
        }
    }
    matrix
}

fn create_moderate_ill_matrix(n: usize) -> Array2<Complex64> {
    let mut matrix = Array2::zeros((n, n));
    for i in 0..n {
        // Exponentially decreasing diagonal
        matrix[[i, i]] = Complex64::new((10.0_f64).powi(-(i as i32) / 10), 0.0);
    }
    matrix
}

fn create_severe_ill_matrix(n: usize) -> Array2<Complex64> {
    // Near-singular matrix
    let mut matrix = Array2::zeros((n, n));
    matrix[[0, 0]] = Complex64::new(1.0, 0.0);
    for i in 1..n {
        matrix[[i, i]] = Complex64::new(1e-12, 0.0);
    }
    matrix
}

fn demo_gw_screening() {
    let naux = 100;

    // Create solver with monitoring enabled
    let config = SolverConfig {
        monitor_condition: true,
        use_cheap_estimate: true,
        adaptive_regularization: true,
        condition_threshold: 1e8,
        ..Default::default()
    };

    let solver = DielectricSolver::with_config(
        naux,
        SolverType::Adaptive,
        SolverBackend::Auto,
        config,
    );

    // Create test polarizability and Coulomb matrices
    let mut p0 = Array2::zeros((naux, naux));
    let mut v_sqrt = Array2::zeros((naux, naux));

    // Setup typical GW matrices
    for i in 0..naux {
        p0[[i, i]] = Complex64::new(-0.1 * (1.0 + i as f64 / naux as f64), 0.0);
        v_sqrt[[i, i]] = (1.0 / (1.0 + i as f64)).sqrt();
    }

    println!("Computing screened interaction W...");
    let start = Instant::now();

    match solver.compute_screened_interaction(&p0, &v_sqrt) {
        Ok(w) => {
            let elapsed = start.elapsed();
            println!("✓ Computation successful in {:.3}ms", elapsed.as_secs_f64() * 1000.0);

            // Check W properties
            let w_diag_mean: f64 = (0..naux).map(|i| w[[i, i]].re).sum::<f64>() / naux as f64;
            println!("  W diagonal mean: {:.6e}", w_diag_mean);

            // Estimate condition of W
            let w_condition = estimate_condition_power_iteration(&w, 20, 1e-4);
            println!("  W condition number: {:.2e}", w_condition);
        },
        Err(e) => {
            println!("✗ Computation failed: {}", e);
            println!("  This would have been a silent failure without monitoring!");
        }
    }

    // Test with ill-conditioned case
    println!("\nTesting ill-conditioned case...");
    for i in 0..naux {
        p0[[i, i]] = Complex64::new(-0.99, 0.0);  // Close to eigenvalue 1
    }

    match solver.compute_screened_interaction(&p0, &v_sqrt) {
        Ok(_) => {
            println!("✓ Adaptive regularization prevented failure");
        },
        Err(e) => {
            println!("✗ Even with regularization, matrix too singular: {}", e);
        }
    }
}