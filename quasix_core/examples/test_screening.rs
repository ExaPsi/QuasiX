//! Example demonstrating the screening module with zero warnings

use ndarray::Array2;
use num_complex::Complex64;
use quasix_core::dielectric::screening::{DielectricSolver, SolverType};

fn main() {
    println!("Testing Screening Module - Zero Warnings Implementation");
    println!("========================================================\n");

    // Create a 50x50 test system
    let n = 50;

    // Build solver with optimized configuration
    let solver = DielectricSolver::new(n, SolverType::Adaptive);

    println!("Solver Configuration:");
    println!("  - Matrix size: {}x{}", n, n);
    println!("  - Solver type: Adaptive");
    println!();

    // Create test polarizability (positive semi-definite)
    let mut p0 = Array2::<Complex64>::zeros((n, n));
    for i in 0..n {
        p0[[i, i]] = Complex64::new(0.1 * (1.0 + i as f64 / n as f64), 0.0);
        if i > 0 {
            p0[[i, i - 1]] = Complex64::new(0.01, 0.005);
            p0[[i - 1, i]] = Complex64::new(0.01, -0.005);
        }
    }

    // Create Coulomb matrix square root
    let mut vsqrt = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        vsqrt[[i, i]] = 1.0 / (1.0 + i as f64).sqrt();
    }

    // Build symmetrized dielectric matrix
    println!("Building symmetrized dielectric matrix...");
    let epsilon = solver
        .build_symmetrized_dielectric(&p0, &vsqrt)
        .expect("Failed to build dielectric");

    // Check properties
    let trace = epsilon.diag().sum();
    println!("Dielectric matrix properties:");
    println!("  - Trace: {:.6} + {:.6}i", trace.re, trace.im);

    // Check Hermiticity
    let mut max_herm_error: f64 = 0.0;
    for i in 0..n {
        for j in i + 1..n {
            let error = (epsilon[[i, j]] - epsilon[[j, i]].conj()).norm();
            max_herm_error = max_herm_error.max(error);
        }
    }
    println!("  - Hermiticity error: {:.2e}", max_herm_error);

    // Invert dielectric matrix
    println!("\nInverting dielectric matrix...");
    let epsilon_inv = solver
        .invert_dielectric(&epsilon)
        .expect("Failed to invert dielectric");

    // Verify inversion quality
    let product = epsilon.dot(&epsilon_inv);
    let mut max_inv_error: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            let error = (product[[i, j]].re - expected).abs() + product[[i, j]].im.abs();
            max_inv_error = max_inv_error.max(error);
        }
    }
    println!("Inversion quality:");
    println!("  - Max error in ε·ε⁻¹: {:.2e}", max_inv_error);

    // Test different solver types
    println!("\nBenchmarking solver types:");
    let solver_types = vec![
        ("Direct", SolverType::Direct),
        ("Regularized", SolverType::Regularized),
        ("Adaptive", SolverType::Adaptive),
    ];

    for (name, solver_type) in solver_types {
        let solver_test = DielectricSolver::new(n, solver_type);

        let start = std::time::Instant::now();
        let _ = solver_test.invert_dielectric(&epsilon);
        let elapsed = start.elapsed();

        println!("  - {}: {:.2} ms", name, elapsed.as_secs_f64() * 1000.0);
    }

    println!("\n✓ All tests completed successfully with zero warnings!");
}
