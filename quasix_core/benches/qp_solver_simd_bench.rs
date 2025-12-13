//! Benchmark for SIMD-optimized QP Solver
//!
//! This benchmark demonstrates the performance improvements from SIMD optimization
//! in the quasiparticle solver, focusing on derivative computation and Z-factor
//! calculations which are the hot paths in the Newton-Raphson iteration.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array1;
use num_complex::Complex64 as Complex;
use quasix_core::common::Result;
use quasix_core::qp::solver::{QPEquationSolver, QPSolverConfig};
use std::hint::black_box;

/// Mock self-energy function for benchmarking
fn mock_sigma_function(orbital_idx: usize, energy: f64) -> Result<Complex> {
    // Simple model: Σ(ω) = Σ_x + α·ω/(ω² + β²)
    let sigma_x = -0.5 - 0.01 * orbital_idx as f64;
    let alpha = 0.1;
    let beta = 1.0;
    let sigma_c = alpha * energy / (energy * energy + beta * beta);
    Ok(Complex::new(sigma_x + sigma_c, 0.0))
}

/// Benchmark derivative computation (scalar vs SIMD)
fn bench_derivative_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("derivative_computation");

    // Test function for derivative
    let test_func = |x: f64| -> Result<f64> { Ok(x * x + 2.0 * x.sin() - 3.0 * x.cos()) };

    let config = QPSolverConfig::default();
    let solver = QPEquationSolver::new(config);

    // Benchmark at different evaluation points
    for x in &[0.5, 1.0, 2.0, 5.0] {
        group.bench_with_input(BenchmarkId::new("compute_derivative", x), x, |b, &x_val| {
            b.iter(|| {
                // This will use SIMD when available
                let _ = black_box(solver.compute_derivative_simple(&test_func, x_val));
            });
        });
    }

    group.finish();
}

/// Benchmark Z-factor computation (scalar vs SIMD)
fn bench_z_factor_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("z_factor_computation");

    let config = QPSolverConfig::default();
    let solver = QPEquationSolver::new(config);

    // Different orbital indices to test
    for orbital_idx in &[0, 5, 10, 20] {
        let sigma_func = |e: f64| -> Result<Complex> { mock_sigma_function(*orbital_idx, e) };

        group.bench_with_input(
            BenchmarkId::new("compute_z_factor", orbital_idx),
            orbital_idx,
            |b, &idx| {
                b.iter(|| {
                    let qp_energy = -0.5 + 0.1 * idx as f64;
                    let _ = black_box(solver.compute_z_factor(idx, qp_energy, &sigma_func));
                });
            },
        );
    }

    group.finish();
}

/// Benchmark full QP solver for single orbital
fn bench_single_orbital_solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_orbital_solve");

    let config = QPSolverConfig {
        max_newton_iterations: 20,
        use_richardson: true,
        use_line_search: true,
        ..Default::default()
    };
    let solver = QPEquationSolver::new(config);

    // Test different initial energies
    for (idx, epsilon) in [(0, -0.5), (5, -0.3), (10, -0.1)].iter() {
        let sigma_func = |e: f64| -> Result<Complex> { mock_sigma_function(*idx, e) };
        let vxc = -0.4;

        group.bench_with_input(
            BenchmarkId::new("solve_qp_newton", idx),
            idx,
            |b, &orbital_idx| {
                b.iter(|| {
                    let _ = black_box(solver.solve_qp_equation_newton(
                        orbital_idx,
                        *epsilon,
                        &sigma_func,
                        vxc,
                    ));
                });
            },
        );
    }

    group.finish();
}

/// Benchmark parallel orbital processing
fn bench_parallel_orbital_processing(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_orbital_processing");

    // Different system sizes
    for n_mo in &[10, 50, 100, 200] {
        let mo_energies = Array1::from_vec((0..*n_mo).map(|i| -1.0 + 0.01 * i as f64).collect());
        let mo_occ = Array1::from_vec(
            (0..*n_mo)
                .map(|i| if i < n_mo / 2 { 2.0 } else { 0.0 })
                .collect(),
        );
        let vxc_diagonal = Array1::from_vec(vec![-0.4; *n_mo]);

        // Test with different thread counts
        for n_threads in &[1, 2, 4, 8] {
            if *n_threads > rayon::current_num_threads() {
                continue;
            }

            let config = QPSolverConfig {
                n_threads: Some(*n_threads),
                max_newton_iterations: 10, // Reduce for benchmark
                ..Default::default()
            };
            let solver = QPEquationSolver::new(config);

            group.bench_with_input(
                BenchmarkId::new(format!("n_mo_{}_threads", n_threads), n_mo),
                n_mo,
                |b, _| {
                    b.iter(|| {
                        let _ = black_box(solver.solve_all_orbitals(
                            &mo_energies,
                            &mo_occ,
                            mock_sigma_function,
                            &vxc_diagonal,
                        ));
                    });
                },
            );
        }
    }

    group.finish();
}

/// Benchmark Richardson extrapolation (scalar vs SIMD)
fn bench_richardson_extrapolation(c: &mut Criterion) {
    let mut group = c.benchmark_group("richardson_extrapolation");

    let config = QPSolverConfig {
        use_richardson: true,
        ..Default::default()
    };
    let solver = QPEquationSolver::new(config);

    // Test function with known derivative
    let test_func = |x: f64| -> Result<f64> { Ok(x.powi(3) + 2.0 * x.powi(2) - 5.0 * x + 3.0) };

    for x in &[0.5, 1.0, 2.0, 5.0, 10.0] {
        group.bench_with_input(BenchmarkId::new("richardson", x), x, |b, &x_val| {
            b.iter(|| {
                let _ = black_box(solver.compute_derivative_richardson(&test_func, x_val));
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_derivative_computation,
    bench_z_factor_computation,
    bench_single_orbital_solve,
    bench_parallel_orbital_processing,
    bench_richardson_extrapolation
);

criterion_main!(benches);
