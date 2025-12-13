//! Benchmarks for the screening module
//!
//! This benchmark suite measures:
//! - Different solver backends (LU, Cholesky, SVD)
//! - SIMD optimization speedups
//! - Memory usage and allocation patterns
//! - Scaling with matrix size

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::{Array2, Array3};
use num_complex::Complex64;
use quasix_core::dielectric::screening::{
    DielectricSolver, ScreenedInteraction, SolverBackend, SolverConfig, SolverType,
};
use std::hint::black_box;

/// Create a test polarizability matrix that is Hermitian
fn create_test_p0(naux: usize) -> Array2<Complex64> {
    let mut p0 = Array2::<Complex64>::zeros((naux, naux));

    // Create a Hermitian matrix with controlled eigenvalue spectrum
    for i in 0..naux {
        // Diagonal elements (real)
        p0[[i, i]] = Complex64::new(0.5 + 0.4 * (i as f64 / naux as f64), 0.0);

        // Off-diagonal elements
        for j in i + 1..naux {
            let val = Complex64::new(
                0.1 * ((i + j) as f64).sin() / (1.0 + (i as i32 - j as i32).abs() as f64),
                0.05 * ((i * j) as f64).cos() / (1.0 + (i as i32 - j as i32).abs() as f64),
            );
            p0[[i, j]] = val;
            p0[[j, i]] = val.conj();
        }
    }

    p0
}

/// Create a test Coulomb metric v^{1/2}
fn create_test_vsqrt(naux: usize) -> Array2<f64> {
    let mut vsqrt = Array2::<f64>::eye(naux);

    // Add some off-diagonal elements to make it more realistic
    for i in 0..naux {
        for j in i + 1..naux.min(i + 3) {
            let val = 0.1 / (1.0 + (i as i32 - j as i32).abs() as f64);
            vsqrt[[i, j]] = val;
            vsqrt[[j, i]] = val;
        }
    }

    // Ensure positive definiteness
    vsqrt = vsqrt.dot(&vsqrt.t());

    // Take square root (simplified - in practice use eigendecomposition)
    vsqrt.mapv(|x| x.sqrt())
}

/// Benchmark different solver backends
fn bench_solver_backends(c: &mut Criterion) {
    let mut group = c.benchmark_group("solver_backends");

    for naux in [50, 100, 200].iter() {
        let p0 = create_test_p0(*naux);
        let vsqrt = create_test_vsqrt(*naux);

        // Benchmark LU backend
        group.bench_with_input(BenchmarkId::new("LU", naux), naux, |b, &size| {
            let solver =
                DielectricSolver::with_backend(size, SolverType::Direct, SolverBackend::LU);
            b.iter(|| {
                let _ = solver.compute_screened_interaction(black_box(&p0), black_box(&vsqrt));
            });
        });

        // Benchmark Cholesky backend
        group.bench_with_input(BenchmarkId::new("Cholesky", naux), naux, |b, &size| {
            let solver =
                DielectricSolver::with_backend(size, SolverType::Direct, SolverBackend::Cholesky);
            b.iter(|| {
                let _ = solver.compute_screened_interaction(black_box(&p0), black_box(&vsqrt));
            });
        });

        // Benchmark SVD backend
        group.bench_with_input(BenchmarkId::new("SVD", naux), naux, |b, &size| {
            let solver =
                DielectricSolver::with_backend(size, SolverType::Direct, SolverBackend::SVD);
            b.iter(|| {
                let _ = solver.compute_screened_interaction(black_box(&p0), black_box(&vsqrt));
            });
        });
    }

    group.finish();
}

/// Benchmark SIMD optimizations
fn bench_simd_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_operations");

    // Test different matrix sizes to see SIMD benefits
    for naux in [64, 128, 256, 512].iter() {
        let p0 = create_test_p0(*naux);
        let vsqrt = create_test_vsqrt(*naux);

        // Benchmark with default block size (optimized for cache)
        group.bench_with_input(BenchmarkId::new("default_block", naux), naux, |b, &size| {
            let solver = DielectricSolver::new(size, SolverType::Direct);
            b.iter(|| {
                let m = solver
                    .build_symmetrized_dielectric(black_box(&p0), black_box(&vsqrt))
                    .unwrap();
                black_box(m);
            });
        });

        // Benchmark with small block size (less cache optimization)
        group.bench_with_input(BenchmarkId::new("small_block", naux), naux, |b, &size| {
            let config = SolverConfig {
                block_size: 32,
                ..Default::default()
            };
            let solver = DielectricSolver::with_config(
                size,
                SolverType::Direct,
                SolverBackend::Auto,
                config,
            );
            b.iter(|| {
                let m = solver
                    .build_symmetrized_dielectric(black_box(&p0), black_box(&vsqrt))
                    .unwrap();
                black_box(m);
            });
        });
    }

    group.finish();
}

/// Benchmark batch operations for multiple frequencies
fn bench_batch_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("batch_operations");

    let naux = 100;
    let vsqrt = create_test_vsqrt(naux);

    for n_freq in [10, 50, 100].iter() {
        // Create batch of P0 matrices
        let mut p0_batch = Array3::<Complex64>::zeros((*n_freq, naux, naux));
        for i in 0..*n_freq {
            let p0 = create_test_p0(naux);
            p0_batch.slice_mut(ndarray::s![i, .., ..]).assign(&p0);
        }

        group.bench_with_input(
            BenchmarkId::new("batch_computation", n_freq),
            n_freq,
            |b, &_| {
                let solver = DielectricSolver::new(naux, SolverType::Direct);
                b.iter(|| {
                    let _ = solver.compute_screened_interaction_batch(
                        black_box(&p0_batch),
                        black_box(&vsqrt),
                    );
                });
            },
        );
    }

    group.finish();
}

/// Benchmark memory usage patterns
fn bench_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_patterns");

    // Benchmark workspace allocation overhead
    for naux in [100, 200, 400].iter() {
        group.bench_with_input(
            BenchmarkId::new("solver_creation", naux),
            naux,
            |b, &size| {
                b.iter(|| {
                    let solver = DielectricSolver::new(size, SolverType::Direct);
                    black_box(solver);
                });
            },
        );
    }

    // Benchmark hermitianization overhead
    for naux in [50, 100, 200].iter() {
        let mut matrix = create_test_p0(*naux);
        // Make it slightly non-Hermitian
        for i in 0..*naux {
            for j in i + 1..*naux {
                matrix[[i, j]] += Complex64::new(0.001, 0.001);
            }
        }

        group.bench_with_input(BenchmarkId::new("hermitianize", naux), naux, |b, &_size| {
            use quasix_core::dielectric::simd_ops::hermitianize_simd;
            b.iter(|| {
                let h = hermitianize_simd(black_box(&matrix));
                black_box(h);
            });
        });
    }

    group.finish();
}

/// Benchmark condition number estimation
fn bench_condition_estimation(c: &mut Criterion) {
    let mut group = c.benchmark_group("condition_estimation");

    for naux in [50, 100, 150].iter() {
        let p0 = create_test_p0(*naux);
        let vsqrt = create_test_vsqrt(*naux);

        // Well-conditioned case
        group.bench_with_input(
            BenchmarkId::new("well_conditioned", naux),
            naux,
            |b, &size| {
                let interaction = ScreenedInteraction::new(size);
                let w = interaction.compute(&p0, &vsqrt).unwrap();
                b.iter(|| {
                    let cond = interaction.condition_number(black_box(&w)).unwrap();
                    black_box(cond);
                });
            },
        );

        // Ill-conditioned case (add near-singular perturbation)
        let mut p0_ill = p0.clone();
        for i in 0..*naux {
            p0_ill[[i, i]] *= 1e-8; // Make nearly singular
        }

        group.bench_with_input(
            BenchmarkId::new("ill_conditioned", naux),
            naux,
            |b, &size| {
                let interaction = ScreenedInteraction::new(size);
                // Don't compute W since it might fail, just test condition estimation
                b.iter(|| {
                    let _ = interaction.condition_number(black_box(&p0_ill));
                });
            },
        );
    }

    group.finish();
}

/// Benchmark self-consistency verification
fn bench_self_consistency(c: &mut Criterion) {
    let mut group = c.benchmark_group("self_consistency");

    for naux in [50, 100, 200].iter() {
        let p0 = create_test_p0(*naux);
        let vsqrt = create_test_vsqrt(*naux);
        let v = vsqrt.dot(&vsqrt.t()); // Reconstruct full V

        let interaction = ScreenedInteraction::new(*naux);
        let w = interaction.compute(&p0, &vsqrt).unwrap();

        group.bench_with_input(BenchmarkId::new("verify", naux), naux, |b, &_| {
            b.iter(|| {
                let error = interaction
                    .check_self_consistency(black_box(&w), black_box(&p0), black_box(&v))
                    .unwrap();
                black_box(error);
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_solver_backends,
    bench_simd_operations,
    bench_batch_operations,
    bench_memory_patterns,
    bench_condition_estimation,
    bench_self_consistency
);
criterion_main!(benches);
