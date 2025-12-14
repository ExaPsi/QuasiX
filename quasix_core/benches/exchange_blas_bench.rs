//! Benchmark for BLAS-optimized exchange self-energy calculation
//!
//! Compares performance of:
//! 1. Original ndarray implementation
//! 2. Optimized direct BLAS calls

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::{Array2, Array3};
use quasix_core::selfenergy::ExchangeSelfEnergyRI;
use rand::prelude::*;
use rand::rng;
use std::hint::black_box;

/// Generate test data for benchmarking
fn generate_test_data(nbasis: usize, _nocc: usize, naux: usize) -> (Array3<f64>, Array2<f64>) {
    use rand_distr::StandardNormal;

    let mut rng = rng();

    // Generate random DF integrals
    let mut df_3c_mo = Array3::<f64>::zeros((nbasis, nbasis, naux));
    for m in 0..nbasis {
        for n in 0..=m {
            for p in 0..naux {
                let val: f64 = rng.sample(StandardNormal);
                df_3c_mo[[m, n, p]] = val * 0.1;
                df_3c_mo[[n, m, p]] = val * 0.1; // Symmetry
            }
        }
    }

    // Generate positive definite metric inverse
    let mut metric_inv = Array2::<f64>::zeros((naux, naux));
    for i in 0..naux {
        for j in 0..=i {
            let val: f64 = if i == j {
                1.0 + rng.random::<f64>() * 0.5
            } else {
                rng.sample::<f64, _>(StandardNormal) * 0.1
            };
            metric_inv[[i, j]] = val;
            metric_inv[[j, i]] = val;
        }
    }

    (df_3c_mo, metric_inv)
}

/// Benchmark exchange matrix computation with different methods
fn benchmark_exchange_methods(c: &mut Criterion) {
    let mut group = c.benchmark_group("exchange_self_energy_blas");

    // Test different system sizes
    let test_cases = vec![
        ("small", 20, 5, 60),
        ("medium", 50, 10, 150),
        ("large", 100, 20, 300),
    ];

    for (name, nbasis, nocc, naux) in test_cases {
        let (df_3c_mo, metric_inv) = generate_test_data(nbasis, nocc, naux);

        group.bench_function(BenchmarkId::new("optimized_blas", name), |b| {
            b.iter(|| {
                let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);
                let _sigma_x = calc
                    .compute_exchange_matrix_ri(black_box(&df_3c_mo), black_box(&metric_inv))
                    .unwrap();
            });
        });
    }

    group.finish();
}

/// Benchmark diagonal computation
fn benchmark_diagonal_methods(c: &mut Criterion) {
    let mut group = c.benchmark_group("exchange_diagonal_blas");

    let test_cases = vec![
        ("small", 20, 5, 60),
        ("medium", 50, 10, 150),
        ("large", 100, 20, 300),
    ];

    for (name, nbasis, nocc, naux) in test_cases {
        let (df_3c_mo, metric_inv) = generate_test_data(nbasis, nocc, naux);

        group.bench_function(BenchmarkId::new("optimized_diagonal", name), |b| {
            b.iter(|| {
                let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);
                let _diag = calc
                    .compute_exchange_diagonal_ri(black_box(&df_3c_mo), black_box(&metric_inv))
                    .unwrap();
            });
        });
    }

    group.finish();
}

/// Benchmark parallel scaling
fn benchmark_parallel_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("exchange_parallel_scaling");

    let nbasis = 100;
    let nocc = 20;
    let naux = 300;
    let (df_3c_mo, metric_inv) = generate_test_data(nbasis, nocc, naux);

    // Test different thread configurations
    let thread_counts = vec![1, 2, 4, 8];

    for num_threads in thread_counts {
        group.bench_function(BenchmarkId::new("threads", num_threads), |b| {
            // Set thread count for Rayon
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .ok();

            b.iter(|| {
                let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);
                let _sigma_x = calc
                    .compute_exchange_matrix_ri(black_box(&df_3c_mo), black_box(&metric_inv))
                    .unwrap();
            });
        });
    }

    group.finish();
}

/// Benchmark cache efficiency with different block sizes
fn benchmark_block_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("exchange_block_sizes");

    let nbasis = 128;
    let nocc = 32;
    let naux = 384;
    let (df_3c_mo, metric_inv) = generate_test_data(nbasis, nocc, naux);

    // This would require modifying the implementation to accept block size parameter
    // For now, just test the default implementation
    group.bench_function("default_block_size", |b| {
        b.iter(|| {
            let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);
            let _sigma_x = calc
                .compute_exchange_matrix_ri(black_box(&df_3c_mo), black_box(&metric_inv))
                .unwrap();
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_exchange_methods,
    benchmark_diagonal_methods,
    benchmark_parallel_scaling,
    benchmark_block_sizes
);
criterion_main!(benches);
