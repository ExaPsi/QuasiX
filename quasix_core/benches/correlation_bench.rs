//! Benchmark for correlation self-energy implementations

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::{Array1, Array3};
use num_complex::Complex64;
use quasix_core::selfenergy::correlation::{ContourDeformationConfig, CorrelationSelfEnergyCD};
use std::hint::black_box;

/// Generate test data for benchmarking
fn generate_test_data(n_mo: usize, n_aux: usize, n_freq: usize) -> TestData {
    // Generate MO energies
    let mut mo_energy = Array1::<f64>::zeros(n_mo);
    for i in 0..n_mo {
        mo_energy[i] = -10.0 + 20.0 * (i as f64) / (n_mo as f64);
    }

    // Generate occupations (half occupied)
    let mut mo_occ = Array1::<f64>::zeros(n_mo);
    for i in 0..n_mo / 2 {
        mo_occ[i] = 2.0;
    }

    // Generate screened interaction
    let mut w_screened = Array3::<Complex64>::zeros((n_freq, n_aux, n_aux));
    for i in 0..n_freq {
        for j in 0..n_aux {
            for k in 0..n_aux {
                let val = 1.0 / (1.0 + (j as f64 - k as f64).abs() + 0.1 * i as f64);
                w_screened[[i, j, k]] = Complex64::new(val, 0.0);
            }
        }
    }

    // Generate frequency grid
    let omega_grid = Array1::linspace(-30.0, 30.0, n_freq);

    // Generate evaluation points
    let eval_points = Array1::linspace(-20.0, 20.0, 20);

    // Generate DF tensor
    let mut df_tensor = Array3::<f64>::zeros((n_mo, n_mo, n_aux));
    for i in 0..n_mo {
        for j in 0..n_mo {
            for k in 0..n_aux {
                df_tensor[[i, j, k]] = 1.0 / (1.0 + (i as f64 - j as f64).abs() + 0.1 * k as f64);
            }
        }
    }

    TestData {
        mo_energy,
        mo_occ,
        w_screened,
        omega_grid,
        eval_points,
        df_tensor,
    }
}

struct TestData {
    mo_energy: Array1<f64>,
    mo_occ: Array1<f64>,
    w_screened: Array3<Complex64>,
    omega_grid: Array1<f64>,
    eval_points: Array1<f64>,
    df_tensor: Array3<f64>,
}

/// Benchmark small system
fn benchmark_small_system(c: &mut Criterion) {
    let data = generate_test_data(10, 50, 32);

    let mut group = c.benchmark_group("correlation_small");

    // Default configuration
    group.bench_function("default_config", |b| {
        let config = ContourDeformationConfig::default();
        let calc = CorrelationSelfEnergyCD::new(10, 50, config);

        b.iter(|| {
            calc.compute_sigma_c(
                black_box(&data.mo_energy),
                black_box(&data.mo_occ),
                black_box(&data.w_screened),
                black_box(&data.omega_grid),
                black_box(&data.eval_points),
                black_box(&data.df_tensor),
            )
        });
    });

    // SIMD enabled
    group.bench_function("simd_enabled", |b| {
        let config = ContourDeformationConfig {
            use_simd: true,
            ..Default::default()
        };
        let calc = CorrelationSelfEnergyCD::new(10, 50, config);

        b.iter(|| {
            calc.compute_sigma_c(
                black_box(&data.mo_energy),
                black_box(&data.mo_occ),
                black_box(&data.w_screened),
                black_box(&data.omega_grid),
                black_box(&data.eval_points),
                black_box(&data.df_tensor),
            )
        });
    });

    // SIMD disabled
    group.bench_function("simd_disabled", |b| {
        let config = ContourDeformationConfig {
            use_simd: false,
            ..Default::default()
        };
        let calc = CorrelationSelfEnergyCD::new(10, 50, config);

        b.iter(|| {
            calc.compute_sigma_c(
                black_box(&data.mo_energy),
                black_box(&data.mo_occ),
                black_box(&data.w_screened),
                black_box(&data.omega_grid),
                black_box(&data.eval_points),
                black_box(&data.df_tensor),
            )
        });
    });

    group.finish();
}

/// Benchmark medium system
fn benchmark_medium_system(c: &mut Criterion) {
    let data = generate_test_data(30, 100, 64);

    let mut group = c.benchmark_group("correlation_medium");

    // Test different thread counts
    for n_threads in [1, 2, 4, 8].iter() {
        group.bench_with_input(
            BenchmarkId::new("threads", n_threads),
            n_threads,
            |b, &n_threads| {
                let config = ContourDeformationConfig {
                    n_threads: Some(n_threads),
                    use_simd: true,
                    ..Default::default()
                };
                let calc = CorrelationSelfEnergyCD::new(30, 100, config);

                b.iter(|| {
                    calc.compute_sigma_c(
                        black_box(&data.mo_energy),
                        black_box(&data.mo_occ),
                        black_box(&data.w_screened),
                        black_box(&data.omega_grid),
                        black_box(&data.eval_points),
                        black_box(&data.df_tensor),
                    )
                });
            },
        );
    }

    group.finish();
}

/// Benchmark large system
fn benchmark_large_system(c: &mut Criterion) {
    let data = generate_test_data(50, 200, 128);

    let mut group = c.benchmark_group("correlation_large");
    group.sample_size(10); // Reduce sample size for large systems

    // Optimized configuration
    group.bench_function("optimized", |b| {
        let config = ContourDeformationConfig {
            use_simd: true,
            n_threads: Some(8),
            n_imag_points: 64,
            ..Default::default()
        };
        let calc = CorrelationSelfEnergyCD::new(50, 200, config);

        b.iter(|| {
            calc.compute_sigma_c(
                black_box(&data.mo_energy),
                black_box(&data.mo_occ),
                black_box(&data.w_screened),
                black_box(&data.omega_grid),
                black_box(&data.eval_points),
                black_box(&data.df_tensor),
            )
        });
    });

    // Baseline configuration
    group.bench_function("baseline", |b| {
        let config = ContourDeformationConfig {
            use_simd: false,
            n_threads: Some(1),
            ..Default::default()
        };
        let calc = CorrelationSelfEnergyCD::new(50, 200, config);

        b.iter(|| {
            calc.compute_sigma_c(
                black_box(&data.mo_energy),
                black_box(&data.mo_occ),
                black_box(&data.w_screened),
                black_box(&data.omega_grid),
                black_box(&data.eval_points),
                black_box(&data.df_tensor),
            )
        });
    });

    group.finish();
}

/// Benchmark integration accuracy vs performance
fn benchmark_integration_accuracy(c: &mut Criterion) {
    let data = generate_test_data(20, 80, 48);

    let mut group = c.benchmark_group("correlation_accuracy");

    // Test different numbers of integration points
    for n_points in [16, 32, 64, 128].iter() {
        group.bench_with_input(
            BenchmarkId::new("gl_points", n_points),
            n_points,
            |b, &n_points| {
                let config = ContourDeformationConfig {
                    n_imag_points: n_points,
                    use_simd: true,
                    ..Default::default()
                };
                let calc = CorrelationSelfEnergyCD::new(20, 80, config);

                b.iter(|| {
                    calc.compute_sigma_c(
                        black_box(&data.mo_energy),
                        black_box(&data.mo_occ),
                        black_box(&data.w_screened),
                        black_box(&data.omega_grid),
                        black_box(&data.eval_points),
                        black_box(&data.df_tensor),
                    )
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_small_system,
    benchmark_medium_system,
    benchmark_large_system,
    benchmark_integration_accuracy
);

criterion_main!(benches);
