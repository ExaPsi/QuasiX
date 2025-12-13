//! Benchmarks for optimized energy-dependent self-energy evaluation
#![warn(clippy::all, clippy::pedantic, clippy::perf)]

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::selfenergy::energy_dependent_optimized::{
    EnergyDependentSelfEnergyOptimized, SigmaCConfig, compute_sigma_c_parallel,
};

/// Generate test data for benchmarking
fn generate_test_data(n_orb: usize, n_freq: usize) -> (Array1<f64>, Array2<Complex64>, Array1<f64>) {
    let sigma_x = Array1::from_vec((0..n_orb).map(|i| -10.0 - i as f64 * 0.5).collect());

    let mut sigma_c_grid = Array2::zeros((n_orb, n_freq));
    for i in 0..n_orb {
        for j in 0..n_freq {
            let omega = -2.0 + 4.0 * j as f64 / (n_freq - 1) as f64;
            sigma_c_grid[[i, j]] = Complex64::new(
                -0.3 - 0.05 * i as f64 + 0.02 * omega,
                -0.01 * (1.0 + omega.abs() / 10.0),
            );
        }
    }

    let omega_grid = Array1::linspace(-2.0, 2.0, n_freq);

    (sigma_x, sigma_c_grid, omega_grid)
}

/// Benchmark single frequency evaluation
fn bench_single_frequency(c: &mut Criterion) {
    let mut group = c.benchmark_group("energy_dependent_single_freq");

    for n_orb in [4, 8, 16, 32].iter() {
        let (sigma_x, sigma_c_grid, omega_grid) = generate_test_data(*n_orb, 101);

        let evaluator = EnergyDependentSelfEnergyOptimized::new(
            sigma_x,
            sigma_c_grid,
            omega_grid,
        ).unwrap();

        group.throughput(Throughput::Elements(*n_orb as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(n_orb),
            n_orb,
            |b, _| {
                b.iter(|| {
                    let omega = black_box(0.5);
                    let result = evaluator.evaluate_all(omega);
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark batch frequency evaluation
fn bench_batch_frequency(c: &mut Criterion) {
    let mut group = c.benchmark_group("energy_dependent_batch_freq");

    let n_orb = 16;
    let (sigma_x, sigma_c_grid, omega_grid) = generate_test_data(n_orb, 201);

    let evaluator = EnergyDependentSelfEnergyOptimized::new(
        sigma_x,
        sigma_c_grid,
        omega_grid.clone(),
    ).unwrap();

    for n_freq_eval in [10, 50, 100, 200].iter() {
        let test_omegas: Vec<f64> = (0..*n_freq_eval)
            .map(|i| -2.0 + 4.0 * i as f64 / (*n_freq_eval - 1) as f64)
            .collect();

        group.throughput(Throughput::Elements((n_orb * n_freq_eval) as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(n_freq_eval),
            n_freq_eval,
            |b, _| {
                b.iter(|| {
                    for &omega in test_omegas.iter() {
                        let result = evaluator.evaluate_all(black_box(omega));
                        black_box(result);
                    }
                });
            },
        );
    }

    group.finish();
}

/// Benchmark interpolation methods
fn bench_interpolation_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("interpolation_methods");

    let n_orb = 16;
    let n_freq = 101;
    let (sigma_x, sigma_c_grid, omega_grid) = generate_test_data(n_orb, n_freq);

    let evaluator = EnergyDependentSelfEnergyOptimized::new(
        sigma_x,
        sigma_c_grid,
        omega_grid,
    ).unwrap();

    // Test interpolation at different points
    let test_points = vec![
        ("grid_point", 0.0),        // Exactly on grid
        ("between_points", 0.123),  // Between grid points
        ("extrapolate_low", -2.5),  // Below grid
        ("extrapolate_high", 2.5),  // Above grid
    ];

    for (name, omega) in test_points {
        group.bench_function(name, |b| {
            b.iter(|| {
                let result = evaluator.evaluate_all(black_box(omega));
                black_box(result)
            });
        });
    }

    group.finish();
}

/// Benchmark derivative computation
fn bench_derivative_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("energy_dependent_derivatives");

    for n_orb in [4, 8, 16, 32].iter() {
        let (sigma_x, sigma_c_grid, omega_grid) = generate_test_data(*n_orb, 101);

        let evaluator = EnergyDependentSelfEnergyOptimized::new(
            sigma_x,
            sigma_c_grid,
            omega_grid,
        ).unwrap();

        group.throughput(Throughput::Elements(*n_orb as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(n_orb),
            n_orb,
            |b, _| {
                b.iter(|| {
                    let omega = black_box(0.5);
                    let result = evaluator.compute_derivatives_all(omega);
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark cache effectiveness
fn bench_cache_effectiveness(c: &mut Criterion) {
    let mut group = c.benchmark_group("cache_effectiveness");

    let n_orb = 16;
    let (sigma_x, sigma_c_grid, omega_grid) = generate_test_data(n_orb, 101);

    let evaluator = EnergyDependentSelfEnergyOptimized::new(
        sigma_x,
        sigma_c_grid,
        omega_grid,
    ).unwrap();

    // Benchmark cache hit
    group.bench_function("cache_hit", |b| {
        let omega = 0.5;
        // Warm up cache
        let _ = evaluator.evaluate_all(omega);

        b.iter(|| {
            let result = evaluator.evaluate_all(black_box(omega));
            black_box(result)
        });
    });

    // Benchmark cache miss
    group.bench_function("cache_miss", |b| {
        let mut omega = 0.0;
        b.iter(|| {
            omega += 0.001; // Different value each time
            evaluator.clear_cache();
            let result = evaluator.evaluate_all(black_box(omega));
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark parallel Sigma_c computation
fn bench_parallel_sigma_c(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_sigma_c");

    for n_mo in [10, 20, 40].iter() {
        let n_occ = n_mo / 2;
        let n_omega = 32;
        let n_aux = 50;

        let mo_energy = Array1::linspace(-20.0, 10.0, *n_mo);
        let mut mo_occ = Array1::zeros(*n_mo);
        mo_occ.slice_mut(ndarray::s![..n_occ]).fill(2.0);

        let df_tensors = ndarray::Array3::ones((*n_mo, *n_mo, n_aux));
        let w_screened = ndarray::Array3::from_elem(
            (n_aux, n_aux, n_omega),
            Complex64::new(0.1, -0.001),
        );
        let omega_grid = Array1::linspace(-5.0, 5.0, n_omega);

        let config = SigmaCConfig::default();

        group.throughput(Throughput::Elements((n_mo * n_omega) as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(n_mo),
            n_mo,
            |b, _| {
                b.iter(|| {
                    let result = compute_sigma_c_parallel(
                        &mo_energy,
                        &mo_occ,
                        &df_tensors,
                        &w_screened,
                        &omega_grid,
                        &config,
                    ).unwrap();
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark SIMD vs scalar performance
fn bench_simd_vs_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_vs_scalar");

    // Test with different orbital counts to show SIMD benefit
    for n_orb in [1, 4, 8, 16, 32, 64].iter() {
        let (sigma_x, sigma_c_grid, omega_grid) = generate_test_data(*n_orb, 101);

        let evaluator = EnergyDependentSelfEnergyOptimized::new(
            sigma_x,
            sigma_c_grid,
            omega_grid,
        ).unwrap();

        group.throughput(Throughput::Elements(*n_orb as u64));
        group.bench_with_input(
            BenchmarkId::new("evaluate_all", n_orb),
            n_orb,
            |b, _| {
                b.iter(|| {
                    let omega = black_box(0.5);
                    let result = evaluator.evaluate_all(omega);
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark memory usage patterns
fn bench_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_patterns");

    // Different grid sizes to test memory scaling
    for n_freq in [51, 101, 201, 401].iter() {
        let n_orb = 16;
        let (sigma_x, sigma_c_grid, omega_grid) = generate_test_data(n_orb, *n_freq);

        let evaluator = EnergyDependentSelfEnergyOptimized::new(
            sigma_x,
            sigma_c_grid,
            omega_grid,
        ).unwrap();

        group.throughput(Throughput::Bytes((n_orb * 16) as u64)); // Complex64 = 16 bytes
        group.bench_with_input(
            BenchmarkId::new("grid_size", n_freq),
            n_freq,
            |b, _| {
                b.iter(|| {
                    // Evaluate at multiple points to test memory access patterns
                    for i in 0..10 {
                        let omega = -1.0 + 0.2 * i as f64;
                        let result = evaluator.evaluate_all(black_box(omega));
                        black_box(result);
                    }
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_single_frequency,
    bench_batch_frequency,
    bench_interpolation_comparison,
    bench_derivative_computation,
    bench_cache_effectiveness,
    bench_parallel_sigma_c,
    bench_simd_vs_scalar,
    bench_memory_patterns
);

criterion_main!(benches);