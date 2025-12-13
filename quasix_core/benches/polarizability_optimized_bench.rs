//! Benchmarks for optimized polarizability computation
#![warn(clippy::all, clippy::pedantic, clippy::perf)]

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;
use quasix_core::dielectric::polarizability_optimized::{PolarizabilityOptimized, RPAOptimized};

/// Generate test data for polarizability benchmarks
fn generate_pol_test_data(
    n_occ: usize,
    n_virt: usize,
    n_aux: usize,
) -> (Array2<f64>, Array1<f64>, Array1<f64>) {
    let n_trans = n_occ * n_virt;

    // Generate DF tensors with realistic structure
    let df_ia = Array2::from_shape_fn((n_trans, n_aux), |(t, p)| {
        let i = t / n_virt;
        let a = t % n_virt;
        ((i + 1) as f64 * (a + 1) as f64 * (p + 1) as f64).sqrt() / 10.0
    });

    let e_occ = Array1::linspace(-20.0, -5.0, n_occ);
    let e_virt = Array1::linspace(1.0, 10.0, n_virt);

    (df_ia, e_occ, e_virt)
}

/// Benchmark single frequency P0 computation
fn bench_p0_single_frequency(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_single_frequency");

    // Test different system sizes
    let configs = vec![
        (2, 4, 20),   // Small: 2 occ, 4 virt, 20 aux
        (5, 10, 50),  // Medium: 5 occ, 10 virt, 50 aux
        (10, 20, 100), // Large: 10 occ, 20 virt, 100 aux
    ];

    for (n_occ, n_virt, n_aux) in configs {
        let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, n_aux);

        let pol_calc = PolarizabilityOptimized::new(
            n_occ,
            n_virt,
            n_aux,
            &df_ia,
            &e_occ,
            &e_virt,
            1e-4,
        ).unwrap();

        let size_label = format!("{}x{}x{}", n_occ, n_virt, n_aux);
        group.throughput(Throughput::Elements((n_aux * n_aux) as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(&size_label),
            &size_label,
            |b, _| {
                let omega = Complex64::new(0.5, 0.01);
                b.iter(|| {
                    let result = pol_calc.compute_p0(black_box(omega)).unwrap();
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark different frequency types
fn bench_p0_frequency_types(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_frequency_types");

    let n_occ = 5;
    let n_virt = 10;
    let n_aux = 50;
    let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, n_aux);

    let pol_calc = PolarizabilityOptimized::new(
        n_occ,
        n_virt,
        n_aux,
        &df_ia,
        &e_occ,
        &e_virt,
        1e-4,
    ).unwrap();

    let frequencies = vec![
        ("static", Complex64::new(0.0, 0.0)),
        ("imaginary", Complex64::new(0.0, 1.0)),
        ("real", Complex64::new(1.0, 0.0)),
        ("complex", Complex64::new(0.5, 0.1)),
    ];

    for (name, omega) in frequencies {
        group.bench_function(name, |b| {
            b.iter(|| {
                let result = pol_calc.compute_p0(black_box(omega)).unwrap();
                black_box(result)
            });
        });
    }

    group.finish();
}

/// Benchmark batch frequency computation
fn bench_p0_batch_frequency(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_batch_frequency");

    let n_occ = 5;
    let n_virt = 10;
    let n_aux = 50;
    let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, n_aux);

    let pol_calc = PolarizabilityOptimized::new(
        n_occ,
        n_virt,
        n_aux,
        &df_ia,
        &e_occ,
        &e_virt,
        1e-4,
    ).unwrap();

    for n_freq in [8, 16, 32, 64].iter() {
        let omega_grid = Array1::from_vec(
            (0..*n_freq)
                .map(|i| Complex64::new(i as f64 * 0.1, 0.01))
                .collect(),
        );

        group.throughput(Throughput::Elements((n_freq * n_aux * n_aux) as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(n_freq),
            n_freq,
            |b, _| {
                b.iter(|| {
                    let result = pol_calc.compute_p0_batch(&omega_grid).unwrap();
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark cache effectiveness
fn bench_p0_cache(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_cache");

    let n_occ = 5;
    let n_virt = 10;
    let n_aux = 50;
    let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, n_aux);

    let pol_calc = PolarizabilityOptimized::new(
        n_occ,
        n_virt,
        n_aux,
        &df_ia,
        &e_occ,
        &e_virt,
        1e-4,
    ).unwrap();

    let omega = Complex64::new(0.5, 0.01);

    // Warm up cache
    let _ = pol_calc.compute_p0(omega).unwrap();

    group.bench_function("cache_hit", |b| {
        b.iter(|| {
            let result = pol_calc.compute_p0(black_box(omega)).unwrap();
            black_box(result)
        });
    });

    group.bench_function("cache_miss", |b| {
        let mut omega_val = 0.0;
        b.iter(|| {
            omega_val += 0.001;
            let omega = Complex64::new(omega_val, 0.01);
            pol_calc.clear_cache();
            let result = pol_calc.compute_p0(black_box(omega)).unwrap();
            black_box(result)
        });
    });

    group.finish();
}

/// Benchmark RPA screened interaction computation
fn bench_rpa_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("rpa_screened_w");

    let configs = vec![
        (3, 6, 30),    // Small
        (5, 10, 50),   // Medium
        (8, 16, 80),   // Large
    ];

    for (n_occ, n_virt, n_aux) in configs {
        let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, n_aux);

        let pol_calc = PolarizabilityOptimized::new(
            n_occ,
            n_virt,
            n_aux,
            &df_ia,
            &e_occ,
            &e_virt,
            1e-4,
        ).unwrap();

        // Create simple Coulomb metric
        let mut v_sqrt = Array2::eye(n_aux);
        for i in 0..n_aux {
            v_sqrt[[i, i]] = ((i + 1) as f64 / n_aux as f64).sqrt();
        }
        let v_sqrt_inv = v_sqrt.mapv(|x| if x != 0.0 { 1.0 / x } else { 0.0 });

        let rpa = RPAOptimized::new(pol_calc, v_sqrt, v_sqrt_inv);

        let size_label = format!("{}x{}x{}", n_occ, n_virt, n_aux);
        group.throughput(Throughput::Elements((n_aux * n_aux) as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(&size_label),
            &size_label,
            |b, _| {
                let omega = Complex64::new(0.5, 0.01);
                b.iter(|| {
                    let result = rpa.compute_screened_interaction(black_box(omega)).unwrap();
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark transition product pre-computation
fn bench_transition_products(c: &mut Criterion) {
    let mut group = c.benchmark_group("transition_products");

    for n_trans in [10, 50, 100, 200].iter() {
        for n_aux in [20, 50, 100].iter() {
            let df_ia = Array2::from_shape_fn((*n_trans, *n_aux), |(t, p)| {
                ((t + 1) as f64 * (p + 1) as f64).sqrt() / 10.0
            });

            let label = format!("{}x{}", n_trans, n_aux);
            group.throughput(Throughput::Elements((n_trans * n_aux) as u64));
            group.bench_with_input(
                BenchmarkId::from_parameter(&label),
                &label,
                |b, _| {
                    b.iter(|| {
                        // This tests the pre-computation in the constructor
                        let pol_calc = PolarizabilityOptimized::new(
                            5,  // n_occ
                            *n_trans / 5,  // n_virt
                            *n_aux,
                            &df_ia,
                            &Array1::ones(5),
                            &Array1::ones(*n_trans / 5),
                            1e-4,
                        ).unwrap();
                        black_box(pol_calc);
                    });
                },
            );
        }
    }

    group.finish();
}

/// Benchmark parallel vs sequential batch computation
fn bench_parallel_vs_sequential(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_vs_sequential");

    let n_occ = 5;
    let n_virt = 10;
    let n_aux = 50;
    let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, n_aux);

    let pol_calc = PolarizabilityOptimized::new(
        n_occ,
        n_virt,
        n_aux,
        &df_ia,
        &e_occ,
        &e_virt,
        1e-4,
    ).unwrap();

    for n_freq in [4, 8, 16, 32].iter() {
        let omega_grid = Array1::from_vec(
            (0..*n_freq)
                .map(|i| Complex64::new(i as f64 * 0.1, 0.01))
                .collect(),
        );

        // Parallel batch computation
        group.bench_with_input(
            BenchmarkId::new("parallel", n_freq),
            n_freq,
            |b, _| {
                b.iter(|| {
                    let result = pol_calc.compute_p0_batch(&omega_grid).unwrap();
                    black_box(result)
                });
            },
        );

        // Sequential computation for comparison
        group.bench_with_input(
            BenchmarkId::new("sequential", n_freq),
            n_freq,
            |b, _| {
                b.iter(|| {
                    let mut results = Vec::with_capacity(*n_freq);
                    for &omega in omega_grid.iter() {
                        let result = pol_calc.compute_p0(omega).unwrap();
                        results.push(result);
                    }
                    black_box(results)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark memory usage scaling
fn bench_memory_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_scaling");

    // Test different auxiliary basis sizes
    for n_aux in [20, 40, 80, 160].iter() {
        let n_occ = 5;
        let n_virt = 10;
        let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, *n_aux);

        let pol_calc = PolarizabilityOptimized::new(
            n_occ,
            n_virt,
            *n_aux,
            &df_ia,
            &e_occ,
            &e_virt,
            1e-4,
        ).unwrap();

        group.throughput(Throughput::Bytes((n_aux * n_aux * 16) as u64)); // Complex64 = 16 bytes
        group.bench_with_input(
            BenchmarkId::from_parameter(n_aux),
            n_aux,
            |b, _| {
                let omega = Complex64::new(0.5, 0.01);
                b.iter(|| {
                    let result = pol_calc.compute_p0(black_box(omega)).unwrap();
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark H2O/cc-pVDZ realistic case
fn bench_h2o_ccpvdz(c: &mut Criterion) {
    let mut group = c.benchmark_group("h2o_ccpvdz");

    // H2O/cc-pVDZ: 5 occupied, 19 virtual, ~40 auxiliary
    let n_occ = 5;
    let n_virt = 19;
    let n_aux = 41;

    let (df_ia, e_occ, e_virt) = generate_pol_test_data(n_occ, n_virt, n_aux);

    let pol_calc = PolarizabilityOptimized::new(
        n_occ,
        n_virt,
        n_aux,
        &df_ia,
        &e_occ,
        &e_virt,
        1e-4,
    ).unwrap();

    // Single frequency computation
    group.bench_function("single_frequency", |b| {
        let omega = Complex64::new(0.5, 0.01);
        b.iter(|| {
            let result = pol_calc.compute_p0(black_box(omega)).unwrap();
            black_box(result)
        });
    });

    // Full frequency grid (typical evGW calculation)
    let n_freq = 32;
    let omega_grid = Array1::from_vec(
        (0..n_freq)
            .map(|i| Complex64::new(0.0, i as f64 * 2.0))
            .collect(),
    );

    group.bench_function("full_grid_32_frequencies", |b| {
        b.iter(|| {
            let result = pol_calc.compute_p0_batch(&omega_grid).unwrap();
            black_box(result)
        });
    });

    // Report timing per frequency
    group.bench_function("per_frequency_average", |b| {
        b.iter(|| {
            for &omega in omega_grid.iter().take(10) {
                let result = pol_calc.compute_p0(black_box(omega)).unwrap();
                black_box(result);
            }
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_p0_single_frequency,
    bench_p0_frequency_types,
    bench_p0_batch_frequency,
    bench_p0_cache,
    bench_rpa_computation,
    bench_transition_products,
    bench_parallel_vs_sequential,
    bench_memory_scaling,
    bench_h2o_ccpvdz
);

criterion_main!(benches);