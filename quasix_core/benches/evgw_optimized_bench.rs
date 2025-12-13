//! Benchmark for optimized evGW implementation
//!
//! Demonstrates performance improvements from SIMD, parallelization,
//! and memory optimizations.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::{Array1, Array2, Array3};
use quasix_core::freq::{FrequencyGrid, GridType};
use quasix_core::gw::evgw::{EvGWConfig, EvGWDriver};
use std::hint::black_box;

type TestData = (
    Array1<f64>,   // mo_energies
    Array1<f64>,   // mo_occ
    Array3<f64>,   // ia_p
    Array3<f64>,   // ij_p
    Array2<f64>,   // chol_v
    Array1<f64>,   // vxc_dft
    FrequencyGrid, // freq_grid
);

/// Generate test data for benchmarking
fn generate_test_data(nbasis: usize, nocc: usize, naux: usize) -> TestData {
    // Generate realistic molecular orbital energies
    let mo_energies = Array1::linspace(-1.0, 1.0, nbasis);

    // Occupation numbers (2 electrons per occupied orbital)
    let mut mo_occ = Array1::zeros(nbasis);
    for i in 0..nocc {
        mo_occ[i] = 2.0;
    }

    // DF tensors with realistic structure
    let nvirt = nbasis - nocc;
    let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| {
        let r = ((i + 1) * (a + 1)) as f64;
        let s = (p + 1) as f64;
        (r * s * 0.01).exp() * 0.1
    });

    let ij_p = Array3::from_shape_fn((nocc, nocc, naux), |(i, j, p)| {
        let r = ((i + 1) * (j + 1)) as f64;
        let s = (p + 1) as f64;
        (r * s * 0.01).exp() * 0.1
    });

    // Cholesky decomposed Coulomb metric
    let mut chol_v = Array2::eye(naux);
    for i in 0..naux {
        for j in i + 1..naux {
            let val = 1.0 / ((i as f64 - j as f64).abs() + 1.0);
            chol_v[[i, j]] = val;
            chol_v[[j, i]] = val;
        }
    }

    // DFT XC potential
    let vxc_dft = Array1::from_shape_fn(nbasis, |i| -0.5 * (1.0 + 0.1 * i as f64));

    // Frequency grid for integration
    let freq_grid = FrequencyGrid::new(16, GridType::GaussLegendre).unwrap();

    (mo_energies, mo_occ, ia_p, ij_p, chol_v, vxc_dft, freq_grid)
}

/// Benchmark polarizability update with different optimizations
fn bench_polarizability_update(c: &mut Criterion) {
    let mut group = c.benchmark_group("polarizability_update");

    for nbasis in [20, 40, 60] {
        let nocc = nbasis / 2;
        let naux = nbasis * 2;

        let (mo_energies, mo_occ, ia_p, ij_p, chol_v, vxc_dft, freq_grid) =
            generate_test_data(nbasis, nocc, naux);

        // Test without SIMD
        group.bench_with_input(BenchmarkId::new("standard", nbasis), &nbasis, |b, _| {
            let config = EvGWConfig {
                n_threads: 1,
                parallel_freq: false,
                max_cycle: 1,
                ..Default::default()
            };

            let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

            b.iter(|| {
                driver
                    .run_evgw_loop(
                        &mo_energies,
                        &mo_occ,
                        &ia_p,
                        &ij_p,
                        &chol_v,
                        &vxc_dft,
                        &freq_grid,
                    )
                    .unwrap()
            });
        });

        // Test with parallel optimization
        group.bench_with_input(BenchmarkId::new("parallel", nbasis), &nbasis, |b, _| {
            let config = EvGWConfig {
                n_threads: 4,
                parallel_freq: true,
                max_cycle: 1,
                ..Default::default()
            };

            let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

            b.iter(|| {
                driver
                    .run_evgw_loop(
                        &mo_energies,
                        &mo_occ,
                        &ia_p,
                        &ij_p,
                        &chol_v,
                        &vxc_dft,
                        &freq_grid,
                    )
                    .unwrap()
            });
        });

        // Test with full optimizations
        group.bench_with_input(BenchmarkId::new("optimized", nbasis), &nbasis, |b, _| {
            let config = EvGWConfig {
                n_threads: num_cpus::get(),
                parallel_freq: true,
                cache_align: true,
                block_size: 64,
                max_cycle: 1,
                ..Default::default()
            };

            let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

            b.iter(|| {
                driver
                    .run_evgw_loop(
                        &mo_energies,
                        &mo_occ,
                        &ia_p,
                        &ij_p,
                        &chol_v,
                        &vxc_dft,
                        &freq_grid,
                    )
                    .unwrap()
            });
        });
    }

    group.finish();
}

/// Benchmark frequency integration with different thread counts
fn bench_frequency_integration(c: &mut Criterion) {
    let mut group = c.benchmark_group("frequency_integration");

    let nbasis = 40;
    let nocc = 20;
    let naux = 80;

    let (mo_energies, mo_occ, ia_p, ij_p, chol_v, vxc_dft, freq_grid) =
        generate_test_data(nbasis, nocc, naux);

    for n_threads in [1, 2, 4, 8] {
        if n_threads > num_cpus::get() {
            continue;
        }

        group.bench_with_input(
            BenchmarkId::new("threads", n_threads),
            &n_threads,
            |b, &threads| {
                let config = EvGWConfig {
                    n_threads: threads,
                    parallel_freq: true,
                    nfreq: 32,
                    max_cycle: 1,
                    ..Default::default()
                };

                let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

                b.iter(|| {
                    driver
                        .run_evgw_loop(
                            &mo_energies,
                            &mo_occ,
                            &ia_p,
                            &ij_p,
                            &chol_v,
                            &vxc_dft,
                            &freq_grid,
                        )
                        .unwrap()
                });
            },
        );
    }

    group.finish();
}

/// Benchmark cache blocking effectiveness
fn bench_cache_blocking(c: &mut Criterion) {
    let mut group = c.benchmark_group("cache_blocking");

    let nbasis = 50;
    let nocc = 25;
    let naux = 100;

    let (mo_energies, mo_occ, ia_p, ij_p, chol_v, vxc_dft, freq_grid) =
        generate_test_data(nbasis, nocc, naux);

    for block_size in [16, 32, 64, 128] {
        group.bench_with_input(
            BenchmarkId::new("block_size", block_size),
            &block_size,
            |b, &bs| {
                let config = EvGWConfig {
                    block_size: bs,
                    n_threads: 4,
                    max_cycle: 1,
                    ..Default::default()
                };

                let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

                b.iter(|| {
                    driver
                        .run_evgw_loop(
                            &mo_energies,
                            &mo_occ,
                            &ia_p,
                            &ij_p,
                            &chol_v,
                            &vxc_dft,
                            &freq_grid,
                        )
                        .unwrap()
                });
            },
        );
    }

    group.finish();
}

/// Benchmark full evGW iteration with realistic molecule sizes
fn bench_evgw_molecules(c: &mut Criterion) {
    let mut group = c.benchmark_group("evgw_molecules");
    group.sample_size(10); // Reduce sample size for longer benchmarks

    // Molecule sizes: H2O, benzene, naphthalene-like
    let molecules = [
        ("H2O", 10, 5, 20),
        ("C6H6", 30, 15, 60),
        ("C10H8", 48, 24, 96),
    ];

    for (name, nbasis, nocc, naux) in molecules {
        let (mo_energies, mo_occ, ia_p, ij_p, chol_v, vxc_dft, freq_grid) =
            generate_test_data(nbasis, nocc, naux);

        group.bench_with_input(BenchmarkId::new("optimized", name), &name, |b, _| {
            let config = EvGWConfig {
                n_threads: num_cpus::get(),
                parallel_freq: true,
                cache_align: true,
                max_cycle: 3, // Run actual iterations
                conv_tol: 1e-3,
                ..Default::default()
            };

            let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

            b.iter(|| {
                driver
                    .run_evgw_loop(
                        black_box(&mo_energies),
                        black_box(&mo_occ),
                        black_box(&ia_p),
                        black_box(&ij_p),
                        black_box(&chol_v),
                        black_box(&vxc_dft),
                        black_box(&freq_grid),
                    )
                    .unwrap()
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_polarizability_update,
    bench_frequency_integration,
    bench_cache_blocking,
    bench_evgw_molecules,
);
criterion_main!(benches);
