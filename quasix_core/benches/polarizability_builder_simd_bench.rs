//! Benchmarks for SIMD-optimized polarizability builder
//!
//! This benchmark measures the performance improvements from SIMD vectorization
//! in the P0 denominator update implementation.

#![warn(clippy::all, clippy::pedantic, clippy::perf)]
#![warn(missing_docs)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::hint::black_box;
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::gw::polarizability_builder::PolarizabilityBuilder;
use rand::{rngs::SmallRng, SeedableRng};
use rand_distr::{Distribution, Normal};

/// Setup test system with given dimensions
fn setup_system(nocc: usize, nvirt: usize, naux: usize) -> (Array1<f64>, Array2<f64>) {
    // Create realistic MO energies
    let mut mo_energies = Array1::zeros(nocc + nvirt);
    for i in 0..nocc {
        mo_energies[i] = -2.0 + 0.15 * i as f64;
    }
    for a in 0..nvirt {
        mo_energies[nocc + a] = 0.1 + 0.08 * a as f64;
    }

    // Create random DF tensor
    let mut rng = SmallRng::seed_from_u64(42);
    let dist = Normal::new(0.0, 0.5).unwrap();
    let n_trans = nocc * nvirt;
    let mut df_ia = Array2::zeros((n_trans, naux));
    for i in 0..n_trans {
        for j in 0..naux {
            df_ia[[i, j]] = dist.sample(&mut rng);
        }
    }

    (mo_energies, df_ia)
}

/// Benchmark gap computation with SIMD
fn bench_gap_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("gap_computation");

    for &(nocc, nvirt) in &[(10, 50), (20, 100), (50, 200), (100, 400)] {
        let naux = nvirt * 2;
        let (mo_energies, df_ia) = setup_system(nocc, nvirt, naux);
        let builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

        let n_gaps = nocc * nvirt;
        group.throughput(Throughput::Elements(n_gaps as u64));

        group.bench_with_input(
            BenchmarkId::new("simd", format!("{nocc}x{nvirt}")),
            &builder,
            |b, builder| {
                b.iter(|| {
                    let gaps = builder.compute_gaps();
                    black_box(gaps)
                });
            },
        );

        // Benchmark scalar version for comparison
        group.bench_with_input(
            BenchmarkId::new("scalar", format!("{nocc}x{nvirt}")),
            &builder,
            |b, builder| {
                b.iter(|| {
                    let mut gaps = vec![0.0; n_gaps];
                    builder.compute_gaps_scalar(&mut gaps);
                    black_box(gaps)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark P0 build for imaginary frequencies
fn bench_p0_imaginary(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_imaginary");

    for &(nocc, nvirt) in &[(10, 50), (20, 100), (30, 150)] {
        let naux = nvirt * 2;
        let (mo_energies, df_ia) = setup_system(nocc, nvirt, naux);
        let builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

        let omega = Complex64::new(0.0, 1.0); // Imaginary frequency

        group.throughput(Throughput::Elements((naux * naux) as u64));

        #[cfg(target_arch = "x86_64")]
        group.bench_with_input(
            BenchmarkId::new("simd", format!("{nocc}x{nvirt}x{naux}")),
            &builder,
            |b, builder| {
                b.iter(|| {
                    let p0 = builder.build_p0_simd_imaginary(black_box(omega)).unwrap();
                    black_box(p0)
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("standard", format!("{nocc}x{nvirt}x{naux}")),
            &builder,
            |b, builder| {
                b.iter(|| {
                    let p0 = builder.build_p0_standard(black_box(omega)).unwrap();
                    black_box(p0)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark energy update with gap statistics
fn bench_energy_update(c: &mut Criterion) {
    let mut group = c.benchmark_group("energy_update");

    for &(nocc, nvirt) in &[(20, 100), (50, 200), (100, 400)] {
        let naux = nvirt;
        let (mo_energies, df_ia) = setup_system(nocc, nvirt, naux);
        let mut builder =
            PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

        // Create updated QP energies
        let mut qp_energies = mo_energies.clone();
        for i in 0..(nocc + nvirt) {
            qp_energies[i] *= 1.05; // 5% shift
        }

        group.throughput(Throughput::Elements((nocc * nvirt) as u64));

        group.bench_with_input(
            BenchmarkId::new("update", format!("{nocc}x{nvirt}")),
            &qp_energies,
            |b, qp_energies| {
                b.iter(|| {
                    let stats = builder.update_energies(black_box(qp_energies)).unwrap();
                    black_box(stats)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark batch P0 computation
fn bench_p0_batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_batch");

    let (nocc, nvirt, naux) = (20, 80, 160);
    let (mo_energies, df_ia) = setup_system(nocc, nvirt, naux);
    let builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

    // Create frequency grids of different sizes
    for n_freq in [10, 20, 50] {
        let frequencies: Vec<Complex64> = (0..n_freq)
            .map(|i| Complex64::new(0.0, 0.1 * i as f64))
            .collect();

        group.throughput(Throughput::Elements(n_freq as u64));

        group.bench_with_input(
            BenchmarkId::new("batch", n_freq),
            &frequencies,
            |b, frequencies| {
                b.iter(|| {
                    let p0_batch = builder.build_p0_batch(black_box(frequencies)).unwrap();
                    black_box(p0_batch)
                });
            },
        );
    }

    group.finish();
}

/// Main benchmark runner with optimized configuration
fn configure_criterion() -> Criterion {
    Criterion::default()
        .sample_size(100)
        .measurement_time(std::time::Duration::from_secs(10))
        .warm_up_time(std::time::Duration::from_secs(2))
}

criterion_group! {
    name = benches;
    config = configure_criterion();
    targets = bench_gap_computation, bench_p0_imaginary, bench_energy_update, bench_p0_batch
}

criterion_main!(benches);