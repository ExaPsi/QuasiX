//! Benchmarks for 3-center integral computation
//!
//! This benchmark suite tests the performance of the optimized
//! 3-center integral implementation against various molecular systems.

#![warn(clippy::all)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use quasix_core::integrals::{
    compute_3center_integrals_optimized, BasisSet, Molecule, OptimizedThreeCenterIntegrals,
};
use std::hint::black_box;

/// Benchmark water molecule (small system)
fn bench_water(c: &mut Criterion) {
    let mol = Molecule::water();
    let basis = BasisSet::for_molecule("H2O", false);
    let aux_basis = BasisSet::for_molecule("H2O", true);

    let mut group = c.benchmark_group("water_3center");

    group.bench_function("optimized", |b| {
        b.iter(|| {
            compute_3center_integrals_optimized(
                black_box(&mol),
                black_box(&basis),
                black_box(&aux_basis),
            )
        });
    });

    group.finish();
}

/// Benchmark ammonia molecule (medium system)
fn bench_ammonia(c: &mut Criterion) {
    let mol = Molecule::ammonia();
    let basis = BasisSet::for_molecule("NH3", false);
    let aux_basis = BasisSet::for_molecule("NH3", true);

    let mut group = c.benchmark_group("ammonia_3center");

    group.bench_function("optimized", |b| {
        b.iter(|| {
            compute_3center_integrals_optimized(
                black_box(&mol),
                black_box(&basis),
                black_box(&aux_basis),
            )
        });
    });

    group.finish();
}

/// Benchmark benzene molecule (target system)
fn bench_benzene(c: &mut Criterion) {
    let mol = Molecule::benzene();
    let basis = BasisSet::for_molecule("C6H6", false);
    let aux_basis = BasisSet::for_molecule("C6H6", true);

    let mut group = c.benchmark_group("benzene_3center");
    group.measurement_time(std::time::Duration::from_secs(10));

    group.bench_function("optimized", |b| {
        b.iter(|| {
            compute_3center_integrals_optimized(
                black_box(&mol),
                black_box(&basis),
                black_box(&aux_basis),
            )
        });
    });

    // Also benchmark with performance metrics
    group.bench_function("optimized_with_metrics", |b| {
        b.iter(|| {
            let evaluator = OptimizedThreeCenterIntegrals::new(36, 72);
            evaluator.compute_optimized(black_box(&mol), black_box(&basis), black_box(&aux_basis))
        });
    });

    group.finish();
}

/// Benchmark scaling with system size
fn bench_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("3center_scaling");

    for nbasis in [8, 16, 32, 64].iter() {
        let naux = nbasis * 2;

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}x{}", nbasis, naux)),
            &(*nbasis, naux),
            |b, &(nb, na)| {
                // Create mock molecule and basis sets
                let mol = Molecule::water(); // Mock, size doesn't matter
                let basis = BasisSet::mock_sto3g(nb);
                let aux_basis = BasisSet::mock_ri(na);

                b.iter(|| {
                    compute_3center_integrals_optimized(
                        black_box(&mol),
                        black_box(&basis),
                        black_box(&aux_basis),
                    )
                });
            },
        );
    }

    group.finish();
}

/// Benchmark parallel efficiency
fn bench_parallel_efficiency(c: &mut Criterion) {
    let mol = Molecule::benzene();
    let basis = BasisSet::for_molecule("C6H6", false);
    let aux_basis = BasisSet::for_molecule("C6H6", true);

    let mut group = c.benchmark_group("parallel_efficiency");

    // Test with different thread counts
    for n_threads in [1, 2, 4, 8].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}_threads", n_threads)),
            n_threads,
            |b, &threads| {
                // Set thread pool size
                let _pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();

                b.iter(|| {
                    compute_3center_integrals_optimized(
                        black_box(&mol),
                        black_box(&basis),
                        black_box(&aux_basis),
                    )
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_water,
    bench_ammonia,
    bench_benzene,
    bench_scaling,
    bench_parallel_efficiency
);
criterion_main!(benches);
