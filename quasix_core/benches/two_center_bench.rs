//! Benchmarks for 2-center integral computation
//!
//! Run with: cargo bench --bench two_center_bench

#![warn(clippy::all)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use quasix_core::integrals::{
    compute_2center_integrals,
    two_center::{compute_2center_integrals_with_config, TwoCenterConfig},
    BasisSet,
};
use std::hint::black_box;
use std::time::Duration;

/// Benchmark different matrix sizes
fn benchmark_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("two_center_sizes");
    group.measurement_time(Duration::from_secs(10));

    // Test different sizes matching the performance targets
    for size in &[50, 100, 200, 300, 400, 500] {
        let aux_basis = BasisSet::mock_ri(*size);

        group.throughput(Throughput::Elements((*size * *size) as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| {
                let result =
                    compute_2center_integrals(&aux_basis).expect("Failed to compute integrals");
                black_box(result)
            });
        });
    }

    group.finish();
}

/// Benchmark serial vs parallel implementations
fn benchmark_parallel(c: &mut Criterion) {
    let mut group = c.benchmark_group("two_center_parallel");
    group.measurement_time(Duration::from_secs(10));

    let sizes = vec![100, 300, 500];

    for size in &sizes {
        let aux_basis = BasisSet::mock_ri(*size);

        // Serial configuration
        let serial_config = TwoCenterConfig {
            parallel: false,
            n_threads: 1,
            cache_blocking: false,
            ..Default::default()
        };

        // Parallel configuration
        let parallel_config = TwoCenterConfig {
            parallel: true,
            n_threads: 0, // Auto-detect
            cache_blocking: false,
            ..Default::default()
        };

        // Cache-blocked configuration
        let blocked_config = TwoCenterConfig {
            parallel: true,
            n_threads: 0,
            cache_blocking: true,
            block_size: 64,
            ..Default::default()
        };

        group.bench_with_input(BenchmarkId::new("serial", size), size, |b, _| {
            b.iter(|| {
                let result = compute_2center_integrals_with_config(&aux_basis, &serial_config)
                    .expect("Failed to compute integrals");
                black_box(result)
            });
        });

        group.bench_with_input(BenchmarkId::new("parallel", size), size, |b, _| {
            b.iter(|| {
                let result = compute_2center_integrals_with_config(&aux_basis, &parallel_config)
                    .expect("Failed to compute integrals");
                black_box(result)
            });
        });

        group.bench_with_input(BenchmarkId::new("blocked", size), size, |b, _| {
            b.iter(|| {
                let result = compute_2center_integrals_with_config(&aux_basis, &blocked_config)
                    .expect("Failed to compute integrals");
                black_box(result)
            });
        });
    }

    group.finish();
}

/// Benchmark cache blocking with different block sizes
fn benchmark_block_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("two_center_blocks");
    group.measurement_time(Duration::from_secs(5));

    let aux_basis = BasisSet::mock_ri(500);

    for block_size in &[16, 32, 64, 128, 256] {
        let config = TwoCenterConfig {
            parallel: true,
            n_threads: 0,
            cache_blocking: true,
            block_size: *block_size,
            ..Default::default()
        };

        group.bench_with_input(
            BenchmarkId::from_parameter(block_size),
            block_size,
            |b, _| {
                b.iter(|| {
                    let result = compute_2center_integrals_with_config(&aux_basis, &config)
                        .expect("Failed to compute integrals");
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark validation overhead
fn benchmark_validation(c: &mut Criterion) {
    let mut group = c.benchmark_group("two_center_validation");

    for size in &[100, 500] {
        let aux_basis = BasisSet::mock_ri(*size);

        // Config with minimal validation
        let fast_config = TwoCenterConfig {
            parallel: true,
            n_threads: 0,
            cache_blocking: true,
            symmetry_tolerance: 1e-10,
            pd_tolerance: 1e-8,
            ..Default::default()
        };

        // Config with strict validation
        let strict_config = TwoCenterConfig {
            parallel: true,
            n_threads: 0,
            cache_blocking: true,
            symmetry_tolerance: 1e-15,
            pd_tolerance: 1e-12,
            ..Default::default()
        };

        group.bench_with_input(BenchmarkId::new("fast", size), size, |b, _| {
            b.iter(|| {
                let result = compute_2center_integrals_with_config(&aux_basis, &fast_config)
                    .expect("Failed to compute integrals");
                black_box(result)
            });
        });

        group.bench_with_input(BenchmarkId::new("strict", size), size, |b, _| {
            b.iter(|| {
                let result = compute_2center_integrals_with_config(&aux_basis, &strict_config)
                    .expect("Failed to compute integrals");
                black_box(result)
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_sizes,
    benchmark_parallel,
    benchmark_block_sizes,
    benchmark_validation
);
criterion_main!(benches);
