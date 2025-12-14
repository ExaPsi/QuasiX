//! Benchmarks for parallel MO transformation optimizations

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array3;
use quasix_core::df::mo_transform::{
    generate_mock_mo_coefficients, transform_mo_3center, transform_mo_3center_with_config,
    TransformConfig,
};
use quasix_core::df::parallel::{
    transform_mo_3center_parallel_optimized, BlockingStrategy, ThreadPoolConfig,
};
use std::hint::black_box;
use std::time::Duration;

/// Generate mock 3-center integrals
fn generate_mock_j3c(n_ao: usize, n_aux: usize) -> Array3<f64> {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);

    let mut j3c = Array3::<f64>::zeros((n_ao, n_ao, n_aux));
    for i in 0..n_ao {
        for j in 0..=i {
            for p in 0..n_aux {
                let val = rng.random::<f64>() * 0.1;
                j3c[[i, j, p]] = val;
                j3c[[j, i, p]] = val; // Symmetry
            }
        }
    }
    j3c
}

fn bench_mo_transform_small(c: &mut Criterion) {
    let mut group = c.benchmark_group("mo_transform_small");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));

    // Small system: 20 AO, 50 aux
    let n_ao = 20;
    let n_aux = 50;
    let n_occ = 8;
    let n_vir = 12;

    let j3c_ao = generate_mock_j3c(n_ao, n_aux);
    let c_occ = generate_mock_mo_coefficients(n_ao, n_occ, 100);
    let c_vir = generate_mock_mo_coefficients(n_ao, n_vir, 101);

    group.bench_function("default", |b| {
        b.iter(|| transform_mo_3center(black_box(&j3c_ao), black_box(&c_occ), black_box(&c_vir)))
    });

    group.bench_function("small_system_config", |b| {
        let config = TransformConfig::small_system();
        b.iter(|| {
            transform_mo_3center_with_config(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                black_box(&config),
            )
        })
    });

    group.bench_function("advanced_parallel", |b| {
        b.iter(|| {
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::default()),
                Some(BlockingStrategy::default()),
            )
        })
    });

    group.finish();
}

fn bench_mo_transform_medium(c: &mut Criterion) {
    let mut group = c.benchmark_group("mo_transform_medium");
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(10));
    group.sample_size(20);

    // Medium system: 100 AO, 200 aux
    let n_ao = 100;
    let n_aux = 200;
    let n_occ = 40;
    let n_vir = 60;

    let j3c_ao = generate_mock_j3c(n_ao, n_aux);
    let c_occ = generate_mock_mo_coefficients(n_ao, n_occ, 102);
    let c_vir = generate_mock_mo_coefficients(n_ao, n_vir, 103);

    group.bench_function("default", |b| {
        b.iter(|| transform_mo_3center(black_box(&j3c_ao), black_box(&c_occ), black_box(&c_vir)))
    });

    group.bench_function("memory_bound", |b| {
        b.iter(|| {
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::memory_bound()),
                Some(BlockingStrategy::default()),
            )
        })
    });

    group.bench_function("compute_bound", |b| {
        b.iter(|| {
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::compute_bound()),
                Some(BlockingStrategy::default()),
            )
        })
    });

    group.finish();
}

fn bench_mo_transform_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("mo_transform_scaling");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(10);

    for n_ao in [20, 40, 60, 80].iter() {
        let n_aux = n_ao * 3;
        let n_occ = n_ao / 3;
        let n_vir = n_ao - n_occ;

        let j3c_ao = generate_mock_j3c(*n_ao, n_aux);
        let c_occ = generate_mock_mo_coefficients(*n_ao, n_occ, 104);
        let c_vir = generate_mock_mo_coefficients(*n_ao, n_vir, 105);

        group.bench_with_input(BenchmarkId::new("default", n_ao), n_ao, |b, _| {
            b.iter(|| {
                transform_mo_3center(black_box(&j3c_ao), black_box(&c_occ), black_box(&c_vir))
            })
        });

        group.bench_with_input(BenchmarkId::new("advanced", n_ao), n_ao, |b, _| {
            b.iter(|| {
                transform_mo_3center_parallel_optimized(
                    black_box(&j3c_ao),
                    black_box(&c_occ),
                    black_box(&c_vir),
                    Some(ThreadPoolConfig::default()),
                    Some(BlockingStrategy::default()),
                )
            })
        });
    }

    group.finish();
}

fn bench_blocking_strategies(c: &mut Criterion) {
    let mut group = c.benchmark_group("blocking_strategies");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(20);

    // Test different blocking strategies
    let n_ao = 80;
    let n_aux = 240;
    let n_occ = 30;
    let n_vir = 50;

    let j3c_ao = generate_mock_j3c(n_ao, n_aux);
    let c_occ = generate_mock_mo_coefficients(n_ao, n_occ, 106);
    let c_vir = generate_mock_mo_coefficients(n_ao, n_vir, 107);

    group.bench_function("default_blocking", |b| {
        b.iter(|| {
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::default()),
                Some(BlockingStrategy::default()),
            )
        })
    });

    group.bench_function("intel_xeon_blocking", |b| {
        b.iter(|| {
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::default()),
                Some(BlockingStrategy::intel_xeon()),
            )
        })
    });

    group.bench_function("amd_epyc_blocking", |b| {
        b.iter(|| {
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::default()),
                Some(BlockingStrategy::amd_epyc()),
            )
        })
    });

    group.finish();
}

fn bench_thread_configs(c: &mut Criterion) {
    let mut group = c.benchmark_group("thread_configs");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(20);

    let n_ao = 60;
    let n_aux = 180;
    let n_occ = 25;
    let n_vir = 35;

    let j3c_ao = generate_mock_j3c(n_ao, n_aux);
    let c_occ = generate_mock_mo_coefficients(n_ao, n_occ, 108);
    let c_vir = generate_mock_mo_coefficients(n_ao, n_vir, 109);

    // Test different thread configurations
    for num_threads in [1, 2, 4, 8].iter() {
        let thread_config = ThreadPoolConfig {
            num_threads: *num_threads,
            ..Default::default()
        };

        group.bench_with_input(
            BenchmarkId::new("threads", num_threads),
            num_threads,
            |b, _| {
                b.iter(|| {
                    transform_mo_3center_parallel_optimized(
                        black_box(&j3c_ao),
                        black_box(&c_occ),
                        black_box(&c_vir),
                        Some(thread_config.clone()),
                        Some(BlockingStrategy::default()),
                    )
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_mo_transform_small,
    bench_mo_transform_medium,
    bench_mo_transform_scaling,
    bench_blocking_strategies,
    bench_thread_configs
);
criterion_main!(benches);
