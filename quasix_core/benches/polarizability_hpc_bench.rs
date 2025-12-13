//! Benchmarks for HPC-optimized polarizability computation
#![allow(clippy::cast_precision_loss)] // Acceptable in benchmarks with known small values

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::dielectric::{PolarizabilityHPC, PolarizabilityRI};
use std::hint::black_box;
use std::time::Duration;

/// Generate test data for benchmarking
fn generate_test_data(
    nocc: usize,
    nvirt: usize,
    naux: usize,
) -> (Array2<f64>, Array1<f64>, Array1<f64>) {
    let n_trans = nocc * nvirt;

    // Generate realistic DF tensor with some structure
    let mut df_ia = Array2::<f64>::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            // Gaussian-like decay
            df_ia[[i, p]] = (-(((i as f64) - (p as f64)).powi(2)) / (10.0 * naux as f64)).exp();
        }
    }

    // Realistic orbital energies
    let e_occ = Array1::linspace(-1.0, -0.3, nocc);
    let e_virt = Array1::linspace(0.05, 2.0, nvirt);

    (df_ia, e_occ, e_virt)
}

/// Benchmark different algorithm variants
fn bench_algorithms(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_algorithms");
    group.measurement_time(Duration::from_secs(10));

    // Test different system sizes
    let test_cases = vec![
        ("small", 10, 20, 100),    // Small molecule
        ("medium", 50, 100, 500),  // Medium molecule
        ("large", 100, 200, 1000), // Large molecule
    ];

    for (name, nocc, nvirt, naux) in test_cases {
        let (df_ia, e_occ, e_virt) = generate_test_data(nocc, nvirt, naux);
        let omega = Complex64::new(0.0, 1.0);

        // Benchmark standard implementation
        group.bench_function(BenchmarkId::new("standard", name), |b| {
            let calc = PolarizabilityRI::new(nocc, nvirt, naux);
            b.iter(|| calc.compute_p0(black_box(omega), &df_ia, &e_occ, &e_virt));
        });

        // Benchmark HPC implementation
        group.bench_function(BenchmarkId::new("hpc_blocked", name), |b| {
            let calc = PolarizabilityHPC::new(nocc, nvirt, naux);
            b.iter(|| calc.compute_p0_blocked(black_box(omega), &df_ia, &e_occ, &e_virt));
        });

        // Benchmark BLAS implementation
        group.bench_function(BenchmarkId::new("hpc_blas", name), |b| {
            let calc = PolarizabilityHPC::new(nocc, nvirt, naux);
            b.iter(|| calc.compute_p0_blas(black_box(omega), &df_ia, &e_occ, &e_virt));
        });
    }

    group.finish();
}

/// Benchmark parallel scaling
fn bench_parallel_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_parallel_scaling");
    group.measurement_time(Duration::from_secs(10));

    let nocc = 50;
    let nvirt = 100;
    let naux = 1000;
    let (df_ia, e_occ, e_virt) = generate_test_data(nocc, nvirt, naux);

    // Test different frequency batch sizes
    let batch_sizes = vec![1, 4, 8, 16, 32];

    for batch_size in batch_sizes {
        let omega_batch: Vec<Complex64> = (0..batch_size)
            .map(|i| Complex64::new(0.0, 0.5 + 0.1 * i as f64))
            .collect();

        group.throughput(Throughput::Elements(batch_size as u64));

        group.bench_function(BenchmarkId::new("batch", batch_size), |b| {
            let calc = PolarizabilityHPC::new(nocc, nvirt, naux);
            b.iter(|| calc.compute_p0_batch(black_box(&omega_batch), &df_ia, &e_occ, &e_virt));
        });
    }

    group.finish();
}

/// Benchmark cache blocking effectiveness
fn bench_cache_blocking(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_cache_blocking");
    group.measurement_time(Duration::from_secs(5));

    let nocc = 50;
    let nvirt = 100;
    let naux = 1000;
    let (df_ia, e_occ, e_virt) = generate_test_data(nocc, nvirt, naux);
    let omega = Complex64::new(0.0, 1.0);

    // Test different block sizes
    let block_sizes = vec![16, 32, 64, 128, 256];

    for block_size in block_sizes {
        group.bench_function(BenchmarkId::new("block_size", block_size), |b| {
            let mut calc = PolarizabilityHPC::new(nocc, nvirt, naux);
            calc.config.cache_block_size = block_size;
            b.iter(|| calc.compute_p0_blocked(black_box(omega), &df_ia, &e_occ, &e_virt));
        });
    }

    group.finish();
}

/// Benchmark SIMD effectiveness
fn bench_simd(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_simd");
    group.measurement_time(Duration::from_secs(5));

    let test_cases = vec![(50, 100, 500), (100, 200, 1000)];

    for (nocc, nvirt, naux) in test_cases {
        let (df_ia, e_occ, e_virt) = generate_test_data(nocc, nvirt, naux);
        let omega = Complex64::new(0.0, 1.0);
        let name = format!("{}x{}x{}", nocc, nvirt, naux);

        // With SIMD
        group.bench_function(BenchmarkId::new("simd_on", &name), |b| {
            let mut calc = PolarizabilityHPC::new(nocc, nvirt, naux);
            calc.config.use_simd = true;
            b.iter(|| calc.compute_p0_blocked(black_box(omega), &df_ia, &e_occ, &e_virt));
        });

        // Without SIMD
        group.bench_function(BenchmarkId::new("simd_off", &name), |b| {
            let mut calc = PolarizabilityHPC::new(nocc, nvirt, naux);
            calc.config.use_simd = false;
            b.iter(|| calc.compute_p0_blocked(black_box(omega), &df_ia, &e_occ, &e_virt));
        });
    }

    group.finish();
}

/// Benchmark memory pool effectiveness
fn bench_memory_pool(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_memory_pool");
    group.measurement_time(Duration::from_secs(5));

    let nocc = 50;
    let nvirt = 100;
    let naux = 500;
    let (df_ia, e_occ, e_virt) = generate_test_data(nocc, nvirt, naux);

    // Multiple frequencies to stress memory allocation
    let omega_batch: Vec<Complex64> = (0..16)
        .map(|i| Complex64::new(0.0, 0.5 + 0.1 * i as f64))
        .collect();

    // With memory pool (default)
    group.bench_function("with_pool", |b| {
        let calc = PolarizabilityHPC::new(nocc, nvirt, naux);
        b.iter(|| calc.compute_p0_batch(black_box(&omega_batch), &df_ia, &e_occ, &e_virt));
    });

    // Without memory pool (pool_size = 0)
    group.bench_function("without_pool", |b| {
        let mut calc = PolarizabilityHPC::new(nocc, nvirt, naux);
        calc.config.pool_size = 0;
        b.iter(|| calc.compute_p0_batch(black_box(&omega_batch), &df_ia, &e_occ, &e_virt));
    });

    group.finish();
}

/// Benchmark real-world GW100 molecules
fn bench_gw100_molecules(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_gw100");
    group.measurement_time(Duration::from_secs(10));
    group.sample_size(10);

    // Realistic GW100 molecule sizes
    let molecules = vec![
        ("H2O", 5, 209, 283),      // cc-pVTZ
        ("NH3", 5, 209, 283),      // cc-pVTZ
        ("CH4", 5, 209, 283),      // cc-pVTZ
        ("CO", 7, 207, 283),       // cc-pVTZ
        ("Benzene", 21, 243, 504), // def2-TZVP
    ];

    for (name, nocc, nvirt, naux) in molecules {
        let (df_ia, e_occ, e_virt) = generate_test_data(nocc, nvirt, naux);

        // Typical evGW frequency grid (48 points)
        let omega_grid: Vec<Complex64> = (0..48)
            .map(|i| {
                let xi = (i as f64 + 0.5) * 0.5; // Imaginary frequency grid
                Complex64::new(0.0, xi)
            })
            .collect();

        group.bench_function(BenchmarkId::new("molecule", name), |b| {
            let calc = PolarizabilityHPC::new(nocc, nvirt, naux);
            b.iter(|| calc.compute_p0_batch(black_box(&omega_grid), &df_ia, &e_occ, &e_virt));
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_algorithms,
    bench_parallel_scaling,
    bench_cache_blocking,
    bench_simd,
    bench_memory_pool,
    bench_gw100_molecules
);

criterion_main!(benches);
