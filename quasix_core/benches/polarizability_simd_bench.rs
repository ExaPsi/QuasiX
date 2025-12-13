//! Benchmarks for SIMD-optimized polarizability calculations
//!
//! This benchmark suite measures the performance improvements from SIMD optimization
//! across different system sizes and frequency configurations.

#![warn(clippy::all, clippy::pedantic, clippy::perf)]
#![allow(clippy::missing_docs_in_private_items)]
#![allow(clippy::cast_precision_loss)] // Acceptable in benchmarks with known small values

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::dielectric::{PolarizabilityRI, PolarizabilitySIMDCore};
use std::hint::black_box;
use std::time::Duration;

/// Generate test system with given dimensions
fn generate_test_system(
    nocc: usize,
    nvirt: usize,
    naux: usize,
) -> (Array2<f64>, Array1<f64>, Array1<f64>) {
    let n_trans = nocc * nvirt;

    // Generate realistic DF tensor
    let mut df_ia = Array2::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            // Gaussian-like decay
            let r = ((i as f64 - n_trans as f64 / 2.0).powi(2)
                + (p as f64 - naux as f64 / 2.0).powi(2))
            .sqrt();
            df_ia[[i, p]] = (-r / (naux as f64 / 4.0)).exp();
        }
    }

    // Realistic orbital energies (Hartree units)
    let e_occ = Array1::linspace(-2.0, -0.3, nocc);
    let e_virt = Array1::linspace(0.05, 2.0, nvirt);

    (df_ia, e_occ, e_virt)
}

/// Benchmark single frequency P0 calculation
fn bench_single_frequency(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_single_frequency");
    group.measurement_time(Duration::from_secs(10));
    group.sample_size(100);

    // Test different system sizes
    let test_cases = [
        ("small", 5, 10, 30),    // Small molecule
        ("medium", 10, 50, 100), // Medium molecule
        ("large", 20, 100, 300), // Large molecule
    ];

    let omega = Complex64::new(0.0, 1.0);

    for (name, nocc, nvirt, naux) in &test_cases {
        let (df_ia, e_occ, e_virt) = generate_test_system(*nocc, *nvirt, *naux);
        let n_elements = naux * naux;

        // Benchmark standard implementation
        group.throughput(Throughput::Elements(n_elements as u64));
        group.bench_with_input(
            BenchmarkId::new("standard", name),
            &(&df_ia, &e_occ, &e_virt),
            |b, (df, e_o, e_v)| {
                let calc = PolarizabilityRI::new(*nocc, *nvirt, *naux);
                b.iter(|| {
                    let result = calc.compute_p0(omega, df, e_o, e_v).unwrap();
                    black_box(result)
                });
            },
        );

        // Benchmark optimized implementation
        group.bench_with_input(
            BenchmarkId::new("optimized", name),
            &(&df_ia, &e_occ, &e_virt),
            |b, (df, e_o, e_v)| {
                let calc = PolarizabilityRI::new(*nocc, *nvirt, *naux);
                b.iter(|| {
                    let result = calc.compute_p0_optimized(omega, df, e_o, e_v).unwrap();
                    black_box(result)
                });
            },
        );

        // Benchmark SIMD implementation
        group.bench_with_input(
            BenchmarkId::new("simd", name),
            &(&df_ia, &e_occ, &e_virt),
            |b, (df, e_o, e_v)| {
                let calc = PolarizabilitySIMDCore::new(*nocc, *nvirt, *naux);
                b.iter(|| {
                    let result = calc.compute_p0_simd(omega, df, e_o, e_v).unwrap();
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark batch frequency P0 calculation
fn bench_batch_frequency(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_batch_frequency");
    group.measurement_time(Duration::from_secs(10));
    group.sample_size(50);

    // Medium system size
    let (nocc, nvirt, naux) = (10, 50, 100);
    let (df_ia, e_occ, e_virt) = generate_test_system(nocc, nvirt, naux);

    // Different batch sizes
    let batch_sizes = [4, 8, 16, 32];

    for n_freq in batch_sizes {
        let frequencies: Vec<Complex64> = (0..n_freq)
            .map(|i| Complex64::new(0.0, 0.1 * (i + 1) as f64))
            .collect();

        let total_elements = (naux * naux * n_freq) as u64;
        group.throughput(Throughput::Elements(total_elements));

        // Benchmark standard batch
        group.bench_with_input(
            BenchmarkId::new("standard_batch", n_freq),
            &(&frequencies, &df_ia, &e_occ, &e_virt),
            |b, (freqs, df, e_o, e_v)| {
                let calc = PolarizabilityRI::new(nocc, nvirt, naux);
                b.iter(|| {
                    let result = calc.compute_p0_batch(freqs, df, e_o, e_v).unwrap();
                    black_box(result)
                });
            },
        );

        // Benchmark SIMD batch
        group.bench_with_input(
            BenchmarkId::new("simd_batch", n_freq),
            &(&frequencies, &df_ia, &e_occ, &e_virt),
            |b, (freqs, df, e_o, e_v)| {
                let calc = PolarizabilitySIMDCore::new(nocc, nvirt, naux);
                b.iter(|| {
                    let result = calc.compute_p0_batch_simd(freqs, df, e_o, e_v).unwrap();
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark denominator computation specifically
fn bench_denominator_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("denominator_computation");

    let test_cases = [("small", 5, 10), ("medium", 10, 50), ("large", 20, 100)];

    for (name, nocc, nvirt) in &test_cases {
        let e_occ = Array1::linspace(-2.0, -0.3, *nocc);
        let e_virt = Array1::linspace(0.05, 2.0, *nvirt);
        let n_trans = nocc * nvirt;

        group.throughput(Throughput::Elements(n_trans as u64));

        // Benchmark scalar denominator computation
        group.bench_with_input(
            BenchmarkId::new("scalar", name),
            &(&e_occ, &e_virt),
            |b, (e_o, e_v)| {
                b.iter(|| {
                    let omega = Complex64::new(0.0, 1.0);
                    let eta = 1e-4;
                    let mut denominators = Vec::with_capacity(n_trans);

                    for i in 0..*nocc {
                        for a in 0..*nvirt {
                            let de = e_v[a] - e_o[i];
                            let denom = 2.0 / (de - omega - Complex64::new(0.0, eta));
                            denominators.push(black_box(denom));
                        }
                    }
                    black_box(denominators)
                });
            },
        );

        // Benchmark SIMD denominator computation (internal method)
        // We'll create a minimal SIMD calculator just for this
        group.bench_with_input(
            BenchmarkId::new("simd", name),
            &(&e_occ, &e_virt),
            |b, (e_o, e_v)| {
                let calc = PolarizabilitySIMDCore::new(*nocc, *nvirt, 10);
                let omega = Complex64::new(0.0, 1.0);

                // Create dummy df_ia just for validation
                let df_ia = Array2::zeros((n_trans, 10));

                b.iter(|| {
                    // The SIMD computation happens internally in compute_p0_simd
                    // We measure it indirectly through a minimal P0 calculation
                    let _ = calc.compute_p0_simd(omega, &df_ia, e_o, e_v).unwrap();
                });
            },
        );
    }

    group.finish();
}

/// Benchmark memory access patterns
fn bench_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_access");

    let (nocc, nvirt, naux) = (10, 50, 100);
    let n_trans = nocc * nvirt;
    let df_ia = Array2::from_elem((n_trans, naux), 0.1_f64);

    // Benchmark row-major access
    group.bench_function("row_major", |b| {
        b.iter(|| {
            let mut sum = 0.0;
            for i in 0..n_trans {
                for p in 0..naux {
                    sum += df_ia[[i, p]];
                }
            }
            black_box(sum)
        });
    });

    // Benchmark column-major access (cache-unfriendly)
    group.bench_function("column_major", |b| {
        b.iter(|| {
            let mut sum = 0.0;
            for p in 0..naux {
                for i in 0..n_trans {
                    sum += df_ia[[i, p]];
                }
            }
            black_box(sum)
        });
    });

    // Benchmark SIMD-aligned access
    group.bench_function("simd_aligned", |b| {
        b.iter(|| {
            let mut sum = 0.0;
            for i in 0..n_trans {
                // Process in chunks of 4 for SIMD
                for p_chunk in (0..naux).step_by(4) {
                    let end = (p_chunk + 4).min(naux);
                    for p in p_chunk..end {
                        sum += df_ia[[i, p]];
                    }
                }
            }
            black_box(sum)
        });
    });

    group.finish();
}

/// Benchmark different frequency types
fn bench_frequency_types(c: &mut Criterion) {
    let mut group = c.benchmark_group("frequency_types");

    let (nocc, nvirt, naux) = (10, 50, 100);
    let (df_ia, e_occ, e_virt) = generate_test_system(nocc, nvirt, naux);

    let test_frequencies = [
        ("static", Complex64::new(0.0, 1e-4)),
        ("imaginary", Complex64::new(0.0, 1.0)),
        ("real", Complex64::new(0.5, 1e-4)),
        ("complex", Complex64::new(0.5, 0.5)),
    ];

    for (freq_type, omega) in &test_frequencies {
        // Standard implementation
        group.bench_with_input(
            BenchmarkId::new("standard", freq_type),
            &(&df_ia, &e_occ, &e_virt, omega),
            |b, (df, e_o, e_v, &omega)| {
                let calc = PolarizabilityRI::new(nocc, nvirt, naux);
                b.iter(|| {
                    let result = calc.compute_p0(omega, df, e_o, e_v).unwrap();
                    black_box(result)
                });
            },
        );

        // SIMD implementation
        group.bench_with_input(
            BenchmarkId::new("simd", freq_type),
            &(&df_ia, &e_occ, &e_virt, omega),
            |b, (df, e_o, e_v, &omega)| {
                let calc = PolarizabilitySIMDCore::new(nocc, nvirt, naux);
                b.iter(|| {
                    let result = calc.compute_p0_simd(omega, df, e_o, e_v).unwrap();
                    black_box(result)
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
    bench_denominator_computation,
    bench_memory_patterns,
    bench_frequency_types
);
criterion_main!(benches);
