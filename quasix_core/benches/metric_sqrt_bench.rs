use criterion::{criterion_group, criterion_main, Criterion};
use ndarray::{Array1, Array2};
use quasix_core::linalg::MetricSqrt;
use std::hint::black_box;

/// Generate a random symmetric positive definite matrix
fn generate_spd_matrix(n: usize) -> Array2<f64> {
    let mut rng = 42u64; // Deterministic for reproducibility
    let mut random = || {
        rng = (rng.wrapping_mul(1_103_515_245).wrapping_add(12345)) & 0x7fff_ffff;
        (rng as f64) / f64::from(0x7fff_ffff_u32) - 0.5
    };

    // Generate random matrix A
    let mut a = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            a[[i, j]] = random();
        }
    }

    // Make it SPD: M = A^T A + I
    a.t().dot(&a) + Array2::<f64>::eye(n)
}

fn bench_metric_sqrt_compute(c: &mut Criterion) {
    let mut group = c.benchmark_group("metric_sqrt_compute");

    for size in [10, 50, 100, 200].iter() {
        let matrix = generate_spd_matrix(*size);
        group.bench_function(format!("size_{}", size), |b| {
            b.iter(|| MetricSqrt::compute(black_box(&matrix), None).unwrap());
        });
    }
    group.finish();
}

fn bench_apply_sqrt_vector(c: &mut Criterion) {
    let mut group = c.benchmark_group("apply_sqrt_vector");

    for size in [50, 100, 200, 500].iter() {
        let matrix = generate_spd_matrix(*size);
        let metric_sqrt = MetricSqrt::compute(&matrix, None).unwrap();
        let x = Array1::ones(*size);

        group.bench_function(format!("size_{}", size), |b| {
            b.iter(|| metric_sqrt.apply_sqrt(black_box(&x.view())).unwrap());
        });
    }
    group.finish();
}

fn bench_apply_sqrt_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("apply_sqrt_matrix");

    for size in [50, 100, 200].iter() {
        let matrix = generate_spd_matrix(*size);
        let metric_sqrt = MetricSqrt::compute(&matrix, None).unwrap();
        let x = Array2::ones((*size, 10)); // Apply to 10 columns at once

        group.bench_function(format!("size_{}_cols_10", size), |b| {
            b.iter(|| metric_sqrt.apply_sqrt_matrix(black_box(&x.view())).unwrap());
        });
    }
    group.finish();
}

fn bench_apply_inv_sqrt_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("apply_inv_sqrt_matrix");

    for size in [50, 100, 200].iter() {
        let matrix = generate_spd_matrix(*size);
        let metric_sqrt = MetricSqrt::compute(&matrix, None).unwrap();
        let x = Array2::ones((*size, 10)); // Apply to 10 columns at once

        group.bench_function(format!("size_{}_cols_10", size), |b| {
            b.iter(|| {
                metric_sqrt
                    .apply_inv_sqrt_matrix(black_box(&x.view()))
                    .unwrap()
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_metric_sqrt_compute,
    bench_apply_sqrt_vector,
    bench_apply_sqrt_matrix,
    bench_apply_inv_sqrt_matrix
);
criterion_main!(benches);
