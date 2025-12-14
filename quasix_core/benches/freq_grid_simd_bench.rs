use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array1;
use quasix_core::freq::grid_optimized::{
    ImaginaryAxisGrid, MinimaxGrid, OptimizedGLGrid, TransformType,
};
use std::hint::black_box;

fn bench_legendre_polynomial_caching(c: &mut Criterion) {
    let mut group = c.benchmark_group("legendre_polynomial_caching");

    for n in [32, 64, 128, 256] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            b.iter(|| {
                let grid = OptimizedGLGrid::new(n, true).unwrap();
                black_box(grid.cached_polynomials);
            });
        });
    }

    group.finish();
}

fn bench_integration_methods(c: &mut Criterion) {
    let mut group = c.benchmark_group("integration");

    for n in [32, 64, 128] {
        let grid = OptimizedGLGrid::new(n, false).unwrap();
        let f = |x: f64| x.sin() * x.cos() + x * x;

        // Benchmark serial integration
        group.bench_with_input(BenchmarkId::new("serial", n), &grid, |b, grid| {
            b.iter(|| grid.integrate_cached(f));
        });

        // Benchmark SIMD integration
        group.bench_with_input(BenchmarkId::new("simd", n), &grid, |b, grid| {
            b.iter(|| grid.integrate_simd(f));
        });

        // Benchmark parallel integration
        group.bench_with_input(BenchmarkId::new("parallel", n), &grid, |b, grid| {
            b.iter(|| grid.integrate_parallel(f));
        });

        // Benchmark dot product integration
        let values = Array1::from_vec(grid.nodes.iter().map(|&x| f(x)).collect());
        group.bench_with_input(BenchmarkId::new("dot_product_simd", n), &grid, |b, grid| {
            b.iter(|| grid.integrate_dot_product(&values));
        });
    }

    group.finish();
}

fn bench_frequency_transformations(c: &mut Criterion) {
    let mut group = c.benchmark_group("frequency_transforms");

    for n in [32, 64, 128] {
        // Benchmark linear transformation
        group.bench_with_input(BenchmarkId::new("linear", n), &n, |b, &n| {
            b.iter(|| {
                ImaginaryAxisGrid::new(n, TransformType::Linear { omega_max: 100.0 }).unwrap()
            });
        });

        // Benchmark tan transformation
        group.bench_with_input(BenchmarkId::new("tan", n), &n, |b, &n| {
            b.iter(|| ImaginaryAxisGrid::new(n, TransformType::Tan { omega_max: 100.0 }).unwrap());
        });

        // Benchmark double exponential transformation
        group.bench_with_input(BenchmarkId::new("double_exp", n), &n, |b, &n| {
            b.iter(|| ImaginaryAxisGrid::new(n, TransformType::DoubleExp { scale: 1.0 }).unwrap());
        });
    }

    group.finish();
}

fn bench_barycentric_interpolation(c: &mut Criterion) {
    let mut group = c.benchmark_group("barycentric_interpolation");

    for n in [16, 32, 64] {
        let grid = MinimaxGrid::chebyshev(n);
        let values = Array1::from_vec(grid.points.iter().map(|&x| x.sin()).collect());
        let test_points: Vec<f64> = (0..100).map(|i| -0.9 + 1.8 * i as f64 / 99.0).collect();

        // Benchmark single point interpolation
        group.bench_with_input(
            BenchmarkId::new("single", n),
            &(&grid, &values),
            |b, (grid, values)| {
                b.iter(|| {
                    for &x in test_points.iter().take(10) {
                        black_box(grid.interpolate(values, x));
                    }
                });
            },
        );

        // Benchmark batch interpolation
        group.bench_with_input(
            BenchmarkId::new("batch", n),
            &(&grid, &values),
            |b, (grid, values)| {
                b.iter(|| {
                    black_box(grid.interpolate_batch(values, &test_points));
                });
            },
        );
    }

    group.finish();
}

fn bench_simd_speedup_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_speedup");

    let n = 128;
    let grid = OptimizedGLGrid::new(n, true).unwrap();

    // Complex function for integration
    let f = |x: f64| {
        let mut sum = 0.0;
        for k in 1..10 {
            sum += (k as f64 * x).sin() / k as f64;
        }
        sum
    };

    // Create values for dot product
    let values = Array1::from_vec(grid.nodes.iter().map(|&x| f(x)).collect());

    // Measure scalar baseline
    group.bench_function("scalar_baseline", |b| {
        b.iter(|| {
            let mut sum = 0.0;
            for i in 0..n {
                sum += grid.weights[i] * f(grid.nodes[i]);
            }
            black_box(sum)
        });
    });

    // Measure SIMD optimized
    group.bench_function("simd_optimized", |b| {
        b.iter(|| black_box(grid.integrate_simd(f)));
    });

    // Measure SIMD dot product
    group.bench_function("simd_dot_product", |b| {
        b.iter(|| black_box(grid.integrate_dot_product(&values)));
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_legendre_polynomial_caching,
    bench_integration_methods,
    bench_frequency_transformations,
    bench_barycentric_interpolation,
    bench_simd_speedup_comparison
);
criterion_main!(benches);
