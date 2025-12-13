//! Benchmarks for condition monitoring optimization
//!
//! This benchmark suite tests the performance of various condition number
//! estimation methods to ensure < 5% overhead for production use.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ndarray::{Array1, Array2};
use ndarray_linalg::{Eigh, Norm, UPLO};
use num_complex::Complex64;
use rand::prelude::*;
use std::time::Duration;

/// Generate a well-conditioned test matrix
fn generate_well_conditioned(n: usize) -> Array2<Complex64> {
    let mut rng = StdRng::seed_from_u64(42);
    let mut matrix = Array2::zeros((n, n));

    // Generate random Hermitian matrix with controlled eigenvalues
    for i in 0..n {
        matrix[[i, i]] = Complex64::new(1.0 + rng.gen::<f64>() * 0.5, 0.0);
        for j in i + 1..n {
            let val = Complex64::new(
                rng.gen::<f64>() * 0.1,
                rng.gen::<f64>() * 0.1,
            );
            matrix[[i, j]] = val;
            matrix[[j, i]] = val.conj();
        }
    }
    matrix
}

/// Generate an ill-conditioned test matrix
fn generate_ill_conditioned(n: usize, condition: f64) -> Array2<Complex64> {
    let mut rng = StdRng::seed_from_u64(42);

    // Use eigenvalue decomposition to control condition number
    let mut eigenvalues = Array1::zeros(n);
    eigenvalues[0] = condition;
    for i in 1..n {
        eigenvalues[i] = 1.0 + rng.gen::<f64>() * (condition - 1.0).log10().exp();
    }
    eigenvalues[n - 1] = 1.0; // Ensure exact condition number

    // Generate random unitary matrix
    let mut q = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            q[[i, j]] = Complex64::new(
                rng.gen::<f64>() - 0.5,
                rng.gen::<f64>() - 0.5,
            );
        }
    }

    // QR decomposition to get orthogonal matrix
    // (simplified - in practice use proper QR)
    for i in 0..n {
        let mut col = q.column_mut(i);
        let norm = col.dot(&col).norm().sqrt();
        col /= Complex64::new(norm, 0.0);
    }

    // Construct matrix with specified eigenvalues
    let lambda = Array2::from_diag(&eigenvalues.mapv(|x| Complex64::new(x, 0.0)));
    let q_h = q.t().mapv(|x| x.conj());
    q.dot(&lambda).dot(&q_h)
}

/// Baseline: Full SVD computation (most accurate but slowest)
fn condition_svd(matrix: &Array2<Complex64>) -> f64 {
    use ndarray_linalg::SVD;

    let (_, s, _) = matrix.svd(false, false).unwrap();
    let max_sv = s[0];
    let min_sv = s[s.len() - 1];

    if min_sv < 1e-15 {
        f64::INFINITY
    } else {
        max_sv / min_sv
    }
}

/// Power iteration for largest eigenvalue
fn power_iteration_largest(matrix: &Array2<Complex64>, iterations: usize) -> f64 {
    let n = matrix.nrows();
    let mut v = Array1::from_elem(n, Complex64::new(1.0, 0.0));
    v /= Complex64::new((n as f64).sqrt(), 0.0);

    let mut lambda = 0.0;

    for _ in 0..iterations {
        let v_new = matrix.dot(&v);
        let norm = v_new.norm_l2();

        if norm < 1e-12 {
            break;
        }

        v = &v_new / Complex64::new(norm, 0.0);

        // Rayleigh quotient
        let mv = matrix.dot(&v);
        lambda = v.mapv(|x| x.conj()).dot(&mv).re;
    }

    lambda
}

/// Inverse power iteration for smallest eigenvalue
fn inverse_power_iteration(matrix: &Array2<Complex64>, iterations: usize) -> f64 {
    use ndarray_linalg::Inverse;

    let n = matrix.nrows();
    let mut v = Array1::from_elem(n, Complex64::new(1.0, 0.0));
    v /= Complex64::new((n as f64).sqrt(), 0.0);

    // Add small regularization for stability
    let mut m_reg = matrix.clone();
    for i in 0..n {
        m_reg[[i, i]] += Complex64::new(1e-10, 0.0);
    }

    let m_inv = match m_reg.inv() {
        Ok(inv) => inv,
        Err(_) => return 0.0, // Matrix is singular
    };

    let mut lambda = 0.0;

    for _ in 0..iterations {
        let v_new = m_inv.dot(&v);
        let norm = v_new.norm_l2();

        if norm < 1e-12 {
            break;
        }

        v = &v_new / Complex64::new(norm, 0.0);

        // Rayleigh quotient for inverse
        let mv = m_inv.dot(&v);
        let lambda_inv = v.mapv(|x| x.conj()).dot(&mv).re;

        if lambda_inv.abs() > 1e-12 {
            lambda = 1.0 / lambda_inv;
        }
    }

    lambda
}

/// Optimized power iteration with BLAS
fn power_iteration_blas(matrix: &Array2<Complex64>, iterations: usize) -> f64 {
    let n = matrix.nrows();
    let mut v = Array1::from_elem(n, Complex64::new(1.0, 0.0));
    let mut v_new = Array1::zeros(n);
    v /= Complex64::new((n as f64).sqrt(), 0.0);

    let mut lambda = 0.0;

    for _ in 0..iterations {
        // Use BLAS ZGEMV for matrix-vector product
        unsafe {
            let alpha = Complex64::new(1.0, 0.0);
            let beta = Complex64::new(0.0, 0.0);

            cblas_sys::cblas_zgemv(
                cblas_sys::CBLAS_LAYOUT::CblasRowMajor,
                cblas_sys::CBLAS_TRANSPOSE::CblasNoTrans,
                n as i32,
                n as i32,
                &alpha as *const _ as *const _,
                matrix.as_ptr() as *const _,
                n as i32,
                v.as_ptr() as *const _,
                1,
                &beta as *const _ as *const _,
                v_new.as_mut_ptr() as *mut _,
                1,
            );
        }

        // Compute norm using BLAS DZNRM2
        let norm = unsafe {
            cblas_sys::cblas_dznrm2(n as i32, v_new.as_ptr() as *const _, 1)
        };

        if norm < 1e-12 {
            break;
        }

        // Scale vector
        v_new /= Complex64::new(norm, 0.0);
        std::mem::swap(&mut v, &mut v_new);

        // Rayleigh quotient
        let mv = matrix.dot(&v);
        lambda = v.mapv(|x| x.conj()).dot(&mv).re;
    }

    lambda
}

/// Fast heuristic based on diagonal dominance
fn quick_condition_check(matrix: &Array2<Complex64>) -> bool {
    let n = matrix.nrows();
    let mut min_dominance = f64::INFINITY;

    for i in 0..n {
        let diag = matrix[[i, i]].norm();
        let off_diag: f64 = (0..n)
            .filter(|&j| j != i)
            .map(|j| matrix[[i, j]].norm())
            .sum();

        if diag <= off_diag {
            return false; // Not diagonally dominant, likely ill-conditioned
        }

        let dominance = diag / (diag + off_diag);
        min_dominance = min_dominance.min(dominance);
    }

    // If strongly diagonally dominant, matrix is well-conditioned
    min_dominance > 0.5
}

/// Gershgorin circle theorem for condition estimation
fn gershgorin_estimate(matrix: &Array2<Complex64>) -> f64 {
    let n = matrix.nrows();
    let mut min_radius = f64::INFINITY;
    let mut max_radius = 0.0;

    for i in 0..n {
        let diag = matrix[[i, i]].norm();
        let off_diag_sum: f64 = (0..n)
            .filter(|&j| j != i)
            .map(|j| matrix[[i, j]].norm())
            .sum();

        let lower = (diag - off_diag_sum).max(1e-12);
        let upper = diag + off_diag_sum;

        min_radius = min_radius.min(lower);
        max_radius = max_radius.max(upper);
    }

    max_radius / min_radius
}

/// Combined fast condition estimator
fn fast_condition_estimate(matrix: &Array2<Complex64>) -> f64 {
    // First, quick check
    if quick_condition_check(matrix) {
        return 1.0; // Well-conditioned, no need for accurate estimate
    }

    // Use Gershgorin for rough estimate
    let gershgorin = gershgorin_estimate(matrix);

    if gershgorin < 1e4 {
        return gershgorin; // Moderately conditioned, Gershgorin is sufficient
    }

    // For potentially ill-conditioned, use power iteration
    let max_eig = power_iteration_blas(matrix, 10);
    let min_eig = inverse_power_iteration(matrix, 10);

    if min_eig.abs() < 1e-12 {
        f64::INFINITY
    } else {
        max_eig / min_eig
    }
}

/// Benchmark different matrix sizes
fn bench_condition_estimation(c: &mut Criterion) {
    let mut group = c.benchmark_group("condition_estimation");
    group.measurement_time(Duration::from_secs(10));

    for n in [10, 50, 100, 200, 500].iter() {
        let well_conditioned = generate_well_conditioned(*n);
        let ill_conditioned = generate_ill_conditioned(*n, 1e8);

        // Benchmark SVD (baseline)
        group.bench_with_input(
            BenchmarkId::new("svd", n),
            &well_conditioned,
            |b, matrix| b.iter(|| condition_svd(black_box(matrix))),
        );

        // Benchmark power iteration
        group.bench_with_input(
            BenchmarkId::new("power_iteration", n),
            &well_conditioned,
            |b, matrix| b.iter(|| {
                let max_eig = power_iteration_largest(black_box(matrix), 20);
                let min_eig = inverse_power_iteration(black_box(matrix), 20);
                max_eig / min_eig.max(1e-12)
            }),
        );

        // Benchmark optimized power iteration
        group.bench_with_input(
            BenchmarkId::new("power_iteration_blas", n),
            &well_conditioned,
            |b, matrix| b.iter(|| {
                let max_eig = power_iteration_blas(black_box(matrix), 20);
                max_eig // Just largest eigenvalue for speed
            }),
        );

        // Benchmark Gershgorin
        group.bench_with_input(
            BenchmarkId::new("gershgorin", n),
            &well_conditioned,
            |b, matrix| b.iter(|| gershgorin_estimate(black_box(matrix))),
        );

        // Benchmark quick check
        group.bench_with_input(
            BenchmarkId::new("quick_check", n),
            &well_conditioned,
            |b, matrix| b.iter(|| quick_condition_check(black_box(matrix))),
        );

        // Benchmark combined fast estimator
        group.bench_with_input(
            BenchmarkId::new("fast_estimate", n),
            &well_conditioned,
            |b, matrix| b.iter(|| fast_condition_estimate(black_box(matrix))),
        );

        // Also test on ill-conditioned matrices
        group.bench_with_input(
            BenchmarkId::new("fast_estimate_ill", n),
            &ill_conditioned,
            |b, matrix| b.iter(|| fast_condition_estimate(black_box(matrix))),
        );
    }

    group.finish();
}

/// Benchmark regularization overhead
fn bench_regularization(c: &mut Criterion) {
    let mut group = c.benchmark_group("regularization");

    for n in [50, 100, 200, 500].iter() {
        let matrix = generate_ill_conditioned(*n, 1e12);

        // Benchmark in-place diagonal modification
        group.bench_with_input(
            BenchmarkId::new("inplace_diagonal", n),
            &matrix,
            |b, matrix| b.iter(|| {
                let mut m = black_box(matrix.clone());
                let alpha = 1e-10;
                for i in 0..*n {
                    m[[i, i]] += Complex64::new(alpha, 0.0);
                }
                m
            }),
        );

        // Benchmark matrix copy + modification
        group.bench_with_input(
            BenchmarkId::new("copy_modify", n),
            &matrix,
            |b, matrix| b.iter(|| {
                let m = black_box(matrix.clone());
                let alpha = 1e-10;
                let identity = Array2::eye(*n).mapv(|x| Complex64::new(x * alpha, 0.0));
                m + identity
            }),
        );
    }

    group.finish();
}

/// Benchmark end-to-end with different strategies
fn bench_end_to_end(c: &mut Criterion) {
    let mut group = c.benchmark_group("end_to_end");
    group.measurement_time(Duration::from_secs(10));

    let n = 200; // Typical size for QuasiX
    let well_conditioned = generate_well_conditioned(n);
    let ill_conditioned = generate_ill_conditioned(n, 1e10);

    // Simulate full computation with condition monitoring
    group.bench_function("with_full_svd", |b| {
        b.iter(|| {
            let matrix = black_box(&well_conditioned);
            let cond = condition_svd(matrix);

            if cond > 1e12 {
                // Apply regularization
                let mut m = matrix.clone();
                for i in 0..n {
                    m[[i, i]] += Complex64::new(1e-10, 0.0);
                }
                m
            } else {
                matrix.clone()
            }
        })
    });

    group.bench_function("with_fast_estimate", |b| {
        b.iter(|| {
            let matrix = black_box(&well_conditioned);
            let cond = fast_condition_estimate(matrix);

            if cond > 1e12 {
                // Apply regularization
                let mut m = matrix.clone();
                for i in 0..n {
                    m[[i, i]] += Complex64::new(1e-10, 0.0);
                }
                m
            } else {
                matrix.clone()
            }
        })
    });

    group.bench_function("with_quick_check_only", |b| {
        b.iter(|| {
            let matrix = black_box(&well_conditioned);

            if !quick_condition_check(matrix) {
                // Apply regularization
                let mut m = matrix.clone();
                for i in 0..n {
                    m[[i, i]] += Complex64::new(1e-10, 0.0);
                }
                m
            } else {
                matrix.clone()
            }
        })
    });

    group.bench_function("no_monitoring", |b| {
        b.iter(|| {
            let matrix = black_box(&well_conditioned);
            matrix.clone()
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_condition_estimation,
    bench_regularization,
    bench_end_to_end
);
criterion_main!(benches);