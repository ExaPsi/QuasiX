//! Performance benchmarks for S5-3 convergence monitoring SIMD evaluation
//!
//! This benchmark evaluates whether SIMD optimization is beneficial for:
//! - RMS computation
//! - Max change computation
//! - Combined metrics computation
//!
//! Target: <100 us per iteration update
//! Decision threshold: SIMD only if >2x speedup

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

// ============================================================================
// SCALAR IMPLEMENTATIONS (Baseline)
// ============================================================================

/// Scalar RMS computation (current approach from convergence.rs)
fn compute_rms_scalar(old: &[f64], new: &[f64]) -> f64 {
    debug_assert_eq!(old.len(), new.len());
    let n = old.len();
    let sum_sq: f64 = old
        .iter()
        .zip(new.iter())
        .map(|(o, n)| {
            let diff = n - o;
            diff * diff
        })
        .sum();
    (sum_sq / n as f64).sqrt()
}

/// Scalar max change computation (current approach from convergence.rs)
fn compute_max_change_scalar(old: &[f64], new: &[f64]) -> f64 {
    old.iter()
        .zip(new.iter())
        .map(|(o, n)| (n - o).abs())
        .fold(0.0, f64::max)
}

/// Scalar combined metrics (RMS + max in single pass)
fn compute_metrics_scalar(old: &[f64], new: &[f64]) -> (f64, f64) {
    let n = old.len();
    let mut sum_sq: f64 = 0.0;
    let mut max_change: f64 = 0.0;

    for (o, n) in old.iter().zip(new.iter()) {
        let diff = (n - o).abs();
        sum_sq += diff * diff;
        max_change = max_change.max(diff);
    }

    ((sum_sq / n as f64).sqrt(), max_change)
}

/// Scalar with 4x loop unrolling (better ILP)
fn compute_metrics_unrolled(old: &[f64], new: &[f64]) -> (f64, f64) {
    debug_assert_eq!(old.len(), new.len());
    let n = old.len();
    let chunks = n / 4;

    let mut sum_sq: f64 = 0.0;
    let mut max_change: f64 = 0.0;

    // Main unrolled loop
    for i in 0..chunks {
        let base = i * 4;
        let d0 = (new[base] - old[base]).abs();
        let d1 = (new[base + 1] - old[base + 1]).abs();
        let d2 = (new[base + 2] - old[base + 2]).abs();
        let d3 = (new[base + 3] - old[base + 3]).abs();

        sum_sq += d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3;
        max_change = max_change.max(d0).max(d1).max(d2).max(d3);
    }

    // Handle remainder
    for i in (chunks * 4)..n {
        let diff = (new[i] - old[i]).abs();
        sum_sq += diff * diff;
        max_change = max_change.max(diff);
    }

    ((sum_sq / n as f64).sqrt(), max_change)
}

// ============================================================================
// SIMD IMPLEMENTATIONS (AVX2)
// ============================================================================

#[cfg(target_arch = "x86_64")]
mod simd_avx2 {
    use std::arch::x86_64::*;

    /// AVX2 RMS computation (4 doubles per iteration)
    #[target_feature(enable = "avx2", enable = "fma")]
    pub unsafe fn compute_rms_avx2(old: &[f64], new: &[f64]) -> f64 {
        debug_assert_eq!(old.len(), new.len());
        let n = old.len();
        let chunks = n / 4;

        let mut sum_vec = _mm256_setzero_pd();

        for i in 0..chunks {
            let idx = i * 4;
            let old_vec = _mm256_loadu_pd(old.as_ptr().add(idx));
            let new_vec = _mm256_loadu_pd(new.as_ptr().add(idx));

            // diff = new - old
            let diff = _mm256_sub_pd(new_vec, old_vec);

            // sum += diff * diff (using FMA: sum = diff * diff + sum)
            sum_vec = _mm256_fmadd_pd(diff, diff, sum_vec);
        }

        // Horizontal sum: [a, b, c, d] -> a + b + c + d
        let sum_high = _mm256_extractf128_pd(sum_vec, 1);
        let sum_low = _mm256_castpd256_pd128(sum_vec);
        let sum_128 = _mm_add_pd(sum_low, sum_high);
        let sum_64 = _mm_add_sd(sum_128, _mm_unpackhi_pd(sum_128, sum_128));
        let sum_scalar = _mm_cvtsd_f64(sum_64);

        // Handle remainder
        let mut remainder_sum = 0.0;
        for i in (chunks * 4)..n {
            let diff = new[i] - old[i];
            remainder_sum += diff * diff;
        }

        ((sum_scalar + remainder_sum) / n as f64).sqrt()
    }

    /// AVX2 max change computation
    #[target_feature(enable = "avx2")]
    pub unsafe fn compute_max_change_avx2(old: &[f64], new: &[f64]) -> f64 {
        let n = old.len();
        let chunks = n / 4;

        // Mask for absolute value (clear sign bit)
        let abs_mask = _mm256_set1_pd(f64::from_bits(0x7FFF_FFFF_FFFF_FFFF));
        let mut max_vec = _mm256_setzero_pd();

        for i in 0..chunks {
            let idx = i * 4;
            let old_vec = _mm256_loadu_pd(old.as_ptr().add(idx));
            let new_vec = _mm256_loadu_pd(new.as_ptr().add(idx));

            // diff = new - old
            let diff = _mm256_sub_pd(new_vec, old_vec);

            // abs_diff = |diff|
            let abs_diff = _mm256_and_pd(diff, abs_mask);

            // Update max
            max_vec = _mm256_max_pd(max_vec, abs_diff);
        }

        // Horizontal max
        let max_high = _mm256_extractf128_pd(max_vec, 1);
        let max_low = _mm256_castpd256_pd128(max_vec);
        let max_128 = _mm_max_pd(max_low, max_high);
        let max_64 = _mm_max_sd(max_128, _mm_unpackhi_pd(max_128, max_128));
        let max_scalar = _mm_cvtsd_f64(max_64);

        // Handle remainder
        let mut remainder_max = max_scalar;
        for i in (chunks * 4)..n {
            remainder_max = remainder_max.max((new[i] - old[i]).abs());
        }

        remainder_max
    }

    /// AVX2 combined metrics (RMS + max in single pass)
    #[target_feature(enable = "avx2", enable = "fma")]
    pub unsafe fn compute_metrics_avx2(old: &[f64], new: &[f64]) -> (f64, f64) {
        let n = old.len();
        let chunks = n / 4;

        let abs_mask = _mm256_set1_pd(f64::from_bits(0x7FFF_FFFF_FFFF_FFFF));
        let mut sum_vec = _mm256_setzero_pd();
        let mut max_vec = _mm256_setzero_pd();

        for i in 0..chunks {
            let idx = i * 4;
            let old_vec = _mm256_loadu_pd(old.as_ptr().add(idx));
            let new_vec = _mm256_loadu_pd(new.as_ptr().add(idx));

            let diff = _mm256_sub_pd(new_vec, old_vec);
            let abs_diff = _mm256_and_pd(diff, abs_mask);

            // For RMS: sum += diff * diff
            sum_vec = _mm256_fmadd_pd(diff, diff, sum_vec);

            // For max: update max
            max_vec = _mm256_max_pd(max_vec, abs_diff);
        }

        // Reduce sum
        let sum_high = _mm256_extractf128_pd(sum_vec, 1);
        let sum_low = _mm256_castpd256_pd128(sum_vec);
        let sum_128 = _mm_add_pd(sum_low, sum_high);
        let sum_64 = _mm_add_sd(sum_128, _mm_unpackhi_pd(sum_128, sum_128));
        let sum_scalar = _mm_cvtsd_f64(sum_64);

        // Reduce max
        let max_high = _mm256_extractf128_pd(max_vec, 1);
        let max_low = _mm256_castpd256_pd128(max_vec);
        let max_128 = _mm_max_pd(max_low, max_high);
        let max_64 = _mm_max_sd(max_128, _mm_unpackhi_pd(max_128, max_128));
        let max_scalar = _mm_cvtsd_f64(max_64);

        // Handle remainder
        let mut remainder_sum = 0.0;
        let mut remainder_max = max_scalar;
        for i in (chunks * 4)..n {
            let diff = (new[i] - old[i]).abs();
            remainder_sum += diff * diff;
            remainder_max = remainder_max.max(diff);
        }

        (
            ((sum_scalar + remainder_sum) / n as f64).sqrt(),
            remainder_max,
        )
    }
}

// ============================================================================
// SIMD IMPLEMENTATIONS (AVX-512)
// ============================================================================

#[cfg(target_arch = "x86_64")]
mod simd_avx512 {
    use std::arch::x86_64::*;

    /// AVX-512 RMS computation (8 doubles per iteration)
    #[target_feature(enable = "avx512f")]
    pub unsafe fn compute_rms_avx512(old: &[f64], new: &[f64]) -> f64 {
        let n = old.len();
        let chunks = n / 8;

        let mut sum_vec = _mm512_setzero_pd();

        for i in 0..chunks {
            let idx = i * 8;
            let old_vec = _mm512_loadu_pd(old.as_ptr().add(idx));
            let new_vec = _mm512_loadu_pd(new.as_ptr().add(idx));

            let diff = _mm512_sub_pd(new_vec, old_vec);
            sum_vec = _mm512_fmadd_pd(diff, diff, sum_vec);
        }

        // Reduce to scalar
        let sum_scalar = _mm512_reduce_add_pd(sum_vec);

        // Handle remainder
        let mut remainder_sum = 0.0;
        for i in (chunks * 8)..n {
            let diff = new[i] - old[i];
            remainder_sum += diff * diff;
        }

        ((sum_scalar + remainder_sum) / n as f64).sqrt()
    }

    /// AVX-512 max change computation
    #[target_feature(enable = "avx512f")]
    pub unsafe fn compute_max_change_avx512(old: &[f64], new: &[f64]) -> f64 {
        let n = old.len();
        let chunks = n / 8;

        let mut max_vec = _mm512_setzero_pd();

        for i in 0..chunks {
            let idx = i * 8;
            let old_vec = _mm512_loadu_pd(old.as_ptr().add(idx));
            let new_vec = _mm512_loadu_pd(new.as_ptr().add(idx));

            let diff = _mm512_sub_pd(new_vec, old_vec);
            let abs_diff = _mm512_abs_pd(diff);

            max_vec = _mm512_max_pd(max_vec, abs_diff);
        }

        // Reduce to scalar
        let max_scalar = _mm512_reduce_max_pd(max_vec);

        // Handle remainder
        let mut remainder_max = max_scalar;
        for i in (chunks * 8)..n {
            remainder_max = remainder_max.max((new[i] - old[i]).abs());
        }

        remainder_max
    }

    /// AVX-512 combined metrics
    #[target_feature(enable = "avx512f")]
    pub unsafe fn compute_metrics_avx512(old: &[f64], new: &[f64]) -> (f64, f64) {
        let n = old.len();
        let chunks = n / 8;

        let mut sum_vec = _mm512_setzero_pd();
        let mut max_vec = _mm512_setzero_pd();

        for i in 0..chunks {
            let idx = i * 8;
            let old_vec = _mm512_loadu_pd(old.as_ptr().add(idx));
            let new_vec = _mm512_loadu_pd(new.as_ptr().add(idx));

            let diff = _mm512_sub_pd(new_vec, old_vec);
            let abs_diff = _mm512_abs_pd(diff);

            sum_vec = _mm512_fmadd_pd(diff, diff, sum_vec);
            max_vec = _mm512_max_pd(max_vec, abs_diff);
        }

        let sum_scalar = _mm512_reduce_add_pd(sum_vec);
        let max_scalar = _mm512_reduce_max_pd(max_vec);

        // Handle remainder
        let mut remainder_sum = 0.0;
        let mut remainder_max = max_scalar;
        for i in (chunks * 8)..n {
            let diff = (new[i] - old[i]).abs();
            remainder_sum += diff * diff;
            remainder_max = remainder_max.max(diff);
        }

        (
            ((sum_scalar + remainder_sum) / n as f64).sqrt(),
            remainder_max,
        )
    }
}


// ============================================================================
// BENCHMARKS
// ============================================================================

fn benchmark_rms(c: &mut Criterion) {
    let mut group = c.benchmark_group("rms_computation");

    for n in [100, 500, 1000, 5000].iter() {
        // Create test data
        let old: Vec<f64> = (0..*n).map(|i| -10.0 + (i as f64) * 0.001).collect();
        let mut new = old.clone();
        new[*n / 2] += 0.01; // One change

        group.bench_with_input(BenchmarkId::new("scalar", n), n, |b, _| {
            b.iter(|| compute_rms_scalar(black_box(&old), black_box(&new)))
        });

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            group.bench_with_input(BenchmarkId::new("avx2", n), n, |b, _| {
                b.iter(|| unsafe {
                    simd_avx2::compute_rms_avx2(black_box(&old), black_box(&new))
                })
            });
        }

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx512f") {
            group.bench_with_input(BenchmarkId::new("avx512", n), n, |b, _| {
                b.iter(|| unsafe {
                    simd_avx512::compute_rms_avx512(black_box(&old), black_box(&new))
                })
            });
        }
    }

    group.finish();
}

fn benchmark_max_change(c: &mut Criterion) {
    let mut group = c.benchmark_group("max_change_computation");

    for n in [100, 500, 1000, 5000].iter() {
        let old: Vec<f64> = (0..*n).map(|i| -10.0 + (i as f64) * 0.001).collect();
        let mut new = old.clone();
        new[*n / 2] += 0.01;

        group.bench_with_input(BenchmarkId::new("scalar", n), n, |b, _| {
            b.iter(|| compute_max_change_scalar(black_box(&old), black_box(&new)))
        });

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            group.bench_with_input(BenchmarkId::new("avx2", n), n, |b, _| {
                b.iter(|| unsafe {
                    simd_avx2::compute_max_change_avx2(black_box(&old), black_box(&new))
                })
            });
        }

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx512f") {
            group.bench_with_input(BenchmarkId::new("avx512", n), n, |b, _| {
                b.iter(|| unsafe {
                    simd_avx512::compute_max_change_avx512(black_box(&old), black_box(&new))
                })
            });
        }
    }

    group.finish();
}

fn benchmark_combined_metrics(c: &mut Criterion) {
    let mut group = c.benchmark_group("combined_metrics");

    for n in [100, 500, 1000, 5000].iter() {
        let old: Vec<f64> = (0..*n).map(|i| -10.0 + (i as f64) * 0.001).collect();
        let mut new = old.clone();
        new[*n / 2] += 0.01;

        group.bench_with_input(BenchmarkId::new("scalar", n), n, |b, _| {
            b.iter(|| compute_metrics_scalar(black_box(&old), black_box(&new)))
        });

        group.bench_with_input(BenchmarkId::new("unrolled_4x", n), n, |b, _| {
            b.iter(|| compute_metrics_unrolled(black_box(&old), black_box(&new)))
        });

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            group.bench_with_input(BenchmarkId::new("avx2", n), n, |b, _| {
                b.iter(|| unsafe {
                    simd_avx2::compute_metrics_avx2(black_box(&old), black_box(&new))
                })
            });
        }

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx512f") {
            group.bench_with_input(BenchmarkId::new("avx512", n), n, |b, _| {
                b.iter(|| unsafe {
                    simd_avx512::compute_metrics_avx512(black_box(&old), black_box(&new))
                })
            });
        }
    }

    group.finish();
}

/// Benchmark full convergence check latency (simulating check_convergence())
fn benchmark_convergence_check_latency(c: &mut Criterion) {
    let mut group = c.benchmark_group("convergence_check_latency");
    group.significance_level(0.01).sample_size(1000);

    for n in [100, 500, 1000, 5000].iter() {
        let old: Vec<f64> = (0..*n).map(|i| -10.0 + (i as f64) * 0.001).collect();
        let mut new = old.clone();
        for i in 0..*n {
            new[i] += 1e-5 * (i as f64);
        }

        // Simulate full convergence check (compute both metrics + compare with threshold)
        group.bench_with_input(
            BenchmarkId::new("full_check_scalar", n),
            n,
            |b, _| {
                let threshold = 1e-6;
                b.iter(|| {
                    let (rms, max_change) =
                        compute_metrics_scalar(black_box(&old), black_box(&new));
                    black_box(rms < threshold && max_change < threshold)
                })
            },
        );

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx512f") {
            group.bench_with_input(
                BenchmarkId::new("full_check_avx512", n),
                n,
                |b, _| {
                    let threshold = 1e-6;
                    b.iter(|| {
                        let (rms, max_change) = unsafe {
                            simd_avx512::compute_metrics_avx512(black_box(&old), black_box(&new))
                        };
                        black_box(rms < threshold && max_change < threshold)
                    })
                },
            );
        }
    }

    group.finish();
}

/// Verify correctness of SIMD implementations
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rms_correctness() {
        let old = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let new = vec![1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1];

        let scalar = compute_rms_scalar(&old, &new);

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            let avx2 = unsafe { simd_avx2::compute_rms_avx2(&old, &new) };
            assert_relative_eq!(scalar, avx2, epsilon = 1e-14);
        }

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx512f") {
            let avx512 = unsafe { simd_avx512::compute_rms_avx512(&old, &new) };
            assert_relative_eq!(scalar, avx512, epsilon = 1e-14);
        }
    }

    #[test]
    fn test_max_change_correctness() {
        let old = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let new = vec![1.1, 2.2, 3.1, 4.1, 5.5, 6.1, 7.1, 8.1, 9.1, 10.1];

        let scalar = compute_max_change_scalar(&old, &new);

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") {
            let avx2 = unsafe { simd_avx2::compute_max_change_avx2(&old, &new) };
            assert_relative_eq!(scalar, avx2, epsilon = 1e-14);
        }

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx512f") {
            let avx512 = unsafe { simd_avx512::compute_max_change_avx512(&old, &new) };
            assert_relative_eq!(scalar, avx512, epsilon = 1e-14);
        }
    }

    #[test]
    fn test_combined_metrics_correctness() {
        let old = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let new = vec![1.1, 2.2, 3.1, 4.1, 5.5, 6.1, 7.1, 8.1, 9.1, 10.1];

        let (scalar_rms, scalar_max) = compute_metrics_scalar(&old, &new);
        let (unrolled_rms, unrolled_max) = compute_metrics_unrolled(&old, &new);

        assert_relative_eq!(scalar_rms, unrolled_rms, epsilon = 1e-14);
        assert_relative_eq!(scalar_max, unrolled_max, epsilon = 1e-14);

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            let (avx2_rms, avx2_max) = unsafe { simd_avx2::compute_metrics_avx2(&old, &new) };
            assert_relative_eq!(scalar_rms, avx2_rms, epsilon = 1e-14);
            assert_relative_eq!(scalar_max, avx2_max, epsilon = 1e-14);
        }

        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx512f") {
            let (avx512_rms, avx512_max) =
                unsafe { simd_avx512::compute_metrics_avx512(&old, &new) };
            assert_relative_eq!(scalar_rms, avx512_rms, epsilon = 1e-14);
            assert_relative_eq!(scalar_max, avx512_max, epsilon = 1e-14);
        }
    }
}

criterion_group!(
    benches,
    benchmark_rms,
    benchmark_max_change,
    benchmark_combined_metrics,
    benchmark_convergence_check_latency,
);

criterion_main!(benches);
