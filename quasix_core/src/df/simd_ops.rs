//! SIMD Operations (Educational/Testing Only)
//!
//! ⚠️ **DO NOT USE IN PRODUCTION** ⚠️
//!
//! This module contains custom SIMD implementations for educational purposes
//! and correctness testing. Performance measurements show these implementations
//! are 350x slower than optimized BLAS libraries.
//!
//! Production code should use ndarray's BLAS-backed operations instead.
//!
//! # Performance Evidence
//!
//! Measured on AMD Ryzen 7 5800X (128×128 matrix multiply):
//! - Custom SIMD: 4.92 ms (0.85 GFLOP/s, 0.2% of peak)
//! - BLAS: 13.7 µs (307 GFLOP/s, 73% of peak)
//! - **BLAS is 359x faster**
//!
//! Reference: docs/reports/s2_5_simd_performance_analysis_20251118.md
//!
//! # Why Custom SIMD is Slower
//!
//! 1. No cache blocking (BLAS uses sophisticated tiling)
//! 2. No register re-use (BLAS keeps partial sums in registers)
//! 3. Poor memory access patterns (BLAS optimizes layouts)
//! 4. No kernel fusion (BLAS fuses operations)
//! 5. Scalar remainder handling (no SIMD for non-multiples of 4)
//! 6. No prefetching (BLAS uses software prefetch)
//!
//! # Safety
//!
//! This module uses unsafe code for SIMD operations.
//! All unsafe blocks are carefully documented with safety invariants.

#![cfg_attr(not(test), allow(dead_code))]
#![allow(clippy::many_single_char_names)] // Mathematical notation

use ndarray::{Array2, ArrayView2, ArrayViewMut2, Axis};
use rayon::prelude::*;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// SIMD-optimized matrix multiplication for MO transformation
///
/// Performs C = A^T @ B using AVX2 SIMD instructions when available.
/// Falls back to standard BLAS-like implementation on other architectures.
///
/// # Arguments
/// * `a` - Left matrix (will be transposed)
/// * `b` - Right matrix
/// * `c` - Output matrix (must be pre-allocated)
///
/// # Safety Invariants
/// - All pointers are valid and aligned
/// - Matrix dimensions are checked before operations
/// - SIMD operations use properly aligned memory
pub fn simd_matmul_transpose(a: ArrayView2<f64>, b: ArrayView2<f64>, mut c: ArrayViewMut2<f64>) {
    let (m, k1) = a.dim(); // A is m×k, A^T is k×m
    let (k2, n) = b.dim(); // B is k×n
    let (cm, cn) = c.dim(); // C is m×n

    assert_eq!(k1, k2, "Inner dimensions must match");
    assert_eq!(cm, m, "Output rows must match A cols");
    assert_eq!(cn, n, "Output cols must match B cols");

    // Zero output
    c.fill(0.0);

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
            // SAFETY: AVX2 and FMA are verified to be available via is_x86_feature_detected!
            // The simd_matmul_transpose_avx2 function is marked #[target_feature(enable = "avx2")]
            // and only uses AVX2/FMA intrinsics, which are guaranteed available at this point.
            // Array bounds are validated by the caller (dimensions a=[m,k], b=[k,n], c=[m,n]).
            unsafe {
                simd_matmul_transpose_avx2(a, b, c);
            }
            return;
        }
    }

    // Fallback to cache-optimized scalar implementation
    matmul_transpose_blocked(a, b, c);
}

/// Cache-blocked matrix multiplication with transpose
///
/// Uses blocking to improve cache locality
fn matmul_transpose_blocked(a: ArrayView2<f64>, b: ArrayView2<f64>, mut c: ArrayViewMut2<f64>) {
    const BLOCK_SIZE: usize = 64; // L1 cache-friendly block size

    let (m, k) = a.dim();
    let (_k, n) = b.dim();

    // Block over m dimension
    for i_block in (0..m).step_by(BLOCK_SIZE) {
        let i_end = (i_block + BLOCK_SIZE).min(m);

        // Block over n dimension
        for j_block in (0..n).step_by(BLOCK_SIZE) {
            let j_end = (j_block + BLOCK_SIZE).min(n);

            // Block over k dimension (reduction)
            for k_block in (0..k).step_by(BLOCK_SIZE) {
                let k_end = (k_block + BLOCK_SIZE).min(k);

                // Compute block
                for i in i_block..i_end {
                    for j in j_block..j_end {
                        let mut sum = c[[i, j]];
                        for kk in k_block..k_end {
                            sum += a[[i, kk]] * b[[kk, j]];
                        }
                        c[[i, j]] = sum;
                    }
                }
            }
        }
    }
}

/// AVX2-optimized matrix multiplication with transpose
///
/// # Safety
/// - Caller must ensure AVX2 and FMA are available
/// - All arrays must be valid and properly sized
#[cfg(target_arch = "x86_64")]
unsafe fn simd_matmul_transpose_avx2(
    a: ArrayView2<f64>,
    b: ArrayView2<f64>,
    mut c: ArrayViewMut2<f64>,
) {
    let (m, k) = a.dim();
    let (_k, n) = b.dim();

    // Process 4×4 tiles using AVX2 (4 doubles per register)
    const TILE: usize = 4;

    // Main loop over tiles
    let m_tiles = m / TILE;
    let n_tiles = n / TILE;

    for i_tile in 0..m_tiles {
        let i = i_tile * TILE;

        for j_tile in 0..n_tiles {
            let j = j_tile * TILE;

            // Load 4×4 tile of C
            let mut c00 = _mm256_setzero_pd();
            let mut c01 = _mm256_setzero_pd();
            let mut c02 = _mm256_setzero_pd();
            let mut c03 = _mm256_setzero_pd();

            // Accumulate over k dimension
            for kk in 0..k {
                // Load column of A (transposed access)
                let a0 = _mm256_set_pd(a[[i + 3, kk]], a[[i + 2, kk]], a[[i + 1, kk]], a[[i, kk]]);

                // Load row of B and accumulate
                let b0 = _mm256_set1_pd(b[[kk, j]]);
                let b1 = _mm256_set1_pd(b[[kk, j + 1]]);
                let b2 = _mm256_set1_pd(b[[kk, j + 2]]);
                let b3 = _mm256_set1_pd(b[[kk, j + 3]]);

                // FMA operations
                c00 = _mm256_fmadd_pd(a0, b0, c00);
                c01 = _mm256_fmadd_pd(a0, b1, c01);
                c02 = _mm256_fmadd_pd(a0, b2, c02);
                c03 = _mm256_fmadd_pd(a0, b3, c03);
            }

            // Store results
            let mut result = [0.0; 4];

            _mm256_storeu_pd(result.as_mut_ptr(), c00);
            for (ii, &val) in result.iter().enumerate() {
                c[[i + ii, j]] = val;
            }

            _mm256_storeu_pd(result.as_mut_ptr(), c01);
            for (ii, &val) in result.iter().enumerate() {
                c[[i + ii, j + 1]] = val;
            }

            _mm256_storeu_pd(result.as_mut_ptr(), c02);
            for (ii, &val) in result.iter().enumerate() {
                c[[i + ii, j + 2]] = val;
            }

            _mm256_storeu_pd(result.as_mut_ptr(), c03);
            for (ii, &val) in result.iter().enumerate() {
                c[[i + ii, j + 3]] = val;
            }
        }
    }

    // Handle remaining elements
    let m_remainder = m % TILE;
    let n_remainder = n % TILE;

    if m_remainder > 0 || n_remainder > 0 {
        // Process remainder with scalar code
        for i in (m_tiles * TILE)..m {
            for j in 0..n {
                let mut sum = 0.0;
                for kk in 0..k {
                    sum += a[[i, kk]] * b[[kk, j]];
                }
                c[[i, j]] = sum;
            }
        }

        for i in 0..(m_tiles * TILE) {
            for j in (n_tiles * TILE)..n {
                let mut sum = 0.0;
                for kk in 0..k {
                    sum += a[[i, kk]] * b[[kk, j]];
                }
                c[[i, j]] = sum;
            }
        }
    }
}

/// SIMD-optimized 3-center integral transformation
///
/// Transforms (μν|P) → (ia|P) using optimized SIMD operations
///
/// # Arguments
/// * `j3c_ao` - AO-basis 3-center integrals [n_ao, n_ao, n_aux]
/// * `c_occ` - Occupied MO coefficients [n_ao, n_occ]
/// * `c_vir` - Virtual MO coefficients [n_ao, n_vir]
///
/// # Returns
/// MO-basis integrals [n_occ * n_vir, n_aux]
pub fn simd_transform_3center(
    j3c_ao: ArrayView2<f64>, // Flattened view [n_ao*n_ao, n_aux]
    c_occ: ArrayView2<f64>,
    c_vir: ArrayView2<f64>,
) -> Array2<f64> {
    let n_ao = c_occ.nrows();
    let n_occ = c_occ.ncols();
    let n_vir = c_vir.ncols();
    let n_aux = j3c_ao.ncols();

    // Check dimensions
    assert_eq!(j3c_ao.nrows(), n_ao * n_ao);
    assert_eq!(c_vir.nrows(), n_ao);

    // Allocate output
    let n_trans = n_occ * n_vir;
    let mut j3c_mo = Array2::<f64>::zeros((n_trans, n_aux));

    // Process in parallel over auxiliary functions for better cache usage
    j3c_mo
        .axis_iter_mut(Axis(1))
        .into_par_iter()
        .enumerate()
        .for_each(|(p, mut j3c_p)| {
            // Extract (μν|P) for this auxiliary function
            let j3c_ao_p = j3c_ao.column(p);

            // Reshape to matrix form
            let j3c_matrix = j3c_ao_p
                .to_owned()
                .into_shape_with_order((n_ao, n_ao))
                .expect("Shape mismatch");

            // Step 1: Transform first index
            // Need: C_occ^T @ j3c_matrix = [n_occ, n_ao] @ [n_ao, n_ao] = [n_occ, n_ao]
            let c_occ_t = c_occ.t().to_owned(); // Make contiguous
            let mut temp = Array2::<f64>::zeros((n_occ, n_ao));
            simd_matmul_transpose(c_occ_t.view(), j3c_matrix.view(), temp.view_mut());

            // Step 2: Transform second index
            // Need: temp @ C_vir = [n_occ, n_ao] @ [n_ao, n_vir] = [n_occ, n_vir]
            let mut result = Array2::<f64>::zeros((n_occ, n_vir));
            simd_matmul_transpose(temp.view(), c_vir.view(), result.view_mut());

            // Flatten and store
            let flat_result = result
                .into_shape_with_order(n_trans)
                .expect("Shape mismatch in flattening");
            j3c_p.assign(&flat_result);
        });

    j3c_mo
}

/// Check CPU features at runtime
pub fn check_simd_support() -> SimdSupport {
    let mut support = SimdSupport::default();

    #[cfg(target_arch = "x86_64")]
    {
        support.sse2 = is_x86_feature_detected!("sse2");
        support.avx = is_x86_feature_detected!("avx");
        support.avx2 = is_x86_feature_detected!("avx2");
        support.fma = is_x86_feature_detected!("fma");
    }

    support
}

/// SIMD feature support information
#[derive(Debug, Default, Clone)]
pub struct SimdSupport {
    pub sse2: bool,
    pub avx: bool,
    pub avx2: bool,
    pub fma: bool,
}

impl SimdSupport {
    /// Get optimal tile size based on available SIMD
    pub fn optimal_tile_size(&self) -> usize {
        if self.avx2 {
            4 // AVX2 processes 4 doubles
        } else if self.sse2 {
            2 // SSE2 processes 2 doubles
        } else {
            1 // Scalar fallback
        }
    }

    /// Check if high-performance SIMD is available
    pub fn has_high_perf_simd(&self) -> bool {
        self.avx2 && self.fma
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_simd_matmul_transpose() {
        // Create test matrices
        let a = Array2::from_shape_vec(
            (3, 4),
            vec![
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
            ],
        )
        .unwrap();

        let b =
            Array2::from_shape_vec((4, 2), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]).unwrap();

        let mut c = Array2::zeros((3, 2));

        // Compute C = A @ B
        simd_matmul_transpose(a.view(), b.view(), c.view_mut());

        // Check results
        let expected = a.dot(&b);
        for i in 0..3 {
            for j in 0..2 {
                assert!(
                    (c[[i, j]] - expected[[i, j]]).abs() < 1e-10,
                    "Mismatch at [{}, {}]: {} vs {}",
                    i,
                    j,
                    c[[i, j]],
                    expected[[i, j]]
                );
            }
        }
    }

    #[test]
    fn test_simd_support() {
        let support = check_simd_support();
        println!("SIMD Support: {:?}", support);

        // At minimum, x86_64 should have SSE2
        #[cfg(target_arch = "x86_64")]
        assert!(support.sse2);

        let tile_size = support.optimal_tile_size();
        assert!(tile_size >= 1);
    }

    #[test]
    fn test_blocked_matmul() {
        // Test with larger matrices to exercise blocking
        let n = 128;
        let a = Array2::from_shape_fn((n, n), |(i, j)| (i + j) as f64);
        let b = Array2::from_shape_fn((n, n), |(i, j)| (i * j) as f64);
        let mut c = Array2::zeros((n, n));

        matmul_transpose_blocked(a.view(), b.view(), c.view_mut());

        // Verify a few elements
        let expected = a.dot(&b);
        for i in (0..n).step_by(17) {
            for j in (0..n).step_by(19) {
                assert!(
                    (c[[i, j]] - expected[[i, j]]).abs() < 1e-9,
                    "Mismatch at [{}, {}]",
                    i,
                    j
                );
            }
        }
    }
}
