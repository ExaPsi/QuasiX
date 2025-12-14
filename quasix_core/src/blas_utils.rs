//! BLAS utilities for high-performance linear algebra operations
//!
//! This module provides optimized BLAS wrappers for performance-critical operations
//! in QuasiX. It integrates with the existing ndarray-linalg BLAS backend.
//!
//! # Performance Considerations
//!
//! - Leverages ndarray's BLAS integration for maximum performance
//! - Properly handles thread-safety with Rayon parallelization
//! - Optimizes memory layout for cache efficiency
//! - Batches small operations into larger BLAS3 calls

use crate::common::Result;
use ndarray::{s, Array1, Array2, ArrayView2, ArrayView3};
use std::sync::atomic::{AtomicUsize, Ordering};

// Global thread count management for OpenBLAS
static OPENBLAS_THREADS: AtomicUsize = AtomicUsize::new(0);

/// Set the number of OpenBLAS threads
///
/// This is critical when using Rayon parallelization to avoid thread oversubscription
pub fn set_num_threads(n: usize) {
    // Use environment variable method which is portable
    std::env::set_var("OPENBLAS_NUM_THREADS", n.to_string());
    std::env::set_var("MKL_NUM_THREADS", n.to_string());
    std::env::set_var("OMP_NUM_THREADS", n.to_string());
    OPENBLAS_THREADS.store(n, Ordering::SeqCst);
}

/// Get the current number of OpenBLAS threads
pub fn get_num_threads() -> usize {
    OPENBLAS_THREADS.load(Ordering::SeqCst).max(1)
}

/// Optimized exchange matrix computation using blocked operations
///
/// This function implements the core computational kernel for exchange self-energy
/// with cache-optimized blocking and BLAS operations through ndarray.
///
/// Formula: Σ_x[m,n] = -Σ_i (mi|P) v^{-1}_{PQ} (ni|Q)
/// where i runs over occupied orbitals
pub fn compute_exchange_matrix_optimized(
    df_occ: &ArrayView3<f64>,     // (n_mo, n_occ, n_aux)
    metric_inv: &ArrayView2<f64>, // (n_aux, n_aux)
    nbasis: usize,
    nocc: usize,
    naux: usize,
) -> Result<Array2<f64>> {
    // Initialize result matrix
    let mut sigma_x = Array2::<f64>::zeros((nbasis, nbasis));

    // Optimal block size for L3 cache (typically 8-32 MB)
    // Assuming 8 bytes per f64 and targeting 2MB working set
    let block_size = ((2 * 1024 * 1024) / (8 * naux)).clamp(16, 64);

    // CRITICAL FIX: Work with 3D tensor to preserve (m,i) indexing
    // Pre-compute (mi|P) @ v^{-1} for all m,i pairs
    use ndarray::Array3;
    let mut df_with_metric = Array3::<f64>::zeros((nbasis, nocc, naux));

    // Apply metric inverse to each (mi|P) vector
    for m in 0..nbasis {
        for i in 0..nocc {
            let mi_p = df_occ.slice(s![m, i, ..]);
            let mi_p_vinv = mi_p.dot(metric_inv);
            df_with_metric.slice_mut(s![m, i, ..]).assign(&mi_p_vinv);
        }
    }

    // Process in blocks for cache efficiency
    use rayon::prelude::*;

    // Set BLAS threads to 1 when using Rayon
    let _guard = BlasThreadGuard::new();

    // Prepare blocks for processing
    let blocks: Vec<(usize, usize, usize, usize)> = (0..nbasis)
        .step_by(block_size)
        .flat_map(|m_start| {
            let m_end = (m_start + block_size).min(nbasis);
            (0..=m_start)
                .step_by(block_size)
                .map(move |n_start| {
                    let n_end = (n_start + block_size).min(nbasis);
                    (m_start, m_end, n_start, n_end)
                })
                .collect::<Vec<_>>()
        })
        .collect();

    // Process blocks in parallel
    let block_results: Vec<_> = blocks
        .par_iter()
        .map(|&(m_start, m_end, n_start, n_end)| {
            let mut block = Array2::<f64>::zeros((m_end - m_start, n_end - n_start));

            // Compute exchange matrix elements for this block
            for (idx_m, m) in (m_start..m_end).enumerate() {
                for (idx_n, n) in (n_start..n_end.min(m + 1)).enumerate() {
                    let mut sigma_mn = 0.0;

                    // Vectorized sum over occupied orbitals using ndarray's dot
                    // Formula: Σ_x[m,n] = -Σ_i (mi|P) v^{-1}_{PQ} (ni|Q)
                    for i in 0..nocc {
                        // Get pre-computed (mi|P) v^{-1}_{PQ}
                        let mi_p_vinv = df_with_metric.slice(s![m, i, ..]);
                        // Get (ni|Q)
                        let ni_q = df_occ.slice(s![n, i, ..]);

                        // Use ndarray's optimized dot product (calls BLAS internally)
                        sigma_mn -= mi_p_vinv.dot(&ni_q);
                    }

                    block[[idx_m, idx_n]] = sigma_mn;
                }
            }

            (m_start, n_start, block)
        })
        .collect();

    // Assemble results
    for (m_start, n_start, block) in block_results {
        let m_end = (m_start + block_size).min(nbasis);
        let n_end = (n_start + block_size).min(nbasis);

        for (idx_m, m) in (m_start..m_end).enumerate() {
            for (idx_n, n) in (n_start..n_end.min(m + 1)).enumerate() {
                let val = block[[idx_m, idx_n]];
                sigma_x[[m, n]] = val;
                if m != n {
                    sigma_x[[n, m]] = val; // Exploit symmetry
                }
            }
        }
    }

    Ok(sigma_x)
}

/// Optimized batched DF tensor contraction for exchange matrix
///
/// Processes multiple occupied orbital blocks efficiently
pub fn df_exchange_contract_batch(
    df_occ: &ArrayView3<f64>,     // (n_mo, n_occ, n_aux)
    metric_inv: &ArrayView2<f64>, // (n_aux, n_aux)
    m_start: usize,
    m_end: usize,
    n_start: usize,
    n_end: usize,
) -> Result<Array2<f64>> {
    let n_occ = df_occ.dim().1;
    let n_aux = df_occ.dim().2;

    let m_size = m_end - m_start;
    let n_size = n_end - n_start;

    let mut result = Array2::<f64>::zeros((m_size, n_size));

    // Pre-allocate workspace for metric multiplication
    let mut workspace = Array2::<f64>::zeros((m_size * n_occ, n_aux));

    // Transform m-block with metric: (mi|P) @ v^(-1)
    for (idx_m, m) in (m_start..m_end).enumerate() {
        for i in 0..n_occ {
            let src_row = df_occ.slice(s![m, i, ..]);
            let dst_idx = idx_m * n_occ + i;

            // Use ndarray's dot for matrix-vector multiplication
            let temp = src_row.dot(metric_inv);
            workspace.row_mut(dst_idx).assign(&temp);
        }
    }

    // Contract with n-block
    for (idx_m, _m) in (m_start..m_end).enumerate() {
        for (idx_n, n) in (n_start..n_end).enumerate() {
            let mut sum = 0.0;

            // Sum over occupied orbitals
            for i in 0..n_occ {
                let m_row = workspace.row(idx_m * n_occ + i);
                let n_row = df_occ.slice(s![n, i, ..]);

                // Use ndarray's dot product
                let contrib = m_row.dot(&n_row);
                sum -= contrib; // No additional factor needed
            }

            result[[idx_m, idx_n]] = sum;
        }
    }

    Ok(result)
}

/// Thread-safe BLAS configuration for use with Rayon
pub struct BlasThreadGuard {
    original_threads: usize,
}

impl BlasThreadGuard {
    /// Create a new thread guard that sets BLAS threads to 1
    pub fn new() -> Self {
        let original = get_num_threads();
        set_num_threads(1);
        Self {
            original_threads: original,
        }
    }

    /// Create a thread guard with specific number of threads
    pub fn with_threads(n: usize) -> Self {
        let original = get_num_threads();
        set_num_threads(n);
        Self {
            original_threads: original,
        }
    }
}

impl Drop for BlasThreadGuard {
    fn drop(&mut self) {
        set_num_threads(self.original_threads);
    }
}

impl Default for BlasThreadGuard {
    fn default() -> Self {
        Self::new()
    }
}

/// Optimized diagonal computation using vectorized operations
pub fn compute_exchange_diagonal_optimized(
    df_occ: &ArrayView3<f64>,     // (n_mo, n_occ, n_aux)
    metric_inv: &ArrayView2<f64>, // (n_aux, n_aux)
    nbasis: usize,
    nocc: usize,
) -> Result<Array1<f64>> {
    use rayon::prelude::*;

    // Set BLAS threads to 1 when using Rayon
    let _guard = BlasThreadGuard::new();

    // Parallel computation of diagonal elements
    let diag: Vec<f64> = (0..nbasis)
        .into_par_iter()
        .map(|m| {
            let mut sigma_mm = 0.0;

            for i in 0..nocc {
                let mi_p = df_occ.slice(s![m, i, ..]);

                // Compute (mi|P) @ v^(-1) @ (mi|P)^T efficiently
                let temp = mi_p.dot(metric_inv);
                let contrib = temp.dot(&mi_p);

                sigma_mm -= contrib; // No additional factor needed
            }

            sigma_mm
        })
        .collect();

    Ok(Array1::from(diag))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array3;

    #[test]
    fn test_thread_guard() {
        let original = get_num_threads();

        {
            let _guard = BlasThreadGuard::new();
            assert_eq!(get_num_threads(), 1);
        }

        // Should be restored after guard is dropped
        assert_eq!(get_num_threads(), original);
    }

    #[test]
    fn test_compute_exchange_matrix_optimized() {
        let nbasis = 4;
        let nocc = 2;
        let naux = 6;

        // Create test DF tensor
        let mut df_occ = Array3::<f64>::zeros((nbasis, nocc, naux));
        for m in 0..nbasis {
            for i in 0..nocc {
                for p in 0..naux {
                    df_occ[[m, i, p]] = (m + i + p) as f64 * 0.1;
                }
            }
        }

        // Create identity metric inverse
        let metric_inv = Array2::<f64>::eye(naux);

        // Compute exchange matrix
        let sigma_x = compute_exchange_matrix_optimized(
            &df_occ.view(),
            &metric_inv.view(),
            nbasis,
            nocc,
            naux,
        )
        .unwrap();

        // Check symmetry
        for i in 0..nbasis {
            for j in 0..i {
                assert_relative_eq!(sigma_x[[i, j]], sigma_x[[j, i]], epsilon = 1e-10);
            }
        }

        // Check that diagonal elements are non-positive (exchange is attractive)
        for i in 0..nbasis {
            assert!(sigma_x[[i, i]] <= 0.0);
        }
    }

    #[test]
    fn test_compute_exchange_diagonal_optimized() {
        let nbasis = 3;
        let nocc = 2;
        let naux = 4;

        // Create test DF tensor
        let mut df_occ = Array3::<f64>::zeros((nbasis, nocc, naux));
        for m in 0..nbasis {
            for i in 0..nocc {
                for p in 0..naux {
                    df_occ[[m, i, p]] = (m * nocc + i) as f64 * 0.1 + p as f64 * 0.01;
                }
            }
        }

        let metric_inv = Array2::<f64>::eye(naux);

        // Compute diagonal
        let diag =
            compute_exchange_diagonal_optimized(&df_occ.view(), &metric_inv.view(), nbasis, nocc)
                .unwrap();

        // Check that all diagonal elements are computed
        assert_eq!(diag.len(), nbasis);

        // Check that diagonal elements are non-positive
        for i in 0..nbasis {
            assert!(diag[i] <= 0.0);
        }
    }

    #[test]
    fn test_df_exchange_contract_batch() {
        let nbasis = 4;
        let nocc = 2;
        let naux = 5;

        // Create test DF tensor
        let mut df_occ = Array3::<f64>::zeros((nbasis, nocc, naux));
        for m in 0..nbasis {
            for i in 0..nocc {
                for p in 0..naux {
                    df_occ[[m, i, p]] = ((m + 1) * (i + 1)) as f64 * 0.1;
                }
            }
        }

        let metric_inv = Array2::<f64>::eye(naux);

        // Test a small batch
        let result = df_exchange_contract_batch(
            &df_occ.view(),
            &metric_inv.view(),
            0,
            2, // m range
            0,
            2, // n range
        )
        .unwrap();

        // Check dimensions
        assert_eq!(result.shape(), &[2, 2]);

        // Check that values are finite
        for val in &result {
            assert!(val.is_finite());
        }
    }
}
