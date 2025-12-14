#![allow(clippy::similar_names)] // n_ao and n_mo are standard in quantum chemistry
//! MO transformation and Cholesky factorization for density fitting
//!
//! This module handles the transformation of 3-center integrals from
//! AO to MO basis and the Cholesky factorization of the metric matrix.
//!
//! Implements high-performance algorithms from S2-2 including:
//! - Two-step MO transformation with BLAS optimization
//! - Blocked Cholesky factorization with pivoting fallback
//! - Memory-efficient blocking strategies
//!
//! # Performance Notes
//!
//! All matrix operations use ndarray's BLAS-backed operations exclusively.
//! Custom SIMD implementations were tested and found to be 350x slower than
//! optimized BLAS libraries (OpenBLAS/MKL).
//!
//! Performance comparison (AMD Ryzen 7 5800X, 128×128 matrices):
//! - Custom SIMD: 4.92 ms (0.85 GFLOP/s, 0.2% of peak)
//! - BLAS (ndarray): 13.7 µs (307 GFLOP/s, 73% of peak)
//!
//! Reference: docs/reports/s2_5_simd_performance_analysis_20251118.md
#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::{QuasixError, Result};
use crate::df::metric::{CholeskyMethod, CholeskyMetric};
use crate::df::parallel::{
    transform_mo_3center_parallel_optimized, BlockingStrategy, ThreadPoolConfig,
};
use ndarray::{s, Array2, Array3};
use rayon::prelude::*;
use std::sync::Arc;
use tracing::{debug, info, instrument, warn};

// CholeskyMetric struct and methods are now imported from the metric module

/// Configuration for MO transformation
///
/// All matrix operations use ndarray's BLAS-backed operations for optimal
/// performance. Custom SIMD paths have been removed after measurements showed
/// BLAS to be 350x faster on typical hardware.
///
/// Reference: docs/reports/s2_5_simd_performance_analysis_20251118.md
#[derive(Debug, Clone)]
pub struct TransformConfig {
    /// Maximum memory to use in GB
    pub max_memory_gb: f64,
    /// Minimum chunk size for parallel processing
    pub min_chunk_size: usize,
    /// Maximum chunk size for parallel processing
    pub max_chunk_size: usize,
    /// Number of threads for BLAS operations
    pub blas_threads: usize,
    /// Use optimized blocked algorithm
    pub use_blocking: bool,
    /// Use advanced parallel optimization
    pub use_advanced_parallel: bool,
    /// Prefetch distance for cache optimization
    pub prefetch_distance: usize,
}

impl Default for TransformConfig {
    fn default() -> Self {
        Self {
            max_memory_gb: 4.0,
            min_chunk_size: 10,
            max_chunk_size: 1000,
            blas_threads: 1, // Set to 1 when using Rayon
            use_blocking: true,
            use_advanced_parallel: true,
            prefetch_distance: 8,
        }
    }
}

impl TransformConfig {
    /// Configuration for small systems (< 100 basis functions)
    pub fn small_system() -> Self {
        Self {
            max_memory_gb: 1.0,
            min_chunk_size: 5,
            max_chunk_size: 50,
            blas_threads: 4,
            use_blocking: false,
            use_advanced_parallel: false,
            prefetch_distance: 4,
        }
    }

    /// Configuration for large systems (> 1000 basis functions)
    pub fn large_system() -> Self {
        Self {
            max_memory_gb: 16.0,
            min_chunk_size: 50,
            max_chunk_size: 2000,
            blas_threads: 1,
            use_blocking: true,
            use_advanced_parallel: true,
            prefetch_distance: 16,
        }
    }

    /// Configuration for memory-constrained environments
    pub fn memory_constrained() -> Self {
        Self {
            max_memory_gb: 2.0,
            min_chunk_size: 5,
            max_chunk_size: 100,
            blas_threads: 2,
            use_blocking: true,
            use_advanced_parallel: true,
            prefetch_distance: 4,
        }
    }
}

/// Transform 3-center integrals from AO to MO basis
///
/// Performs the transformation (μν|P) → (ia|P) where:
/// - μ,ν are AO basis functions
/// - i,a are occupied and virtual MO indices
/// - P is an auxiliary basis function
///
/// The transformation is done in two steps:
/// 1. (μν|P) → (iν|P) using C_occ
/// 2. (iν|P) → (ia|P) using C_vir
///
/// # Arguments
/// * `j3c_ao` - 3-center AO integrals of shape [n_ao, n_ao, n_aux]
/// * `c_occ` - Occupied MO coefficients of shape [n_ao, n_occ]
/// * `c_vir` - Virtual MO coefficients of shape [n_ao, n_vir]
///
/// # Returns
/// MO-transformed integrals of shape [n_occ*n_vir, n_aux] for efficiency
#[instrument(skip(j3c_ao, c_occ, c_vir))]
pub fn transform_mo_3center(
    j3c_ao: &Array3<f64>,
    c_occ: &Array2<f64>,
    c_vir: &Array2<f64>,
) -> Result<Array2<f64>> {
    transform_mo_3center_with_config(j3c_ao, c_occ, c_vir, &TransformConfig::default())
}

/// Transform 3-center integrals with custom configuration
#[instrument(skip(j3c_ao, c_occ, c_vir, config))]
pub fn transform_mo_3center_with_config(
    j3c_ao: &Array3<f64>,
    c_occ: &Array2<f64>,
    c_vir: &Array2<f64>,
    config: &TransformConfig,
) -> Result<Array2<f64>> {
    let (n_ao_1, n_ao_2, n_aux) = j3c_ao.dim();
    let (n_ao_occ, n_occ) = c_occ.dim();
    let (n_ao_vir, n_vir) = c_vir.dim();

    // Validate dimensions
    if n_ao_1 != n_ao_2 {
        return Err(QuasixError::InvalidInput(format!(
            "3-center integrals must have square AO dimensions, got [{}, {}, {}]",
            n_ao_1, n_ao_2, n_aux
        )));
    }
    if n_ao_1 != n_ao_occ || n_ao_1 != n_ao_vir {
        return Err(QuasixError::InvalidInput(format!(
            "MO coefficient dimensions don't match AO basis: AO={}, C_occ={}, C_vir={}",
            n_ao_1, n_ao_occ, n_ao_vir
        )));
    }

    info!(
        "Transforming 3-center integrals: ({}, {}, {}) -> ({}, {})",
        n_ao_1,
        n_ao_1,
        n_aux,
        n_occ * n_vir,
        n_aux
    );

    // All matrix operations use BLAS-backed operations (350x faster than custom SIMD)
    // Evidence: docs/reports/s2_5_simd_performance_analysis_20251118.md
    info!("Using BLAS-optimized matrix operations (307 GFLOP/s typical)");

    // Determine optimal algorithm based on system size and configuration
    let memory_required = estimate_memory_requirement(n_ao_1, n_occ, n_vir, n_aux);
    info!(
        "Memory requirement: {:.2} GB (available: {:.2} GB)",
        memory_required, config.max_memory_gb
    );

    let j3c_ia = if config.use_advanced_parallel && n_aux > 50 {
        // Use advanced parallel optimization for medium to large systems
        let thread_config = if n_aux * n_ao_1 * n_ao_1 > 100_000_000 {
            ThreadPoolConfig::compute_bound()
        } else {
            ThreadPoolConfig::memory_bound()
        };

        let blocking = if n_aux > 1000 {
            BlockingStrategy::intel_xeon()
        } else {
            BlockingStrategy::default()
        };

        transform_mo_3center_parallel_optimized(
            j3c_ao,
            c_occ,
            c_vir,
            Some(thread_config),
            Some(blocking),
        )?
    } else if config.use_blocking && n_aux > 100 {
        // Use standard blocked parallel algorithm
        transform_mo_3center_blocked(j3c_ao, c_occ, c_vir, config)?
    } else {
        // Use optimized sequential algorithm for small systems
        transform_mo_3center_sequential_optimized(j3c_ao, c_occ, c_vir)?
    };

    // Validate output
    validate_mo_3center(&j3c_ia, n_occ, n_vir, n_aux)?;

    info!(
        "MO transformation complete: min={:.6}, max={:.6}, mean={:.6}",
        j3c_ia.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
        j3c_ia.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),
        j3c_ia.iter().sum::<f64>() / j3c_ia.len() as f64
    );

    Ok(j3c_ia)
}

/// Validate MO-transformed 3-center integrals
fn validate_mo_3center(
    j3c_ia: &Array2<f64>,
    n_occ: usize,
    n_vir: usize,
    n_aux: usize,
) -> Result<()> {
    let (n_trans, n_aux_actual) = j3c_ia.dim();

    if n_trans != n_occ * n_vir {
        return Err(QuasixError::InvalidInput(format!(
            "Wrong number of transitions: expected {}, got {}",
            n_occ * n_vir,
            n_trans
        )));
    }

    if n_aux_actual != n_aux {
        return Err(QuasixError::InvalidInput(format!(
            "Wrong number of auxiliary functions: expected {}, got {}",
            n_aux, n_aux_actual
        )));
    }

    // Check for NaN or infinite values
    for value in j3c_ia {
        if !value.is_finite() {
            return Err(QuasixError::NumericalError(
                "MO-transformed integrals contain non-finite values".to_string(),
            ));
        }
    }

    Ok(())
}

/// Build Cholesky factorization of the DF metric matrix (deprecated)
///
/// This function is deprecated. Use `compute_cholesky_v` from the metric module instead,
/// which provides more robust factorization with eigenvalue-based condition number
/// computation and automatic fallback to pivoted Cholesky.
///
/// # Arguments
/// * `metric` - 2-center metric matrix of shape [n_aux, n_aux]
///
/// # Returns
/// CholeskyMetric struct containing L matrix and condition number
#[deprecated(
    since = "0.2.0",
    note = "Use compute_cholesky_v from metric module instead"
)]
#[instrument(skip(metric))]
pub fn build_cholesky_metric(metric: &Array2<f64>) -> Result<CholeskyMetric> {
    let (n_aux_1, n_aux_2) = metric.dim();

    if n_aux_1 != n_aux_2 {
        return Err(QuasixError::InvalidInput(format!(
            "Metric matrix must be square, got shape [{}, {}]",
            n_aux_1, n_aux_2
        )));
    }

    let n_aux = n_aux_1;
    info!(
        "Computing Cholesky factorization of {}x{} metric matrix",
        n_aux, n_aux
    );

    // Check symmetry
    let tol = 1e-10;
    for i in 0..n_aux {
        for j in i + 1..n_aux {
            if (metric[[i, j]] - metric[[j, i]]).abs() > tol {
                return Err(QuasixError::NumericalError(format!(
                    "Metric matrix not symmetric at [{},{}]: {:.6} != {:.6}",
                    i,
                    j,
                    metric[[i, j]],
                    metric[[j, i]]
                )));
            }
        }
    }

    // Manual Cholesky decomposition (since we don't have ndarray-linalg)
    // Algorithm: L L^T = A
    // L_ij = (A_ij - Σ_{k<j} L_ik L_jk) / L_jj for i > j
    // L_ii = sqrt(A_ii - Σ_{k<i} L_ik^2)

    let mut l_matrix = Array2::<f64>::zeros((n_aux, n_aux));

    for i in 0..n_aux {
        for j in 0..=i {
            if i == j {
                // Diagonal element
                let sum_sq: f64 = (0..j).map(|k| l_matrix[[i, k]].powi(2)).sum();
                let diag_val = metric[[i, i]] - sum_sq;

                if diag_val <= 0.0 {
                    return Err(QuasixError::NumericalError(format!(
                        "Metric matrix not positive definite: diagonal element {} gives {:.6}",
                        i, diag_val
                    )));
                }

                l_matrix[[i, j]] = diag_val.sqrt();
            } else {
                // Off-diagonal element
                let sum_prod: f64 = (0..j).map(|k| l_matrix[[i, k]] * l_matrix[[j, k]]).sum();
                l_matrix[[i, j]] = (metric[[i, j]] - sum_prod) / l_matrix[[j, j]];
            }
        }
    }

    // Estimate condition number using diagonal elements
    // This is a rough approximation; full eigenvalue computation would be more accurate
    let diag_min = (0..n_aux)
        .map(|i| l_matrix[[i, i]])
        .fold(f64::INFINITY, f64::min);
    let diag_max = (0..n_aux)
        .map(|i| l_matrix[[i, i]])
        .fold(f64::NEG_INFINITY, f64::max);

    // For a Cholesky factor L, if λ_min and λ_max are eigenvalues of original matrix,
    // then sqrt(λ_min) and sqrt(λ_max) are roughly the min/max diagonal elements of L
    let condition_number = (diag_max / diag_min).powi(2);

    if condition_number > 1e8 {
        warn!(
            "Metric matrix may be ill-conditioned: condition number = {:.2e}",
            condition_number
        );
    } else {
        info!("Metric matrix condition number: {:.2e}", condition_number);
    }

    // Validate reconstruction
    let reconstructed = l_matrix.dot(&l_matrix.t());
    let max_error = metric
        .iter()
        .zip(reconstructed.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(f64::NEG_INFINITY, f64::max);

    if max_error > 1e-8 {
        warn!("Cholesky reconstruction error: {:.2e}", max_error);
        if max_error > 1e-6 {
            return Err(QuasixError::NumericalError(format!(
                "Cholesky factorization failed: reconstruction error {:.2e}",
                max_error
            )));
        }
    } else {
        debug!("Cholesky reconstruction error: {:.2e}", max_error);
    }

    let result = CholeskyMetric {
        l_matrix,
        condition_number,
        naux: n_aux,
        method: CholeskyMethod::Standard,
        rank: None,
        permutation: None,
        eigenvalues: None, // Not computed in deprecated function
    };

    result.validate()?;

    Ok(result)
}

/// Generate mock orthonormal MO coefficients for testing
///
/// Creates physically reasonable MO coefficients that are orthonormal
/// and suitable for testing transformations.
///
/// # Arguments
/// * `n_ao` - Number of atomic orbitals
/// * `n_mo` - Number of molecular orbitals to generate
/// * `seed` - Random seed for reproducibility
///
/// # Returns
/// Orthonormal MO coefficient matrix of shape [n_ao, n_mo]
pub fn generate_mock_mo_coefficients(n_ao: usize, n_mo: usize, seed: u64) -> Array2<f64> {
    use rand::{Rng, SeedableRng};
    use rand_distr::StandardNormal;

    assert!(
        n_mo <= n_ao,
        "Cannot have more MOs than AOs: {} > {}",
        n_mo,
        n_ao
    );

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    // Generate random matrix
    let mut coeffs = Array2::<f64>::zeros((n_ao, n_mo));
    for i in 0..n_ao {
        for j in 0..n_mo {
            coeffs[[i, j]] = rng.sample(StandardNormal);
        }
    }

    // Gram-Schmidt orthogonalization
    for j in 0..n_mo {
        // Extract column j
        let mut v = coeffs.column(j).to_owned();

        // Subtract projections onto previous columns
        for k in 0..j {
            let u_k = coeffs.column(k);
            let proj = v.dot(&u_k);
            v = v - proj * &u_k;
        }

        // Normalize
        let norm = v.dot(&v).sqrt();
        if norm > 1e-10 {
            v /= norm;
        }

        // Store back
        coeffs.column_mut(j).assign(&v);
    }

    // Verify orthonormality
    let overlap = coeffs.t().dot(&coeffs);
    for i in 0..n_mo {
        for j in 0..n_mo {
            let expected = if i == j { 1.0 } else { 0.0 };
            let actual = overlap[[i, j]];
            if (actual - expected).abs() > 1e-10 {
                warn!(
                    "MO coefficients not perfectly orthonormal at [{},{}]: {:.6}",
                    i, j, actual
                );
            }
        }
    }

    coeffs
}

// s! macro already imported above

/// Optimized sequential MO transformation using BLAS operations
///
/// Always uses ndarray's BLAS-backed operations (350x faster than custom SIMD).
/// Evidence: docs/reports/s2_5_simd_performance_analysis_20251118.md
/// BLAS achieves 307 GFLOP/s vs custom SIMD 0.85 GFLOP/s on typical hardware.
fn transform_mo_3center_sequential_optimized(
    j3c_ao: &Array3<f64>,
    c_occ: &Array2<f64>,
    c_vir: &Array2<f64>,
) -> Result<Array2<f64>> {
    let (n_ao, _, n_aux) = j3c_ao.dim();
    let n_occ = c_occ.ncols();
    let n_vir = c_vir.ncols();

    debug!(
        "Using BLAS-optimized sequential transformation (n_ao={}, n_occ={}, n_vir={}, n_aux={})",
        n_ao, n_occ, n_vir, n_aux
    );

    // Step 1: Transform first index (μν|P) → (iν|P)
    // Use BLAS matrix multiplication for each P
    let mut j3c_half = Array3::<f64>::zeros((n_occ, n_ao, n_aux));

    for p in 0..n_aux {
        // Extract (μν|P) slice
        let j3c_p = j3c_ao.slice(s![.., .., p]);

        // Always use ndarray's BLAS-backed operations (350x faster than custom SIMD)
        // Evidence: docs/reports/s2_5_simd_performance_analysis_20251118.md
        // BLAS achieves 307 GFLOP/s vs custom SIMD 0.85 GFLOP/s
        let j3c_transformed = c_occ.t().dot(&j3c_p);
        j3c_half.slice_mut(s![.., .., p]).assign(&j3c_transformed);
    }

    // Step 2: Transform second index (iν|P) → (ia|P)
    let n_trans = n_occ * n_vir;
    let mut j3c_ia = Array2::<f64>::zeros((n_trans, n_aux));

    for p in 0..n_aux {
        for i in 0..n_occ {
            // Extract (iν|P) for this occupied orbital
            let j3c_iv = j3c_half.slice(s![i, .., p]);

            // Transform: (ia|P) = (iν|P) @ C_vir
            // Uses BLAS for optimal performance
            let j3c_ia_p = j3c_iv.dot(c_vir);

            // Store in composite index format
            let start_idx = i * n_vir;
            for (a, &val) in j3c_ia_p.iter().enumerate() {
                j3c_ia[[start_idx + a, p]] = val;
            }
        }
    }

    Ok(j3c_ia)
}

/// Blocked parallel MO transformation for large systems
fn transform_mo_3center_blocked(
    j3c_ao: &Array3<f64>,
    c_occ: &Array2<f64>,
    c_vir: &Array2<f64>,
    config: &TransformConfig,
) -> Result<Array2<f64>> {
    let (n_ao, _, n_aux) = j3c_ao.dim();
    let n_occ = c_occ.ncols();
    let n_vir = c_vir.ncols();
    let n_trans = n_occ * n_vir;

    debug!("Using blocked parallel transformation");

    // Calculate optimal chunk size based on memory constraints
    let chunk_size = calculate_optimal_chunk_size(n_aux, n_ao, n_occ, n_vir, config);
    info!(
        "Using chunk size: {} for {} auxiliary functions",
        chunk_size, n_aux
    );

    // Pre-allocate output
    let j3c_ia = Arc::new(parking_lot::RwLock::new(Array2::<f64>::zeros((
        n_trans, n_aux,
    ))));

    // Process in parallel chunks
    (0..n_aux)
        .into_par_iter()
        .chunks(chunk_size)
        .enumerate()
        .for_each(|(chunk_idx, chunk)| {
            let chunk_vec: Vec<usize> = chunk.into_iter().collect();
            let chunk_start = chunk_vec[0];
            let chunk_end = *chunk_vec.last().unwrap() + 1;
            let chunk_size_actual = chunk_end - chunk_start;

            debug!(
                "Processing chunk {} [{}, {})",
                chunk_idx, chunk_start, chunk_end
            );

            // Extract chunk of j3c_ao
            let _j3c_chunk = j3c_ao.slice(s![.., .., chunk_start..chunk_end]);

            // Step 1: Transform first index for this chunk
            let mut j3c_half_chunk = Array3::<f64>::zeros((n_occ, n_ao, chunk_size_actual));

            for (p_local, &p_global) in chunk_vec.iter().enumerate() {
                let j3c_p = j3c_ao.slice(s![.., .., p_global]);
                let j3c_transformed = c_occ.t().dot(&j3c_p);
                j3c_half_chunk
                    .slice_mut(s![.., .., p_local])
                    .assign(&j3c_transformed);
            }

            // Step 2: Transform second index for this chunk
            let mut j3c_ia_chunk = Array2::<f64>::zeros((n_trans, chunk_size_actual));

            for (p_local, &_p_global) in chunk_vec.iter().enumerate() {
                for i in 0..n_occ {
                    let j3c_iv = j3c_half_chunk.slice(s![i, .., p_local]);
                    // Transform: (ia|P) = (iν|P) @ C_vir
                    let j3c_ia_p = j3c_iv.to_owned().dot(c_vir);

                    let start_idx = i * n_vir;
                    for (a, &val) in j3c_ia_p.iter().enumerate() {
                        j3c_ia_chunk[[start_idx + a, p_local]] = val;
                    }
                }
            }

            // Write chunk to output
            let mut j3c_ia_write = j3c_ia.write();
            j3c_ia_write
                .slice_mut(s![.., chunk_start..chunk_end])
                .assign(&j3c_ia_chunk);
        });

    // Extract the result
    let j3c_ia_final = Arc::try_unwrap(j3c_ia)
        .map_or_else(|arc| arc.read().clone(), parking_lot::RwLock::into_inner);

    Ok(j3c_ia_final)
}

/// Calculate optimal chunk size for blocked operations
fn calculate_optimal_chunk_size(
    n_aux: usize,
    n_ao: usize,
    n_occ: usize,
    n_vir: usize,
    config: &TransformConfig,
) -> usize {
    const BYTES_PER_F64: usize = 8;

    // Estimate memory per auxiliary function
    let mem_per_aux = BYTES_PER_F64
        * (
            n_ao * n_ao +           // Input slice
        n_occ * n_ao +          // Half-transformed
        n_occ * n_vir
            // Fully transformed
        );

    // Safe conversion: max_memory_gb is always positive in practice
    let available_bytes = (config.max_memory_gb.abs() * 1e9) as usize;
    let optimal_chunk = available_bytes / mem_per_aux;

    // Clamp to configured bounds
    optimal_chunk
        .min(config.max_chunk_size)
        .max(config.min_chunk_size)
        .min(n_aux)
}

/// Estimate memory requirement in GB
fn estimate_memory_requirement(n_ao: usize, n_occ: usize, n_vir: usize, n_aux: usize) -> f64 {
    const BYTES_PER_F64: usize = 8;

    let bytes = BYTES_PER_F64
        * (
            n_ao * n_ao * n_aux +       // Input j3c_ao
        n_occ * n_ao * n_aux +      // Half-transformed
        n_occ * n_vir * n_aux
            // Output j3c_ia
        );

    bytes as f64 / 1e9
}

use parking_lot;
use rayon::iter::IndexedParallelIterator;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrals::{BasisSet, IntegralEngine, Molecule};

    #[test]
    fn test_mo_transform_h2o() {
        // Setup H2O molecule
        let mol = Molecule::water();
        let basis = BasisSet::mock_sto3g(8);
        let aux_basis = BasisSet::mock_ri(18);

        // Create integral engine
        let engine = IntegralEngine::new(mol, basis, aux_basis);

        // Get integrals
        let j3c_ao = engine.compute_3center().unwrap();
        let (n_ao, _, n_aux) = j3c_ao.dim();

        // Generate MO coefficients (5 occupied, 3 virtual for H2O)
        let c_occ = generate_mock_mo_coefficients(n_ao, 5, 42);
        let c_vir = generate_mock_mo_coefficients(n_ao, 3, 43);

        // Transform to MO basis
        let j3c_ia = transform_mo_3center(&j3c_ao, &c_occ, &c_vir).unwrap();

        // Check dimensions
        assert_eq!(j3c_ia.dim(), (15, n_aux)); // 5*3 = 15 transitions

        // Check finite values
        for val in &j3c_ia {
            assert!(val.is_finite());
        }

        // Test with custom config
        let config = TransformConfig {
            max_memory_gb: 1.0,
            use_blocking: true,
            ..Default::default()
        };

        let j3c_ia_blocked =
            transform_mo_3center_with_config(&j3c_ao, &c_occ, &c_vir, &config).unwrap();
        assert_eq!(j3c_ia_blocked.dim(), j3c_ia.dim());
    }

    #[test]
    fn test_mo_transform_nh3() {
        // Setup NH3 molecule
        let mol = Molecule::ammonia();
        let basis = BasisSet::mock_sto3g(11);
        let aux_basis = BasisSet::mock_ri(23);

        // Create integral engine
        let engine = IntegralEngine::new(mol, basis, aux_basis);

        // Get integrals
        let j3c_ao = engine.compute_3center().unwrap();
        let (n_ao, _, n_aux) = j3c_ao.dim();

        // Generate MO coefficients (5 occupied, 6 virtual for NH3)
        let c_occ = generate_mock_mo_coefficients(n_ao, 5, 44);
        let c_vir = generate_mock_mo_coefficients(n_ao, 6, 45);

        // Transform to MO basis
        let j3c_ia = transform_mo_3center(&j3c_ao, &c_occ, &c_vir).unwrap();

        // Check dimensions
        assert_eq!(j3c_ia.dim(), (30, n_aux)); // 5*6 = 30 transitions
    }

    #[test]
    fn test_cholesky_metric_h2o() {
        // Use a known positive definite matrix for testing
        // since the mock IntegralEngine might not produce positive definite metrics
        let n = 10;
        let mut metric = Array2::<f64>::zeros((n, n));

        // Create a positive definite matrix: A = B^T B + I
        for i in 0..n {
            for j in 0..=i {
                metric[[i, j]] = if i == j {
                    2.0 + (i as f64) * 0.1 // Diagonal dominance
                } else {
                    0.1 * ((i + j) as f64).sin() // Small off-diagonal
                };
                if i != j {
                    metric[[j, i]] = metric[[i, j]]; // Symmetry
                }
            }
        }

        let cholesky = crate::df::metric::compute_cholesky_v(&metric, None).unwrap();

        assert_eq!(cholesky.naux, n);
        assert_eq!(cholesky.l_matrix.dim(), (n, n));

        // Check lower triangular
        for i in 0..n {
            for j in i + 1..n {
                assert!(cholesky.l_matrix[[i, j]].abs() < 1e-14);
            }
        }

        // Check reconstruction
        let reconstructed = cholesky.reconstruct_metric();
        let max_error = metric
            .iter()
            .zip(reconstructed.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(f64::NEG_INFINITY, f64::max);

        assert!(
            max_error < 1e-10,
            "Reconstruction error too large: {:.2e}",
            max_error
        );

        // Check condition number is reasonable
        assert!(cholesky.condition_number > 1.0);
        assert!(cholesky.condition_number < 1e10);
    }

    #[test]
    fn test_cholesky_metric_nh3() {
        // Use a known positive definite matrix for testing
        let n = 15;
        let mut metric = Array2::<f64>::eye(n);

        // Add some structure while maintaining positive definiteness
        for i in 0..n {
            for j in 0..i {
                let val = 0.05 * ((i - j) as f64).exp().recip();
                metric[[i, j]] = val;
                metric[[j, i]] = val;
            }
            metric[[i, i]] += 1.0; // Ensure diagonal dominance
        }

        let cholesky = crate::df::metric::compute_cholesky_v(&metric, None).unwrap();

        assert_eq!(cholesky.naux, n);

        // Check reconstruction
        let reconstructed = cholesky.reconstruct_metric();
        let max_error = metric
            .iter()
            .zip(reconstructed.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(f64::NEG_INFINITY, f64::max);

        assert!(max_error < 1e-10);
    }

    #[test]
    fn test_mock_mo_coefficients() {
        let coeffs = generate_mock_mo_coefficients(10, 6, 42);

        assert_eq!(coeffs.dim(), (10, 6));

        // Check orthonormality
        let overlap = coeffs.t().dot(&coeffs);
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((overlap[[i, j]] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_cholesky_validation() {
        let mut cholesky = CholeskyMetric {
            l_matrix: Array2::eye(5),
            condition_number: 1.0,
            naux: 5,
            method: CholeskyMethod::Standard,
            rank: None,
            permutation: None,
            eigenvalues: None,
        };

        // Valid Cholesky
        assert!(cholesky.validate().is_ok());

        // Break lower triangular structure
        cholesky.l_matrix[[0, 1]] = 1.0;
        assert!(cholesky.validate().is_err());

        // Fix and test wrong dimensions
        cholesky.l_matrix = Array2::eye(5);
        cholesky.naux = 6;
        cholesky.eigenvalues = None;
        assert!(cholesky.validate().is_err());
    }

    #[test]
    #[should_panic(expected = "Cannot have more MOs than AOs")]
    fn test_invalid_mo_dimensions() {
        // Try to create more MOs than AOs
        generate_mock_mo_coefficients(5, 10, 42);
    }

    #[test]
    fn test_transform_dimension_mismatch() {
        let j3c_ao = Array3::zeros((8, 8, 18));
        let c_occ = Array2::zeros((7, 5)); // Wrong AO dimension
        let c_vir = Array2::zeros((8, 3));

        let result = transform_mo_3center(&j3c_ao, &c_occ, &c_vir);
        assert!(result.is_err());
    }

    #[test]
    fn test_transform_config() {
        // Test configuration options
        let config = TransformConfig::default();
        assert!((config.max_memory_gb - 4.0).abs() < f64::EPSILON);
        assert_eq!(config.blas_threads, 1); // Changed to 1 for Rayon compatibility
        assert!(config.use_blocking);

        // Test memory estimation
        let mem = estimate_memory_requirement(100, 20, 80, 200);
        assert!(mem > 0.0);

        // Test chunk size calculation
        let chunk = calculate_optimal_chunk_size(1000, 100, 20, 80, &config);
        assert!(chunk > 0);
        assert!(chunk <= 1000);
    }

    #[test]
    fn test_performance_comparison() {
        // Compare sequential vs blocked algorithms for a medium-sized system
        let j3c_ao = Array3::zeros((50, 50, 100));
        let c_occ = generate_mock_mo_coefficients(50, 10, 100);
        let c_vir = generate_mock_mo_coefficients(50, 40, 101);

        // Sequential
        let start = std::time::Instant::now();
        let j3c_seq = transform_mo_3center_sequential_optimized(&j3c_ao, &c_occ, &c_vir).unwrap();
        let seq_time = start.elapsed();

        // Blocked (if available)
        let config = TransformConfig {
            use_blocking: true,
            max_memory_gb: 1.0,
            ..Default::default()
        };
        let start = std::time::Instant::now();
        let j3c_blocked = transform_mo_3center_blocked(&j3c_ao, &c_occ, &c_vir, &config).unwrap();
        let blocked_time = start.elapsed();

        // Results should be identical
        let max_diff = j3c_seq
            .iter()
            .zip(j3c_blocked.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(max_diff < 1e-12, "Results differ: max_diff = {}", max_diff);

        println!("Sequential: {:?}, Blocked: {:?}", seq_time, blocked_time);
    }
}
