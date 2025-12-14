//! Exchange self-energy Σˣ (Story S3-4)
//!
//! **Status**: Clean placeholder - awaiting re-implementation
//! **Reference**: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-4)
// Note: into_shape is deprecated but we use it for stable array reshaping
// TODO: Migrate to into_shape_with_order when API stabilizes
#![allow(deprecated)]
//!
//! # Formula (GW Convention)
//!
//! Σˣₙₙ = -Σᵢ Σ_PQ (ni|P) v⁻¹_PQ (Q|ni)
//!
//! where:
//! - n is the orbital index
//! - i runs over occupied orbitals
//! - P,Q are auxiliary basis indices
//! - v⁻¹ is the inverse DF metric
//!
//! **CRITICAL**: This is NOT the same as Hartree-Fock exchange!
//! GW exchange includes correlation effects through the density fitting.
//!
//! # PySCF Reference
//!
//! File: `pyscf/gw/gw_cd.py`, line 77
//! ```python
//! def get_sigma_x(gw, mo_energy, mo_coeff):
//!     # Exchange self-energy in GW approximation
//!     # Uses DF metric v⁻¹ (not v as in HF)
//! ```
//!
//! # Implementation Plan (TDD Workflow)
//!
//! 1. **Extract PySCF baseline** (S3-4.1):
//!    - Run: `tests/validation/extract_pyscf_sigma_x.py`
//!    - Save: `tests/validation/pyscf_baselines/sigma_x_h2_sto3g.npy`
//!
//! 2. **Write RED test** (S3-4.2):
//!    - Test: `tests/validation/test_sigma_x_h2.py`
//!    - Expected: FAIL with unimplemented!()
//!
//! 3. **Implement GREEN code** (S3-4.3):
//!    - Implement `compute_exchange_diagonal()` below
//!    - Minimal code to pass test
//!
//! 4. **Validate vs PySCF** (S3-4.4):
//!    - Tolerance: max_diff < 1e-8 Ha
//!    - Test must PASS
//!
//! 5. **REFACTOR** (S3-4.5):
//!    - Clean up code
//!    - Add documentation
//!    - Optimize if needed (AFTER validation!)

use crate::common::Result;
use ndarray::{Array1, Array2, Array3, Axis};
use rayon::prelude::*;

/// Compute exchange self-energy diagonal (GW convention)
///
/// # Type Signature
///
/// From `docs/G0W0/04_TYPE_SIGNATURES.md`:
/// ```text
/// pub fn compute_exchange_diagonal(
///     ia_P: &Array3<f64>,      // [n_occ, n_virt, n_aux]
///     v_aux_inv: &Array2<f64>, // [n_aux, n_aux]
///     n_mo: usize,
/// ) -> Result<Array1<f64>>     // [n_mo]
/// ```
///
/// # Arguments
///
/// * `ia_P` - Density fitting tensor (occupied-virtual blocks), shape [n_occ, n_virt, n_aux]
/// * `v_aux_inv` - Inverse DF metric matrix, shape [n_aux, n_aux]
/// * `n_mo` - Total number of molecular orbitals
///
/// # Returns
///
/// Exchange self-energy diagonal, shape [n_mo]
///
/// # PySCF Reference
///
/// ```python
/// # pyscf/gw/gw_cd.py line 77
/// sigma_x = -np.einsum('Ppq,Pqr,Rrs->ps', Lpq, Lpq, metric_inv)
/// ```
///
/// # Errors
///
/// Returns error if:
/// - Array dimensions are inconsistent
/// - BLAS operations fail
///
/// # Examples
///
/// ```ignore
/// use quasix_core::selfenergy::exchange::compute_exchange_diagonal;
/// use ndarray::Array1;
///
/// let sigma_x = compute_exchange_diagonal(&ia_P, &v_aux_inv, n_mo)?;
/// assert_eq!(sigma_x.len(), n_mo);
/// ```
pub fn compute_exchange_diagonal(
    ia_p: &Array3<f64>,
    v_aux_inv: &Array2<f64>,
    n_mo: usize,
) -> Result<Array1<f64>> {
    // Validate inputs
    let (n_occ, n_virt, n_aux) = ia_p.dim();

    if v_aux_inv.dim() != (n_aux, n_aux) {
        return Err(crate::common::QuasixError::InvalidInput(format!(
            "DF metric dimension mismatch: expected ({}, {}), got ({}, {})",
            n_aux,
            n_aux,
            v_aux_inv.dim().0,
            v_aux_inv.dim().1
        )));
    }

    if n_occ + n_virt > n_mo {
        return Err(crate::common::QuasixError::InvalidInput(format!(
            "Inconsistent MO dimensions: n_occ={}, n_virt={}, n_mo={}",
            n_occ, n_virt, n_mo
        )));
    }

    // Check for NaN/Inf in inputs
    if ia_p.iter().any(|x| !x.is_finite()) {
        return Err(crate::common::QuasixError::InvalidInput(
            "ia_p contains NaN or Inf values".to_string(),
        ));
    }
    if v_aux_inv.iter().any(|x| !x.is_finite()) {
        return Err(crate::common::QuasixError::InvalidInput(
            "v_aux_inv contains NaN or Inf values".to_string(),
        ));
    }

    // Compute full exchange matrix using BLAS-3 operations
    let sigma_x_full = compute_exchange_matrix(ia_p, v_aux_inv)?;

    // Extract diagonal
    let mut sigma_x_diag = Array1::<f64>::zeros(n_mo);

    // Occupied-occupied block (diagonal only)
    for i in 0..n_occ {
        sigma_x_diag[i] = sigma_x_full[[i, i]];
    }

    // Virtual-virtual block is zero in G₀W₀ (no virtual-virtual exchange)
    // This is already initialized to zero above

    Ok(sigma_x_diag)
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Compute full exchange self-energy matrix
///
/// # Formula (GW Convention)
///
/// Σˣ_nn' = -Σᵢ Σ_PQ (ni|P) v⁻¹_PQ (Q|n'i)
///
/// For occupied-occupied block (i,j both occupied):
/// Σˣ_ij = -Σ_PQ Σₐ (ia|P) v⁻¹_PQ (Q|ja)
///
/// where:
/// - i,j are occupied orbital indices
/// - a runs over virtual orbitals
/// - P,Q are auxiliary basis indices
/// - v⁻¹ is the inverse DF metric
///
/// # Implementation Strategy
///
/// Use BLAS-3 operations for performance:
/// 1. Reshape ia_P to [n_occ, n_virt * n_aux]
/// 2. Compute intermediate: tmp_P = ia_P @ v_aux_inv (contract over aux)
/// 3. Compute Σˣ = -tmp @ ia_P.T (contract over virt * aux)
///
/// # Arguments
///
/// * `ia_P` - DF tensor [n_occ, n_virt, n_aux]
/// * `v_aux_inv` - Inverse DF metric [n_aux, n_aux]
///
/// # Returns
///
/// Exchange self-energy matrix [n_occ, n_occ] (occupied-occupied block only)
///
/// # Errors
///
/// Returns error if:
/// - Dimensions incompatible
/// - Result contains NaN/Inf
pub fn compute_exchange_matrix(ia_p: &Array3<f64>, v_aux_inv: &Array2<f64>) -> Result<Array2<f64>> {
    use ndarray::s;

    let (n_occ, n_virt, n_aux) = ia_p.dim();

    // Validate dimensions
    if v_aux_inv.dim() != (n_aux, n_aux) {
        return Err(crate::common::QuasixError::InvalidInput(format!(
            "DF metric dimension mismatch: expected ({}, {}), got ({}, {})",
            n_aux,
            n_aux,
            v_aux_inv.dim().0,
            v_aux_inv.dim().1
        )));
    }

    // Formula: Σˣ_ij = -Σ_PQ Σₐ (ia|P) v⁻¹_PQ (Q|ja)
    //
    // BLAS-3 Implementation:
    // 1. Reshape ia_P[i,a,P] to matrix [i, a*P]
    // 2. For each occupied i:
    //    - Extract ia_P[i,:,:] → [n_virt, n_aux]
    //    - Compute tmp[i,a,Q] = Σ_P ia_P[i,a,P] * v_inv[P,Q] via GEMM
    // 3. Compute Σˣ[i,j] = -Σₐ Σ_Q tmp[i,a,Q] * ia_P[j,a,Q] via GEMM

    // Allocate intermediate storage: tmp[i,a,Q] = ia_P[i,a,P] @ v_inv[P,Q]
    //
    // Step 1: Contract ia_p with v_aux_inv over auxiliary index (PARALLELIZED)
    //
    // Each occupied orbital i is independent:
    //   tmp[i,:,:] = ia_p[i,:,:] @ v_aux_inv
    //   Shape: [n_virt, n_aux] @ [n_aux, n_aux] -> [n_virt, n_aux]
    //
    // Performance: Expected 3-5x speedup on 16+ core systems (Xeon Silver 4314)
    // since each DGEMM is independent and memory-bound operations benefit from
    // parallel memory bandwidth utilization.

    // Compute all tmp[i] slices in parallel using Rayon
    // Each iteration performs: tmp[i,:,:] = ia_p[i,:,:] @ v_aux_inv
    let tmp_slices: Vec<Array2<f64>> = (0..n_occ)
        .into_par_iter()
        .map(|i| {
            // Extract ia_p[i,:,:] -> [n_virt, n_aux]
            let ia_slice = ia_p.slice(s![i, .., ..]);

            // Compute tmp[i,:,:] = ia_p[i,:,:] @ v_aux_inv via BLAS DGEMM
            // Shape: [n_virt, n_aux] @ [n_aux, n_aux] -> [n_virt, n_aux]
            ia_slice.dot(v_aux_inv)
        })
        .collect();

    // Stack the parallel results into Array3
    // Each tmp_slices[i] has shape [n_virt, n_aux], stack along axis 0 -> [n_occ, n_virt, n_aux]
    let tmp_views: Vec<_> = tmp_slices.iter().map(|arr| arr.view()).collect();
    let tmp = ndarray::stack(Axis(0), &tmp_views).map_err(|e| {
        crate::common::QuasixError::InvalidInput(format!("Failed to stack parallel results: {}", e))
    })?;

    // Step 2: Compute exchange matrix Σˣ[i,j] = -Σₐ Σ_Q tmp[i,a,Q] * ia_p[j,a,Q]
    //
    // Reshape to enable BLAS-3:
    // tmp[i,a,Q] → [i, a*Q] and ia_p[j,a,Q] → [j, a*Q]
    // Then: Σˣ = -tmp.reshape([i, a*Q]) @ ia_p.reshape([j, a*Q]).T
    //         = -tmp.reshape([i, a*Q]) @ ia_p.reshape([a*Q, j])

    let tmp_flat = tmp.into_shape((n_occ, n_virt * n_aux)).map_err(|e| {
        crate::common::QuasixError::InvalidInput(format!("Failed to reshape tmp: {}", e))
    })?;

    let ia_p_flat = ia_p
        .to_owned()
        .into_shape((n_occ, n_virt * n_aux))
        .map_err(|e| {
            crate::common::QuasixError::InvalidInput(format!("Failed to reshape ia_p: {}", e))
        })?;

    // Compute Σˣ = -tmp_flat @ ia_p_flat.T
    // Shape: [n_occ, n_virt*n_aux] @ [n_virt*n_aux, n_occ] → [n_occ, n_occ]
    let sigma_x = -tmp_flat.dot(&ia_p_flat.t());

    // Validate result
    if sigma_x.iter().any(|x| !x.is_finite()) {
        return Err(crate::common::QuasixError::NumericalError(
            "Exchange self-energy contains NaN or Inf".to_string(),
        ));
    }

    Ok(sigma_x)
}

// ============================================================================
// Legacy Type Stubs (for compilation compatibility)
// ============================================================================

/// LEGACY: Exchange self-energy calculator (RI basis)
/// Will be removed during S3-4 re-implementation
pub struct ExchangeSelfEnergyRI {
    pub nocc: usize,
    pub nbasis: usize,
    pub naux: usize,
}

impl ExchangeSelfEnergyRI {
    pub fn new(nbasis: usize, nocc: usize, naux: usize) -> Self {
        Self { nocc, nbasis, naux }
    }

    pub fn compute_exchange_matrix_ri(
        &mut self,
        _df_tensor: &Array3<f64>,
        _df_metric_inv: &Array2<f64>,
    ) -> Result<Array2<f64>> {
        unimplemented!("Legacy ExchangeSelfEnergyRI - use S3-4 implementation")
    }

    pub fn compute_exchange_diagonal_ri(
        &mut self,
        _df_tensor: &Array3<f64>,
        _df_metric_inv: &Array2<f64>,
    ) -> Result<Array1<f64>> {
        unimplemented!("Legacy ExchangeSelfEnergyRI - use S3-4 implementation")
    }
}

/// LEGACY: Exchange metadata
pub struct ExchangeMetadata;

// ============================================================================
// Legacy Code Removed
// ============================================================================
//
// Previous implementation (1,811 LOC) removed on 2025-11-23.
// Reasons for removal:
// - Used incorrect formulas (HF exchange instead of GW exchange)
// - Complex numerical issues
// - Failed PySCF validation
//
// See: docs/reports/2025-11-23/SPRINT3_COMPLETE_CLEANUP.md
// ============================================================================
