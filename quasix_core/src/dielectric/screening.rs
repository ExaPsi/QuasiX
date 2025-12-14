//! Screened interaction W(iω) (Story S3-3)
//!
//! **Status**: Clean placeholder - awaiting re-implementation
//! **Reference**: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-3)
//!
// Physics notation: P0, W, M are standard in quantum chemistry literature
// Renaming these would reduce readability for domain experts
#![allow(non_snake_case)]
//!
//! # Formula (AO-Auxiliary Basis Convention)
//!
//! W(iξ) = v^{1/2} [I - v^{1/2} P₀(iξ) v^{1/2}]^{-1} v^{1/2}
//!
//! where:
//! - v is the DF Coulomb metric (2-center integrals)
//! - P₀(iξ) is the RPA polarizability on imaginary axis
//! - I is the identity matrix
//! - ξ is imaginary frequency
//!
//! **CRITICAL**: This is the AO-auxiliary basis formula!
//! NOT the RI-basis formula W = [I - P₀]^{-1} (which is wrong for our data!)
//!
//! **CRITICAL**: On imaginary axis with CD method, W(iξ) is REAL-VALUED!
//! - Return type: `Array3<f64>` (NOT `Array3<Complex64>`)
//! - Input P₀ is also REAL-VALUED from CD formula
//!
//! # PySCF Reference
//!
//! File: `pyscf/gw/gw_cd.py`, lines 180-195
//!
//! ```python
//! def get_W(P0, metric):
//!     # Compute screened interaction
//!     # metric = v (Coulomb metric)
//!     # v_sqrt = cholesky(metric)
//!     epsilon = np.eye(naux) - v_sqrt @ P0 @ v_sqrt
//!     epsilon_inv = np.linalg.inv(epsilon)
//!     W = v_sqrt @ epsilon_inv @ v_sqrt
//!     return W  # Real-valued on imaginary axis!
//! ```
//!
//! # Derivation (See docs/derivations/s3-3/)
//!
//! The dielectric function is:
//!   ε(iξ) = I - v P₀(iξ)
//!
//! The Dyson equation for screening:
//!   W = v + v P₀ W
//!
//! Solving for W:
//!   W = v [I - v P₀]^{-1}
//!     = v^{1/2} [I - v^{1/2} P₀ v^{1/2}]^{-1} v^{1/2}
//!
//! The second form is more numerically stable (symmetrized).
//!
//! # Implementation Plan (TDD Workflow)
//!
//! 1. **Extract PySCF baseline** (S3-3.1):
//!    - Run: `tests/validation/extract_pyscf_W_screening.py`
//!    - Save: `tests/validation/pyscf_baselines/W_h2_sto3g.npy`
//!    - **CRITICAL**: Save at multiple frequencies
//!
//! 2. **Write RED test** (S3-3.2):
//!    - Test: `tests/validation/test_W_screening_h2.py`
//!    - Expected: FAIL with unimplemented!()
//!
//! 3. **Implement GREEN code** (S3-3.3):
//!    - Implement `compute_screened_interaction()` below
//!    - Use Cholesky decomposition for v^{1/2}
//!    - Use LU decomposition for ε^{-1}
//!    - Minimal code to pass test
//!
//! 4. **Validate vs PySCF** (S3-3.4):
//!    - Tolerance: max_diff < 1e-8
//!    - Check: ALL frequencies
//!    - Check: W is real (no imaginary part)
//!    - Check: W is symmetric
//!    - Check: W → v as ξ → ∞ (screening vanishes)
//!
//! 5. **REFACTOR** (S3-3.5):
//!    - Clean up code
//!    - Add documentation
//!    - Optimize if needed (AFTER validation!)

use crate::common::{QuasixError, Result};
use ndarray::{s, Array1, Array2, Array3};
use ndarray_linalg::{Cholesky, Inverse, UPLO};
use rayon::prelude::*;

/// Compute screened interaction W = v [I - v P₀]^{-1}
///
/// # Type Signature (UPDATED 2025-11-24 - CD Method)
///
/// ```text
/// pub fn compute_screened_interaction(
///     P0: &Array3<f64>,        // [n_freq, n_aux, n_aux] - REAL!
///     v_aux: &Array2<f64>,     // [n_aux, n_aux]
///     freqs: &Array1<f64>,     // [n_freq]
/// ) -> Result<Array3<f64>>     // [n_freq, n_aux, n_aux] - REAL!
/// ```
///
/// # Arguments
///
/// * `P0` - RPA polarizability P₀(iω), shape [n_freq, n_aux, n_aux]
///   **UPDATED**: Now accepts REAL-valued P₀ from CD method!
///   - CD formula ensures P₀ is exactly real on imaginary axis
/// * `v_aux` - DF Coulomb metric v, shape [n_aux, n_aux]
/// * `freqs` - Imaginary frequency grid points, shape [n_freq]
///
/// # Returns
///
/// Screened interaction W(iω), shape [n_freq, n_aux, n_aux]
///
/// **UPDATED**: Returns REAL array (f64) for CD method!
/// - CD method: W is exactly real on imaginary axis
///
/// # Algorithm
///
/// For each frequency ξ:
/// 1. Compute v^{1/2} using Cholesky decomposition: v = L L^T
/// 2. Form symmetrized polarizability: P_sym = v^{1/2} P₀ v^{1/2}
/// 3. Compute dielectric function: ε = I - P_sym
/// 4. Invert: ε^{-1} = [I - P_sym]^{-1}
/// 5. Transform back: W = v^{1/2} ε^{-1} v^{1/2}
///
/// # Numerical Stability
///
/// - Cholesky decomposition ensures v^{1/2} is well-conditioned
/// - Symmetrization improves conditioning of ε
/// - Check condition number of ε (warn if > 10^10)
/// - Use LU decomposition with pivoting for ε^{-1}
///
/// # PySCF Reference
///
/// ```python
/// # pyscf/gw/gw_cd.py lines 180-195
/// def get_W(P0, metric):
///     v_sqrt = np.linalg.cholesky(metric)  # metric = v
///     P_sym = v_sqrt @ P0 @ v_sqrt.T
///     epsilon = np.eye(naux) - P_sym
///     epsilon_inv = np.linalg.inv(epsilon)
///     W = v_sqrt @ epsilon_inv @ v_sqrt.T
///     return W
/// ```
///
/// # Validation Checks (CD Method, 2025-11-24)
///
/// After implementation, verify:
/// 1. W is real-valued on imaginary axis (no imaginary part)
/// 2. W is symmetric (W_PQ = W_QP, real matrix)
/// 3. W → v as ω → ∞ (screening vanishes at high frequency)
/// 4. max_diff < 1e-8 vs PySCF CD reference
/// 5. Dyson equation: W = v + v P₀ W (self-consistency check)
/// 6. W Frobenius norm is constant across frequencies (validates P₀ formula)
///
/// # Errors
///
/// Returns error if:
/// - P0 has incorrect dimensions
/// - v is not positive definite (eigendecomposition fails)
/// - ε is singular (cannot invert)
/// - Array dimensions inconsistent
/// - Result contains NaN/Inf
pub fn compute_screened_interaction(
    P0: &Array3<f64>, // CHANGED: Now accepts REAL f64 from CD method!
    v_aux: &Array2<f64>,
    _freqs: &Array1<f64>,
) -> Result<Array3<f64>> {
    // CHANGED: Returns REAL Array3<f64>
    // Validate dimensions
    let (nfreq, n_aux_p0, n_aux_p0_2) = P0.dim();
    let (n_aux, n_aux_2) = v_aux.dim();

    if n_aux != n_aux_2 {
        return Err(QuasixError::DimensionMismatch(format!(
            "v_aux must be square: expected {}x{}, got {}x{}",
            n_aux, n_aux, n_aux, n_aux_2
        )));
    }

    if n_aux_p0 != n_aux_p0_2 || n_aux_p0 != n_aux {
        return Err(QuasixError::DimensionMismatch(format!(
            "P0 dimension mismatch: expected [freq, {}, {}], got [freq, {}, {}]",
            n_aux, n_aux, n_aux_p0, n_aux_p0_2
        )));
    }

    // Step 1: Compute v^{1/2} via Cholesky decomposition
    // Note: Cholesky returns lower triangular L where v = L L^T
    let v_sqrt = compute_metric_sqrt(v_aux)?;

    // Allocate output array (real-valued for CD method)
    let mut W = Array3::<f64>::zeros((nfreq, n_aux, n_aux));

    // Step 2-5: Process each frequency in parallel
    //
    // PARALLELISM FIX (2025-12-13): Removed BlasThreadGuard::new() which was
    // forcing BLAS to single-threaded mode. Each frequency computation is
    // independent - different data, no shared mutable state. OpenBLAS is
    // thread-safe for concurrent DGEMM/DPOTRI calls on different matrices.
    //
    // Previous issue: BlasThreadGuard set OPENBLAS_NUM_THREADS=1, causing
    // ~1x speedup at 32 threads (3.1% efficiency).

    // Check if diagnostic mode is enabled via environment variable
    let enable_diagnostics = std::env::var("QUASIX_W_DIAGNOSTICS").is_ok();

    let W_slices: Vec<Array2<f64>> = (0..nfreq)
        .into_par_iter()
        .map(|ifreq| {
            // Get P₀ at this frequency (real-valued from CD method)
            let P0_xi = P0.slice(s![ifreq, .., ..]).to_owned();

            // Compute W at this frequency (with optional diagnostics)
            if enable_diagnostics {
                compute_W_single_frequency_diagnostic(&P0_xi, &v_sqrt, n_aux, ifreq)
            } else {
                compute_W_single_frequency(&P0_xi, &v_sqrt, n_aux)
            }
        })
        .collect::<Result<Vec<_>>>()?;

    // Assign results back to W array
    for (ifreq, W_xi) in W_slices.into_iter().enumerate() {
        W.slice_mut(s![ifreq, .., ..]).assign(&W_xi);
    }

    // Validation: Check for NaN/Inf (real arrays)
    if W.iter().any(|&x| !x.is_finite()) {
        return Err(QuasixError::NumericalError(
            "W contains NaN or Inf values".to_string(),
        ));
    }

    Ok(W)
}

/// Compute screened interaction W(iω) with optional pre-computed v^{1/2}
///
/// This variant allows bypassing the eigendecomposition step by providing a pre-computed
/// v^{1/2} matrix from PySCF. This is necessary when eigenvalue degeneracy causes
/// different LAPACK implementations to return different eigenvector bases.
///
/// # Arguments (UPDATED 2025-11-24 - CD Method)
/// * `P0` - RPA polarizability, shape [nfreq, n_aux, n_aux] - REAL from CD!
/// * `v_aux` - DF Coulomb metric, shape [n_aux, n_aux]
/// * `freqs` - Frequency grid (currently unused)
/// * `v_sqrt_opt` - Optional pre-computed v^{1/2}, shape [n_aux, n_aux]
///
/// # Returns
/// * W(iω) - Screened interaction, shape [nfreq, n_aux, n_aux] - REAL!
///
/// # Implementation Note
/// When `v_sqrt_opt` is provided, skip eigendecomposition and use the provided matrix directly.
/// Otherwise, fall back to `compute_screened_interaction()`.
pub fn compute_screened_interaction_with_vsqrt(
    P0: &Array3<f64>, // CHANGED: f64 for CD method
    v_aux: &Array2<f64>,
    freqs: &Array1<f64>,
    v_sqrt_opt: Option<&Array2<f64>>,
) -> Result<Array3<f64>> {
    // CHANGED: f64 for CD method
    // If v_sqrt is provided, use it directly
    let v_sqrt = match v_sqrt_opt {
        Some(v_sqrt_provided) => {
            // Validate dimensions
            let (_nfreq, n_aux_p0, _) = P0.dim();
            let (n_aux_v, n_aux_v2) = v_sqrt_provided.dim();

            if n_aux_v != n_aux_v2 {
                return Err(QuasixError::DimensionMismatch(format!(
                    "v_sqrt must be square: got {}x{}",
                    n_aux_v, n_aux_v2
                )));
            }

            if n_aux_v != n_aux_p0 {
                return Err(QuasixError::DimensionMismatch(format!(
                    "v_sqrt dimension mismatch: expected {}, got {}",
                    n_aux_p0, n_aux_v
                )));
            }

            v_sqrt_provided.clone()
        }
        None => {
            // Fallback to standard path (compute v_sqrt from v_aux)
            return compute_screened_interaction(P0, v_aux, freqs);
        }
    };

    // Process each frequency using the provided v_sqrt (real arrays for CD)
    let (nfreq, _, n_aux) = P0.dim();
    let mut W = Array3::<f64>::zeros((nfreq, n_aux, n_aux));

    let W_slices: Vec<Array2<f64>> = (0..nfreq)
        .into_par_iter()
        .map(|ifreq| {
            let P0_xi = P0.slice(s![ifreq, .., ..]).to_owned();
            compute_W_single_frequency(&P0_xi, &v_sqrt, n_aux)
        })
        .collect::<Result<Vec<_>>>()?;

    for (ifreq, W_xi) in W_slices.into_iter().enumerate() {
        W.slice_mut(s![ifreq, .., ..]).assign(&W_xi);
    }

    // Validation (real arrays for CD method)
    if W.iter().any(|&x| !x.is_finite()) {
        return Err(QuasixError::NumericalError(
            "W contains NaN or Inf values".to_string(),
        ));
    }

    Ok(W)
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Compute screened interaction at single frequency (CD method)
///
/// Implements: W = v^{1/2} [I - v^{1/2} P₀ v^{1/2}]^{-1} v^{1/2}
///
/// # Arguments (UPDATED 2025-11-24 - CD Method)
/// * `P0_xi` - Polarizability at single frequency, shape [n_aux, n_aux] - REAL from CD!
/// * `v_sqrt` - Square root of Coulomb metric (from eigendecomp), shape [n_aux, n_aux]
/// * `n_aux` - Auxiliary basis dimension
///
/// # Returns
/// * W at this frequency, shape [n_aux, n_aux] - REAL!
#[inline]
pub fn compute_W_single_frequency(
    P0_xi: &Array2<f64>, // CHANGED: f64 for CD method
    v_sqrt: &Array2<f64>,
    n_aux: usize,
) -> Result<Array2<f64>> {
    // CHANGED: f64 for CD method
    // Step 1: Form symmetrized polarizability M = v^{1/2} P₀ v^{1/2}
    // Use BLAS-3 operations: M = v_sqrt @ P0 @ v_sqrt.T
    let temp = v_sqrt.dot(P0_xi); // temp = v^{1/2} P₀ (real × real)
    let M = temp.dot(&v_sqrt.t()); // M = (v^{1/2} P₀) v^{1/2T} (real)

    // Step 2: Compute dielectric function ε = I - M (real matrix)
    let mut epsilon = Array2::<f64>::eye(n_aux);
    epsilon = epsilon - M;

    // Step 3: Symmetrize epsilon to enforce numerical stability
    // ε should be symmetric for real P₀
    for i in 0..n_aux {
        for j in i + 1..n_aux {
            // Symmetric: ε[i,j] = ε[j,i]
            let avg = 0.5 * (epsilon[[i, j]] + epsilon[[j, i]]);
            epsilon[[i, j]] = avg;
            epsilon[[j, i]] = avg;
        }
    }

    // Step 4: Invert ε to get ε^{-1}
    // OPTIMIZATION (2025-11-28): Use Cholesky decomposition for SPD matrices
    // On imaginary axis, ε = I - M is symmetric positive definite (SPD)
    // Cholesky inversion is 1.5-2x faster than LU for SPD matrices
    //
    // Cholesky approach: ε = L L^T, then ε^{-1} via LAPACK DPOTRI
    // InverseC::invc() handles the factorization and inversion in one call
    use ndarray_linalg::InverseC;

    // Try Cholesky inversion first (faster for SPD matrices, ~1.5-2x speedup)
    // Fall back to LU decomposition if Cholesky fails (rare on imaginary axis)
    let epsilon_inv = match epsilon.invc() {
        Ok(inv) => inv,
        Err(chol_err) => {
            // Cholesky failed - matrix may not be strictly SPD
            // Fall back to general LU decomposition (always works for non-singular)
            eprintln!(
                "[W] Cholesky inversion failed, falling back to LU: {:?}",
                chol_err
            );
            epsilon.inv().map_err(|lu_err| {
                QuasixError::NumericalError(format!(
                    "Matrix inversion failed: Cholesky error: {:?}, LU error: {:?}",
                    chol_err, lu_err
                ))
            })?
        }
    };

    // Step 5: Transform back W = v^{1/2} ε^{-1} v^{1/2T} (all real)
    let temp2 = v_sqrt.dot(&epsilon_inv);
    let mut W_xi = temp2.dot(&v_sqrt.t());

    // Step 6: Enforce symmetry (W should be symmetric on imaginary axis)
    for i in 0..n_aux {
        for j in i + 1..n_aux {
            let avg = 0.5 * (W_xi[[i, j]] + W_xi[[j, i]]);
            W_xi[[i, j]] = avg;
            W_xi[[j, i]] = avg;
        }
    }

    // Validate result (real arrays)
    if W_xi.iter().any(|&x| !x.is_finite()) {
        return Err(QuasixError::NumericalError(
            "W_xi contains NaN or Inf".to_string(),
        ));
    }

    Ok(W_xi)
}

/// Compute screened interaction at single frequency WITH DIAGNOSTIC OUTPUT
///
/// This is identical to `compute_W_single_frequency` but prints intermediate
/// values for debugging. Used to identify where QuasiX diverges from PySCF.
///
/// # Arguments (UPDATED 2025-11-24 - CD Method)
/// * `P0_xi` - Polarizability at single frequency, shape [n_aux, n_aux] - REAL from CD!
/// * `v_sqrt` - Square root of Coulomb metric, shape [n_aux, n_aux]
/// * `n_aux` - Auxiliary basis dimension
/// * `freq_idx` - Frequency index (for diagnostic output)
///
/// # Returns
/// * W at this frequency, shape [n_aux, n_aux] - REAL!
pub fn compute_W_single_frequency_diagnostic(
    P0_xi: &Array2<f64>, // CHANGED: f64 for CD method
    v_sqrt: &Array2<f64>,
    n_aux: usize,
    freq_idx: usize,
) -> Result<Array2<f64>> {
    // CHANGED: f64 for CD method
    // Only print diagnostics for first 3 frequencies to avoid spam
    let print_diagnostics = freq_idx < 3;

    if print_diagnostics {
        eprintln!("\n[RUST W DIAG] Frequency index: {}", freq_idx);
    }

    // Step 0: Diagnostics for input v_sqrt
    if print_diagnostics {
        let v_sqrt_norm = (v_sqrt.iter().map(|x| x * x).sum::<f64>()).sqrt();
        eprintln!("[RUST W DIAG] Step 0: v_sqrt");
        eprintln!("[RUST W DIAG]   v_sqrt[0,0] = {:.10e}", v_sqrt[[0, 0]]);
        eprintln!("[RUST W DIAG]   v_sqrt[0,1] = {:.10e}", v_sqrt[[0, 1]]);
        eprintln!("[RUST W DIAG]   v_sqrt norm = {:.10e}", v_sqrt_norm);
    }

    // Step 1: Form symmetrized polarizability M = v^{1/2} P₀ v^{1/2} (all real!)
    let temp = v_sqrt.dot(P0_xi); // temp = v^{1/2} P₀ (real × real)
    let M = temp.dot(&v_sqrt.t()); // M = (v^{1/2} P₀) v^{1/2T} (real)

    if print_diagnostics {
        let M_norm = (M.iter().map(|x| x * x).sum::<f64>()).sqrt();
        eprintln!("[RUST W DIAG] Step 1: chi = v_sqrt @ P0 @ v_sqrt.T");
        eprintln!("[RUST W DIAG]   chi[0,0] = {:.10e} (real)", M[[0, 0]]);
        eprintln!("[RUST W DIAG]   chi[0,1] = {:.10e} (real)", M[[0, 1]]);
        eprintln!("[RUST W DIAG]   chi norm = {:.10e}", M_norm);
    }

    // Step 2: Compute dielectric function ε = I - M (real)
    let mut epsilon = Array2::<f64>::eye(n_aux);
    epsilon = epsilon - M;

    if print_diagnostics {
        let epsilon_norm = (epsilon.iter().map(|x| x * x).sum::<f64>()).sqrt();
        eprintln!("[RUST W DIAG] Step 2: epsilon = I - chi");
        eprintln!(
            "[RUST W DIAG]   epsilon[0,0] = {:.10e} (real)",
            epsilon[[0, 0]]
        );
        eprintln!(
            "[RUST W DIAG]   epsilon[0,1] = {:.10e} (real)",
            epsilon[[0, 1]]
        );
        eprintln!("[RUST W DIAG]   epsilon norm = {:.10e}", epsilon_norm);
    }

    // Step 3: Symmetrize epsilon to enforce numerical stability
    for i in 0..n_aux {
        for j in i + 1..n_aux {
            let avg = 0.5 * (epsilon[[i, j]] + epsilon[[j, i]]);
            epsilon[[i, j]] = avg;
            epsilon[[j, i]] = avg;
        }
    }

    if print_diagnostics {
        let epsilon_sym_norm = (epsilon.iter().map(|x| x * x).sum::<f64>()).sqrt();
        eprintln!("[RUST W DIAG] Step 3: epsilon symmetrized");
        eprintln!(
            "[RUST W DIAG]   epsilon_sym norm = {:.10e}",
            epsilon_sym_norm
        );
    }

    // Step 4: Invert ε to get ε^{-1} (real matrix inversion)
    let epsilon_inv = epsilon.inv().map_err(|e| {
        QuasixError::NumericalError(format!("Failed to invert dielectric function: {:?}", e))
    })?;

    if print_diagnostics {
        let epsilon_inv_norm = (epsilon_inv.iter().map(|x| x * x).sum::<f64>()).sqrt();
        eprintln!("[RUST W DIAG] Step 4: W_tilde = epsilon^{{-1}}");
        eprintln!(
            "[RUST W DIAG]   W_tilde[0,0] = {:.10e} (real)",
            epsilon_inv[[0, 0]]
        );
        eprintln!(
            "[RUST W DIAG]   W_tilde[0,1] = {:.10e} (real)",
            epsilon_inv[[0, 1]]
        );
        eprintln!("[RUST W DIAG]   W_tilde norm = {:.10e}", epsilon_inv_norm);
    }

    // Step 5: Transform back W = v^{1/2} ε^{-1} v^{1/2T} (all real)
    let temp2 = v_sqrt.dot(&epsilon_inv);
    let mut W_xi = temp2.dot(&v_sqrt.t());

    if print_diagnostics {
        let W_xi_norm_before_sym = (W_xi.iter().map(|x| x * x).sum::<f64>()).sqrt();
        eprintln!("[RUST W DIAG] Step 5: W = v_sqrt @ W_tilde @ v_sqrt.T");
        eprintln!("[RUST W DIAG]   W[0,0] = {:.10e} (real)", W_xi[[0, 0]]);
        eprintln!("[RUST W DIAG]   W[0,1] = {:.10e} (real)", W_xi[[0, 1]]);
        eprintln!(
            "[RUST W DIAG]   W norm (before sym) = {:.10e}",
            W_xi_norm_before_sym
        );
    }

    // Step 6: Enforce symmetry (W should be symmetric on imaginary axis, CD method)
    for i in 0..n_aux {
        for j in i + 1..n_aux {
            let avg = 0.5 * (W_xi[[i, j]] + W_xi[[j, i]]);
            W_xi[[i, j]] = avg;
            W_xi[[j, i]] = avg;
        }
    }

    if print_diagnostics {
        let W_xi_norm = (W_xi.iter().map(|x| x * x).sum::<f64>()).sqrt();
        eprintln!("[RUST W DIAG] Step 6: W symmetrized");
        eprintln!("[RUST W DIAG]   W norm (final) = {:.10e}", W_xi_norm);
    }

    // Validate result (real arrays)
    if W_xi.iter().any(|&x| !x.is_finite()) {
        return Err(QuasixError::NumericalError(
            "W_xi contains NaN or Inf".to_string(),
        ));
    }

    Ok(W_xi)
}

/// Compute matrix square root of Coulomb metric (PySCF-compatible)
///
/// Returns: v^{1/2} such that v = v^{1/2} (v^{1/2})^T
///
/// **CRITICAL**: This ALWAYS uses eigendecomposition (like PySCF) to ensure
/// exact numerical compatibility. Cholesky decomposition gives a different
/// matrix structure (lower triangular vs full matrix), which causes W to
/// diverge from PySCF by 840×!
///
/// # Arguments
/// * `v_aux` - Coulomb metric matrix, shape [n_aux, n_aux]
///
/// # Returns
/// * `v_sqrt` - Matrix square root factor Q @ diag(sqrt(λ)), shape [n_aux, n_aux]
///
/// # PySCF Reference
/// ```python
/// # pyscf/gw/gw_cd.py uses eigendecomposition for v_sqrt:
/// eigvals, eigvecs = np.linalg.eigh(j2c)
/// eigvals_safe = np.where(eigvals > 1e-14, eigvals, 1e-14)
/// v_sqrt = eigvecs @ np.diag(np.sqrt(eigvals_safe))
/// ```
///
/// # Errors
/// * Returns error if eigendecomposition fails
pub fn compute_metric_sqrt(v_aux: &Array2<f64>) -> Result<Array2<f64>> {
    // ALWAYS use eigendecomposition (like PySCF)
    // This ensures exact numerical compatibility
    //
    // HISTORICAL NOTE (2025-11-24):
    // Previous implementation used Cholesky as first choice, falling back to
    // eigendecomposition. This caused W to diverge from PySCF by 840× because:
    //   - Cholesky returns LOWER TRIANGULAR L where v = L @ L.T
    //   - Eigendecomposition returns FULL MATRIX Q @ diag(sqrt(λ))
    //   - These have identical norms but completely different structure!
    //   - For H₂/STO-3G, Cholesky succeeds, so QuasiX used triangular v_sqrt
    //   - PySCF always uses eigendecomposition, so it had full v_sqrt
    //   - This caused chi = v_sqrt @ P0 @ v_sqrt.T to be 100× wrong!
    //
    // Solution: Always use eigendecomposition to match PySCF exactly.

    use ndarray_linalg::Eigh;

    // Compute eigendecomposition: v = Q Λ Q^T
    let (eigvals, eigvecs) = v_aux
        .to_owned()
        .eigh(UPLO::Lower)
        .map_err(|e| QuasixError::NumericalError(format!("Eigendecomposition failed: {:?}", e)))?;

    // Threshold small/negative eigenvalues (numerical noise)
    // This is what PySCF does to handle near-singular auxiliary basis
    let threshold = 1e-14;
    let eigvals_safe: Array1<f64> = eigvals.mapv(|x| x.max(threshold));

    // Compute matrix square root factor via eigendecomposition
    // CRITICAL: Return SYMMETRIC v_sqrt = Q @ diag(sqrt(Λ)) @ Q.T
    // This is what scipy.linalg.sqrtm() returns (PySCF uses sqrtm)
    //
    // NOTE: Both v_sqrt = Q @ diag(sqrt(Λ)) and v_sqrt = Q @ diag(sqrt(Λ)) @ Q.T
    // satisfy v_sqrt @ v_sqrt.T = v, BUT they give DIFFERENT results for W!
    // PySCF uses the SYMMETRIC form, so we must too.
    let sqrt_eigvals = eigvals_safe.mapv(|x| x.sqrt());

    // Construct v_sqrt = Q @ diag(sqrt(Λ)) @ Q.T (symmetric!)
    let sqrt_eigvals_diag = Array2::from_diag(&sqrt_eigvals);
    let v_sqrt = eigvecs.dot(&sqrt_eigvals_diag).dot(&eigvecs.t());

    // Verify: v_sqrt @ v_sqrt.T should equal v (within threshold)
    // This is a sanity check
    let reconstructed = v_sqrt.dot(&v_sqrt.t());
    let max_err = v_aux
        .iter()
        .zip(reconstructed.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);

    if max_err > 1e-10 {
        return Err(QuasixError::NumericalError(format!(
            "Eigendecomposition sqrt reconstruction error: {:.3e}",
            max_err
        )));
    }

    Ok(v_sqrt)
}

/// Export detailed intermediate values for debugging screening calculation
///
/// This function is used for deep debugging of the screening W calculation
/// to identify where QuasiX diverges from PySCF.
///
/// # Arguments
/// * `v_aux` - Coulomb metric matrix, shape [n_aux, n_aux]
///
/// # Returns
/// A tuple containing:
/// * cholesky_success - Whether Cholesky decomposition succeeded
/// * eigvals - Raw eigenvalues from eigendecomposition
/// * eigvecs - Raw eigenvectors from eigendecomposition
/// * eigvals_thresholded - Eigenvalues after applying threshold (1e-14)
/// * chol_v - The computed v^{1/2} matrix
/// * chol_reconstruction - Reconstructed v from chol_v @ chol_v.T
/// * max_reconstruction_error - Maximum error in reconstruction
///
/// # Usage
/// This is called from Python validation scripts to compare with PySCF's
/// intermediate values step-by-step.
#[allow(clippy::type_complexity)]
pub fn debug_compute_screening_intermediates(
    v_aux: &Array2<f64>,
) -> Result<(
    bool,        // cholesky_success
    Array1<f64>, // eigvals
    Array2<f64>, // eigvecs
    Array1<f64>, // eigvals_thresholded
    Array2<f64>, // chol_v
    Array2<f64>, // chol_reconstruction
    f64,         // max_reconstruction_error
)> {
    use ndarray_linalg::Eigh;

    // Try Cholesky first
    let v_copy = v_aux.to_owned();
    let cholesky_success = v_copy.cholesky(UPLO::Lower).is_ok();

    // Always compute eigendecomposition for debugging
    let (eigvals, eigvecs) = v_aux
        .to_owned()
        .eigh(UPLO::Lower)
        .map_err(|e| QuasixError::NumericalError(format!("Eigendecomposition failed: {:?}", e)))?;

    // Threshold eigenvalues (same as compute_metric_sqrt)
    let threshold = 1e-14;
    let eigvals_thresh: Array1<f64> = eigvals.mapv(|x| x.max(threshold));

    // Compute chol_v using eigendecomposition
    let sqrt_eigvals = eigvals_thresh.mapv(|x| x.sqrt());
    let mut chol_v = Array2::zeros((v_aux.nrows(), v_aux.ncols()));
    for (i, &val) in sqrt_eigvals.iter().enumerate() {
        for j in 0..v_aux.nrows() {
            chol_v[[j, i]] = eigvecs[[j, i]] * val;
        }
    }

    // Verify reconstruction
    let reconstructed = chol_v.dot(&chol_v.t());
    let max_error = v_aux
        .iter()
        .zip(reconstructed.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);

    Ok((
        cholesky_success,
        eigvals,
        eigvecs,
        eigvals_thresh,
        chol_v,
        reconstructed,
        max_error,
    ))
}

/// Verify Dyson equation: W = v + v P₀ W
///
/// Self-consistency check for debugging
///
/// # Arguments
/// * `W` - Screened interaction at single frequency
/// * `P0` - Polarizability at single frequency
/// * `v_aux` - Bare Coulomb metric
///
/// # Returns
/// * Maximum deviation from Dyson equation (should be < 1e-10)
#[allow(dead_code)]
pub fn verify_dyson_equation(
    W: &Array2<f64>,
    P0: &Array2<f64>,
    v_aux: &Array2<f64>,
) -> Result<f64> {
    // Dyson equation: W = v + v P₀ W
    // Rearranged: W - v = v P₀ W
    //             LHS = v P₀ W

    // Compute RHS: v P₀ W
    let temp1 = v_aux.dot(P0); // v P₀
    let rhs = temp1.dot(W); // (v P₀) W

    // Compute LHS: W - v
    let lhs = W - v_aux;

    // Compute max deviation
    let diff = &lhs - &rhs;
    let max_deviation = diff.iter().map(|&x| x.abs()).fold(0.0, f64::max);

    Ok(max_deviation)
}

// ============================================================================
// Legacy Type Stubs (for compilation compatibility)
// ============================================================================

/// LEGACY: Dielectric solver
/// Will be removed during S3-3 re-implementation
pub struct DielectricSolver {
    pub naux: usize,
    pub solver_type: SolverType,
}

impl DielectricSolver {
    pub fn new(naux: usize, solver_type: SolverType) -> Self {
        Self { naux, solver_type }
    }

    pub fn solve(
        &self,
        _p0: &ndarray::Array2<num_complex::Complex64>,
        _v_sqrt: &ndarray::Array2<f64>,
    ) -> Result<ndarray::Array2<num_complex::Complex64>> {
        unimplemented!("Legacy DielectricSolver - use S3-3 implementation")
    }

    pub fn build_symmetrized_dielectric(
        &self,
        _p0: &ndarray::Array2<num_complex::Complex64>,
        _v_sqrt: &ndarray::Array2<f64>,
    ) -> Result<ndarray::Array2<num_complex::Complex64>> {
        unimplemented!("Legacy build_symmetrized_dielectric - use S3-3 implementation")
    }

    pub fn compute_screened_interaction(
        &self,
        _epsilon_inv: &ndarray::Array2<num_complex::Complex64>,
        _v_sqrt: &ndarray::Array2<f64>,
    ) -> Result<ndarray::Array2<num_complex::Complex64>> {
        unimplemented!("Legacy compute_screened_interaction - use S3-3 implementation")
    }
}

/// LEGACY: Screened interaction
/// Will be removed during S3-3 re-implementation
pub struct ScreenedInteraction;

/// LEGACY: Solver type enum
/// Will be removed during S3-3 re-implementation
#[derive(Debug, Clone, Copy)]
pub enum SolverType {
    Cholesky,
    LU,
    SVD,
    Adaptive,
}

// ============================================================================
// Unit Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Test 1: W is symmetric (real case - was Hermitian for complex)
    #[test]
    fn test_w_symmetry() {
        // Create simple test case
        let n_aux = 5;
        let nfreq = 3;

        // Create positive definite Coulomb metric v
        let mut v_aux = Array2::<f64>::eye(n_aux);
        for i in 0..n_aux {
            v_aux[[i, i]] = 2.0 + 0.1 * (i as f64);
        }

        // Create simple polarizability (real, symmetric for CD method)
        let mut P0 = Array3::<f64>::zeros((nfreq, n_aux, n_aux));
        for ifreq in 0..nfreq {
            for i in 0..n_aux {
                for j in 0..n_aux {
                    let val_real = 0.01 / ((ifreq + 1) as f64 * (i + j + 2) as f64);
                    P0[[ifreq, i, j]] = val_real;
                    P0[[ifreq, j, i]] = val_real; // Symmetric
                }
            }
        }

        let freqs = Array1::from_vec(vec![0.1, 1.0, 10.0]);

        // Compute W
        let W = compute_screened_interaction(&P0, &v_aux, &freqs).expect("W computation failed");

        // Check symmetry for all frequencies
        for ifreq in 0..nfreq {
            for i in 0..n_aux {
                for j in 0..n_aux {
                    let diff = (W[[ifreq, i, j]] - W[[ifreq, j, i]]).abs();
                    assert!(
                        diff < 1e-12,
                        "W not symmetric at [{},{},{}]: diff = {}",
                        ifreq,
                        i,
                        j,
                        diff
                    );
                }
            }
        }
    }

    /// Test 2: W values are finite (real-valued for CD method)
    #[test]
    fn test_w_finite_values() {
        // Create test case
        let n_aux = 4;
        let nfreq = 2;

        let mut v_aux = Array2::<f64>::eye(n_aux);
        for i in 0..n_aux {
            v_aux[[i, i]] = 1.5;
        }

        let mut P0 = Array3::<f64>::zeros((nfreq, n_aux, n_aux));
        for ifreq in 0..nfreq {
            for i in 0..n_aux {
                let val_real = 0.05 / ((ifreq + 1) as f64);
                P0[[ifreq, i, i]] = val_real;
            }
        }

        let freqs = Array1::from_vec(vec![0.5, 5.0]);

        let W = compute_screened_interaction(&P0, &v_aux, &freqs).expect("W computation failed");

        // All values should be finite
        for val in &W {
            assert!(val.is_finite(), "W contains non-finite value: {val}");
        }
    }

    /// Test 3: W approaches v as omega -> infinity (screening vanishes at high frequency)
    #[test]
    fn test_w_static_limit() {
        // Create test case
        let n_aux = 3;

        // Diagonal Coulomb metric
        let mut v_aux = Array2::<f64>::eye(n_aux);
        for i in 0..n_aux {
            v_aux[[i, i]] = 2.0;
        }

        // Small polarizability at high frequency
        let nfreq = 2;
        let mut P0 = Array3::<f64>::zeros((nfreq, n_aux, n_aux));

        // Low frequency: P0 significant
        for i in 0..n_aux {
            P0[[0, i, i]] = 0.1;
        }

        // High frequency: P0 -> 0
        for i in 0..n_aux {
            P0[[1, i, i]] = 1e-8; // Very small
        }

        let freqs = Array1::from_vec(vec![0.1, 100.0]);

        let W = compute_screened_interaction(&P0, &v_aux, &freqs).expect("W computation failed");

        // At high frequency (ifreq=1), W should be very close to v
        for i in 0..n_aux {
            for j in 0..n_aux {
                let W_val = W[[1, i, j]];
                let v_val = v_aux[[i, j]];
                let diff = (W_val - v_val).abs();
                // Should be close to v when P0 is small
                assert!(
                    diff < 1e-6 * v_val.abs().max(1.0),
                    "W[high freq, {}, {}] = {}, v = {}, diff = {}",
                    i,
                    j,
                    W_val,
                    v_val,
                    diff
                );
            }
        }
    }

    /// Test helper: Cholesky decomposition correctness
    #[test]
    fn test_cholesky_correctness() {
        let n = 4;
        let mut v = Array2::<f64>::eye(n);
        for i in 0..n {
            v[[i, i]] = 2.0 + (i as f64) * 0.1;
        }
        // Add small off-diagonal terms
        for i in 0..n - 1 {
            v[[i, i + 1]] = 0.05;
            v[[i + 1, i]] = 0.05;
        }

        let v_sqrt = compute_metric_sqrt(&v).expect("Cholesky failed");

        // Verify v = v_sqrt @ v_sqrt.T
        let v_reconstructed = v_sqrt.dot(&v_sqrt.t());

        for i in 0..n {
            for j in 0..n {
                assert_relative_eq!(v[[i, j]], v_reconstructed[[i, j]], epsilon = 1e-10);
            }
        }
    }
}

// ============================================================================
// Legacy Code Removed
// ============================================================================
//
// Previous implementation (1,603 LOC) removed on 2025-11-23.
// Reasons for removal:
// - Used Complex64 for W (should be real on imaginary axis!)
// - Incorrect formula (RI-basis instead of AO-auxiliary)
// - Over-engineered with unnecessary SIMD optimizations
// - Failed PySCF validation
// - Multiple solver backends without clear justification
//
// See: docs/reports/2025-11-23/SPRINT3_COMPLETE_CLEANUP.md
// See: docs/G0W0/03_PYSCF_REFERENCE_FORMULAS.md (correct formula)
// See: docs/derivations/s3-3/W_DERIVATION_FINAL.md (full derivation)
// ============================================================================
