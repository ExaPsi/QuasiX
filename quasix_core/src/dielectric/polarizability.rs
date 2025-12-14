//! RPA Polarizability P₀(iω) (Story S3-2)
//!
//! **Status**: FIXED - Changed to Contour Deformation (CD) formula (2025-11-24)
//! **Reference**: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-2)
//!
// Physics notation: ia_P represents MO-transformed 3-center integrals (ia|P)
// This is standard notation in quantum chemistry literature
#![allow(non_snake_case)]
//!
//! # Formula (CD Method)
//!
//! P⁰_PQ(iω) = 4 Σᵢₐ (ia|P)(ia|Q) * εᵢₐ / (ω² + εᵢₐ²)
//!
//! where:
//! - i runs over occupied orbitals
//! - a runs over virtual orbitals
//! - P,Q are auxiliary basis indices
//! - εᵢₐ = εᵢ - εₐ (orbital energy difference, NEGATIVE for excitations)
//! - ω is imaginary frequency (REAL positive value, NOT iω)
//!
//! **CRITICAL FORMULA CHANGE (2025-11-24)**:
//! - OLD (AC): weight = 4 / (iω + εᵢₐ)  → Complex-valued P₀  ← WRONG for CD!
//! - NEW (CD): weight = 4 * εᵢₐ / (ω² + εᵢₐ²)  → REAL-valued P₀  ← CORRECT!
//!
//! **CRITICAL**: P₀ is now REAL-VALUED on imaginary axis
//! - Return type: `Array3<f64>` (changed from Complex64)
//! - CD formula ensures P₀ is symmetric and real
//! - No imaginary part on imaginary axis
//!
//! # PySCF Reference (CD Method)
//!
//! File: `pyscf/gw/gw_cd.py`, function `get_rho_response()`, lines 135-145
//!
//! ```python
//! def get_rho_response(omega, mo_energy, Lpq):
//!     '''Compute density response function in auxiliary basis at freq iw'''
//!     naux, nocc, nvir = Lpq.shape
//!     eia = mo_energy[:nocc,None] - mo_energy[None,nocc:]
//!     eia = eia/(omega**2+eia*eia)  # CD formula: ε_ia / (ω² + ε_ia²)
//!     Pia = einsum('Pia,ia->Pia',Lpq,eia)
//!     # Response from both spin-up and spin-down density
//!     Pi = 4. * einsum('Pia,Qia->PQ',Pia,Lpq)
//!     return Pi
//! ```
//!
//! # Implementation Plan (TDD Workflow)
//!
//! 1. **Extract PySCF baseline** (S3-2.1):
//!    - Run: `tests/validation/extract_pyscf_polarizability.py`
//!    - Save: `tests/validation/pyscf_baselines/P0_h2_sto3g.npy`
//!    - **CRITICAL**: Save at multiple frequencies to validate ξ dependence
//!
//! 2. **Write RED test** (S3-2.2):
//!    - Test: `tests/validation/test_polarizability_h2.py`
//!    - Expected: FAIL with unimplemented!()
//!
//! 3. **Implement GREEN code** (S3-2.3):
//!    - Implement `compute_polarizability_p0()` below
//!    - Use BLAS-3 operations (DGEMM)
//!    - Minimal code to pass test
//!
//! 4. **Validate vs PySCF** (S3-2.4):
//!    - Tolerance: max_diff < 1e-8
//!    - Check: ALL frequencies
//!    - Check: P₀ is real (no imaginary part)
//!    - Check: P₀ is symmetric (Hermitian)
//!
//! 5. **REFACTOR** (S3-2.5):
//!    - Clean up code
//!    - Add documentation
//!    - Optimize if needed (AFTER validation!)

use crate::common::{QuasixError, Result};
use ndarray::{s, Array1, Array2, Array3};
use rayon::prelude::*;

/// Compute RPA polarizability on imaginary axis (CD Method)
///
/// # Type Signature (UPDATED 2025-11-24 - CD Formula)
///
/// ```text
/// pub fn compute_polarizability_p0(
///     ia_P: &Array3<f64>,       // [n_occ, n_virt, n_aux]
///     mo_energy: &Array1<f64>,  // [n_mo]
///     mo_occ: &Array1<f64>,     // [n_mo]
///     freqs: &Array1<f64>,      // [n_freq]
/// ) -> Result<Array3<f64>>  // [n_freq, n_aux, n_aux] - REAL-valued!
/// ```
///
/// # Arguments
///
/// * `ia_P` - DF tensor (occupied-virtual blocks), shape [n_occ, n_virt, n_aux]
/// * `mo_energy` - Molecular orbital energies, shape [n_mo]
/// * `mo_occ` - Occupation numbers (2.0 for occupied, 0.0 for virtual), shape [n_mo]
/// * `freqs` - Imaginary frequency grid points (ω, treated as imaginary), shape [n_freq]
///
/// # Returns
///
/// Polarizability P₀(iω), shape [n_freq, n_aux, n_aux]
///
/// **CRITICAL**: Returns REAL array (f64) for CD method!
/// - CD formula ensures P₀ is real-valued on imaginary axis
/// - P₀ is symmetric (not Hermitian, since it's real)
/// - No imaginary part: max|Im(P₀)| = 0 exactly
///
/// # Formula Breakdown (CD Method, 2025-11-24)
///
/// For each frequency ω (positive real value):
/// 1. Compute transition energies: εᵢₐ = εₐ - εᵢ
/// 2. Compute real weight: 4 * εᵢₐ / (ω² + εᵢₐ²)
/// 3. Weight DF products: weight * (ia|P)(ia|Q)
/// 4. Sum over ia pairs
///
/// **CRITICAL FORMULA (CD vs AC)**:
/// - AC: weight = 4 / (iω + εᵢₐ)  → Complex P₀  ← WRONG for CD!
/// - CD: weight = 4 * εᵢₐ / (ω² + εᵢₐ²)  → Real P₀  ← CORRECT!
///
/// # BLAS-3 Implementation
///
/// ```ignore
/// // Efficient DGEMM-based implementation
/// for (f_idx, &omega) in freqs.iter().enumerate() {
///     let omega_imag = Complex64::new(0.0, omega);  // iω
///
///     // Step 1: Compute complex weighted DF tensor
///     let mut weighted = Array3::<Complex64>::zeros((n_occ, n_virt, n_aux));
///     for i in 0..n_occ {
///         for a in 0..n_virt {
///             let eps_ia = mo_energy[nocc + a] - mo_energy[i];
///             let denom = omega_imag + Complex64::new(eps_ia, 0.0);
///             let weight = Complex64::new(4.0, 0.0) / denom;
///             // Apply weight to DF tensor
///             for p in 0..n_aux {
///                 weighted[[i, a, p]] = Complex64::new(ia_P[[i, a, p]], 0.0) * weight;
///             }
///         }
///     }
///
///     // Step 2: Contract: P₀_PQ = Σᵢₐ weighted[i,a,P] * ia_P[i,a,Q]
///     // (Similar BLAS-3 contraction, but with complex numbers)
/// }
/// ```
///
/// # PySCF Reference
///
/// ```python
/// # pyscf/gw/gw_ac.py line 67
/// def get_rho_response(omega, mo_energy, mo_occ, Lpq):
///     eia = mo_energy[occ_idx] - mo_energy[virt_idx]
///     p0_pq = -2.0 * einsum('iaP,iaQ->PQ', Lpq, Lpq) / (omega**2 + eia**2)
///     # Effective -4 for RHF (includes spin)
/// ```
///
/// # Validation Checks (CD Method, 2025-11-24)
///
/// After implementation, verify:
/// 1. P₀ is real-valued (no imaginary part by construction)
/// 2. P₀ is symmetric (P₀_PQ = P₀_QP)
/// 3. P₀ → 0 as ω → ∞
/// 4. max_diff < 1e-8 vs PySCF CD reference
/// 5. W computed from P₀ is also real-valued
///
/// # Errors
///
/// Returns error if:
/// - Array dimensions inconsistent
/// - BLAS operations fail
/// - Result contains NaN/Inf
pub fn compute_polarizability_p0(
    ia_P: &Array3<f64>,
    mo_energy: &Array1<f64>,
    mo_occ: &Array1<f64>,
    freqs: &Array1<f64>,
) -> Result<Array3<f64>> {
    // CHANGED: Returns REAL Array3<f64> for CD method!
    let (n_occ, n_virt, n_aux) = ia_P.dim();
    let nfreq = freqs.len();

    // Identify occupied/virtual orbitals
    let occ_indices: Vec<usize> = mo_occ
        .iter()
        .enumerate()
        .filter(|(_, &occ)| occ > 0.5)
        .map(|(i, _)| i)
        .collect();
    let vir_indices: Vec<usize> = mo_occ
        .iter()
        .enumerate()
        .filter(|(_, &occ)| occ < 0.5)
        .map(|(i, _)| i)
        .collect();

    if occ_indices.len() != n_occ || vir_indices.len() != n_virt {
        return Err(QuasixError::InvalidInput(format!(
            "Occupation mismatch: expected {}/{}, got {}/{}",
            n_occ,
            n_virt,
            occ_indices.len(),
            vir_indices.len()
        )));
    }

    // Parallel computation over frequencies (CD method)
    //
    // PARALLELISM FIX (2025-12-13): Removed BlasThreadGuard::new() which was
    // forcing BLAS to single-threaded mode. Each frequency computation is
    // independent (no shared mutable state), so BLAS operations in different
    // Rayon tasks don't conflict. OpenBLAS is thread-safe for this pattern.
    //
    // Previous issue: BlasThreadGuard set OPENBLAS_NUM_THREADS=1, causing
    // ~1x speedup at 32 threads (3.1% efficiency).

    let p0_vec: Vec<Array2<f64>> = (0..nfreq)
        .into_par_iter()
        .map(|f_idx| {
            let xi = freqs[f_idx]; // ω (positive real value, NOT iω!)

            let mut p0_xi = Array2::<f64>::zeros((n_aux, n_aux));

            // Sum over occupied-virtual pairs
            for (i_idx, &i) in occ_indices.iter().enumerate() {
                for (a_idx, &a) in vir_indices.iter().enumerate() {
                    // CRITICAL FIX (2025-11-24): Use εᵢₐ = εᵢ - εₐ (NEGATIVE for excitations)
                    // This matches PySCF convention: eia = mo_energy[:nocc] - mo_energy[nocc:]
                    // Previous bug: Used εₐ - εᵢ (positive), causing SIGN ERROR in Σᶜ!
                    let eps_ia = mo_energy[i] - mo_energy[a]; // CORRECTED SIGN

                    // CD FORMULA (2025-11-24): weight = 4 * ε_ia / (ω² + ε_ia²)
                    // This is REAL-VALUED (both ω and ε_ia are real)
                    let denominator = xi * xi + eps_ia * eps_ia; // ω² + ε_ia²
                    let weight = 4.0 * eps_ia / denominator; // REAL weight!

                    // Contract: P⁰_PQ += weight * (ia|P)(ia|Q)
                    // Outer product: ia_P[i,a,:] ⊗ ia_P[i,a,:]
                    for p in 0..n_aux {
                        let ia_p = ia_P[[i_idx, a_idx, p]];
                        for q in 0..n_aux {
                            let ia_q = ia_P[[i_idx, a_idx, q]];
                            p0_xi[[p, q]] += weight * ia_p * ia_q; // Real arithmetic!
                        }
                    }
                }
            }

            // Symmetrize to enforce symmetry (symmetric, not Hermitian)
            // CD formula with real weight creates a symmetric matrix:
            //   P₀[p,q] = weight * (ia|p) * (ia|q) = P₀[q,p]
            for p in 0..n_aux {
                for q in p + 1..n_aux {
                    // Symmetric: P[p,q] = P[q,p] (real-valued)
                    let avg = 0.5 * (p0_xi[[p, q]] + p0_xi[[q, p]]);
                    p0_xi[[p, q]] = avg;
                    p0_xi[[q, p]] = avg;
                }
            }

            p0_xi
        })
        .collect();

    // Stack into 3D array
    let mut p0 = Array3::<f64>::zeros((nfreq, n_aux, n_aux));
    for (f_idx, p0_xi) in p0_vec.into_iter().enumerate() {
        p0.slice_mut(s![f_idx, .., ..]).assign(&p0_xi);
    }

    // VALIDATION: Check for NaN/Inf
    if p0.iter().any(|&x| !x.is_finite()) {
        return Err(QuasixError::NumericalError(
            "P₀ contains NaN or Inf values".to_string(),
        ));
    }

    Ok(p0)
}

// ============================================================================
// Helper Functions (to be implemented)
// ============================================================================

/// Compute polarizability at single frequency (helper for batching)
///
/// UPDATED 2025-11-24: Now returns f64 for CD method
pub fn compute_p0_single_frequency(
    _ia_P: &Array3<f64>,
    _mo_energy: &Array1<f64>,
    _mo_occ: &Array1<f64>,
    _omega: f64,
) -> Result<ndarray::Array2<f64>> {
    // CHANGED: f64 for CD method
    // TODO: Implement in S3-2
    unimplemented!("Single frequency P₀ awaiting implementation")
}

// ============================================================================
// Legacy Type Stubs (for compilation compatibility)
// ============================================================================

/// LEGACY: Polarizability configuration
/// Will be removed during S3-2 re-implementation
#[derive(Debug, Clone)]
pub struct PolarizabilityConfig {
    pub eta: f64,
}

impl Default for PolarizabilityConfig {
    fn default() -> Self {
        Self { eta: 1e-4 }
    }
}

/// LEGACY: Polarizability calculator (RI basis)
/// Will be removed during S3-2 re-implementation
pub struct PolarizabilityRI {
    pub nocc: usize,
    pub nvirt: usize,
    pub naux: usize,
}

impl PolarizabilityRI {
    pub fn new(nocc: usize, nvirt: usize, naux: usize) -> Self {
        Self { nocc, nvirt, naux }
    }

    pub fn compute_p0(
        &self,
        _omega: num_complex::Complex64,
        _df_ia: &ndarray::Array2<f64>,
        _e_occ: &ndarray::Array1<f64>,
        _e_virt: &ndarray::Array1<f64>,
    ) -> Result<ndarray::Array2<num_complex::Complex64>> {
        unimplemented!("Legacy PolarizabilityRI - use S3-2 implementation")
    }

    pub fn symmetrize_p0(
        &self,
        _p0: &ndarray::Array2<num_complex::Complex64>,
        _vsqrt: &ndarray::Array2<f64>,
    ) -> Result<ndarray::Array2<num_complex::Complex64>> {
        unimplemented!("Legacy symmetrize_p0 - use S3-2 implementation")
    }
}

/// LEGACY: Static polarizability
/// Will be removed during S3-2 re-implementation
pub struct StaticPolarizability;

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr1;

    #[test]
    fn test_p0_symmetry() {
        // Test that P₀ is SYMMETRIC on imaginary axis (CD method)
        // CD formula creates REAL-valued symmetric matrix
        let n_occ = 1;
        let n_virt = 1;
        let n_aux = 3;
        let nfreq = 2;

        let ia_p = Array3::<f64>::from_shape_fn((n_occ, n_virt, n_aux), |(i, a, p)| {
            (i + a + p) as f64 * 0.1
        });

        let mo_energy = arr1(&[-0.5, 0.5]);
        let mo_occ = arr1(&[2.0, 0.0]);
        let freqs = arr1(&[1.0, 2.0]);

        let p0 = compute_polarizability_p0(&ia_p, &mo_energy, &mo_occ, &freqs).unwrap();

        // Check SYMMETRY (P[p,q] = P[q,p] for real matrices)
        for f in 0..nfreq {
            for p in 0..n_aux {
                for q in 0..n_aux {
                    let diff = (p0[[f, p, q]] - p0[[f, q, p]]).abs();
                    assert!(
                        diff < 1e-14,
                        "P0[{},{},{}] not symmetric: diff = {}",
                        f,
                        p,
                        q,
                        diff
                    );
                }
            }
        }

        // P₀ is REAL-valued for CD method (no imaginary part)
        for f in 0..nfreq {
            for p in 0..n_aux {
                for q in 0..n_aux {
                    assert!(
                        p0[[f, p, q]].is_finite(),
                        "P0[{},{},{}] is not finite: {}",
                        f,
                        p,
                        q,
                        p0[[f, p, q]]
                    );
                }
            }
        }
    }

    #[test]
    fn test_p0_dimensions() {
        let n_occ = 2;
        let n_virt = 3;
        let n_aux = 5;
        let nfreq = 4;

        let ia_p = Array3::<f64>::from_shape_fn((n_occ, n_virt, n_aux), |(i, a, p)| {
            ((i + 1) * (a + 1) * (p + 1)) as f64 * 0.01
        });

        let mo_energy = arr1(&[-1.0, -0.5, 0.2, 0.5, 0.8]);
        let mo_occ = arr1(&[2.0, 2.0, 0.0, 0.0, 0.0]);
        let freqs = arr1(&[0.5, 1.0, 1.5, 2.0]);

        let p0 = compute_polarizability_p0(&ia_p, &mo_energy, &mo_occ, &freqs).unwrap();

        // Check dimensions
        assert_eq!(p0.dim(), (nfreq, n_aux, n_aux));
    }

    #[test]
    fn test_p0_prefactor_sign() {
        // Test that CD formula gives correct sign and magnitude
        let n_occ = 1;
        let n_virt = 1;
        let n_aux = 2;

        // Use constant DF values for predictability
        let ia_p = Array3::<f64>::from_elem((n_occ, n_virt, n_aux), 0.1);

        // IMPORTANT: eps_ia = mo_energy[occ] - mo_energy[virt]
        // For occupied orbital below virtual: eps_ia < 0 (normal case)
        // The CD formula weight = 4 * eps_ia / (omega^2 + eps_ia^2)
        // So for negative eps_ia, weight is NEGATIVE, and P0 is NEGATIVE
        //
        // To test positive P0, we need eps_ia > 0, which means occ energy > virt energy
        // This is unphysical but valid for testing the formula
        let mo_energy = arr1(&[0.5, -0.5]); // occ=0.5, virt=-0.5 => eps_ia = 1.0 > 0
        let mo_occ = arr1(&[2.0, 0.0]);
        let freqs = arr1(&[1.0]);

        let p0 = compute_polarizability_p0(&ia_p, &mo_energy, &mo_occ, &freqs).unwrap();

        // CD FORMULA:
        // eps_ia = mo_energy[0] - mo_energy[1] = 0.5 - (-0.5) = 1.0
        // omega = 1.0 (real)
        // denominator = omega² + eps_ia² = 1.0 + 1.0 = 2.0
        // weight = 4 * eps_ia / denominator = 4 * 1.0 / 2.0 = 2.0 (POSITIVE!)
        // contribution per element: 2.0 * 0.1 * 0.1 = 0.02 (POSITIVE!)
        // P₀ should be POSITIVE and REAL

        // Check actual values
        eprintln!("P0[0,0,0] = {:.6}", p0[[0, 0, 0]]);

        for p in 0..n_aux {
            for q in 0..n_aux {
                assert!(
                    p0[[0, p, q]] > 0.0,
                    "P0[{},{}] should be positive: {}",
                    p,
                    q,
                    p0[[0, p, q]]
                );
                // CD method: P₀ is exactly real (no imaginary part)
                assert!(
                    p0[[0, p, q]].is_finite(),
                    "P0[{},{}] not finite: {}",
                    p,
                    q,
                    p0[[0, p, q]]
                );
            }
        }
    }
}

// ============================================================================
// Legacy Code Removed
// ============================================================================
//
// Previous implementation (821 LOC) removed on 2025-11-23.
// Reasons for removal:
// - Used Complex64 for P₀ (should be real on imaginary axis!)
// - Incorrect prefactor (-2 instead of -4)
// - SIMD code without validation
// - Failed PySCF comparison
//
// See: docs/reports/2025-11-23/SPRINT3_COMPLETE_CLEANUP.md
// See: docs/G0W0/03_PYSCF_REFERENCE_FORMULAS.md (correct -4 prefactor)
// ============================================================================
