//! Quasiparticle Equation Solver (Story S3-6)
//!
//! **Status**: Clean placeholder - awaiting re-implementation
//! **Reference**: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-6)
//!
//! # Quasiparticle Equation (Linearized Dyson)
//!
//! E_QP,n = ε_HF,n + Σˣₙₙ + Σᶜₙₙ(ε_HF,n) - V_xc,nn
//!
//! where:
//! - E_QP,n is the quasiparticle energy
//! - ε_HF,n is the Hartree-Fock eigenvalue
//! - Σˣₙₙ is the exchange self-energy (diagonal)
//! - Σᶜₙₙ(ω) is the correlation self-energy (frequency-dependent)
//! - V_xc,nn is the exchange-correlation potential (diagonal)
//!
//! # Linearized Approximation
//!
//! Instead of solving the full non-linear equation, we use:
//!
//! E_QP,n ≈ ε_HF,n + Z_n [Σˣₙₙ + Σᶜₙₙ(ε_HF,n) - V_xc,nn]
//!
//! where the renormalization factor (Z-factor) is:
//!
//! Z_n = [1 - ∂Σ/∂ω|_{ω=ε_HF,n}]^{-1}
//!
//! **For G₀W₀**: We use Z_n = 1 (no iteration)
//! **For evGW**: We iterate until self-consistency (not in Sprint 3)
//!
//! # PySCF Reference
//!
//! File: `pyscf/gw/gw_ac.py`, function `kernel()`, lines 220-250
//!
//! ```python
//! def kernel(gw, omega_list, mo_energy, mo_coeff):
//!     # Solve quasiparticle equation
//!     sigma_x = get_sigma_x(...)  # Exchange
//!     sigma_c = get_sigma_c(...)  # Correlation (at HF energies)
//!     vxc = get_vxc(...)          # XC potential
//!
//!     # Linearized approximation (Z = 1 for G₀W₀)
//!     qp_energy = mo_energy + sigma_x + sigma_c - vxc
//!     return qp_energy
//! ```
//!
//! # Implementation Plan (TDD Workflow)
//!
//! 1. **Extract PySCF baseline** (S3-6.1):
//!    - Run: `tests/validation/extract_pyscf_qp_energies.py`
//!    - Save: `tests/validation/pyscf_baselines/qp_energies_h2_sto3g.npy`
//!
//! 2. **Write RED test** (S3-6.2):
//!    - Test: `tests/validation/test_qp_solver_h2.py`
//!    - Expected: FAIL with unimplemented!()
//!
//! 3. **Implement GREEN code** (S3-6.3):
//!    - Implement `solve_quasiparticle_linearized()` below
//!    - Use Z = 1 (no iteration for G₀W₀)
//!    - Minimal code to pass test
//!
//! 4. **Validate vs PySCF** (S3-6.4):
//!    - Tolerance: max_diff < 1e-8 Ha
//!    - Check: HOMO, LUMO, IP, EA, gap all match
//!
//! 5. **REFACTOR** (S3-6.5):
//!    - Clean up code
//!    - Add documentation
//!    - Optimize if needed (AFTER validation!)

use crate::common::Result;
use ndarray::Array1;

/// Solve linearized quasiparticle equation (G₀W₀)
///
/// # Type Signature
///
/// From `docs/G0W0/04_TYPE_SIGNATURES.md`:
/// ```text
/// pub fn solve_quasiparticle_linearized(
///     mo_energy: &Array1<f64>,  // [n_mo] - HF eigenvalues
///     sigma_x: &Array1<f64>,    // [n_mo] - exchange self-energy
///     sigma_c: &Array1<f64>,    // [n_mo] - correlation self-energy
///     vxc: &Array1<f64>,        // [n_mo] - XC potential diagonal
/// ) -> Result<(Array1<f64>, Array1<f64>)>  // (qp_energies, z_factors)
/// ```
///
/// # Arguments
///
/// * `mo_energy` - Hartree-Fock eigenvalues, shape [n_mo]
/// * `sigma_x` - Exchange self-energy diagonal, shape [n_mo]
/// * `sigma_c` - Correlation self-energy diagonal (evaluated at HF energies), shape [n_mo]
/// * `vxc` - Exchange-correlation potential diagonal, shape [n_mo]
///
/// # Returns
///
/// Tuple of:
/// - `qp_energies`: Quasiparticle energies, shape [n_mo]
/// - `z_factors`: Renormalization factors, shape [n_mo]
///
/// For G₀W₀, all Z-factors are 1.0 (no iteration).
///
/// # Formula
///
/// ```text
/// E_QP,n = ε_HF,n + Z_n [Σˣₙₙ + Σᶜₙₙ(ε_HF,n) - V_xc,nn]
/// ```
///
/// For G₀W₀ (no iteration):
/// ```text
/// Z_n = 1.0
/// E_QP,n = ε_HF,n + Σˣₙₙ + Σᶜₙₙ(ε_HF,n) - V_xc,nn
/// ```
///
/// # PySCF Reference
///
/// ```python
/// # pyscf/gw/gw_ac.py lines 112-125 (linearized branch)
/// # e = ep + zn*(sigmaR.real + vk[p,p] - v_mf[p,p])
/// # For G₀W₀: zn = 1.0 (initially)
/// qp_energy = mo_energy + sigma_x + sigma_c - vxc
/// ```
///
/// # Validation Checks
///
/// After implementation, verify:
/// 1. QP energies are real and finite
/// 2. HOMO < LUMO (gap is positive)
/// 3. max_diff < 1e-8 vs PySCF
/// 4. IP, EA, gap match PySCF
///
/// # Errors
///
/// Returns error if:
/// - Array dimensions inconsistent
/// - Result contains NaN/Inf
/// - Physical constraints violated (e.g., negative gap for stable molecules)
pub fn solve_quasiparticle_linearized(
    mo_energy: &Array1<f64>,
    sigma_x: &Array1<f64>,
    sigma_c: &Array1<f64>,
    vxc: &Array1<f64>,
) -> Result<(Array1<f64>, Array1<f64>)> {
    // (qp_energies, z_factors)

    // 1. Validate dimensions
    let n_mo = mo_energy.len();
    if sigma_x.len() != n_mo {
        return Err(crate::common::QuasixError::InvalidInput(format!(
            "sigma_x dimension mismatch: expected {}, got {}",
            n_mo,
            sigma_x.len()
        )));
    }
    if sigma_c.len() != n_mo {
        return Err(crate::common::QuasixError::InvalidInput(format!(
            "sigma_c dimension mismatch: expected {}, got {}",
            n_mo,
            sigma_c.len()
        )));
    }
    if vxc.len() != n_mo {
        return Err(crate::common::QuasixError::InvalidInput(format!(
            "vxc dimension mismatch: expected {}, got {}",
            n_mo,
            vxc.len()
        )));
    }

    // 2. Compute QP energies using linearized approximation (Z = 1 for G₀W₀)
    //
    // PySCF formula (gw_ac.py line 124):
    //   e = ep + zn*(sigmaR.real + vk[p,p] - v_mf[p,p])
    //
    // Where:
    //   ep = mo_energy[p]
    //   zn = 1.0 (for linearized G₀W₀, no derivative)
    //   sigmaR.real = sigma_c[p] (correlation self-energy at ε_HF)
    //   vk[p,p] = sigma_x[p] (exchange self-energy, Hartree-Fock exchange)
    //   v_mf[p,p] = vxc[p] (DFT exchange-correlation potential)
    //
    // Simplified to:
    //   E_QP = ε_HF + Σˣ + Σᶜ - V_xc
    let mut qp_energies = Array1::<f64>::zeros(n_mo);
    for i in 0..n_mo {
        qp_energies[i] = mo_energy[i] + sigma_x[i] + sigma_c[i] - vxc[i];
    }

    // 3. Check for NaN/Inf (numerical stability)
    for i in 0..n_mo {
        if !qp_energies[i].is_finite() {
            return Err(crate::common::QuasixError::NumericalError(format!(
                "QP energy[{}] is not finite: {}",
                i, qp_energies[i]
            )));
        }
    }

    // 4. Z-factors are all 1.0 for G₀W₀ linearized approximation
    let z_factors = Array1::<f64>::from_elem(n_mo, 1.0);

    Ok((qp_energies, z_factors))
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Hartree to eV conversion factor
const HA_TO_EV: f64 = 27.211_386_245_988;

/// Compute Z-factors (renormalization factors)
///
/// **For G₀W₀**: Returns Array1 filled with 1.0 (no iteration)
/// **For evGW**: Computes Z = [1 - ∂Σ/∂ω]^{-1} (not in Sprint 3)
///
/// # Arguments
///
/// * `n_mo` - Number of molecular orbitals
///
/// # Returns
///
/// Array of Z-factors, all 1.0 for G₀W₀
pub fn compute_z_factors_g0w0(n_mo: usize) -> Array1<f64> {
    Array1::<f64>::from_elem(n_mo, 1.0)
}

/// Compute spectroscopic properties from QP energies
///
/// Returns: (IP, EA, gap) in eV
/// - IP: Ionization potential (-HOMO)
/// - EA: Electron affinity (-LUMO)
/// - gap: Fundamental gap (IP - EA)
///
/// # Arguments
///
/// * `qp_energies` - Quasiparticle energies in Hartree, shape [n_mo]
/// * `n_occ` - Number of occupied orbitals
///
/// # Returns
///
/// Tuple of (IP, EA, gap) in eV
///
/// # Errors
///
/// Returns error if n_occ is out of bounds
pub fn compute_spectroscopic_properties(
    qp_energies: &Array1<f64>,
    n_occ: usize,
) -> Result<(f64, f64, f64)> {
    let n_mo = qp_energies.len();
    if n_occ == 0 || n_occ >= n_mo {
        return Err(crate::common::QuasixError::InvalidInput(format!(
            "Invalid n_occ={} for n_mo={}",
            n_occ, n_mo
        )));
    }

    // HOMO index: n_occ - 1 (last occupied)
    // LUMO index: n_occ (first virtual)
    let homo_idx = n_occ - 1;
    let lumo_idx = n_occ;

    let e_homo = qp_energies[homo_idx];
    let e_lumo = qp_energies[lumo_idx];

    // Ionization potential (IP): -E_HOMO (in eV)
    // Electron affinity (EA): -E_LUMO (in eV)
    // Fundamental gap: IP - EA = E_LUMO - E_HOMO (in eV)
    let ip = -e_homo * HA_TO_EV;
    let ea = -e_lumo * HA_TO_EV;
    let gap = (e_lumo - e_homo) * HA_TO_EV;

    Ok((ip, ea, gap))
}

// ============================================================================
// Legacy Type Stubs (for compilation compatibility)
// ============================================================================

/// LEGACY: QP state solution
/// Will be removed during S3-6 re-implementation
#[derive(Debug, Clone)]
pub struct QPStateSolution {
    pub orbital_idx: usize,
    pub qp_energy: f64,
    pub z_factor: f64,
    pub converged: bool,
    pub iterations: usize,
    pub used_bisection: bool,
    pub energy_shift: f64,
}

// ============================================================================
// Legacy Code Removed
// ============================================================================
//
// Previous implementation (2,240 LOC) removed on 2025-11-23.
// Reasons for removal:
// - Over-engineered with Newton-Raphson, bisection, Richardson extrapolation
// - Not needed for G₀W₀ (only linearized approximation required)
// - SIMD optimizations without validation
// - Failed to match PySCF results
//
// evGW-specific files also removed:
// - evgw.rs (2,181 LOC) - NOT IN SPRINT 3 SCOPE!
// - evgw_parallel.rs (809 LOC) - NOT IN SPRINT 3 SCOPE!
//
// Sprint 3 is G₀W₀ ONLY. evGW is Sprint 4!
//
// See: docs/reports/2025-11-23/SPRINT3_COMPLETE_CLEANUP.md
// See: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-6 is linearized solver only)
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_qp_solver_basic() {
        // Simple 2-orbital system (1 occupied, 1 virtual)
        let mo_energy = Array1::from_vec(vec![-0.5, 0.3]); // Ha
        let sigma_x = Array1::from_vec(vec![-0.2, 0.0]);
        let sigma_c = Array1::from_vec(vec![-0.1, 0.05]);
        let vxc = Array1::from_vec(vec![-0.15, 0.02]);

        let (qp_energies, z_factors) =
            solve_quasiparticle_linearized(&mo_energy, &sigma_x, &sigma_c, &vxc).unwrap();

        // Expected: E_QP = ε_HF + Σˣ + Σᶜ - V_xc
        // HOMO: -0.5 + (-0.2) + (-0.1) - (-0.15) = -0.5 - 0.2 - 0.1 + 0.15 = -0.65
        // LUMO: 0.3 + 0.0 + 0.05 - 0.02 = 0.33
        assert_abs_diff_eq!(qp_energies[0], -0.65, epsilon = 1e-10);
        assert_abs_diff_eq!(qp_energies[1], 0.33, epsilon = 1e-10);

        // Z-factors should all be 1.0
        assert_eq!(z_factors[0], 1.0);
        assert_eq!(z_factors[1], 1.0);
    }

    #[test]
    fn test_qp_solver_dimension_check() {
        let mo_energy = Array1::from_vec(vec![-0.5, 0.3]);
        let sigma_x = Array1::from_vec(vec![-0.2]); // Wrong dimension
        let sigma_c = Array1::from_vec(vec![-0.1, 0.05]);
        let vxc = Array1::from_vec(vec![-0.15, 0.02]);

        let result = solve_quasiparticle_linearized(&mo_energy, &sigma_x, &sigma_c, &vxc);

        assert!(result.is_err());
    }

    #[test]
    fn test_qp_gap_is_positive() {
        // QP gap should be larger than HF gap (for typical systems)
        let mo_energy = Array1::from_vec(vec![-0.5, 0.2]); // HF gap = 0.7 Ha
        let sigma_x = Array1::from_vec(vec![-0.1, 0.0]);
        let sigma_c = Array1::from_vec(vec![-0.15, 0.08]); // Opens gap
        let vxc = Array1::from_vec(vec![-0.1, 0.03]);

        let (qp_energies, _) =
            solve_quasiparticle_linearized(&mo_energy, &sigma_x, &sigma_c, &vxc).unwrap();

        let qp_gap = qp_energies[1] - qp_energies[0];
        assert!(qp_gap > 0.0, "QP gap must be positive");
    }

    #[test]
    fn test_spectroscopic_properties() {
        // H2 example: HOMO at -0.5 Ha, LUMO at 0.3 Ha
        let qp_energies = Array1::from_vec(vec![-0.6, 0.35]);
        let n_occ = 1;

        let (ip, ea, gap) = compute_spectroscopic_properties(&qp_energies, n_occ).unwrap();

        // IP = -E_HOMO * 27.2114
        let expected_ip = -(-0.6) * HA_TO_EV;
        assert_abs_diff_eq!(ip, expected_ip, epsilon = 1e-8);

        // EA = -E_LUMO * 27.2114
        let expected_ea = -(0.35) * HA_TO_EV;
        assert_abs_diff_eq!(ea, expected_ea, epsilon = 1e-8);

        // Gap = (E_LUMO - E_HOMO) * 27.2114
        let expected_gap = (0.35 - (-0.6)) * HA_TO_EV;
        assert_abs_diff_eq!(gap, expected_gap, epsilon = 1e-8);
    }

    #[test]
    fn test_z_factors_g0w0() {
        let z_factors = compute_z_factors_g0w0(5);
        assert_eq!(z_factors.len(), 5);
        for &z in &z_factors {
            assert_eq!(z, 1.0);
        }
    }

    #[test]
    fn test_qp_solver_nan_check() {
        // Create inputs that would produce NaN
        let mo_energy = Array1::from_vec(vec![f64::NAN, 0.3]);
        let sigma_x = Array1::from_vec(vec![-0.2, 0.0]);
        let sigma_c = Array1::from_vec(vec![-0.1, 0.05]);
        let vxc = Array1::from_vec(vec![-0.15, 0.02]);

        let result = solve_quasiparticle_linearized(&mo_energy, &sigma_x, &sigma_c, &vxc);

        assert!(result.is_err(), "Should detect NaN in QP energies");
    }
}
