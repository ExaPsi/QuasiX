//! Green's Function G(iω) for GW Self-Energy (Story S3-5)
//!
//! **Status**: Clean implementation following TDD workflow
//! **Reference**: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-5.1)
//!
//! # Formula
//!
//! G_n(iω) = 1 / (iω + μ - ε_n)
//!
//! where:
//! - iω = imaginary frequency (purely imaginary)
//! - μ = chemical potential (Fermi level)
//! - ε_n = HF orbital energy
//!
//! # Physical Properties
//!
//! - Analytic in upper half-plane (iω > 0)
//! - Poles at iω = ε_n - μ
//! - Decay as 1/ω for large ω
//! - Complex-valued on imaginary axis
//!
//! # PySCF Reference
//!
//! File: `pyscf/gw/gw_ac.py`, lines 156-212
//!
//! ```python
//! def get_sigma_element(gw, omega, eia, evi):
//!     for w, wt in zip(freqs, weights):
//!         # Green's function: 1/(omega - e_p ± iη)
//!         # On imaginary axis: G_p(iω) = 1/(iω + μ - e_p)
//! ```

use crate::common::{QuasixError, Result};
use ndarray::{Array1, Array2};
use num_complex::Complex;

/// Compute Green's function kernel at complex frequencies
///
/// # Formula
///
/// G_n(ω) = 1 / (ω + μ - ε_n)
///
/// # Arguments
///
/// * `omega` - Complex frequencies [nfreq] (can be ±iω for occupied/virtual)
/// * `mo_energy` - HF orbital energies [nmo] (Hartree)
/// * `chemical_potential` - Fermi level μ (Hartree)
///
/// # Returns
///
/// Green's function [nfreq, nmo] (complex-valued)
///
/// Each element G[i, n] represents G_n(ω_i) for:
/// - i-th frequency point
/// - n-th molecular orbital
///
/// # Physical Validation
///
/// The function checks:
/// - All frequencies must be finite (no NaN/Inf)
/// - No NaN or Inf in result
///
/// # PySCF Convention
///
/// PySCF uses different signs for occupied and virtual orbitals:
/// - Occupied: ω = -iω (negative imaginary part)
/// - Virtual: ω = +iω (positive imaginary part)
///
/// # Examples
///
/// ```ignore
/// use quasix_core::selfenergy::green_function::compute_green_function;
/// use ndarray::array;
/// use num_complex::Complex;
///
/// let omega = array![Complex::new(0.0, 0.5), Complex::new(0.0, -0.5)];
/// let mo_energy = array![-0.5, -0.3, 0.2, 0.5];
/// let mu = -0.25;
///
/// let green = compute_green_function(&omega, &mo_energy, mu)?;
/// assert_eq!(green.dim(), (2, 4));
/// ```
///
/// # Errors
///
/// Returns error if:
/// - Any frequency is not finite (NaN or Inf)
/// - Result contains NaN or Inf
pub fn compute_green_function(
    omega: &Array1<Complex<f64>>,
    mo_energy: &Array1<f64>,
    chemical_potential: f64,
) -> Result<Array2<Complex<f64>>> {
    let nfreq = omega.len();
    let nmo = mo_energy.len();

    // Validate inputs
    if omega.iter().any(|z| !z.is_finite()) {
        return Err(QuasixError::InvalidInput(
            "Omega must be finite (no NaN/Inf)".to_string(),
        ));
    }

    if !chemical_potential.is_finite() {
        return Err(QuasixError::InvalidInput(format!(
            "Chemical potential must be finite, got: {}",
            chemical_potential
        )));
    }

    if mo_energy.iter().any(|x| !x.is_finite()) {
        return Err(QuasixError::InvalidInput(
            "MO energies contain NaN or Inf values".to_string(),
        ));
    }

    // Allocate result
    let mut green = Array2::<Complex<f64>>::zeros((nfreq, nmo));

    // Compute G_n(ω) = 1 / (ω + μ - ε_n)
    // where ω can be ±iω for occupied/virtual orbitals (PySCF convention)
    for (i, &w) in omega.iter().enumerate() {
        for (n, &en) in mo_energy.iter().enumerate() {
            // Denominator: ω + μ - ε_n
            let denominator = w + chemical_potential - en;

            // Green's function: G = 1 / denominator
            green[[i, n]] = Complex::new(1.0, 0.0) / denominator;
        }
    }

    // Validate output: check for NaN/Inf
    if green.iter().any(|z| !z.is_finite()) {
        return Err(QuasixError::NumericalError(
            "Green's function contains NaN or Inf".to_string(),
        ));
    }

    Ok(green)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::array;

    #[test]
    fn test_green_function_basic() {
        // Simple test case with known values
        let omega = array![Complex::new(0.0, 1.0)]; // Single frequency ω = +i
        let mo_energy = array![0.0]; // Single orbital at E=0
        let mu = 0.0; // Fermi level at 0

        let green = compute_green_function(&omega, &mo_energy, mu)
            .expect("Green's function computation failed");

        // Expected: G(i) = 1 / (i + 0 - 0) = 1/i = -i
        assert_eq!(green.dim(), (1, 1));

        let expected = Complex::new(0.0, -1.0);
        assert_abs_diff_eq!(green[[0, 0]].re, expected.re, epsilon = 1e-12);
        assert_abs_diff_eq!(green[[0, 0]].im, expected.im, epsilon = 1e-12);
    }

    #[test]
    fn test_green_function_multiple_orbitals() {
        // Test with HOMO-LUMO gap
        let omega = array![
            Complex::new(0.0, 0.5),
            Complex::new(0.0, 1.0),
            Complex::new(0.0, 2.0)
        ];
        let mo_energy = array![-0.5, -0.3, 0.2, 0.5]; // HOMO=-0.3, LUMO=0.2
        let mu = -0.05; // Fermi level between HOMO and LUMO

        let green = compute_green_function(&omega, &mo_energy, mu)
            .expect("Green's function computation failed");

        assert_eq!(green.dim(), (3, 4));

        // All values should be finite
        assert!(green.iter().all(|z| z.is_finite()));

        // Check that magnitude decreases with frequency for fixed orbital
        // |G(iω)| ~ 1/ω for large ω
        let mag_low = green[[0, 0]].norm();
        let mag_high = green[[2, 0]].norm();
        assert!(
            mag_high < mag_low,
            "Green's function should decay with frequency"
        );
    }

    #[test]
    fn test_green_function_pyscf_convention() {
        // Test PySCF convention: different signs for occupied/virtual
        let omega_occ = array![Complex::new(0.0, -0.5)]; // Occupied: -iω
        let omega_vir = array![Complex::new(0.0, 0.5)]; // Virtual: +iω
        let mo_energy = array![-0.5, 0.0, 0.5];
        let mu = 0.0;

        // Both should succeed
        let result_occ = compute_green_function(&omega_occ, &mo_energy, mu);
        assert!(
            result_occ.is_ok(),
            "Should accept negative imaginary frequency"
        );

        let result_vir = compute_green_function(&omega_vir, &mo_energy, mu);
        assert!(
            result_vir.is_ok(),
            "Should accept positive imaginary frequency"
        );
    }

    #[test]
    fn test_green_function_validation() {
        // Test error handling for invalid inputs
        let omega = array![Complex::new(0.0, 0.5), Complex::new(0.0, 1.0)];
        let mo_energy = array![-0.5, 0.0, 0.5];
        let mu = 0.0;

        // Test NaN omega (should fail)
        let omega_invalid = array![Complex::new(f64::NAN, 0.5), Complex::new(0.0, 1.0)];
        let result = compute_green_function(&omega_invalid, &mo_energy, mu);
        assert!(result.is_err(), "Should reject NaN in omega");

        // Test valid inputs (should succeed)
        let result = compute_green_function(&omega, &mo_energy, mu);
        assert!(result.is_ok(), "Should accept valid complex frequencies");
    }

    #[test]
    fn test_green_function_complex_structure() {
        // Verify that Green's function has correct complex structure
        let omega = array![Complex::new(0.0, 1.0)];
        let mo_energy = array![0.5];
        let mu = 0.0;

        let green = compute_green_function(&omega, &mo_energy, mu).unwrap();

        // G(iω) = 1/(iω + μ - ε) = 1/(i*1.0 + 0 - 0.5) = 1/(−0.5 + i)
        // Expected: 1/(−0.5 + i) = (−0.5 − i)/((−0.5)² + 1²) = (−0.5 − i)/1.25
        let expected_re = -0.5 / 1.25;
        let expected_im = -1.0 / 1.25;

        assert_abs_diff_eq!(green[[0, 0]].re, expected_re, epsilon = 1e-12);
        assert_abs_diff_eq!(green[[0, 0]].im, expected_im, epsilon = 1e-12);
    }
}
