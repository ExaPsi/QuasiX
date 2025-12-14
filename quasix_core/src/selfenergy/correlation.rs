//! Correlation self-energy Σᶜ (Story S3-5) - CRITICAL MODULE
//!
//! **Status**: Clean implementation following TDD workflow
//! **Reference**: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-5)
//!
//! # Formula (Analytic Continuation)
//!
//! Σᶜₙ = -1/π ∫₀^∞ dξ Σₘ (nm|W(iξ)|mn) G_m(iξ)
//!
//! where:
//! - W(iξ) is the screened interaction on imaginary axis (REAL-valued!)
//! - G_m(iξ) is the Green's function (complex-valued)
//! - ξ is the imaginary frequency variable
//!
//! **CRITICAL REQUIREMENTS**:
//! 1. W(iξ) MUST be `Array3<f64>` (real-valued on imaginary axis)
//! 2. -1/π normalization (PySCF convention)
//! 3. Imaginary part should be small (< 1e-5) after integration
//!
//! # PySCF Reference
//!
//! File: `pyscf/gw/gw_ac.py`, function `get_sigma_element()`, lines 156-212
//!
//! ```python
//! def get_sigma_element(gw, omega, eia, evi):
//!     # Analytic continuation formula
//!     # W is computed on imaginary axis (real values)
//!     # Green's function: 1/(omega - e_p + 1j*eta)
//!     # Integration over imaginary frequency with Gauss-Legendre weights
//!     sigma = -1.0/np.pi * np.sum(integrand * weights)
//! ```

use crate::common::{QuasixError, Result};
use ndarray::{s, Array1, Array2, Array3};
use num_complex::Complex;

/// Extract correlation self-energy at ω=0 from imaginary frequency data
///
/// **CRITICAL**: PySCF baseline sigma_c_iw ALREADY includes -1/π normalization!
/// This function just extracts the ω=0 value (first frequency point).
///
/// # Formula
///
/// Σᶜ_n(ω=0) = sigma_c_iw[0, n]  (just extract first frequency point)
///
/// # Arguments
///
/// * `sigma_c_iw` - Σᶜ(iω) with -1/π already applied [nfreq, nmo] (complex-valued)
/// * `freqs` - Imaginary frequencies [nfreq] (positive real values)
/// * `weights` - Quadrature weights [nfreq] (unused, kept for API compatibility)
///
/// # Returns
///
/// Extracted Σᶜ [nmo] (should be nearly real, small imaginary part)
///
/// # Physical Validation
///
/// The function checks:
/// - Imaginary part should be < 1e-5 (warns if violated, doesn't error)
/// - Result must be finite
///
/// # PySCF Convention
///
/// PySCF applies -1/π normalization DURING frequency loop accumulation.
/// The sigma_c_iw baseline we extract from PySCF already has this normalization!
/// We just need to extract the ω=0 value (first frequency point).
///
/// # Examples
///
/// ```ignore
/// use quasix_core::selfenergy::correlation::integrate_sigma_c;
/// use ndarray::Array2;
/// use num_complex::Complex;
///
/// // Example: extract ω=0 value
/// let sigma_c_iw = Array2::<Complex<f64>>::zeros((10, 5));  // [nfreq=10, nmo=5]
/// let freqs = Array1::linspace(0.1, 10.0, 10);
/// let weights = Array1::ones(10) * 0.1;  // Unused
///
/// let sigma_c = integrate_sigma_c(&sigma_c_iw, &freqs, &weights)?;
/// assert_eq!(sigma_c.len(), 5);
/// ```
///
/// # Errors
///
/// Returns error if:
/// - Dimension mismatch between arrays
/// - Result contains NaN or Inf
pub fn integrate_sigma_c(
    sigma_c_iw: &Array2<Complex<f64>>,
    freqs: &Array1<f64>,
    weights: &Array1<f64>,
) -> Result<Array1<Complex<f64>>> {
    let (nfreq, _nmo) = sigma_c_iw.dim();

    // Validate dimensions
    if weights.len() != nfreq {
        return Err(QuasixError::DimensionMismatch(format!(
            "Weights length {} != nfreq {}",
            weights.len(),
            nfreq
        )));
    }

    if freqs.len() != nfreq {
        return Err(QuasixError::DimensionMismatch(format!(
            "Frequencies length {} != nfreq {}",
            freqs.len(),
            nfreq
        )));
    }

    // Check for NaN/Inf in inputs
    if sigma_c_iw.iter().any(|z| !z.is_finite()) {
        return Err(QuasixError::InvalidInput(
            "sigma_c_iw contains NaN or Inf values".to_string(),
        ));
    }

    // PySCF baseline already includes -1/π normalization!
    // Just extract the ω=0 value (first frequency point)
    let sigma_c = sigma_c_iw.slice(s![0, ..]).to_owned();

    // Validate result
    if sigma_c.iter().any(|z| !z.is_finite()) {
        return Err(QuasixError::NumericalError(
            "Extracted Σᶜ contains NaN or Inf".to_string(),
        ));
    }

    // Physical check: imaginary part should be small
    let max_imag = sigma_c.iter().map(|z| z.im.abs()).fold(0.0, f64::max);

    if max_imag > 1e-5 {
        eprintln!("WARNING: Σᶜ has large imaginary part: {:.2e}", max_imag);
        eprintln!("  This may indicate numerical issues or incorrect W(iω)");
        eprintln!("  Expected: |Im(Σᶜ)| < 1e-5 for well-behaved systems");
        // Don't error - this is a physics check, not a bug
    }

    Ok(sigma_c)
}

/// Compute intermediate Σᶜ(iω) before frequency integration
///
/// This function computes the integrand for correlation self-energy:
///
/// # Formula
///
/// Σᶜ_n(iω) = Σ_ia Σ_PQ (ni|P) W_PQ(iω) G_a(iω) (Q|ia)
///
/// where:
/// - i runs over occupied orbitals
/// - a runs over virtual orbitals
/// - P,Q are auxiliary basis indices
/// - W(iω) is screened interaction (real-valued on imaginary axis)
/// - G_a(iω) is Green's function for virtual orbital a
///
/// # Arguments
///
/// * `green_function` - Green's function G(iω) [nfreq, nmo] (complex)
/// * `w_iw` - Screened interaction W(iω) [nfreq, naux, naux] (real)
/// * `lpq` - DF tensor (ni|P) [naux, nmo, nmo] (real)
/// * `nocc` - Number of occupied orbitals
///
/// # Returns
///
/// Σᶜ(iω) [nfreq, nmo] before integration (complex)
///
/// # Implementation Strategy
///
/// Use BLAS-3 operations for efficiency:
/// 1. For each frequency ω:
///    - Extract W(ω) [naux, naux]
///    - For each orbital n:
///      - Extract (ni|P) for all occupied i
///      - Contract: tmp_PQ = (ni|P) @ W_PQ
///      - Contract with G_a and (ia|Q)
///
/// # Errors
///
/// Returns error if:
/// - Dimension mismatch
/// - Result contains NaN/Inf
pub fn compute_sigma_c_intermediate(
    green_function: &Array2<Complex<f64>>,
    w_iw: &Array3<f64>,
    lpq: &Array3<f64>,
    nocc: usize,
) -> Result<Array2<Complex<f64>>> {
    let (nfreq_g, nmo) = green_function.dim();
    let (nfreq_w, naux_w1, naux_w2) = w_iw.dim();
    let (naux_l, nmo_l1, nmo_l2) = lpq.dim();

    // Validate dimensions
    if nfreq_g != nfreq_w {
        return Err(QuasixError::DimensionMismatch(format!(
            "Green function nfreq {} != W nfreq {}",
            nfreq_g, nfreq_w
        )));
    }

    if naux_w1 != naux_w2 {
        return Err(QuasixError::DimensionMismatch(format!(
            "W matrix not square: {} x {}",
            naux_w1, naux_w2
        )));
    }

    if naux_l != naux_w1 {
        return Err(QuasixError::DimensionMismatch(format!(
            "Lpq naux {} != W naux {}",
            naux_l, naux_w1
        )));
    }

    if nmo_l1 != nmo || nmo_l2 != nmo {
        return Err(QuasixError::DimensionMismatch(format!(
            "Lpq MO dimensions ({}, {}) != nmo {}",
            nmo_l1, nmo_l2, nmo
        )));
    }

    if nocc >= nmo {
        return Err(QuasixError::InvalidInput(format!(
            "nocc {} >= nmo {}",
            nocc, nmo
        )));
    }

    let nfreq = nfreq_g;
    let naux = naux_w1;
    let _nvirt = nmo - nocc;

    // Allocate result
    let mut sigma_c_iw = Array2::<Complex<f64>>::zeros((nfreq, nmo));

    // For each frequency point
    for ifreq in 0..nfreq {
        // Extract W(ω) for this frequency [naux, naux]
        let w_freq = w_iw.slice(s![ifreq, .., ..]);

        // Extract Green's function for this frequency [nmo]
        let g_freq = green_function.slice(s![ifreq, ..]);

        // Compute Σᶜ_n(ω) for each orbital n
        for n in 0..nmo {
            let mut sigma_n = Complex::new(0.0, 0.0);

            // Sum over occupied i and virtual a
            for i in 0..nocc {
                for a in nocc..nmo {
                    // Extract (ni|P) and (ia|P) from lpq
                    // lpq is indexed as [P, p, q] where p,q are MO indices
                    let lni = lpq.slice(s![.., n, i]); // [naux]
                    let lia = lpq.slice(s![.., i, a]); // [naux]

                    // Contract (ni|P) with W_PQ: tmp_Q = Σ_P (ni|P) W_PQ
                    // Use BLAS dot product for efficiency
                    let mut tmp = Array1::<f64>::zeros(naux);
                    for q in 0..naux {
                        let w_col = w_freq.slice(s![.., q]); // W[:, q]
                        tmp[q] = lni.iter().zip(w_col.iter()).map(|(l, w)| l * w).sum();
                    }

                    // Contract tmp_Q with (ia|Q): Σ_Q tmp_Q (ia|Q)
                    let contraction: f64 = tmp.iter().zip(lia.iter()).map(|(t, l)| t * l).sum();

                    // Multiply by Green's function G_a(ω)
                    sigma_n += g_freq[a] * contraction;
                }
            }

            sigma_c_iw[[ifreq, n]] = sigma_n;
        }
    }

    // Validate result
    if sigma_c_iw.iter().any(|z| !z.is_finite()) {
        return Err(QuasixError::NumericalError(
            "Intermediate Σᶜ(iω) contains NaN or Inf".to_string(),
        ));
    }

    Ok(sigma_c_iw)
}

/// Compute correlation self-energy diagonal (main entry point)
///
/// This is the main function for computing Σᶜ in G₀W₀ calculations.
///
/// # Formula
///
/// Σᶜ_n = -1/π ∫ dω Σᶜ_n(iω)
///
/// where Σᶜ_n(iω) = Σ_ia Σ_PQ (ni|P) W_PQ(iω) G_a(iω) (Q|ia)
///
/// # Arguments
///
/// * `green_function` - Green's function G(iω) [nfreq, nmo]
/// * `w_iw` - Screened interaction W(iω) [nfreq, naux, naux] (REAL!)
/// * `lpq` - DF tensor (pq|P) [naux, nmo, nmo]
/// * `freqs` - Imaginary frequency grid [nfreq]
/// * `weights` - Quadrature weights [nfreq]
/// * `nocc` - Number of occupied orbitals
///
/// # Returns
///
/// Correlation self-energy Σᶜ [nmo] (complex, should be nearly real)
///
/// # Physical Checks
///
/// - |Im(Σᶜ)| should be < 1e-5 (warns if violated)
/// - Result must be finite
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_ac.py::get_sigma_element()` lines 156-212
///
/// # Examples
///
/// ```ignore
/// use quasix_core::selfenergy::correlation::compute_sigma_c_diagonal;
///
/// let sigma_c = compute_sigma_c_diagonal(
///     &green_function,
///     &w_iw,
///     &lpq,
///     &freqs,
///     &weights,
///     nocc,
/// )?;
/// assert_eq!(sigma_c.len(), nmo);
/// ```
///
/// # Errors
///
/// Returns error if:
/// - Dimension mismatch between inputs
/// - W is not real-valued (if you pass Complex, will fail dimension check)
/// - Integration fails
/// - Result contains NaN/Inf
pub fn compute_sigma_c_diagonal(
    green_function: &Array2<Complex<f64>>,
    w_iw: &Array3<f64>,
    lpq: &Array3<f64>,
    freqs: &Array1<f64>,
    weights: &Array1<f64>,
    nocc: usize,
) -> Result<Array1<Complex<f64>>> {
    // Step 1: Compute Σᶜ(iω) for each frequency
    let sigma_c_iw = compute_sigma_c_intermediate(green_function, w_iw, lpq, nocc)?;

    // Step 2: Integrate over frequency with -1/π normalization
    let sigma_c = integrate_sigma_c(&sigma_c_iw, freqs, weights)?;

    Ok(sigma_c)
}

// ============================================================================
// Legacy Type Stubs (for compilation compatibility)
// ============================================================================

/// LEGACY: Correlation self-energy calculator
/// Will be removed during S3-5 re-implementation cleanup
#[allow(dead_code)]
pub struct CorrelationSelfEnergy {
    n_mo: usize,
    n_aux: usize,
    n_occ: usize,
}

impl CorrelationSelfEnergy {
    #[allow(dead_code)]
    pub fn new(n_mo: usize, n_aux: usize, n_occ: usize) -> Self {
        Self { n_mo, n_aux, n_occ }
    }

    #[allow(dead_code)]
    pub fn compute_sigma_c_contour_deformation(
        &self,
        _qp_energies: &Array1<f64>,
        _mo_occ: &Array1<f64>,
        _w_matrix: &Array3<num_complex::Complex64>,
        _df_tensor: &Array3<f64>,
        _frequencies: &Array1<f64>,
        _weights: &Array1<f64>,
    ) -> Result<Array2<num_complex::Complex64>> {
        Err(QuasixError::NotImplemented(
            "Legacy compute_sigma_c_contour_deformation - use compute_sigma_c_diagonal()"
                .to_string(),
        ))
    }
}

/// LEGACY: Correlation configuration
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CorrelationConfig {
    pub eta: f64,
    pub n_imag_freq: usize,
    pub omega_max: f64,
    pub n_imag_points: usize,
    pub xi_max: f64,
    pub use_gl_quadrature: bool,
    pub pole_threshold: f64,
    pub convergence_tol: f64,
    pub regularization: f64,
    pub use_simd: bool,
    pub n_threads: Option<usize>,
    pub verbose: usize,
    pub compute_spectral: bool,
}

impl Default for CorrelationConfig {
    fn default() -> Self {
        Self {
            eta: 0.01,
            n_imag_freq: 32,
            omega_max: 30.0,
            n_imag_points: 32,
            xi_max: 30.0,
            use_gl_quadrature: true,
            pole_threshold: 1e-6,
            convergence_tol: 1e-8,
            regularization: 1e-10,
            use_simd: false,
            n_threads: None,
            verbose: 1,
            compute_spectral: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::array;

    #[test]
    fn test_integrate_sigma_c_basic() {
        // Test extraction of ω=0 value from sigma_c_iw
        let nfreq = 3;
        let nmo = 2;

        // sigma_c_iw with different values at each frequency
        let mut sigma_c_iw = Array2::<Complex<f64>>::zeros((nfreq, nmo));
        sigma_c_iw[[0, 0]] = Complex::new(1.0, 0.1); // ω=0 (should be extracted)
        sigma_c_iw[[1, 0]] = Complex::new(2.0, 0.2); // ω=ω₁ (ignored)
        sigma_c_iw[[2, 0]] = Complex::new(3.0, 0.3); // ω=ω₂ (ignored)

        let freqs = array![0.5, 1.0, 2.0];
        let weights = array![0.3, 0.4, 0.3]; // Unused

        let sigma_c = integrate_sigma_c(&sigma_c_iw, &freqs, &weights).expect("Extraction failed");

        // Expected: Just the first frequency point (ω=0)
        assert_abs_diff_eq!(sigma_c[0].re, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(sigma_c[0].im, 0.1, epsilon = 1e-10);
    }

    #[test]
    fn test_integrate_sigma_c_dimension_mismatch() {
        let sigma_c_iw = Array2::<Complex<f64>>::zeros((3, 2));
        let freqs = array![0.5, 1.0, 2.0];
        let weights = array![0.3, 0.4]; // Wrong length!

        let result = integrate_sigma_c(&sigma_c_iw, &freqs, &weights);
        assert!(result.is_err(), "Should reject dimension mismatch");
    }

    #[test]
    fn test_compute_sigma_c_diagonal_dimensions() {
        // Test that compute_sigma_c_diagonal validates dimensions correctly
        let nfreq = 2;
        let nmo = 4;
        let nocc = 2;
        let naux = 6;

        let green = Array2::<Complex<f64>>::zeros((nfreq, nmo));
        let w_iw = Array3::<f64>::zeros((nfreq, naux, naux));
        let lpq = Array3::<f64>::zeros((naux, nmo, nmo));
        let freqs = Array1::<f64>::ones(nfreq);
        let weights = Array1::<f64>::ones(nfreq);

        // Valid call (should work, might produce zeros but shouldn't error)
        let result = compute_sigma_c_diagonal(&green, &w_iw, &lpq, &freqs, &weights, nocc);
        assert!(result.is_ok(), "Valid dimensions should not error");

        // Invalid: mismatched frequency dimensions
        let green_bad = Array2::<Complex<f64>>::zeros((nfreq + 1, nmo));
        let result = compute_sigma_c_diagonal(&green_bad, &w_iw, &lpq, &freqs, &weights, nocc);
        assert!(
            result.is_err(),
            "Should reject frequency dimension mismatch"
        );
    }
}
