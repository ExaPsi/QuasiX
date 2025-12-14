//! Optical properties from BSE-TDA calculations
//!
//! Computes oscillator strengths and absorption spectra from
//! BSE-TDA eigenvalues and eigenvectors.
//!
//! # Oscillator Strength Formula
//!
//! ```text
//! f_n = (2/3) * E_n * sum_alpha |d_n^alpha|^2
//! d_n^alpha = sum_ia X_n(ia) * mu_{ia}^alpha
//! ```
//!
//! where:
//! - `E_n` is the excitation energy of state n (Hartree)
//! - `X_n(ia)` is the BSE eigenvector component for transition i->a
//! - `mu_{ia}^alpha` is the transition dipole integral <i|r_alpha|a>
//!
//! # Absorption Spectrum
//!
//! The broadened absorption spectrum is computed as:
//!
//! **Lorentzian** (natural linewidth):
//! ```text
//! sigma(omega) = sum_n f_n * gamma / ((omega - E_n)^2 + gamma^2)
//! ```
//!
//! **Gaussian** (inhomogeneous broadening):
//! ```text
//! sigma(omega) = sum_n f_n * exp(-(omega - E_n)^2 / (2*sigma^2))
//! ```
//!
//! # References
//!
//! - Rohlfing-Louie (2000): Optical excitations in conjugated polymers
//! - Blase et al. (2011): GW calculations for molecules
//! - docs/derivations/s6-1/theory.md (Section 6.2: Oscillator Strengths)

use crate::common::{QuasixError, Result};
use ndarray::{Array1, Array2, ArrayView1, ArrayView2};

/// Module version for compatibility tracking
pub const MODULE_VERSION: &str = "0.6.3";

/// Broadening type for absorption spectra
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BroadeningType {
    /// Lorentzian broadening (natural linewidth)
    /// L(x) = gamma / (x^2 + gamma^2)
    #[default]
    Lorentzian,

    /// Gaussian broadening (inhomogeneous)
    /// G(x) = exp(-x^2 / (2*sigma^2))
    Gaussian,
}

/// Configuration for optical properties calculation
#[derive(Debug, Clone)]
pub struct OpticalConfig {
    /// Broadening type for spectrum
    pub broadening_type: BroadeningType,

    /// Broadening parameter (FWHM in Hartree)
    /// Default: 0.1 eV (~0.0037 Ha)
    pub broadening: f64,

    /// Energy range for spectrum [min, max] in Hartree
    /// If None, auto-detect from eigenvalues with padding
    pub energy_range: Option<(f64, f64)>,

    /// Number of points in spectrum
    pub n_points: usize,
}

impl Default for OpticalConfig {
    fn default() -> Self {
        Self {
            broadening_type: BroadeningType::Lorentzian,
            broadening: 0.1 / 27.2114, // 0.1 eV in Hartree
            energy_range: None,        // Auto-detect from eigenvalues
            n_points: 1000,
        }
    }
}

/// Result of oscillator strength calculation
#[derive(Debug, Clone)]
pub struct OpticalResult {
    /// Oscillator strengths for each root [n_roots]
    pub oscillator_strengths: Array1<f64>,

    /// Transition dipole moments [n_roots, 3] for x, y, z components
    pub transition_dipoles: Array2<f64>,

    /// Sum of oscillator strengths (for Thomas-Reiche-Kuhn sum rule check)
    pub f_sum: f64,
}

/// Result of absorption spectrum calculation
#[derive(Debug, Clone)]
pub struct SpectrumResult {
    /// Energy grid (Hartree) [n_points]
    pub energies: Array1<f64>,

    /// Absorption intensity (arbitrary units) [n_points]
    pub intensities: Array1<f64>,

    /// Broadening parameter used (Hartree)
    pub broadening: f64,

    /// Broadening type used
    pub broadening_type: BroadeningType,
}

/// Lorentzian lineshape function
///
/// L(x; gamma) = gamma / (x^2 + gamma^2)
///
/// This gives a peak of height 1/gamma at x=0 with FWHM = 2*gamma.
#[inline]
fn lorentzian(x: f64, gamma: f64) -> f64 {
    gamma / (x * x + gamma * gamma)
}

/// Gaussian lineshape function
///
/// G(x; sigma) = exp(-x^2 / (2*sigma^2))
///
/// Normalized Gaussian with peak at x=0.
#[inline]
fn gaussian(x: f64, sigma: f64) -> f64 {
    (-0.5 * x * x / (sigma * sigma)).exp()
}

/// Compute oscillator strengths from BSE-TDA eigenvectors
///
/// # Formula
///
/// ```text
/// f_n = (2/3) * E_n * sum_alpha |d_n^alpha|^2
/// d_n^alpha = sum_ia X_n(ia) * mu_{ia}^alpha
/// ```
///
/// # Arguments
///
/// * `eigenvalues` - Excitation energies from Davidson solver [n_roots] (Hartree)
/// * `eigenvectors` - BSE eigenvectors X_n(ia), shape [n_trans, n_roots]
/// * `transition_dipoles` - <i|r_alpha|a> integrals, shape [n_trans, 3] for x, y, z
///
/// # Returns
///
/// `OpticalResult` containing oscillator strengths and transition dipoles
///
/// # Errors
///
/// Returns error if:
/// - Dimensions are inconsistent (eigenvectors rows != transition_dipoles rows)
/// - Transition dipoles does not have 3 columns
/// - Eigenvalues contain non-finite or negative values
///
/// # Example
///
/// ```ignore
/// use ndarray::{arr1, arr2, Array2};
/// use quasix_core::bse::optical::compute_oscillator_strengths;
///
/// let eigenvalues = arr1(&[0.3, 0.5]);  // Ha
/// let eigenvectors = Array2::eye(2);    // Identity (simple case)
/// let transition_dipoles = arr2(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
///
/// let result = compute_oscillator_strengths(
///     eigenvalues.view(),
///     eigenvectors.view(),
///     transition_dipoles.view()
/// ).unwrap();
///
/// assert!(result.oscillator_strengths[0] > 0.0);
/// ```
pub fn compute_oscillator_strengths(
    eigenvalues: ArrayView1<f64>,
    eigenvectors: ArrayView2<f64>,
    transition_dipoles: ArrayView2<f64>,
) -> Result<OpticalResult> {
    // Validate dimensions
    let n_roots = eigenvalues.len();
    let (n_trans_ev, n_roots_ev) = eigenvectors.dim();
    let (n_trans_dip, n_comp) = transition_dipoles.dim();

    if n_roots_ev != n_roots {
        return Err(QuasixError::DimensionMismatch(format!(
            "Eigenvalues length {} != eigenvectors columns {}",
            n_roots, n_roots_ev
        )));
    }

    if n_trans_ev != n_trans_dip {
        return Err(QuasixError::DimensionMismatch(format!(
            "Eigenvectors rows {} != transition_dipoles rows {}",
            n_trans_ev, n_trans_dip
        )));
    }

    if n_comp != 3 {
        return Err(QuasixError::DimensionMismatch(format!(
            "Transition dipoles must have 3 components (x,y,z), got {}",
            n_comp
        )));
    }

    // Validate eigenvalues are finite and non-negative
    for (i, &e) in eigenvalues.iter().enumerate() {
        if !e.is_finite() {
            return Err(QuasixError::NumericalError(format!(
                "Non-finite eigenvalue at index {}: {}",
                i, e
            )));
        }
        if e < 0.0 {
            return Err(QuasixError::NumericalError(format!(
                "Negative eigenvalue at index {}: {} (BSE-TDA should give positive energies)",
                i, e
            )));
        }
    }

    // Compute transition dipoles for each root
    // d[n, alpha] = sum_ia X[ia, n] * mu[ia, alpha]
    // eigenvectors: [n_trans, n_roots]
    // transition_dipoles: [n_trans, 3]
    // Result: [n_roots, 3]
    let d = eigenvectors.t().dot(&transition_dipoles);

    // Compute oscillator strengths
    // f_n = (2/3) * E_n * sum_alpha |d_n[alpha]|^2
    let mut oscillator_strengths = Array1::zeros(n_roots);

    for n in 0..n_roots {
        let d_squared: f64 = d.row(n).iter().map(|&x| x * x).sum();
        oscillator_strengths[n] = (2.0 / 3.0) * eigenvalues[n] * d_squared;
    }

    // Compute sum for Thomas-Reiche-Kuhn sum rule check
    // Sum should approximately equal N_electrons for well-converged calculations
    let f_sum = oscillator_strengths.sum();

    Ok(OpticalResult {
        oscillator_strengths,
        transition_dipoles: d,
        f_sum,
    })
}

/// Compute absorption spectrum with broadening
///
/// # Formula
///
/// **Lorentzian**:
/// ```text
/// sigma(omega) = sum_n f_n * gamma / ((omega - E_n)^2 + gamma^2)
/// ```
///
/// **Gaussian**:
/// ```text
/// sigma(omega) = sum_n f_n * exp(-(omega - E_n)^2 / (2*sigma^2))
/// ```
///
/// # Arguments
///
/// * `eigenvalues` - Excitation energies [n_roots] (Hartree)
/// * `oscillator_strengths` - Oscillator strengths [n_roots]
/// * `config` - Spectrum configuration (broadening type, width, energy range, grid size)
///
/// # Returns
///
/// `SpectrumResult` containing energy grid and absorption intensities
///
/// # Errors
///
/// Returns error if eigenvalues and oscillator_strengths have different lengths.
///
/// # Example
///
/// ```ignore
/// use ndarray::arr1;
/// use quasix_core::bse::optical::{compute_absorption_spectrum, OpticalConfig, BroadeningType};
///
/// let eigenvalues = arr1(&[0.2, 0.4]);  // Ha
/// let oscillator_strengths = arr1(&[0.5, 0.3]);
///
/// let config = OpticalConfig {
///     broadening_type: BroadeningType::Lorentzian,
///     broadening: 0.01,
///     energy_range: Some((0.1, 0.5)),
///     n_points: 100,
/// };
///
/// let result = compute_absorption_spectrum(
///     eigenvalues.view(),
///     oscillator_strengths.view(),
///     &config
/// ).unwrap();
///
/// assert_eq!(result.energies.len(), 100);
/// ```
pub fn compute_absorption_spectrum(
    eigenvalues: ArrayView1<f64>,
    oscillator_strengths: ArrayView1<f64>,
    config: &OpticalConfig,
) -> Result<SpectrumResult> {
    let n_roots = eigenvalues.len();

    if oscillator_strengths.len() != n_roots {
        return Err(QuasixError::DimensionMismatch(format!(
            "Eigenvalues length {} != oscillator_strengths length {}",
            n_roots,
            oscillator_strengths.len()
        )));
    }

    // Determine energy range
    let (e_min, e_max) = match config.energy_range {
        Some((min, max)) => (min, max),
        None => {
            // Auto-detect: extend beyond eigenvalue range
            let eig_min = eigenvalues.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let eig_max = eigenvalues.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            let range = eig_max - eig_min;
            let padding = (range * 0.2).max(config.broadening * 5.0);
            ((eig_min - padding).max(0.0), eig_max + padding)
        }
    };

    // Generate energy grid
    let n_points = config.n_points;
    let de = if n_points > 1 {
        (e_max - e_min) / (n_points - 1) as f64
    } else {
        0.0
    };

    let energies: Array1<f64> = (0..n_points).map(|i| e_min + i as f64 * de).collect();

    // Compute spectrum
    let mut intensities = Array1::zeros(n_points);
    let gamma = config.broadening;

    for i in 0..n_points {
        let e = energies[i];
        let mut intensity = 0.0;

        for n in 0..n_roots {
            let x = e - eigenvalues[n];
            let lineshape = match config.broadening_type {
                BroadeningType::Lorentzian => lorentzian(x, gamma),
                BroadeningType::Gaussian => gaussian(x, gamma),
            };
            intensity += oscillator_strengths[n] * lineshape;
        }

        intensities[i] = intensity;
    }

    Ok(SpectrumResult {
        energies,
        intensities,
        broadening: config.broadening,
        broadening_type: config.broadening_type,
    })
}

/// Convenience function to compute full optical properties
///
/// Combines [`compute_oscillator_strengths`] and [`compute_absorption_spectrum`]
/// into a single call for convenience.
///
/// # Arguments
///
/// * `eigenvalues` - Excitation energies [n_roots] (Hartree)
/// * `eigenvectors` - BSE eigenvectors [n_trans, n_roots]
/// * `transition_dipoles` - <i|r|a> integrals [n_trans, 3]
/// * `config` - Spectrum configuration
///
/// # Returns
///
/// Tuple of (OpticalResult, SpectrumResult)
///
/// # Errors
///
/// Returns error if any validation fails in the underlying functions.
pub fn compute_optical_properties(
    eigenvalues: ArrayView1<f64>,
    eigenvectors: ArrayView2<f64>,
    transition_dipoles: ArrayView2<f64>,
    config: &OpticalConfig,
) -> Result<(OpticalResult, SpectrumResult)> {
    let optical = compute_oscillator_strengths(eigenvalues, eigenvectors, transition_dipoles)?;
    let spectrum =
        compute_absorption_spectrum(eigenvalues, optical.oscillator_strengths.view(), config)?;

    Ok((optical, spectrum))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::{arr1, arr2, Array2};

    #[test]
    fn test_oscillator_strength_simple() {
        // Simple 2-transition system
        // eigenvalues: [0.3, 0.5] Ha
        // eigenvectors: identity (each root = one transition)
        // transition_dipoles: [[1, 0, 0], [0, 1, 0]]

        let eigenvalues = arr1(&[0.3, 0.5]);
        let eigenvectors = Array2::eye(2);
        let transition_dipoles = arr2(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);

        let result = compute_oscillator_strengths(
            eigenvalues.view(),
            eigenvectors.view(),
            transition_dipoles.view(),
        )
        .unwrap();

        // f_0 = (2/3) * 0.3 * 1^2 = 0.2
        // f_1 = (2/3) * 0.5 * 1^2 = 0.333...
        assert_relative_eq!(result.oscillator_strengths[0], 0.2, epsilon = 1e-10);
        assert_relative_eq!(result.oscillator_strengths[1], 1.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_oscillator_strength_normalization() {
        // Test with non-unit eigenvector
        let eigenvalues = arr1(&[0.4]);
        let eigenvectors = arr2(&[[0.6], [0.8]]); // normalized: |v|^2 = 1
        let transition_dipoles = arr2(&[[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);

        let result = compute_oscillator_strengths(
            eigenvalues.view(),
            eigenvectors.view(),
            transition_dipoles.view(),
        )
        .unwrap();

        // d_x = 0.6*1.0 + 0.8*0.0 = 0.6
        // d_z = 0.6*0.0 + 0.8*1.0 = 0.8
        // |d|^2 = 0.36 + 0.64 = 1.0
        // f = (2/3) * 0.4 * 1.0 = 0.2667
        assert_relative_eq!(
            result.oscillator_strengths[0],
            0.4 * 2.0 / 3.0,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_oscillator_strengths_non_negative() {
        // Oscillator strengths must be >= 0 for positive eigenvalues
        let eigenvalues = arr1(&[0.1, 0.2, 0.3]);
        let eigenvectors = Array2::from_shape_fn((3, 3), |(i, j)| if i == j { 0.8 } else { 0.1 });
        let transition_dipoles =
            Array2::from_shape_fn((3, 3), |(i, j)| ((i + j) as f64 * 0.1).sin());

        let result = compute_oscillator_strengths(
            eigenvalues.view(),
            eigenvectors.view(),
            transition_dipoles.view(),
        )
        .unwrap();

        for &f in &result.oscillator_strengths {
            assert!(f >= 0.0, "Oscillator strength must be non-negative");
        }
    }

    #[test]
    fn test_absorption_spectrum_lorentzian() {
        let eigenvalues = arr1(&[0.2, 0.4]);
        let oscillator_strengths = arr1(&[0.5, 0.3]);

        let config = OpticalConfig {
            broadening_type: BroadeningType::Lorentzian,
            broadening: 0.01, // Small broadening for distinct peaks
            energy_range: Some((0.1, 0.5)),
            n_points: 100,
        };

        let result =
            compute_absorption_spectrum(eigenvalues.view(), oscillator_strengths.view(), &config)
                .unwrap();

        assert_eq!(result.energies.len(), 100);
        assert_eq!(result.intensities.len(), 100);

        // Check peaks near eigenvalues
        let idx_peak1 = result
            .energies
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                ((**a - 0.2).abs())
                    .partial_cmp(&((**b - 0.2).abs()))
                    .unwrap()
            })
            .map(|(i, _)| i)
            .unwrap();

        // Peak should be at maximum near eigenvalue
        assert!(result.intensities[idx_peak1] > result.intensities[idx_peak1.saturating_sub(5)]);
    }

    #[test]
    fn test_absorption_spectrum_gaussian() {
        let eigenvalues = arr1(&[0.3]);
        let oscillator_strengths = arr1(&[1.0]);

        let config = OpticalConfig {
            broadening_type: BroadeningType::Gaussian,
            broadening: 0.02,
            energy_range: Some((0.1, 0.5)),
            n_points: 200,
        };

        let result =
            compute_absorption_spectrum(eigenvalues.view(), oscillator_strengths.view(), &config)
                .unwrap();

        // Peak should be at E = 0.3
        let max_idx = result
            .intensities
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();

        assert_relative_eq!(result.energies[max_idx], 0.3, epsilon = 0.01);
    }

    #[test]
    fn test_dimension_validation() {
        let eigenvalues = arr1(&[0.3, 0.4]);
        let eigenvectors = Array2::eye(3); // Wrong dimension: 3x3 instead of 2x2
        let transition_dipoles = arr2(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);

        let result = compute_oscillator_strengths(
            eigenvalues.view(),
            eigenvectors.view(),
            transition_dipoles.view(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_transition_dipole_shape() {
        let eigenvalues = arr1(&[0.3]);
        let eigenvectors = arr2(&[[1.0]]);
        let bad_dipoles = arr2(&[[1.0, 0.0]]); // Should be [n_trans, 3], not [n_trans, 2]

        let result = compute_oscillator_strengths(
            eigenvalues.view(),
            eigenvectors.view(),
            bad_dipoles.view(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_negative_eigenvalue_error() {
        let eigenvalues = arr1(&[-0.3]); // Negative - should error
        let eigenvectors = arr2(&[[1.0]]);
        let transition_dipoles = arr2(&[[1.0, 0.0, 0.0]]);

        let result = compute_oscillator_strengths(
            eigenvalues.view(),
            eigenvectors.view(),
            transition_dipoles.view(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_f_sum_positive() {
        // Sum of oscillator strengths should be positive for physical systems
        let eigenvalues = arr1(&[0.2, 0.3, 0.4]);
        let eigenvectors = Array2::eye(3);
        let transition_dipoles = arr2(&[[1.0, 0.5, 0.0], [0.5, 1.0, 0.3], [0.3, 0.2, 1.0]]);

        let result = compute_oscillator_strengths(
            eigenvalues.view(),
            eigenvectors.view(),
            transition_dipoles.view(),
        )
        .unwrap();

        assert!(result.f_sum > 0.0);
        assert!(result.f_sum.is_finite());
    }

    #[test]
    fn test_optical_config_default() {
        let config = OpticalConfig::default();
        assert_eq!(config.broadening_type, BroadeningType::Lorentzian);
        assert_eq!(config.n_points, 1000);
        assert!(config.energy_range.is_none());
        // 0.1 eV in Hartree
        assert_relative_eq!(config.broadening, 0.1 / 27.2114, epsilon = 1e-10);
    }

    #[test]
    fn test_compute_optical_properties() {
        let eigenvalues = arr1(&[0.3, 0.4]);
        let eigenvectors = Array2::eye(2);
        let transition_dipoles = arr2(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let config = OpticalConfig {
            n_points: 50,
            ..OpticalConfig::default()
        };

        let (optical, spectrum) = compute_optical_properties(
            eigenvalues.view(),
            eigenvectors.view(),
            transition_dipoles.view(),
            &config,
        )
        .unwrap();

        // Check optical result
        assert_eq!(optical.oscillator_strengths.len(), 2);
        assert_eq!(optical.transition_dipoles.dim(), (2, 3));

        // Check spectrum result
        assert_eq!(spectrum.energies.len(), 50);
        assert_eq!(spectrum.intensities.len(), 50);
    }

    #[test]
    fn test_lorentzian_function() {
        // At peak (x=0), L(0; gamma) = 1/gamma
        assert_relative_eq!(lorentzian(0.0, 0.1), 10.0, epsilon = 1e-10);

        // At x = gamma, L(gamma; gamma) = gamma / (2*gamma^2) = 1/(2*gamma)
        assert_relative_eq!(lorentzian(0.1, 0.1), 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gaussian_function() {
        // At peak (x=0), G(0; sigma) = 1
        assert_relative_eq!(gaussian(0.0, 0.1), 1.0, epsilon = 1e-10);

        // At x = sigma, G(sigma; sigma) = exp(-0.5)
        assert_relative_eq!(gaussian(0.1, 0.1), (-0.5_f64).exp(), epsilon = 1e-10);
    }

    #[test]
    fn test_spectrum_dimension_mismatch() {
        let eigenvalues = arr1(&[0.3, 0.4]);
        let oscillator_strengths = arr1(&[0.5]); // Wrong length

        let config = OpticalConfig::default();

        let result =
            compute_absorption_spectrum(eigenvalues.view(), oscillator_strengths.view(), &config);

        assert!(result.is_err());
    }

    #[test]
    fn test_broadening_type_equality() {
        assert_eq!(BroadeningType::Lorentzian, BroadeningType::Lorentzian);
        assert_eq!(BroadeningType::Gaussian, BroadeningType::Gaussian);
        assert_ne!(BroadeningType::Lorentzian, BroadeningType::Gaussian);
    }
}
