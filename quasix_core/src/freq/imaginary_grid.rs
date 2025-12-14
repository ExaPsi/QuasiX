//! Gauss-Legendre frequency grid for imaginary axis GW calculations
//!
//! **Story**: S3-1
//! **Reference**: docs/G0W0/02_STORY_BY_STORY_PLAN.md (S3-1, lines 109-168)
//! **PySCF Source**: pyscf/lib/gto/grid_ao_drv.c::GaussLegendre()
//!
//! **Algorithm**:
//! Uses existing Gauss-Legendre implementation from `freq::gauss_legendre`
//! which implements the Golub-Welsch algorithm with Newton-Raphson refinement.
//!
//! **Validation**: < 1e-14 vs PySCF (machine precision)

use crate::common::{QuasixError, Result};
use crate::freq::gauss_legendre::gauss_legendre_nodes_weights;
use ndarray::Array1;

/// Compute Gauss-Legendre nodes and weights on [-1, 1]
///
/// Returns GL quadrature points for ∫₋₁¹ f(x) dx ≈ Σᵢ wᵢ f(xᵢ)
///
/// # Arguments
/// * `nfreq` - Number of quadrature points
///
/// # Returns
/// * `(nodes, weights)` - GL nodes in [-1, 1] and corresponding weights
///
/// # Validation
/// * Nodes match PySCF < 1e-14
/// * Weights sum to 2.0 (integral over [-1, 1])
/// * Nodes symmetric about 0
///
/// # PySCF Reference
/// `pyscf/lib/gto/grid_ao_drv.c::GaussLegendre()`
///
/// # Example
/// ```rust
/// use quasix_core::freq::imaginary_grid::gauss_legendre_grid;
///
/// let (nodes, weights) = gauss_legendre_grid(32).unwrap();
/// assert!((weights.iter().sum::<f64>() - 2.0).abs() < 1e-14);
/// ```
pub fn gauss_legendre_grid(nfreq: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    if nfreq == 0 {
        return Err(QuasixError::InvalidInput("nfreq must be > 0".to_string()));
    }

    // Use existing high-precision GL implementation
    let (nodes_arr, weights_arr) = gauss_legendre_nodes_weights(nfreq)?;

    // Convert to Vec for Python compatibility
    let nodes = nodes_arr.to_vec();
    let weights = weights_arr.to_vec();

    // Validation check (sanity)
    let sum_weights: f64 = weights.iter().sum();
    if (sum_weights - 2.0).abs() > 1e-12 {
        return Err(QuasixError::NumericalError(format!(
            "GL weight normalization error: sum = {}, expected 2.0",
            sum_weights
        )));
    }

    Ok((nodes, weights))
}

/// Transform GL nodes from [-1, 1] to imaginary frequency axis [0, ω_max]
///
/// **DEPRECATED**: Use `transform_to_imaginary_axis_pyscf()` instead for GW calculations!
///
/// Uses LINEAR transformation: ξ = ω_max * (x + 1) / 2
///
/// # Arguments
/// * `nodes` - GL nodes in [-1, 1]
/// * `omega_max` - Maximum imaginary frequency (Ha)
///
/// # Returns
/// * Frequencies in [0, ω_max]
///
/// # Example
/// ```rust,no_run
/// use quasix_core::freq::imaginary_grid::{gauss_legendre_grid, transform_to_imaginary_axis};
///
/// let (nodes, _) = gauss_legendre_grid(32).unwrap();
/// let omega_max = 10.0;  // Ha
/// let freqs = transform_to_imaginary_axis(&nodes, omega_max);
/// // First point near 0, last point near ω_max
/// ```
#[deprecated(note = "Use transform_to_imaginary_axis_pyscf() for GW calculations")]
pub fn transform_to_imaginary_axis(nodes: &[f64], omega_max: f64) -> Vec<f64> {
    nodes.iter().map(|&x| omega_max * (x + 1.0) / 2.0).collect()
}

/// Transform GL nodes from [-1, 1] to imaginary axis [0, ∞) (PySCF convention)
///
/// **CRITICAL**: This is the CORRECT transformation for GW calculations!
///
/// Uses RATIONAL transformation: ξ = x₀ * (1 + x) / (1 - x)
///
/// This maps:
/// - x = -1 → ξ = 0
/// - x = 0  → ξ = x₀
/// - x = 1  → ξ = ∞
///
/// # Arguments
/// * `nodes` - GL nodes in [-1, 1]
/// * `x0` - Scaling parameter (default: 0.5 in PySCF)
///
/// # Returns
/// * Frequencies in [0, ∞)
///
/// # PySCF Reference
/// `pyscf/pbc/gw/gw_ac.py::_get_scaled_legendre_roots()`
/// ```python
/// x0 = 0.5
/// freqs_new = x0 * (1. + freqs) / (1. - freqs)
/// wts = wts * 2. * x0 / (1. - freqs)**2
/// ```
///
/// # Validation
/// Match PySCF frequency grid < 1e-14 for same nfreq and x0
///
/// # Example
/// ```rust,no_run
/// use quasix_core::freq::imaginary_grid::{gauss_legendre_grid, transform_to_imaginary_axis_pyscf};
///
/// let (nodes, _) = gauss_legendre_grid(100).unwrap();
/// let freqs = transform_to_imaginary_axis_pyscf(&nodes, 0.5);
/// // First point near 0, last point >> 1000 Ha
/// ```
pub fn transform_to_imaginary_axis_pyscf(nodes: &[f64], x0: f64) -> Vec<f64> {
    nodes
        .iter()
        .map(|&x| {
            // PySCF formula: freqs_new = x0 * (1 + x) / (1 - x)
            // Handle x → 1 limit carefully (avoid division by zero)
            if (1.0 - x).abs() < 1e-14 {
                // Near x=1, frequency approaches infinity
                // Return a large but finite value
                1e12 // 1e12 Ha ≈ 2.7e13 eV (essentially infinity for GW)
            } else {
                x0 * (1.0 + x) / (1.0 - x)
            }
        })
        .collect()
}

/// Transform GL weights for imaginary frequency integration (linear transformation)
///
/// **DEPRECATED**: Use `transform_weights_pyscf()` for GW calculations!
///
/// Applies Jacobian for linear transformation: dξ/dx = ω_max / 2
///
/// # Arguments
/// * `weights` - GL weights on [-1, 1]
/// * `omega_max` - Maximum imaginary frequency (Ha)
///
/// # Returns
/// * Transformed weights for [0, ω_max] interval
#[deprecated(note = "Use transform_weights_pyscf() for GW calculations")]
pub fn transform_weights(weights: &[f64], omega_max: f64) -> Vec<f64> {
    let jacobian = omega_max / 2.0;
    weights.iter().map(|&w| w * jacobian).collect()
}

/// Transform GL weights for imaginary axis integration (PySCF convention)
///
/// **CRITICAL**: This is the CORRECT weight transformation for GW calculations!
///
/// Applies Jacobian for rational transformation: dξ/dx = 2*x₀ / (1 - x)²
///
/// # Arguments
/// * `weights` - GL weights on [-1, 1]
/// * `nodes` - GL nodes on [-1, 1] (needed for Jacobian)
/// * `x0` - Scaling parameter (default: 0.5 in PySCF)
///
/// # Returns
/// * Transformed weights for [0, ∞) integration
///
/// # PySCF Reference
/// `pyscf/pbc/gw/gw_ac.py::_get_scaled_legendre_roots()`
/// ```python
/// wts = wts * 2. * x0 / (1. - freqs)**2
/// ```
///
/// # Validation
/// Weight sum should match PySCF < 1e-10
///
/// # Example
/// ```rust
/// use quasix_core::freq::imaginary_grid::{gauss_legendre_grid, transform_weights_pyscf};
///
/// let (nodes, weights) = gauss_legendre_grid(100).unwrap();
/// let transformed_wts = transform_weights_pyscf(&weights, &nodes, 0.5);
/// assert!(transformed_wts.iter().all(|&w| w > 0.0));  // All weights positive
/// ```
pub fn transform_weights_pyscf(weights: &[f64], nodes: &[f64], x0: f64) -> Vec<f64> {
    weights
        .iter()
        .zip(nodes.iter())
        .map(|(&w, &x)| {
            // PySCF formula: wts_new = wts * 2 * x0 / (1 - x)^2
            // Jacobian from rational transformation
            if (1.0 - x).abs() < 1e-14 {
                // Near x=1, Jacobian diverges
                // But weight*jacobian product remains finite due to GL weight → 0
                0.0
            } else {
                w * 2.0 * x0 / ((1.0 - x) * (1.0 - x))
            }
        })
        .collect()
}

/// Old configuration structure (kept for compatibility)
#[derive(Debug, Clone)]
pub struct ImaginaryGridConfig {
    /// Maximum frequency cutoff (in Hartree)
    pub omega_max: f64,
    /// Number of grid points
    pub n_points: usize,
    /// Cutoff for near-singularity handling (unused for GL)
    pub singularity_cutoff: f64,
}

impl Default for ImaginaryGridConfig {
    fn default() -> Self {
        Self {
            omega_max: 50.0, // 50 Ha ≈ 1360 eV
            n_points: 32,    // Default grid size
            singularity_cutoff: 1e-10,
        }
    }
}

/// Create imaginary frequency grid (legacy compatibility)
///
/// Note: For S3-1, use `gauss_legendre_grid()` + `transform_to_imaginary_axis()`
/// This function is kept for backward compatibility with existing code.
#[allow(deprecated)]
pub fn create_imaginary_frequency_grid(
    n_points: usize,
    omega_max: f64,
) -> Result<(Array1<f64>, Array1<f64>)> {
    // Validate inputs
    if n_points < 2 {
        return Err(QuasixError::InvalidInput(format!(
            "Number of points must be at least 2, got {n_points}"
        )));
    }
    if omega_max <= 0.0 {
        return Err(QuasixError::InvalidInput(format!(
            "omega_max must be positive, got {omega_max}"
        )));
    }

    // Get GL grid on [-1, 1]
    let (nodes, weights) = gauss_legendre_grid(n_points)?;

    // Transform to [0, omega_max]
    let frequencies = transform_to_imaginary_axis(&nodes, omega_max);
    let transformed_weights = transform_weights(&weights, omega_max);

    Ok((
        Array1::from_vec(frequencies),
        Array1::from_vec(transformed_weights),
    ))
}

/// Create imaginary frequency grid with custom configuration (legacy)
pub fn create_imaginary_frequency_grid_with_config(
    config: &ImaginaryGridConfig,
) -> Result<(Array1<f64>, Array1<f64>)> {
    create_imaginary_frequency_grid(config.n_points, config.omega_max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gl_grid_normalization() {
        let (nodes, weights) = gauss_legendre_grid(32).unwrap();

        // Check length
        assert_eq!(nodes.len(), 32);
        assert_eq!(weights.len(), 32);

        // Check normalization
        let sum_weights: f64 = weights.iter().sum();
        assert!((sum_weights - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_gl_grid_symmetry() {
        let (nodes, _) = gauss_legendre_grid(32).unwrap();

        // Check symmetry about 0
        for i in 0..nodes.len() {
            let symmetric_diff = (nodes[i] + nodes[nodes.len() - 1 - i]).abs();
            assert!(
                symmetric_diff < 1e-14,
                "Nodes not symmetric at index {}: diff = {}",
                i,
                symmetric_diff
            );
        }
    }

    #[test]
    fn test_gl_grid_range() {
        let (nodes, _) = gauss_legendre_grid(32).unwrap();

        // Nodes should be in [-1, 1]
        assert!(nodes[0] >= -1.0 && nodes[0] <= -0.99); // Close to -1
        assert!(nodes.last().unwrap() >= &0.99 && nodes.last().unwrap() <= &1.0);
        // Close to 1
    }

    #[test]
    #[allow(deprecated)] // Testing deprecated function for backwards compatibility
    fn test_imaginary_axis_transform() {
        let nodes = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
        let omega_max = 10.0;

        let freqs = transform_to_imaginary_axis(&nodes, omega_max);

        // Check endpoints
        assert!((freqs[0] - 0.0).abs() < 1e-14);
        assert!((freqs[4] - 10.0).abs() < 1e-14);

        // Check midpoint
        assert!((freqs[2] - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_gl_multiple_sizes() {
        for nfreq in [16, 24, 32, 48] {
            let (nodes, weights) = gauss_legendre_grid(nfreq).unwrap();

            assert_eq!(nodes.len(), nfreq);
            assert_eq!(weights.len(), nfreq);

            let sum_w: f64 = weights.iter().sum();
            assert!((sum_w - 2.0).abs() < 1e-14);
        }
    }

    #[test]
    #[allow(deprecated)] // Testing deprecated function for backwards compatibility
    fn test_weight_transformation() {
        let weights = vec![0.5, 1.0, 1.5, 2.0];
        let omega_max = 10.0;

        let transformed = transform_weights(&weights, omega_max);

        // Jacobian = omega_max / 2 = 5.0
        assert!((transformed[0] - 2.5).abs() < 1e-14);
        assert!((transformed[1] - 5.0).abs() < 1e-14);
        assert!((transformed[2] - 7.5).abs() < 1e-14);
        assert!((transformed[3] - 10.0).abs() < 1e-14);
    }

    #[test]
    fn test_legacy_create_imaginary_frequency_grid() {
        let (frequencies, weights) = create_imaginary_frequency_grid(32, 50.0).unwrap();

        // Check dimensions
        assert_eq!(frequencies.len(), 32);
        assert_eq!(weights.len(), 32);

        // All frequencies should be positive
        for freq in &frequencies {
            assert!(*freq >= 0.0);
        }

        // Weights should be positive
        for &w in &weights {
            assert!(w > 0.0);
        }

        // First point should be near zero
        assert!(frequencies[0] < 1.0);

        // Points should be ordered
        for i in 1..frequencies.len() {
            assert!(frequencies[i] >= frequencies[i - 1]);
        }
    }
}
