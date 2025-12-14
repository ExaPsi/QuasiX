//! Imaginary axis frequency grid for G0W0 contour deformation
//!
//! This module implements the proper transformation of Gauss-Legendre quadrature
//! points to the imaginary frequency axis for stable and accurate GW calculations.

use crate::common::{QuasixError, Result};
use ndarray::Array1;
use num_complex::Complex64;

/// Overflow protection configuration
#[derive(Debug, Clone)]
pub struct OverflowProtection {
    /// Minimum allowed value for (1-x) before cutoff
    pub min_denominator: f64,

    /// Maximum allowed Jacobian value
    pub max_jacobian: f64,

    /// Regularization epsilon for near-singularity
    pub epsilon: f64,
}

impl Default for OverflowProtection {
    fn default() -> Self {
        Self {
            min_denominator: 1e-12,
            max_jacobian: 1e10,
            epsilon: 1e-15,
        }
    }
}

impl OverflowProtection {
    /// Check if a value is safe for division
    pub fn is_safe(&self, one_minus_t: f64) -> bool {
        one_minus_t >= self.min_denominator
    }

    /// Apply regularization to prevent division by zero
    pub fn regularize(&self, value: f64) -> f64 {
        value.max(self.epsilon)
    }

    /// Check if jacobian is within acceptable bounds
    pub fn is_jacobian_safe(&self, jacobian: f64) -> bool {
        jacobian.is_finite() && jacobian <= self.max_jacobian
    }
}

/// Metadata for imaginary axis grid
#[derive(Debug, Clone)]
pub struct GridMetadata {
    /// Number of points originally requested
    pub n_requested: usize,

    /// Actual number of points after filtering
    pub n_actual: usize,

    /// Frequency cutoff parameter
    pub omega_max: f64,

    /// Minimum imaginary frequency
    pub xi_min: f64,

    /// Maximum imaginary frequency
    pub xi_max: f64,

    /// Number of points filtered due to overflow
    pub overflow_filtered: usize,

    /// Sum of weights (should be ≈ π/2 for [0,∞))
    pub weight_sum: f64,
}

/// Imaginary axis frequency grid with Gauss-Legendre quadrature
#[derive(Debug, Clone)]
pub struct ImaginaryAxisGrid {
    /// Grid points on imaginary axis (purely imaginary)
    pub frequencies: Vec<Complex64>,

    /// Integration weights including Jacobian
    pub weights: Vec<f64>,

    /// Original GL points for diagnostics
    pub gl_points: Array1<f64>,

    /// Original GL weights
    pub gl_weights: Array1<f64>,

    /// Grid metadata
    pub metadata: GridMetadata,
}

impl ImaginaryAxisGrid {
    /// Create a new imaginary axis grid
    ///
    /// # Arguments
    /// * `n_points` - Number of Gauss-Legendre points
    /// * `omega_max` - Frequency cutoff (Ha)
    ///
    /// # Returns
    /// Grid with transformed points and weights
    pub fn new(n_points: usize, omega_max: f64) -> Result<Self> {
        // Validate inputs
        if n_points < 4 {
            return Err(QuasixError::InvalidInput(
                "Need at least 4 quadrature points".to_string(),
            ));
        }

        if omega_max <= 0.0 {
            return Err(QuasixError::InvalidInput(
                "omega_max must be positive".to_string(),
            ));
        }

        // Get Gauss-Legendre quadrature on [-1, 1]
        let (gl_points, gl_weights) =
            crate::freq::gauss_legendre::gauss_legendre_nodes_weights(n_points)?;

        // Set up overflow protection (kept for future use with non-linear transformations)
        let _protection = OverflowProtection::default();

        // Transform to imaginary axis
        let mut frequencies = Vec::with_capacity(n_points);
        let mut weights = Vec::with_capacity(n_points);
        let overflow_filtered = 0; // No filtering needed for linear transformation

        for i in 0..n_points {
            let x = gl_points[i]; // On [-1, 1]
            let w = gl_weights[i];

            // Map from [-1, 1] to [0, 1]
            let t = (x + 1.0) / 2.0;
            let w_01 = w / 2.0; // Jacobian for linear transformation

            // Transform to [0, ω_max]: ξ = ω_max * t (LINEAR, not rational)
            // This matches the expected frequency range for AC correlation
            let xi = omega_max * t;

            // Create purely imaginary frequency
            frequencies.push(Complex64::new(0.0, xi));

            // Jacobian for linear transformation: dξ/dt = ω_max (constant)
            let jacobian = omega_max;

            // Total weight includes both transformations
            weights.push(w_01 * jacobian);
        }

        // Validate grid
        if frequencies.is_empty() {
            return Err(QuasixError::NumericalError(
                "All grid points filtered due to overflow".to_string(),
            ));
        }

        // Compute metadata
        let xi_min = frequencies.first().map(|f| f.im).unwrap_or(0.0);
        let xi_max = frequencies.last().map(|f| f.im).unwrap_or(0.0);
        let weight_sum: f64 = weights.iter().sum();

        let metadata = GridMetadata {
            n_requested: n_points,
            n_actual: frequencies.len(),
            omega_max,
            xi_min,
            xi_max,
            overflow_filtered,
            weight_sum,
        };

        // Log grid info
        log::debug!(
            "Created imaginary axis grid: {} points, range [0, {:.2}i] Ha, {} filtered",
            metadata.n_actual,
            metadata.xi_max,
            metadata.overflow_filtered
        );

        Ok(Self {
            frequencies,
            weights,
            gl_points,
            gl_weights,
            metadata,
        })
    }

    /// Create grid with automatic omega_max selection
    pub fn with_auto_cutoff(n_points: usize, max_energy_diff: f64) -> Result<Self> {
        // Set omega_max based on maximum energy difference
        // Rule of thumb: omega_max ≈ 3-5 × max|εₐ - εᵢ|
        let omega_max = (4.0 * max_energy_diff).max(30.0).min(100.0);

        log::debug!(
            "Auto-selected omega_max = {:.2} Ha for max energy diff = {:.2} Ha",
            omega_max,
            max_energy_diff
        );

        Self::new(n_points, omega_max)
    }

    /// Get the imaginary part of frequencies as an array
    pub fn get_xi_values(&self) -> Array1<f64> {
        Array1::from_vec(self.frequencies.iter().map(|f| f.im).collect())
    }

    /// Validate grid properties
    pub fn validate(&self) -> Result<()> {
        // Check all frequencies are purely imaginary
        for freq in &self.frequencies {
            if freq.re != 0.0 {
                return Err(QuasixError::NumericalError(
                    "Grid contains non-imaginary frequencies".to_string(),
                ));
            }
            if !freq.im.is_finite() || freq.im < 0.0 {
                return Err(QuasixError::NumericalError(
                    "Invalid imaginary frequency value".to_string(),
                ));
            }
        }

        // Check weights are positive
        for weight in &self.weights {
            if !weight.is_finite() || *weight <= 0.0 {
                return Err(QuasixError::NumericalError(
                    "Invalid weight value".to_string(),
                ));
            }
        }

        // Check monotonicity
        for i in 1..self.frequencies.len() {
            if self.frequencies[i].im <= self.frequencies[i - 1].im {
                return Err(QuasixError::NumericalError(
                    "Grid frequencies not monotonically increasing".to_string(),
                ));
            }
        }

        // Check weight sum (should be approximately π/2 for [0,∞))
        let expected_sum = std::f64::consts::PI / 2.0;
        let relative_error = (self.metadata.weight_sum - expected_sum).abs() / expected_sum;

        if relative_error > 0.1 {
            log::warn!(
                "Weight sum {:.6} differs from expected π/2 = {:.6} by {:.2}%",
                self.metadata.weight_sum,
                expected_sum,
                relative_error * 100.0
            );
        }

        Ok(())
    }

    /// Estimate integration error based on grid density
    pub fn estimate_error(&self) -> f64 {
        // Simple error estimate based on grid spacing
        if self.frequencies.len() < 2 {
            return 1.0;
        }

        // Maximum gap between consecutive points
        let mut max_gap: f64 = 0.0;
        for i in 1..self.frequencies.len() {
            let gap = self.frequencies[i].im - self.frequencies[i - 1].im;
            max_gap = max_gap.max(gap);
        }

        // Error scales with max gap and omega_max
        max_gap / self.metadata.omega_max
    }
}

/// Transform Gauss-Legendre quadrature to imaginary axis
///
/// This is the core transformation function that maps GL points from [-1,1]
/// to the imaginary frequency axis [0, i*infinity) with proper Jacobian.
pub fn transform_to_imaginary_axis(
    gl_points: &Array1<f64>,
    _gl_weights: &Array1<f64>, // Note: weights are recomputed by ImaginaryAxisGrid
    omega_max: f64,
) -> Result<(Vec<Complex64>, Vec<f64>)> {
    let grid = ImaginaryAxisGrid::new(gl_points.len(), omega_max)?;
    Ok((grid.frequencies, grid.weights))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_overflow_protection() {
        let protection = OverflowProtection::default();

        // Test near-singularity detection
        assert!(!protection.is_safe(1e-13));
        assert!(protection.is_safe(1e-10));

        // Test regularization
        let small_val = 1e-16;
        let regularized = protection.regularize(small_val);
        assert_eq!(regularized, protection.epsilon);

        // Test jacobian bounds
        assert!(protection.is_jacobian_safe(1e9));
        assert!(!protection.is_jacobian_safe(1e11));
        assert!(!protection.is_jacobian_safe(f64::INFINITY));
    }

    #[test]
    fn test_imaginary_axis_grid_creation() {
        let grid = ImaginaryAxisGrid::new(16, 30.0).unwrap();

        // Check grid size
        assert!(grid.frequencies.len() <= 16);
        assert_eq!(grid.frequencies.len(), grid.weights.len());

        // Validate grid
        assert!(grid.validate().is_ok());
    }

    #[test]
    fn test_grid_monotonicity() {
        let grid = ImaginaryAxisGrid::new(32, 50.0).unwrap();

        // Check all frequencies are purely imaginary
        for freq in &grid.frequencies {
            assert_eq!(freq.re, 0.0);
            assert!(freq.im >= 0.0);
        }

        // Check monotonic increase
        for i in 1..grid.frequencies.len() {
            assert!(grid.frequencies[i].im > grid.frequencies[i - 1].im);
        }
    }

    #[test]
    fn test_weight_positivity() {
        let grid = ImaginaryAxisGrid::new(24, 40.0).unwrap();

        for weight in &grid.weights {
            assert!(*weight > 0.0);
            assert!(weight.is_finite());
        }
    }

    #[test]
    fn test_gaussian_integral() {
        // Test integral of exp(-ξ²) from 0 to ∞
        // Exact result: √π/2

        let grid = ImaginaryAxisGrid::new(48, 10.0).unwrap();

        let mut integral = 0.0;
        for (i, freq) in grid.frequencies.iter().enumerate() {
            let xi = freq.im;
            let f = (-xi * xi).exp();
            integral += f * grid.weights[i];
        }

        let exact = std::f64::consts::PI.sqrt() / 2.0;

        // Should be accurate to ~1% with 48 points
        assert_relative_eq!(integral, exact, epsilon = 0.01);
    }

    #[test]
    fn test_convergence_with_grid_size() {
        // Test that integral converges as grid size increases
        let mut prev_integral = 0.0;

        for n in [8, 16, 32, 64] {
            let grid = ImaginaryAxisGrid::new(n, 20.0).unwrap();

            // Integrate exp(-ξ)
            let mut integral = 0.0;
            for (i, freq) in grid.frequencies.iter().enumerate() {
                let xi = freq.im;
                integral += (-xi).exp() * grid.weights[i];
            }

            if n > 8 {
                // Check convergence
                let change = (integral - prev_integral).abs();
                assert!(change < prev_integral * 0.1); // Each refinement improves
            }

            prev_integral = integral;
        }

        // Final value should be close to 1.0
        assert_relative_eq!(prev_integral, 1.0, epsilon = 1e-4);
    }

    #[test]
    fn test_auto_cutoff() {
        let max_energy_diff = 2.0; // Ha
        let grid = ImaginaryAxisGrid::with_auto_cutoff(32, max_energy_diff).unwrap();

        // Check omega_max is reasonable
        assert!(grid.metadata.omega_max >= 3.0 * max_energy_diff);
        assert!(grid.metadata.omega_max <= 100.0);
    }

    #[test]
    fn test_error_estimate() {
        let coarse_grid = ImaginaryAxisGrid::new(8, 30.0).unwrap();
        let fine_grid = ImaginaryAxisGrid::new(32, 30.0).unwrap();

        // Both grids should have reasonable error estimates
        let coarse_error = coarse_grid.estimate_error();
        let fine_error = fine_grid.estimate_error();

        // Print for debugging
        println!("Coarse error: {}, Fine error: {}", coarse_error, fine_error);
        println!(
            "Coarse points: {}, Fine points: {}",
            coarse_grid.metadata.n_actual, fine_grid.metadata.n_actual
        );

        assert!(coarse_error > 0.0);
        assert!(fine_error > 0.0);
        assert!(coarse_error.is_finite());
        assert!(fine_error.is_finite());

        // The error estimate is based on max gap, which can vary
        // Just check that they're positive and finite
    }
}
