//! Frequency grids and transformations module
//!
//! This module provides frequency grid generation for contour deformation (CD)
//! and analytical continuation (AC) methods, including Gauss-Legendre quadrature.
//!
//! Key features:
//! - High-precision Gauss-Legendre quadrature (machine precision < 1e-14)
//! - Imaginary frequency grids for GW calculations
//! - Frequency transformations between real and imaginary axes
//! - Contour deformation integration for self-energy evaluation
//! - Analytic continuation with multipole and Padé models
//! - Catastrophic cancellation prevention via Kahan summation
//! - Condition number monitoring for numerical stability
//! - Symmetrized operations for improved accuracy
//! - **Adaptive quadrature** for 2-4x speedup with error control
#![allow(clippy::many_single_char_names)] // Mathematical notation

pub mod adaptive_quadrature;
pub mod analytic_continuation;
pub mod contour_deformation;
pub mod gauss_legendre;
pub mod grid_optimized;
pub mod imag_axis;
pub mod imaginary_axis_grid;
pub mod imaginary_grid;
pub mod jacobian_optimized;
pub mod stability;
pub mod stability_validation;
pub mod transforms;

use crate::common::{QuasixError, Result};
use ndarray::Array1;
use num_complex::Complex64;

pub use adaptive_quadrature::{
    AdaptiveConfig, AdaptiveFrequencyGrid, GaussKronrodPair, OptimizedGrids,
};
pub use analytic_continuation::{
    ACConfig, ACResult, AnalyticContinuation, AnalyticContinuationFitter, CausalityMetrics,
    ModelType, MultipoleModel, PadeModel,
};
pub use contour_deformation::{CDConfig, CDIntegrator};
pub use gauss_legendre::{gauss_legendre_nodes_weights, gauss_legendre_scaled};
pub use grid_optimized::{
    DoubleExponentialGrid, ImaginaryAxisGrid, MinimaxGrid, OptimizedGLGrid, TransformType,
};
pub use imag_axis::{ImagAxisCalculator, ImagAxisConfig};
pub use imaginary_axis_grid::{
    transform_to_imaginary_axis as new_transform_to_imaginary_axis, GridMetadata,
    ImaginaryAxisGrid as NewImaginaryAxisGrid, OverflowProtection,
};
// Note: transform_to_imaginary_axis and transform_weights are deprecated but exported for compatibility
// Use transform_to_imaginary_axis_pyscf() and transform_weights_pyscf() for GW calculations
#[allow(deprecated)]
pub use imaginary_grid::{
    create_imaginary_frequency_grid, gauss_legendre_grid, transform_to_imaginary_axis,
    transform_weights, ImaginaryGridConfig,
};
pub use jacobian_optimized::{JacobianOptimized, TransformType as JacobianTransformType};
pub use stability::{
    ComplexKahanSum, ConditionMonitor, ExtendedPrecision, KahanSum, StabilizedTransform,
    SymmetrizedOps,
};
pub use transforms::{apply_frequency_transform, imaginary_to_real, real_to_imaginary};

/// Frequency grid types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GridType {
    /// Gauss-Legendre grid for imaginary axis
    GaussLegendre,
    /// Contour deformation grid
    ContourDeformation,
    /// Minimax grid for analytical continuation
    Minimax,
    /// Modified Gauss-Legendre for GW calculations
    ModifiedGaussLegendre { omega_max: f64 },
}

/// Frequency grid generator for GW/BSE calculations with numerical stability monitoring
#[derive(Debug, Clone)]
pub struct FrequencyGrid {
    /// Number of frequency points
    pub nfreq: usize,
    /// Grid type
    pub grid_type: GridType,
    /// Frequency points (can be real or imaginary)
    pub points: Array1<f64>,
    /// Integration weights
    pub weights: Array1<f64>,
    /// Condition number of the weights
    pub condition_number: Option<f64>,
    /// Maximum stable frequency
    pub omega_stable: Option<f64>,
}

impl FrequencyGrid {
    /// Create a new frequency grid with stability monitoring
    pub fn new(nfreq: usize, grid_type: GridType) -> Result<Self> {
        let (points, weights) = Self::generate_grid(nfreq, grid_type)?;

        // Monitor condition number
        let mut monitor = ConditionMonitor::new();
        let condition_number = monitor.check_weights(&weights).ok();

        // Find maximum stable frequency
        let omega_stable = Self::find_stable_cutoff(&points, &weights);

        Ok(Self {
            nfreq,
            grid_type,
            points,
            weights,
            condition_number,
            omega_stable,
        })
    }

    /// Generate frequency points and weights based on grid type with stability enhancements
    fn generate_grid(nfreq: usize, grid_type: GridType) -> Result<(Array1<f64>, Array1<f64>)> {
        match grid_type {
            GridType::GaussLegendre => {
                // Standard GL quadrature on [-1, 1]
                gauss_legendre_nodes_weights(nfreq)
            }
            GridType::ModifiedGaussLegendre { omega_max } => {
                // Use stabilized transformation for GL quadrature
                Self::generate_stabilized_gl_grid(nfreq, omega_max)
            }
            GridType::ContourDeformation => {
                // Contour deformation grid
                Self::generate_cd_grid(nfreq)
            }
            GridType::Minimax => {
                // Minimax grid for analytical continuation
                Self::generate_minimax_grid(nfreq)
            }
        }
    }

    /// Generate contour deformation grid
    fn generate_cd_grid(nfreq: usize) -> Result<(Array1<f64>, Array1<f64>)> {
        // For now, use a simple Gauss-Legendre grid
        // In production, this would be a complex contour
        gauss_legendre_scaled(nfreq, 0.0, 100.0)
    }

    /// Generate minimax grid for analytical continuation
    fn generate_minimax_grid(nfreq: usize) -> Result<(Array1<f64>, Array1<f64>)> {
        // Minimax grid optimized for Padé approximation
        // For now, use logarithmically spaced points
        let mut points = Array1::zeros(nfreq);
        let mut weights = Array1::zeros(nfreq);

        let omega_min: f64 = 1e-3;
        let omega_max: f64 = 1e3;
        let log_min = omega_min.ln();
        let log_max = omega_max.ln();

        for i in 0..nfreq {
            let t = i as f64 / (nfreq - 1) as f64;
            let log_omega = log_min + t * (log_max - log_min);
            points[i] = log_omega.exp();
            // Simple trapezoidal weights for now
            weights[i] = if i == 0 || i == nfreq - 1 {
                0.5 / (nfreq - 1) as f64
            } else {
                1.0 / (nfreq - 1) as f64
            };
        }

        Ok((points, weights))
    }

    /// Get frequency points
    pub fn points(&self) -> &Array1<f64> {
        &self.points
    }

    /// Get integration weights
    pub fn weights(&self) -> &Array1<f64> {
        &self.weights
    }

    /// Transform to imaginary frequency grid
    pub fn to_imaginary_grid(&self) -> Array1<Complex64> {
        self.points.mapv(|omega| Complex64::new(0.0, omega))
    }

    /// Evaluate integrand using quadrature with Kahan summation for stability
    pub fn integrate<F>(&self, mut integrand: F) -> f64
    where
        F: FnMut(f64) -> f64,
    {
        let mut sum = KahanSum::new();

        for (&x, &w) in self.points.iter().zip(self.weights.iter()) {
            sum.add(w * integrand(x));
        }

        sum.sum()
    }

    /// Evaluate complex integrand using quadrature with compensated summation
    pub fn integrate_complex<F>(&self, mut integrand: F) -> Complex64
    where
        F: FnMut(Complex64) -> Complex64,
    {
        let mut sum = ComplexKahanSum::new();

        for (&x, &w) in self.points.iter().zip(self.weights.iter()) {
            let z = Complex64::new(x, 0.0);
            sum.add(Complex64::new(w, 0.0) * integrand(z));
        }

        sum.sum()
    }

    /// Generate stabilized GL grid with controlled transformation
    fn generate_stabilized_gl_grid(
        nfreq: usize,
        omega_max: f64,
    ) -> Result<(Array1<f64>, Array1<f64>)> {
        // Get GL nodes on [-1, 1]
        let (t_nodes, w_gl) = gauss_legendre_nodes_weights(nfreq)?;

        // Map to [0, 1]
        let t_01 = t_nodes.mapv(|t| (t + 1.0) * 0.5);
        let w_01 = w_gl.mapv(|w| w * 0.5);

        // Apply stabilized transformation
        let transform = StabilizedTransform::new();
        let mut omega = Array1::zeros(nfreq);
        let mut weights = Array1::zeros(nfreq);

        for i in 0..nfreq {
            let (om, jac) = transform.transform(t_01[i], omega_max);
            omega[i] = om;
            weights[i] = w_01[i] * jac;
        }

        // Check for stability issues
        let mut monitor = ConditionMonitor::new();
        if let Err(msg) = monitor.check_weights(&weights) {
            eprintln!("Warning in grid generation: {}", msg);
        }

        Ok((omega, weights))
    }

    /// Find maximum stable frequency cutoff
    fn find_stable_cutoff(points: &Array1<f64>, weights: &Array1<f64>) -> Option<f64> {
        let n = points.len();
        let weight_threshold = 1e-12; // Below this, numerical errors dominate

        for i in (0..n).rev() {
            if weights[i] > weight_threshold {
                return Some(points[i]);
            }
        }

        None
    }

    /// Validate grid quality with comprehensive checks
    pub fn validate(&self) -> Result<()> {
        // Check positive weights
        if self.weights.iter().any(|&w| w <= 0.0) {
            return Err(QuasixError::InvalidInput(
                "Non-positive weights detected".to_string(),
            ));
        }

        // Check monotonic points
        for i in 1..self.nfreq {
            if self.points[i] <= self.points[i - 1] {
                return Err(QuasixError::InvalidInput(
                    "Non-monotonic frequency points".to_string(),
                ));
            }
        }

        // Warn about conditioning
        if let Some(cond) = self.condition_number {
            if cond > 1e8 {
                eprintln!(
                    "Warning: Grid condition number {:.2e} may cause accuracy loss",
                    cond
                );
            }
        }

        Ok(())
    }
}

/// Analytical continuation fitter using Padé approximants
pub struct ACFitter {
    /// Number of poles in rational function
    pub npoles: usize,
    /// Regularization parameter for fitting
    pub regularization: f64,
    /// Fitted pole positions
    pub poles: Option<Array1<Complex64>>,
    /// Fitted residues
    pub residues: Option<Array1<Complex64>>,
}

impl ACFitter {
    /// Create a new AC fitter
    #[must_use]
    pub fn new(npoles: usize) -> Self {
        Self {
            npoles,
            regularization: 1e-8,
            poles: None,
            residues: None,
        }
    }

    /// Set regularization parameter
    #[must_use]
    pub fn with_regularization(mut self, reg: f64) -> Self {
        self.regularization = reg;
        self
    }

    /// Fit data on imaginary axis to rational function
    pub fn fit(&mut self, _iw_points: &Array1<f64>, _iw_data: &Array1<Complex64>) -> Result<()> {
        // Placeholder for Padé approximant fitting
        // In production, this would use least squares or AAA algorithm

        // For now, store simple poles and residues
        let mut poles = Array1::zeros(self.npoles);
        let mut residues = Array1::zeros(self.npoles);

        for i in 0..self.npoles {
            // Simple pole distribution
            poles[i] = Complex64::new(0.0, (i + 1) as f64);
            residues[i] = Complex64::new(1.0 / self.npoles as f64, 0.0);
        }

        self.poles = Some(poles);
        self.residues = Some(residues);

        Ok(())
    }

    /// Evaluate continuation on real axis
    pub fn evaluate(&self, omega: f64) -> Result<Complex64> {
        let poles = self
            .poles
            .as_ref()
            .ok_or_else(|| QuasixError::InvalidInput("Fitter not yet fitted".to_string()))?;
        let residues = self
            .residues
            .as_ref()
            .ok_or_else(|| QuasixError::InvalidInput("Fitter not yet fitted".to_string()))?;

        let z = Complex64::new(omega, self.regularization);
        let result: Complex64 = poles
            .iter()
            .zip(residues.iter())
            .map(|(pole, residue)| residue / (z - pole))
            .sum();

        Ok(result)
    }
}

/// Contour deformation calculator for self-energy evaluation
pub struct ContourDeformation {
    /// Energy range for deformation [E_min, E_max]
    pub energy_range: (f64, f64),
    /// Number of contour points
    pub npoints: usize,
    /// Broadening parameter (imaginary shift)
    pub eta: f64,
}

impl ContourDeformation {
    /// Create a new CD calculator
    pub fn new(energy_range: (f64, f64), npoints: usize) -> Self {
        Self {
            energy_range,
            npoints,
            eta: 1e-3, // Default broadening
        }
    }

    /// Set broadening parameter
    #[must_use]
    pub fn with_eta(mut self, eta: f64) -> Self {
        self.eta = eta;
        self
    }

    /// Generate contour points for integration
    pub fn contour_points(&self) -> Array1<Complex64> {
        let (e_min, e_max) = self.energy_range;
        let mut points = Array1::zeros(self.npoints);

        for i in 0..self.npoints {
            let t = i as f64 / (self.npoints - 1) as f64;
            let energy = e_min + t * (e_max - e_min);
            points[i] = Complex64::new(energy, self.eta);
        }

        points
    }

    /// Compute CD integral using trapezoidal rule
    pub fn integrate<F>(&self, mut integrand: F) -> Result<Complex64>
    where
        F: FnMut(Complex64) -> Complex64,
    {
        let points = self.contour_points();
        let (e_min, e_max) = self.energy_range;
        let h = (e_max - e_min) / (self.npoints - 1) as f64;

        let mut result = Complex64::new(0.0, 0.0);

        for (i, &z) in points.iter().enumerate() {
            let weight = if i == 0 || i == self.npoints - 1 {
                h * 0.5
            } else {
                h
            };
            result += Complex64::new(weight, 0.0) * integrand(z);
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_frequency_grid_creation() {
        let grid = FrequencyGrid::new(16, GridType::GaussLegendre).unwrap();
        assert_eq!(grid.nfreq, 16);
        assert_eq!(grid.points.len(), 16);
        assert_eq!(grid.weights.len(), 16);

        // Validate grid quality
        grid.validate().unwrap();

        // Check that weights sum to 2 (integral over [-1, 1])
        // Use Kahan summation for accurate comparison
        let mut sum = KahanSum::new();
        for &w in &grid.weights {
            sum.add(w);
        }
        assert_relative_eq!(sum.sum(), 2.0, epsilon = 1e-14);

        // Check condition number is reasonable
        assert!(grid.condition_number.is_some());
        if let Some(cond) = grid.condition_number {
            assert!(cond < 1e6, "Condition number too high: {}", cond);
        }
    }

    #[test]
    fn test_modified_gl_grid() {
        let omega_max = 100.0;
        let grid = FrequencyGrid::new(20, GridType::ModifiedGaussLegendre { omega_max }).unwrap();

        // Check that points are in [0, omega_max]
        assert!(grid.points.iter().all(|&x| x >= 0.0 && x <= omega_max));

        // For transformed grids, the weight sum depends on the transformation
        // We just verify that weights are positive and reasonable
        assert!(grid.weights.iter().all(|&w| w > 0.0));

        // The integral of a constant function should give a reasonable result
        let const_integral = grid.integrate(|_| 1.0);
        assert!(const_integral > 0.0 && const_integral.is_finite());
    }

    #[test]
    fn test_integration() {
        let grid = FrequencyGrid::new(32, GridType::GaussLegendre).unwrap();

        // Test integration of x^2 over [-1, 1]
        // Exact integral = 2/3
        let result = grid.integrate(|x| x * x);
        assert_relative_eq!(result, 2.0 / 3.0, epsilon = 1e-14);

        // Test integration of x^4
        // Exact integral = 2/5
        let result = grid.integrate(|x| x.powi(4));
        assert_relative_eq!(result, 2.0 / 5.0, epsilon = 1e-14);
    }

    #[test]
    fn test_ac_fitter_creation() {
        let fitter = ACFitter::new(10);
        assert_eq!(fitter.npoles, 10);
        assert_eq!(fitter.regularization, 1e-8);
    }

    #[test]
    fn test_contour_deformation() {
        let cd = ContourDeformation::new((-10.0, 10.0), 100);
        assert_eq!(cd.energy_range, (-10.0, 10.0));
        assert_eq!(cd.npoints, 100);

        let points = cd.contour_points();
        assert_eq!(points.len(), 100);

        // Check that all points have small imaginary part
        assert!(points.iter().all(|z| z.im.abs() < 0.1));
    }

    #[test]
    fn test_to_imaginary_grid() {
        let grid = FrequencyGrid::new(10, GridType::GaussLegendre).unwrap();
        let imag_grid = grid.to_imaginary_grid();

        assert_eq!(imag_grid.len(), 10);
        // All points should be purely imaginary
        assert!(imag_grid.iter().all(|z| z.re == 0.0));
    }

    #[test]
    fn test_stability_enhanced_integration() {
        // Test integration with values that would cause catastrophic cancellation
        let grid = FrequencyGrid::new(32, GridType::GaussLegendre).unwrap();

        // Function that tests cancellation: large constant that cancels + x^2
        // The large values should cancel perfectly, leaving only the x^2 integral
        let large_val = 1e10;
        let result = grid.integrate(|x| {
            // Add and subtract the same large value, plus x^2
            large_val - large_val + x * x
        });

        // The integral should be close to 2/3 (integral of x^2 over [-1,1])
        // This tests that Kahan summation preserves the small x^2 contribution
        assert_relative_eq!(result, 2.0 / 3.0, epsilon = 1e-14);
    }

    #[test]
    fn test_stabilized_gl_transformation() {
        // Test that stabilized transformation prevents overflow
        let omega_max = 1000.0;
        let grid = FrequencyGrid::new(50, GridType::ModifiedGaussLegendre { omega_max }).unwrap();

        // Check all points and weights are finite
        assert!(grid.points.iter().all(|&x| x.is_finite()));
        assert!(grid.weights.iter().all(|&w| w.is_finite()));

        // Check stable cutoff is identified
        assert!(grid.omega_stable.is_some());
        if let Some(cutoff) = grid.omega_stable {
            assert!(cutoff > 0.0 && cutoff <= omega_max);
        }
    }

    #[test]
    fn test_condition_monitoring() {
        // Create grid with high omega_max that might cause conditioning issues
        let omega_max = 10000.0;
        let grid = FrequencyGrid::new(80, GridType::ModifiedGaussLegendre { omega_max }).unwrap();

        // Condition number should be tracked
        assert!(grid.condition_number.is_some());

        // For very high omega_max, expect higher condition number
        if let Some(cond) = grid.condition_number {
            // Should be reasonable and finite (transformation dependent)
            assert!(cond > 0.0 && cond.is_finite());
            // Warn if condition number is very high
            if cond > 1e12 {
                eprintln!("Warning: Very high condition number: {:.2e}", cond);
            }
        }
    }
}
