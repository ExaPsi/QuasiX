//! Gauss-Legendre quadrature implementation with numerical stability
//!
//! This module provides high-precision Gauss-Legendre quadrature nodes and weights
//! computation for frequency integration in GW calculations.
#![allow(clippy::many_single_char_names)] // Mathematical notation

use ndarray::Array1;
use std::f64::consts::PI;

/// Tolerance for Newton-Raphson convergence
const NEWTON_TOL: f64 = 1e-15;

/// Maximum Newton iterations
const MAX_NEWTON_ITER: usize = 100;

/// Machine epsilon safety factor
const EPS_SAFETY: f64 = 10.0;

/// Gauss-Legendre quadrature rule
#[derive(Debug, Clone)]
pub struct GaussLegendre {
    /// Quadrature nodes in [-1, 1]
    pub nodes: Array1<f64>,
    /// Quadrature weights
    pub weights: Array1<f64>,
    /// Number of quadrature points
    pub n: usize,
}

impl GaussLegendre {
    /// Create a new Gauss-Legendre quadrature rule with n points
    ///
    /// # Arguments
    /// * `n` - Number of quadrature points
    ///
    /// # Returns
    /// * Gauss-Legendre quadrature rule with nodes and weights
    #[must_use]
    pub fn new(n: usize) -> Self {
        assert!(n > 0, "Number of quadrature points must be positive");

        let mut nodes = Array1::zeros(n);
        let mut weights = Array1::zeros(n);

        // The roots of P_n are symmetric about 0
        // We compute roots from largest to smallest (right to left)
        for i in 0..n {
            // Get initial guess for (i+1)-th root counting from the right
            let x0 = Self::initial_guess(n, i + 1);

            // Find root using Newton-Raphson with Halley correction
            let (xi, dpn) = Self::find_root_newton_halley(n, x0);

            // Compute weight using stable formula
            let wi = Self::compute_weight(n, xi, dpn);

            // Store root and weight
            nodes[i] = xi;
            weights[i] = wi;
        }

        // Sort nodes in ascending order
        let mut indices: Vec<usize> = (0..n).collect();
        indices.sort_by(|&i, &j| nodes[i].partial_cmp(&nodes[j]).unwrap());

        let sorted_nodes = indices.iter().map(|&i| nodes[i]).collect();
        let sorted_weights = indices.iter().map(|&i| weights[i]).collect();

        Self {
            nodes: sorted_nodes,
            weights: sorted_weights,
            n,
        }
    }

    /// Compute initial guess for k-th root using Tricomi approximation
    /// k ranges from 1 to n, where k=1 is the largest (rightmost) root
    fn initial_guess(n: usize, k: usize) -> f64 {
        let n_f64 = n as f64;
        let k_f64 = k as f64;

        // Use standard Chebyshev approximation
        // The k-th root (counting from right) is approximately:
        let theta = PI * (k_f64 - 0.25) / (n_f64 + 0.5);
        let x0 = theta.cos();

        // Apply Tricomi correction for better accuracy
        if k == 1 {
            // Special case for root near 1
            let h = PI / (2.0 * n_f64 + 1.0);
            1.0 - h * h / 3.0 * (1.0 + h * h / 15.0)
        } else if k == n {
            // Special case for root near -1
            let h = PI / (2.0 * n_f64 + 1.0);
            -1.0 + h * h / 3.0 * (1.0 + h * h / 15.0)
        } else {
            // Interior roots with first-order correction
            x0
        }
    }

    /// Find root using Newton-Raphson with Halley correction
    /// Returns (root, derivative at root)
    fn find_root_newton_halley(n: usize, x0: f64) -> (f64, f64) {
        let mut x = x0;
        let n_f64 = n as f64;

        for _ in 0..MAX_NEWTON_ITER {
            // Compute P_n(x) and P'_n(x) using stable recurrence
            let (pn, dpn) = Self::legendre_pair(n, x);

            // Check for convergence
            if pn.abs() < NEWTON_TOL {
                return (x, dpn);
            }

            // Newton correction
            let dx = pn / dpn;

            // Halley correction for cubic convergence
            // P''_n(x) = (2x P'_n(x) - n(n+1) P_n(x)) / (1 - x²)
            let x2 = x * x;
            if (1.0 - x2).abs() > EPS_SAFETY * f64::EPSILON {
                let d2pn = (2.0 * x * dpn - n_f64 * (n_f64 + 1.0) * pn) / (1.0 - x2);
                let halley_factor = 1.0 / (1.0 - 0.5 * dx * d2pn / dpn);
                x -= dx * halley_factor;
            } else {
                // Near boundaries, use simple Newton
                x -= dx;
            }

            // Safeguard: keep x in (-1, 1)
            x = x
                .max(-1.0 + EPS_SAFETY * f64::EPSILON)
                .min(1.0 - EPS_SAFETY * f64::EPSILON);
        }

        // If not converged, use bisection as fallback
        let (root, dpn) = Self::bisection_fallback(n, x0);
        (root, dpn)
    }

    /// Compute P_n(x) and P'_n(x) simultaneously using stable recurrence
    fn legendre_pair(n: usize, x: f64) -> (f64, f64) {
        if n == 0 {
            return (1.0, 0.0);
        }
        if n == 1 {
            return (x, 1.0);
        }

        // Three-term recurrence for P_n
        let mut p0 = 1.0;
        let mut p1 = x;

        // Recurrence for P'_n
        let mut dp0 = 0.0;
        let mut dp1 = 1.0;

        for k in 2..=n {
            let k_f64 = k as f64;

            // P_k = ((2k-1) x P_{k-1} - (k-1) P_{k-2}) / k
            let pk = ((2.0 * k_f64 - 1.0) * x * p1 - (k_f64 - 1.0) * p0) / k_f64;

            // P'_k = ((2k-1) (P_{k-1} + x P'_{k-1}) - (k-1) P'_{k-2}) / k
            let dpk = ((2.0 * k_f64 - 1.0) * (p1 + x * dp1) - (k_f64 - 1.0) * dp0) / k_f64;

            p0 = p1;
            p1 = pk;
            dp0 = dp1;
            dp1 = dpk;
        }

        (p1, dp1)
    }

    /// Compute weight using numerically stable formula
    fn compute_weight(n: usize, x: f64, dpn: f64) -> f64 {
        let n_f64 = n as f64;
        let x2 = x * x;

        // Check if near endpoints
        if (x2 - 1.0).abs() < 100.0 * f64::EPSILON {
            // Special formula for endpoints
            2.0 / (n_f64 * (n_f64 + 1.0))
        } else {
            // Standard formula: w = 2 / ((1 - x²) [P'_n(x)]²)
            2.0 / ((1.0 - x2) * dpn * dpn)
        }
    }

    /// Bisection fallback for difficult roots
    fn bisection_fallback(n: usize, x0: f64) -> (f64, f64) {
        // Define search interval around initial guess
        let mut a = (x0 - 0.1).max(-1.0 + f64::EPSILON);
        let mut b = (x0 + 0.1).min(1.0 - f64::EPSILON);

        let (mut pa, _) = Self::legendre_pair(n, a);
        let (mut pb, _) = Self::legendre_pair(n, b);

        // Ensure opposite signs
        if pa * pb > 0.0 {
            // Expand search interval
            a = -1.0 + f64::EPSILON;
            b = 1.0 - f64::EPSILON;
            let (pa_new, _) = Self::legendre_pair(n, a);
            let (pb_new, _) = Self::legendre_pair(n, b);
            pa = pa_new;
            pb = pb_new;
        }

        // Bisection iteration
        for _ in 0..100 {
            let c = 0.5 * (a + b);
            let (pc, dpc) = Self::legendre_pair(n, c);

            if pc.abs() < NEWTON_TOL || (b - a) < NEWTON_TOL {
                return (c, dpc);
            }

            if pa * pc < 0.0 {
                b = c;
                pb = pc;
            } else {
                a = c;
                pa = pc;
            }
        }

        let c = 0.5 * (a + b);
        let (_, dpc) = Self::legendre_pair(n, c);
        (c, dpc)
    }

    /// Transform nodes and weights to [0, ω_max] for frequency integration
    pub fn transform_to_frequency(&self, omega_max: f64) -> (Array1<f64>, Array1<f64>) {
        let mut freq_nodes = Array1::zeros(self.n);
        let mut freq_weights = Array1::zeros(self.n);

        for i in 0..self.n {
            // Linear transformation: [-1, 1] → [0, ω_max]
            freq_nodes[i] = omega_max * (self.nodes[i] + 1.0) / 2.0;
            freq_weights[i] = self.weights[i] * omega_max / 2.0;
        }

        (freq_nodes, freq_weights)
    }

    /// Transform nodes and weights to [0, ∞) using t/(1-t) mapping
    pub fn transform_to_semi_infinite(&self, scale: f64) -> (Array1<f64>, Array1<f64>) {
        let mut freq_nodes = Array1::zeros(self.n);
        let mut freq_weights = Array1::zeros(self.n);

        for i in 0..self.n {
            let t = self.nodes[i];

            // Transformation: ω = scale * (1+t)/(1-t)
            // This maps [-1, 1] → [0, ∞)
            if (1.0 - t).abs() > f64::EPSILON {
                freq_nodes[i] = scale * (1.0 + t) / (1.0 - t);
                freq_weights[i] = self.weights[i] * 2.0 * scale / ((1.0 - t) * (1.0 - t));
            } else {
                // Handle near t = 1 carefully
                freq_nodes[i] = 1e10 * scale; // Large but finite
                freq_weights[i] = 0.0; // Negligible weight
            }
        }

        (freq_nodes, freq_weights)
    }

    /// Compute integral using Kahan summation for reduced round-off error
    pub fn integrate<F>(&self, f: F) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let mut sum = 0.0;
        let mut c = 0.0; // Compensation term

        for i in 0..self.n {
            let y = self.weights[i] * f(self.nodes[i]) - c;
            let t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }

        sum
    }

    /// Validate quadrature rule by checking weight sum and polynomial exactness
    pub fn validate(&self) -> Result<(), String> {
        // Check weight sum (should be 2 for [-1, 1])
        let weight_sum: f64 = self.weights.iter().sum();
        if (weight_sum - 2.0).abs() > 1e-14 {
            return Err(format!("Weight sum = {}, expected 2.0", weight_sum));
        }

        // Check nodes are in [-1, 1] and sorted
        for i in 0..self.n {
            if self.nodes[i] < -1.0 || self.nodes[i] > 1.0 {
                return Err(format!("Node {} = {} is outside [-1, 1]", i, self.nodes[i]));
            }
            if i > 0 && self.nodes[i] <= self.nodes[i - 1] {
                return Err(format!("Nodes not sorted at index {i}"));
            }
        }

        // Check exactness for x^{2n-1}
        let degree = 2 * self.n - 1;
        let exact = if degree % 2 == 0 {
            2.0 / (degree as f64 + 1.0)
        } else {
            0.0 // Odd powers integrate to 0
        };

        let computed = self.integrate(|x| x.powi(degree as i32));
        if (computed - exact).abs() > 1e-13 {
            return Err(format!(
                "Failed exactness test for x^{}: computed = {}, exact = {}",
                degree, computed, exact
            ));
        }

        Ok(())
    }
}

/// Adaptive Gauss-Legendre quadrature for automatic convergence
pub struct AdaptiveGaussLegendre {
    min_points: usize,
    max_points: usize,
    tolerance: f64,
}

impl AdaptiveGaussLegendre {
    /// Create adaptive quadrature with specified parameters
    #[must_use]
    pub fn new(min_points: usize, max_points: usize, tolerance: f64) -> Self {
        assert!(min_points > 0 && min_points <= max_points);
        assert!(tolerance > 0.0);

        Self {
            min_points,
            max_points,
            tolerance,
        }
    }

    /// Perform adaptive integration
    pub fn integrate<F>(&self, f: F) -> Result<f64, String>
    where
        F: Fn(f64) -> f64,
    {
        let mut n = self.min_points;
        let mut prev_result = 0.0;
        let mut converged = false;
        let mut result = 0.0;

        while n <= self.max_points {
            let gl = GaussLegendre::new(n);
            result = gl.integrate(&f);

            if n > self.min_points {
                let error = (result - prev_result).abs();
                if error < self.tolerance {
                    converged = true;
                    break;
                }
            }

            prev_result = result;

            // Increase points (golden ratio for smooth increase)
            n = ((n as f64 * 1.618).ceil() as usize).min(self.max_points);
            if n == self.max_points && !converged {
                break;
            }
        }

        if !converged {
            Err(format!(
                "Failed to converge with {} points, error = {}",
                self.max_points,
                (result - prev_result).abs()
            ))
        } else {
            Ok(result)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gauss_legendre_creation() {
        for n in [1, 2, 5, 10, 20, 40] {
            let gl = GaussLegendre::new(n);
            assert_eq!(gl.n, n);
            assert_eq!(gl.nodes.len(), n);
            assert_eq!(gl.weights.len(), n);

            // Validate the rule
            gl.validate()
                .expect(&format!("Validation failed for n={n}"));
        }
    }

    #[test]
    fn test_weight_sum() {
        for n in [5, 10, 20, 30, 50] {
            let gl = GaussLegendre::new(n);
            let sum: f64 = gl.weights.iter().sum();
            assert_relative_eq!(sum, 2.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn test_polynomial_exactness() {
        let n = 10;
        let gl = GaussLegendre::new(n);

        // Test exactness for polynomials up to degree 2n-1
        for degree in 0..2 * n {
            let exact = if degree % 2 == 0 {
                2.0 / (degree as f64 + 1.0)
            } else {
                0.0
            };

            let computed = gl.integrate(|x| x.powi(degree as i32));
            assert_relative_eq!(computed, exact, epsilon = 1e-13);
        }
    }

    #[test]
    fn test_gaussian_integration() {
        let gl = GaussLegendre::new(30);

        // Integrate exp(-x²) over [-1, 1]
        let computed = gl.integrate(|x| (-x * x).exp());
        let exact = 2.0 * std::f64::consts::PI.sqrt() * erf::erf(1.0); // Using error function

        assert_relative_eq!(computed, exact, epsilon = 1e-12);
    }

    #[test]
    fn test_frequency_transform() {
        let gl = GaussLegendre::new(10);
        let omega_max = 100.0;

        let (freq_nodes, freq_weights) = gl.transform_to_frequency(omega_max);

        // Check range
        assert!(freq_nodes.iter().all(|&x| x >= 0.0 && x <= omega_max));

        // Check weight sum
        let sum: f64 = freq_weights.iter().sum();
        assert_relative_eq!(sum, omega_max, epsilon = 1e-14);
    }

    #[test]
    fn test_adaptive_integration() {
        let adaptive = AdaptiveGaussLegendre::new(10, 100, 1e-10);

        // Integrate sin(x) over [-1, 1]
        let result = adaptive.integrate(|x| x.sin()).unwrap();
        let exact = 0.0; // sin is odd, integral over symmetric interval is 0

        assert_relative_eq!(result, exact, epsilon = 1e-10);
    }

    #[test]
    fn test_symmetry() {
        let gl = GaussLegendre::new(20);

        // Check symmetry of nodes and weights
        let n = gl.n;
        for i in 0..n / 2 {
            assert_relative_eq!(gl.nodes[i], -gl.nodes[n - 1 - i], epsilon = 1e-15);
            assert_relative_eq!(gl.weights[i], gl.weights[n - 1 - i], epsilon = 1e-15);
        }
    }
}

// Helper module for error function (if not available in dependencies)
#[cfg(test)]
mod erf {
    pub fn erf(x: f64) -> f64 {
        // Simple approximation for testing
        let a1 = 0.254829592;
        let a2 = -0.284496736;
        let a3 = 1.421413741;
        let a4 = -1.453152027;
        let a5 = 1.061405429;
        let p = 0.3275911;

        let sign = if x < 0.0 { -1.0 } else { 1.0 };
        let x = x.abs();

        let t = 1.0 / (1.0 + p * x);
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;
        let t5 = t4 * t;

        let y = 1.0 - (((((a5 * t5 + a4 * t4) + a3 * t3) + a2 * t2) + a1 * t) * (-x * x).exp());

        sign * y
    }
}
