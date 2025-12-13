// Richardson extrapolation for accurate numerical derivatives
// Specifically designed for Z-factor computation in GW calculations

use ndarray::{Array1, Array2};

/// Richardson extrapolation for computing derivatives
/// Achieves O(h^(2n)) accuracy where n is the number of extrapolation levels
pub struct RichardsonDerivative {
    /// Initial step size (typically 1e-3 to 1e-4 for self-energy)
    initial_delta: f64,
    /// Refinement factor (typically 2.0 for optimal extrapolation)
    refinement_factor: f64,
    /// Number of extrapolation levels (3-4 optimal for GW)
    n_levels: usize,
    /// Tolerance for convergence check
    tolerance: f64,
}

impl RichardsonDerivative {
    /// Create new Richardson extrapolator with optimal GW parameters
    pub fn new_for_gw() -> Self {
        Self {
            initial_delta: 1e-3,      // Optimal for eV-scale energies
            refinement_factor: 2.0,    // Power-of-2 for exact extrapolation
            n_levels: 3,               // Balance accuracy vs cost
            tolerance: 1e-10,          // High precision for Z-factors
        }
    }

    /// Create custom Richardson extrapolator
    pub fn new(initial_delta: f64, refinement_factor: f64, n_levels: usize) -> Self {
        Self {
            initial_delta,
            refinement_factor,
            n_levels,
            tolerance: 1e-10,
        }
    }

    /// Compute derivative using Richardson extrapolation
    /// Returns (derivative, error_estimate)
    pub fn compute_derivative<F>(&self, f: F, x: f64) -> (f64, f64)
    where
        F: Fn(f64) -> f64,
    {
        // Build Richardson tableau
        let mut tableau = Array2::<f64>::zeros((self.n_levels, self.n_levels));
        let mut h = self.initial_delta;

        // First column: finite differences at different step sizes
        for i in 0..self.n_levels {
            tableau[[i, 0]] = self.central_difference(&f, x, h);
            h /= self.refinement_factor;
        }

        // Richardson extrapolation: build remaining columns
        let mut power = 1.0;
        for j in 1..self.n_levels {
            power *= self.refinement_factor.powi(2); // For central differences: O(h²)
            let factor = power / (power - 1.0);

            for i in j..self.n_levels {
                tableau[[i, j]] = factor * tableau[[i, j-1]]
                                 - tableau[[i-1, j-1]] / (power - 1.0);
            }
        }

        // Best estimate is bottom-right corner
        let derivative = tableau[[self.n_levels-1, self.n_levels-1]];

        // Error estimate from convergence pattern
        let error = if self.n_levels >= 2 {
            (tableau[[self.n_levels-1, self.n_levels-1]]
             - tableau[[self.n_levels-1, self.n_levels-2]]).abs()
        } else {
            self.initial_delta.powi(2) // Fallback error estimate
        };

        (derivative, error)
    }

    /// Adaptive Richardson with automatic level selection
    pub fn compute_derivative_adaptive<F>(&self, f: F, x: f64, max_levels: usize) -> (f64, f64)
    where
        F: Fn(f64) -> f64,
    {
        let mut tableau = Vec::with_capacity(max_levels);
        let mut h = self.initial_delta;
        let mut converged = false;
        let mut best_derivative = 0.0;
        let mut best_error = f64::INFINITY;

        for level in 0..max_levels {
            // Add new row to tableau
            let mut row = Vec::with_capacity(level + 1);

            // First entry: finite difference
            row.push(self.central_difference(&f, x, h));

            // Richardson extrapolation
            let mut power = 1.0;
            for j in 1..=level {
                power *= self.refinement_factor.powi(2);
                let factor = power / (power - 1.0);

                let extrap = factor * row[j-1] - tableau[level-1][j-1] / (power - 1.0);
                row.push(extrap);
            }

            // Check convergence
            if level > 0 {
                let current = row[level];
                let previous = tableau[level-1][level-1];
                let error = (current - previous).abs();

                if error < best_error {
                    best_derivative = current;
                    best_error = error;
                }

                // Stop if converged or error is increasing (noise dominating)
                if error < self.tolerance || (level > 2 && error > 2.0 * best_error) {
                    converged = true;
                    break;
                }
            } else {
                best_derivative = row[0];
            }

            tableau.push(row);
            h /= self.refinement_factor;
        }

        (best_derivative, best_error)
    }

    /// Central difference formula
    #[inline]
    fn central_difference<F>(&self, f: &F, x: f64, h: f64) -> f64
    where
        F: Fn(f64) -> f64,
    {
        (f(x + h) - f(x - h)) / (2.0 * h)
    }

    /// Higher-order central difference (5-point stencil)
    #[inline]
    fn central_difference_5point<F>(&self, f: &F, x: f64, h: f64) -> f64
    where
        F: Fn(f64) -> f64,
    {
        (-f(x + 2.0*h) + 8.0*f(x + h) - 8.0*f(x - h) + f(x - 2.0*h)) / (12.0 * h)
    }
}

/// Specialized implementation for Z-factor computation
pub struct ZFactorDerivative {
    richardson: RichardsonDerivative,
}

impl ZFactorDerivative {
    /// Create with optimal parameters for self-energy derivatives
    pub fn new() -> Self {
        Self {
            richardson: RichardsonDerivative {
                initial_delta: 5e-4,     // Smaller for steep self-energy
                refinement_factor: 2.0,   // Standard refinement
                n_levels: 3,              // 3 levels optimal for GW
                tolerance: 1e-12,         // High precision needed
            },
        }
    }

    /// Compute Z-factor from self-energy function
    /// Z = 1 / (1 - ∂Σᶜ/∂ω)
    pub fn compute_z_factor<F>(&self, sigma_c: F, omega: f64) -> Result<f64, String>
    where
        F: Fn(f64) -> f64,
    {
        let (derivative, error) = self.richardson.compute_derivative_adaptive(sigma_c, omega, 4);

        // Check for reasonable derivative (causality: ∂Σ/∂ω < 0 typically)
        if derivative.abs() > 10.0 {
            return Err(format!("Derivative too large: {:.6}, possible numerical instability", derivative));
        }

        let z = 1.0 / (1.0 - derivative);

        // Physical bounds check
        if z < 0.0 || z > 1.0 {
            return Err(format!("Unphysical Z-factor: {:.6} (derivative: {:.6})", z, derivative));
        }

        // Warn if close to bounds
        if z < 0.01 {
            eprintln!("Warning: Z-factor very small: {:.6} (strong quasiparticle decay)", z);
        } else if z > 0.99 {
            eprintln!("Warning: Z-factor very large: {:.6} (weak correlation)", z);
        }

        Ok(z)
    }

    /// Compute Z-factor with fallback strategies
    pub fn compute_z_factor_robust<F>(&self, sigma_c: F, omega: f64) -> f64
    where
        F: Fn(f64) -> f64,
    {
        // Try Richardson first
        if let Ok(z) = self.compute_z_factor(&sigma_c, omega) {
            return z;
        }

        // Fallback 1: Try larger step size
        let richardson_large = RichardsonDerivative::new(1e-2, 2.0, 3);
        let (deriv_large, _) = richardson_large.compute_derivative(&sigma_c, omega);
        let z_large = 1.0 / (1.0 - deriv_large);

        if z_large > 0.0 && z_large < 1.0 {
            eprintln!("Warning: Using larger step size for Z-factor: {:.6}", z_large);
            return z_large;
        }

        // Fallback 2: Simple finite difference
        let h = 1e-3;
        let deriv_simple = (sigma_c(omega + h) - sigma_c(omega - h)) / (2.0 * h);
        let z_simple = 1.0 / (1.0 - deriv_simple);

        // Clamp to physical bounds
        let z_clamped = z_simple.max(0.01).min(0.99);
        eprintln!("Warning: Clamping Z-factor to physical range: {:.6} -> {:.6}", z_simple, z_clamped);

        z_clamped
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_richardson_polynomial() {
        // Test on polynomial where derivative is exact
        let f = |x: f64| x.powi(3) - 2.0 * x;
        let df_exact = |x: f64| 3.0 * x.powi(2) - 2.0;

        let richardson = RichardsonDerivative::new_for_gw();

        for x in &[0.0, 1.0, -1.0, 2.5] {
            let (deriv, error) = richardson.compute_derivative(f, *x);
            let exact = df_exact(*x);

            assert!((deriv - exact).abs() < 1e-10,
                    "x={}: got {}, expected {}", x, deriv, exact);
            assert!(error < 1e-10, "Error estimate too large: {}", error);
        }
    }

    #[test]
    fn test_richardson_exponential() {
        // Test on exponential (sensitive to step size)
        let f = |x: f64| (-x.powi(2)).exp();
        let df_exact = |x: f64| -2.0 * x * (-x.powi(2)).exp();

        let richardson = RichardsonDerivative::new_for_gw();

        for x in &[0.0, 0.5, 1.0, 1.5] {
            let (deriv, _) = richardson.compute_derivative(f, *x);
            let exact = df_exact(*x);
            let rel_error = ((deriv - exact) / exact.max(1e-10)).abs();

            assert!(rel_error < 1e-8,
                    "x={}: relative error {} too large", x, rel_error);
        }
    }

    #[test]
    fn test_richardson_rational() {
        // Test on rational function (like self-energy)
        let f = |x: f64| 1.0 / (1.0 + x.powi(2));
        let df_exact = |x: f64| -2.0 * x / (1.0 + x.powi(2)).powi(2);

        let richardson = RichardsonDerivative::new_for_gw();

        for x in &[0.0, 0.5, 1.0, 2.0] {
            let (deriv, _) = richardson.compute_derivative_adaptive(f, *x, 4);
            let exact = df_exact(*x);
            let rel_error = ((deriv - exact) / exact.max(1e-10)).abs();

            assert!(rel_error < 1e-9,
                    "x={}: relative error {} too large", x, rel_error);
        }
    }

    #[test]
    fn test_z_factor_physical() {
        // Mock self-energy with known derivative
        let sigma_c = |omega: f64| -0.3 * omega - 0.1 * omega.powi(2);
        // ∂Σ/∂ω = -0.3 - 0.2*ω

        let z_calc = ZFactorDerivative::new();

        // At ω = 0: ∂Σ/∂ω = -0.3, Z = 1/(1-(-0.3)) = 1/1.3 ≈ 0.769
        let z = z_calc.compute_z_factor(sigma_c, 0.0).unwrap();
        assert!((z - 0.769).abs() < 0.001, "Z-factor at ω=0: {}", z);

        // At ω = 1: ∂Σ/∂ω = -0.5, Z = 1/(1-(-0.5)) = 1/1.5 ≈ 0.667
        let z = z_calc.compute_z_factor(sigma_c, 1.0).unwrap();
        assert!((z - 0.667).abs() < 0.001, "Z-factor at ω=1: {}", z);
    }

    #[test]
    fn test_richardson_vs_simple() {
        // Compare accuracy: Richardson vs simple finite difference
        let f = |x: f64| (2.0 * x).sin() * (-x).exp();
        let df_exact = |x: f64| {
            2.0 * (2.0 * x).cos() * (-x).exp() - (2.0 * x).sin() * (-x).exp()
        };

        let x = 1.0;
        let exact = df_exact(x);

        // Simple finite difference
        let h = 1e-3;
        let simple = (f(x + h) - f(x - h)) / (2.0 * h);
        let error_simple = (simple - exact).abs();

        // Richardson extrapolation
        let richardson = RichardsonDerivative::new_for_gw();
        let (deriv_rich, _) = richardson.compute_derivative(f, x);
        let error_rich = (deriv_rich - exact).abs();

        // Richardson should be >100x more accurate
        assert!(error_rich < error_simple / 100.0,
                "Richardson not accurate enough: {} vs {}", error_rich, error_simple);

        println!("Accuracy improvement: {:.2}x", error_simple / error_rich);
    }
}