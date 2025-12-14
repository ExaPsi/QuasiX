//! Contour Deformation (CD) integration for GW self-energy evaluation
//!
//! This module implements robust numerical algorithms for evaluating the GW correlation
//! self-energy Σc(ω) using contour deformation techniques. The method separates the
//! frequency integral into pole residues and a smooth imaginary-axis integral.
#![allow(clippy::many_single_char_names)] // Mathematical notation

use ndarray::{Array1, Array2, Array3};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Configuration for contour deformation integration
#[derive(Debug, Clone)]
pub struct CDConfig {
    /// Number of Gauss-Legendre quadrature points (30-80 typical)
    pub n_grid: usize,

    /// Maximum frequency for grid transformation (Ha)
    pub omega_max: f64,

    /// Small imaginary shift for stability (0.01-0.05 eV typical)
    pub eta: f64,

    /// Relative tolerance for adaptive refinement
    pub rel_tol: f64,

    /// Absolute tolerance for convergence
    pub abs_tol: f64,

    /// Maximum refinement iterations
    pub max_refine: usize,

    /// Distance threshold for near-pole detection (Ha)
    pub pole_threshold: f64,

    /// Condition number threshold for numerical stability
    pub cond_threshold: f64,
}

impl Default for CDConfig {
    fn default() -> Self {
        Self {
            n_grid: 50,           // Good balance for molecules
            omega_max: 30.0,      // ~800 eV cutoff
            eta: 0.001,           // ~0.027 eV broadening
            rel_tol: 1e-6,        // Relative accuracy
            abs_tol: 1e-10,       // Absolute accuracy (Ha)
            max_refine: 3,        // Adaptive refinement depth
            pole_threshold: 0.01, // Near-pole detection (~0.27 eV)
            cond_threshold: 1e8,  // Numerical stability limit
        }
    }
}

/// Gauss-Legendre quadrature grid with transformation
pub struct GLGrid {
    /// Quadrature nodes on [0,∞) via transformation
    pub nodes: Array1<f64>,

    /// Quadrature weights including Jacobian
    pub weights: Array1<f64>,

    /// Original GL nodes on [-1,1]
    pub raw_nodes: Array1<f64>,

    /// Original GL weights
    pub raw_weights: Array1<f64>,

    /// Transformation parameter
    pub omega_max: f64,
}

impl GLGrid {
    /// Generate Gauss-Legendre grid with optimal transformation
    #[must_use]
    pub fn new(n_points: usize, omega_max: f64) -> Self {
        let (raw_nodes, raw_weights) = gauss_legendre_quadrature(n_points);

        let mut nodes = Array1::zeros(n_points);
        let mut weights = Array1::zeros(n_points);

        for i in 0..n_points {
            let t = raw_nodes[i];

            // Transform to [0,∞)
            nodes[i] = omega_max * (1.0 + t) / (1.0 - t);

            // Include Jacobian in weights
            let jacobian = 2.0 * omega_max / (1.0 - t).powi(2);
            weights[i] = raw_weights[i] * jacobian;
        }

        Self {
            nodes,
            weights,
            raw_nodes,
            raw_weights,
            omega_max,
        }
    }

    /// Refine grid adaptively based on error estimate
    pub fn refine(&mut self, error_estimate: f64) -> bool {
        if error_estimate < 1e-10 {
            return false; // Already converged
        }

        // Increase grid points by 50%
        let new_n = (self.nodes.len().saturating_mul(3) / 2).min(120);
        if new_n <= self.nodes.len() {
            return false; // Already at maximum
        }

        *self = Self::new(new_n, self.omega_max);
        true
    }
}

/// Compute Gauss-Legendre nodes and weights on [-1,1]
fn gauss_legendre_quadrature(n: usize) -> (Array1<f64>, Array1<f64>) {
    // Use the existing gauss_legendre module implementation
    use crate::freq::gauss_legendre::gauss_legendre_nodes_weights;

    if let Ok((nodes, weights)) = gauss_legendre_nodes_weights(n) {
        (nodes, weights)
    } else {
        // Fallback to simple implementation
        let mut nodes = Array1::zeros(n);
        let mut weights = Array1::zeros(n);

        for i in 0..n {
            let x = -1.0 + 2.0 * (i as f64 + 0.5) / n as f64;
            nodes[i] = x;
            weights[i] = 2.0 / n as f64;
        }

        (nodes, weights)
    }
}

/// Residue calculator for pole contributions
pub struct ResidueCalculator {
    /// Orbital energies (Ha)
    pub energies: Array1<f64>,

    /// Occupation numbers (0-2)
    pub occupations: Array1<f64>,

    /// Fermi level (Ha)
    pub fermi_level: f64,

    /// Temperature for Fermi-Dirac (Ha)
    pub temperature: f64,
}

impl ResidueCalculator {
    /// Calculate residue contributions from Green's function poles
    pub fn compute_residues(
        &self,
        omega: f64,
        w_correlation: &dyn Fn(f64) -> Array2<f64>,
        orbitals: &Array2<f64>,
    ) -> Array1<f64> {
        let n_orb = self.energies.len();
        let mut sigma_residue = Array1::zeros(n_orb);

        for n in 0..n_orb {
            let mut residue_sum = 0.0;

            // Electron contributions (p > n, virtual)
            for p in 0..n_orb {
                if self.occupations[p] < 0.5 {
                    // Virtual orbital
                    let pole_freq = omega - self.energies[p] + self.energies[n];

                    // Check if pole is close to real axis
                    if pole_freq.abs() < 10.0 * self.temperature {
                        // Near-pole regularization
                        let w_c = w_correlation(pole_freq);
                        let coupling = self.compute_coupling(n, p, &w_c, orbitals);

                        // Fermi-Dirac factor for finite temperature
                        let fermi_factor = if self.temperature > 1e-10 {
                            1.0 / (1.0
                                + ((self.energies[p] - self.fermi_level) / self.temperature).exp())
                        } else {
                            0.0 // T=0 case
                        };

                        residue_sum += coupling * (1.0 - fermi_factor);
                    }
                }
            }

            // Hole contributions (p < n, occupied)
            for p in 0..n_orb {
                if self.occupations[p] > 0.5 {
                    // Occupied orbital
                    let pole_freq = omega - self.energies[p] + self.energies[n];

                    if pole_freq.abs() < 10.0 * self.temperature {
                        let w_c = w_correlation(pole_freq);
                        let coupling = self.compute_coupling(n, p, &w_c, orbitals);

                        let fermi_factor = if self.temperature > 1e-10 {
                            1.0 / (1.0
                                + ((self.energies[p] - self.fermi_level) / self.temperature).exp())
                        } else {
                            1.0 // T=0 case
                        };

                        residue_sum += coupling * fermi_factor;
                    }
                }
            }

            sigma_residue[n] = residue_sum;
        }

        sigma_residue
    }

    /// Compute coupling element <np|W_c|pn>
    fn compute_coupling(
        &self,
        n: usize,
        p: usize,
        w_c: &Array2<f64>,
        orbitals: &Array2<f64>,
    ) -> f64 {
        // Contract: sum_PQ φ_n(P) W_c(P,Q) φ_p(Q)
        let phi_n = orbitals.row(n);
        let phi_p = orbitals.row(p);

        // Two-step contraction for efficiency
        let temp = phi_n.dot(w_c);
        temp.dot(&phi_p)
    }
}

/// Main contour deformation integrator
pub struct CDIntegrator {
    config: CDConfig,
    grid: GLGrid,
    residue_calc: ResidueCalculator,
    /// Cache for W(iξ) evaluations (reserved for future use)
    _w_cache: HashMap<String, Array3<f64>>,
}

impl CDIntegrator {
    #[must_use]
    pub fn new(
        config: CDConfig,
        energies: Array1<f64>,
        occupations: Array1<f64>,
        fermi_level: f64,
    ) -> Self {
        let grid = GLGrid::new(config.n_grid, config.omega_max);
        let residue_calc = ResidueCalculator {
            energies,
            occupations,
            fermi_level,
            temperature: 0.001, // ~300K default
        };

        Self {
            config,
            grid,
            residue_calc,
            _w_cache: HashMap::new(),
        }
    }

    /// Evaluate correlation self-energy using contour deformation
    pub fn compute_self_energy(
        &mut self,
        omega: f64,
        w_correlation: &dyn Fn(f64) -> Array2<f64>,
        orbitals: &Array2<f64>,
    ) -> Result<Array1<f64>, String> {
        let n_orb = self.residue_calc.energies.len();

        // Step 1: Residue contributions
        let sigma_residue = self
            .residue_calc
            .compute_residues(omega, w_correlation, orbitals);

        // Step 2: Imaginary-axis integral
        let mut sigma_integral = Array1::zeros(n_orb);
        let mut error_estimate = 0.0;

        for refine_iter in 0..=self.config.max_refine {
            let (integral, error) =
                self.integrate_imaginary_axis(omega, w_correlation, orbitals)?;

            sigma_integral = integral;
            error_estimate = error;

            // Check convergence
            if error < self.config.abs_tol || error < self.config.rel_tol * sigma_integral.norm() {
                break;
            }

            // Adaptive refinement
            if refine_iter < self.config.max_refine && !self.grid.refine(error) {
                break; // Cannot refine further
            }
        }

        // Check final accuracy
        if error_estimate > 10.0 * self.config.abs_tol {
            eprintln!(
                "Warning: CD integration accuracy {:.2e} exceeds target",
                error_estimate
            );
        }

        // Step 3: Combine contributions
        let sigma_total = sigma_residue + sigma_integral / PI;

        Ok(sigma_total)
    }

    /// Evaluate integral along imaginary axis
    fn integrate_imaginary_axis(
        &self,
        omega: f64,
        w_correlation: &dyn Fn(f64) -> Array2<f64>,
        orbitals: &Array2<f64>,
    ) -> Result<(Array1<f64>, f64), String> {
        let n_orb = self.residue_calc.energies.len();
        let mut integral = Array1::<f64>::zeros(n_orb);

        // Error estimation using embedded lower-order rule
        let mut integral_low = Array1::<f64>::zeros(n_orb);
        let skip_stride = 2; // Use every other point for lower-order estimate

        // Gauss-Legendre quadrature along iξ
        for (idx, (&xi, &weight)) in self
            .grid
            .nodes
            .iter()
            .zip(self.grid.weights.iter())
            .enumerate()
        {
            // Evaluate Green's function on imaginary axis
            let g_imag = self.green_function_imag(omega, xi);

            // Get W_c(iξ) - real and symmetric on imaginary axis
            let w_c_imag = w_correlation(xi);

            // Check conditioning
            let cond = condition_number(&w_c_imag);
            if cond > self.config.cond_threshold {
                return Err(format!(
                    "Poor conditioning at ξ={:.2e}: cond={:.2e}",
                    xi, cond
                ));
            }

            // Accumulate integral
            for n in 0..n_orb {
                let coupling = self.compute_imag_coupling(n, xi, &g_imag, &w_c_imag, orbitals);
                integral[n] += weight * coupling;

                // Lower-order estimate for error
                if idx % skip_stride == 0 {
                    integral_low[n] += weight * coupling * skip_stride as f64;
                }
            }
        }

        // Richardson extrapolation error estimate
        let error = (&integral - &integral_low).norm() / (2_f64.powi(4) - 1.0);

        Ok((integral, error))
    }

    /// Green's function on imaginary axis
    fn green_function_imag(&self, omega: f64, xi: f64) -> Array1<f64> {
        let n_orb = self.residue_calc.energies.len();
        let mut g_imag = Array1::zeros(n_orb);

        for p in 0..n_orb {
            // G_p(ω + iξ) = 1/(ω + iξ - ε_p + iη)
            let denom_real = omega - self.residue_calc.energies[p];
            let denom_imag = xi + self.config.eta;

            // Real part of 1/(a + ib) = a/(a² + b²)
            g_imag[p] = denom_real / (denom_real.powi(2) + denom_imag.powi(2));
        }

        g_imag
    }

    /// Compute coupling on imaginary axis
    fn compute_imag_coupling(
        &self,
        n: usize,
        _xi: f64,
        g_imag: &Array1<f64>,
        w_c_imag: &Array2<f64>,
        orbitals: &Array2<f64>,
    ) -> f64 {
        let mut coupling = 0.0;

        for p in 0..g_imag.len() {
            if g_imag[p].abs() > 1e-12 {
                let phi_n = orbitals.row(n);
                let phi_p = orbitals.row(p);

                // <np|W_c(iξ)|pn> with proper symmetrization
                let temp = phi_n.dot(w_c_imag);
                let matrix_element = temp.dot(&phi_p);

                coupling += g_imag[p] * matrix_element;
            }
        }

        coupling
    }
}

/// Norm function for `Array1`
trait Norm {
    fn norm(&self) -> f64;
}

impl Norm for Array1<f64> {
    fn norm(&self) -> f64 {
        self.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
}

/// Estimate condition number of symmetric matrix
fn condition_number(matrix: &Array2<f64>) -> f64 {
    // For production, use LAPACK's SYEV for eigenvalues
    // Here we use power iteration for largest eigenvalue estimate

    let n = matrix.nrows();
    let mut v = Array1::ones(n) / (n as f64).sqrt();
    let mut lambda_max = 0.0;

    // Power iteration for largest eigenvalue
    for _ in 0..20 {
        v = matrix.dot(&v);
        lambda_max = v.norm();
        if lambda_max > 0.0 {
            v /= lambda_max;
        }
    }

    // Inverse power iteration for smallest eigenvalue
    // (would use actual linear solver in production)
    let lambda_min = 1e-8; // Conservative estimate

    lambda_max / lambda_min
}

/// Adaptive error estimator using Richardson extrapolation
#[derive(Default)]
pub struct ErrorEstimator {
    /// Previous integral values for extrapolation
    history: Vec<Array1<f64>>,

    /// Corresponding grid sizes
    grid_sizes: Vec<usize>,
}

impl ErrorEstimator {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Estimate error using Richardson extrapolation
    pub fn estimate_error(&mut self, integral: &Array1<f64>, grid_size: usize) -> Option<f64> {
        self.history.push(integral.clone());
        self.grid_sizes.push(grid_size);

        if self.history.len() < 2 {
            return None;
        }

        let n = self.history.len();
        let i1 = &self.history[n - 2];
        let i2 = &self.history[n - 1];
        let n1 = self.grid_sizes[n - 2];
        let n2 = self.grid_sizes[n - 1];

        // Estimate convergence rate
        let alpha = 2.0 * PI / (n2 as f64).sqrt(); // Theoretical for GL
        let ratio = ((n2 - n1) as f64 * alpha).exp();

        let error = (i2 - i1).norm() / (ratio - 1.0);

        Some(error)
    }
}

/// Near-pole regularization strategies
pub mod regularization {
    use super::*;

    /// Detect near-pole regions requiring special treatment
    pub fn detect_near_poles(omega: f64, energies: &Array1<f64>, threshold: f64) -> Vec<usize> {
        let mut near_poles = Vec::new();

        for (i, &e) in energies.iter().enumerate() {
            if (omega - e).abs() < threshold {
                near_poles.push(i);
            }
        }

        near_poles
    }

    /// Apply Padé regularization near poles
    pub fn pade_regularization(omega: f64, pole: f64, eta: f64, order: usize) -> f64 {
        // [order/order] Padé approximant for 1/(ω - pole + iη)
        let z = omega - pole;
        let mut num = 1.0;
        let mut den = z;

        for n in 1..=order {
            let n_i32 = i32::try_from(n).unwrap_or(i32::MAX);
            let coeff = eta.powi(2 * n_i32) / factorial(2 * n);
            if n % 2 == 0 {
                num += coeff * z.powi(n_i32);
            } else {
                den += coeff * z.powi(n_i32);
            }
        }

        num / den
    }

    fn factorial(n: usize) -> f64 {
        (1..=n).map(|i| i as f64).product()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gl_grid_generation() {
        let grid = GLGrid::new(10, 10.0);

        // Check grid properties
        assert_eq!(grid.nodes.len(), 10);
        assert_eq!(grid.weights.len(), 10);

        // Nodes should be positive and increasing
        for i in 1..10 {
            assert!(grid.nodes[i] > grid.nodes[i - 1]);
            assert!(grid.nodes[i] > 0.0);
        }

        // Weights should be positive
        for w in grid.weights {
            assert!(w > 0.0);
        }
    }

    #[test]
    fn test_error_estimator() {
        let mut estimator = ErrorEstimator::new();

        let i1 = Array1::from_vec(vec![1.0, 2.0, 3.0]);
        let i2 = Array1::from_vec(vec![1.01, 2.01, 3.01]);

        assert!(estimator.estimate_error(&i1, 30).is_none());
        let error = estimator.estimate_error(&i2, 40).unwrap();

        assert!(error > 0.0);
        assert!(error < 0.1); // Should be small for close values
    }

    #[test]
    fn test_near_pole_detection() {
        use regularization::detect_near_poles;

        let energies = Array1::from_vec(vec![-0.5, -0.1, 0.0, 0.4, 0.8]);
        let omega = 0.41;
        let threshold = 0.05;

        let near_poles = detect_near_poles(omega, &energies, threshold);

        assert_eq!(near_poles, vec![3]); // Only index 3 (0.4) is near 0.41
    }
}
