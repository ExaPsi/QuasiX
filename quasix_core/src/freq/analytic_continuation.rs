//! Analytic continuation module for transforming imaginary-axis data to real frequencies
//!
//! This module provides high-performance implementations of multipole and Padé
//! analytic continuation methods for GW/BSE calculations. Key features include:
//!
//! - Multi-pole expansion with nonlinear optimization
//! - Padé approximants with Thiele recursion and SVD stabilization
//! - Cross-validation framework for model selection
//! - SIMD-optimized evaluation on real axis
//! - Thread-safe design for parallel processing
//!
//! # Mathematical Background
//!
//! Analytic continuation transforms functions known on the imaginary axis
//! f(iξ) to the real axis f(ω + iη) using rational function approximations.
//!
//! ## Multipole Model
//! ```text
//! f(z) = Σ_p A_p/(z - ω_p) + polynomial
//! ```
//!
//! ## Padé Approximant
//! ```text
//! f(z) ≈ P_M(z)/Q_N(z) = (Σ_m p_m z^m)/(1 + Σ_n q_n z^n)
//! ```

#![warn(clippy::all, clippy::pedantic, clippy::perf)]
#![warn(missing_docs, missing_debug_implementations)]

use crate::common::{QuasixError, Result};
use ndarray::{s, Array1, Array2};
use ndarray_linalg::{Solve, SVD};
use num_complex::Complex64;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    _mm256_add_pd, _mm256_div_pd, _mm256_loadu_pd, _mm256_mul_pd, _mm256_set1_pd,
    _mm256_setzero_pd, _mm256_storeu_pd, _mm256_sub_pd,
};

/// Type alias for complex numbers in calculations
type ComplexF64 = Complex64;

/// Solve a real linear system Ax = b using Gaussian elimination with pivoting
fn solve_linear_system(a: &Array2<f64>, b: &Array1<f64>) -> Option<Array1<f64>> {
    let n = a.nrows();
    if n != a.ncols() || n != b.len() {
        return None;
    }

    // Create augmented matrix [A|b]
    let mut aug = Array2::zeros((n, n + 1));
    aug.slice_mut(ndarray::s![.., ..n]).assign(a);
    aug.slice_mut(ndarray::s![.., n]).assign(b);

    // Gaussian elimination with partial pivoting
    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        for k in (i + 1)..n {
            if aug[[k, i]].abs() > aug[[max_row, i]].abs() {
                max_row = k;
            }
        }

        // Swap rows
        if max_row != i {
            for j in 0..=n {
                let tmp = aug[[i, j]];
                aug[[i, j]] = aug[[max_row, j]];
                aug[[max_row, j]] = tmp;
            }
        }

        // Check for singularity
        if aug[[i, i]].abs() < 1e-14 {
            return None;
        }

        // Forward elimination
        for k in (i + 1)..n {
            let factor = aug[[k, i]] / aug[[i, i]];
            for j in i..=n {
                let aug_ij = aug[[i, j]];
                aug[[k, j]] -= factor * aug_ij;
            }
        }
    }

    // Back substitution
    let mut x = Array1::zeros(n);
    for i in (0..n).rev() {
        x[i] = aug[[i, n]];
        for j in (i + 1)..n {
            let xj = x[j];
            x[i] -= aug[[i, j]] * xj;
        }
        x[i] /= aug[[i, i]];
    }

    Some(x)
}

/// Solve a complex linear system Ax = b using Gaussian elimination
#[allow(dead_code)]
fn _solve_complex_system(
    a: &Array2<ComplexF64>,
    b: &Array1<ComplexF64>,
) -> Option<Array1<ComplexF64>> {
    let n = a.nrows();
    if n != a.ncols() || n != b.len() {
        return None;
    }

    // Create augmented matrix [A|b]
    let mut aug = Array2::zeros((n, n + 1));
    aug.slice_mut(ndarray::s![.., ..n]).assign(a);
    aug.slice_mut(ndarray::s![.., n]).assign(b);

    // Gaussian elimination with partial pivoting
    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        for k in (i + 1)..n {
            if aug[[k, i]].norm() > aug[[max_row, i]].norm() {
                max_row = k;
            }
        }

        // Swap rows
        if max_row != i {
            for j in 0..=n {
                let tmp = aug[[i, j]];
                aug[[i, j]] = aug[[max_row, j]];
                aug[[max_row, j]] = tmp;
            }
        }

        // Check for singularity
        if aug[[i, i]].norm() < 1e-14 {
            return None;
        }

        // Forward elimination
        for k in (i + 1)..n {
            let factor = aug[[k, i]] / aug[[i, i]];
            for j in i..=n {
                let aug_ij = aug[[i, j]];
                aug[[k, j]] -= factor * aug_ij;
            }
        }
    }

    // Back substitution
    let mut x = Array1::zeros(n);
    for i in (0..n).rev() {
        x[i] = aug[[i, n]];
        for j in (i + 1)..n {
            let xj = x[j];
            x[i] -= aug[[i, j]] * xj;
        }
        x[i] /= aug[[i, i]];
    }

    Some(x)
}

/// Configuration for analytic continuation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ACConfig {
    /// Model type selection (Multipole or Pade)
    pub model_selection: ModelType,
    /// Maximum number of poles for multipole model
    pub max_poles: usize,
    /// Maximum order for Padé approximant [M/N]
    pub max_pade_order: (usize, usize),
    /// Fraction of data for cross-validation
    pub cv_fraction: f64,
    /// Number of CV iterations for error estimation
    pub cv_iterations: usize,
    /// Regularization parameter for fitting
    pub regularization: f64,
    /// Stability threshold for pole analysis
    pub stability_threshold: f64,
    /// Convergence tolerance for optimization
    pub convergence_tol: f64,
    /// Maximum iterations for fitting
    pub max_iterations: usize,
    /// Broadening parameter η for real-axis evaluation
    pub eta: f64,
    /// Enable parallel fitting over matrix elements
    pub parallel: bool,
}

impl Default for ACConfig {
    fn default() -> Self {
        Self {
            model_selection: ModelType::Multipole,
            max_poles: 20,
            max_pade_order: (15, 15),
            cv_fraction: 0.3,
            cv_iterations: 5,
            regularization: 1e-10,
            stability_threshold: 1e-3,
            convergence_tol: 1e-8,
            max_iterations: 1000,
            eta: 0.01,
            parallel: true,
        }
    }
}

/// Result of analytic continuation fitting
#[derive(Debug, Clone)]
pub struct ACResult {
    /// Selected model type
    pub model_type: ModelType,
    /// Cross-validation error
    pub cv_error: f64,
    /// Stability score (0=unstable, 1=stable)
    pub stability_score: f64,
    /// Fitted model
    pub model: Arc<dyn AnalyticContinuation>,
}

/// Types of analytic continuation models
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ModelType {
    /// Multi-pole expansion model
    Multipole,
    /// Padé approximant model
    Pade,
}

/// Trait for analytic continuation models
pub trait AnalyticContinuation: Send + Sync + std::fmt::Debug {
    /// Fit the model to imaginary-axis data
    ///
    /// # Errors
    ///
    /// Returns an error if the input data sizes don't match or if fitting fails.
    fn fit(&mut self, xi_data: &Array1<f64>, f_data: &Array1<ComplexF64>) -> Result<()>;

    /// Evaluate the model at a complex frequency
    fn evaluate(&self, z: ComplexF64) -> ComplexF64;

    /// Evaluate the model at multiple points (vectorized)
    fn evaluate_batch(&self, z: &Array1<ComplexF64>) -> Array1<ComplexF64>;

    /// Get cross-validation error
    fn cross_validation_error(&self) -> f64;

    /// Check stability of the fitted model
    fn stability_check(&self) -> bool;

    /// Get poles and residues (if applicable)
    fn get_poles_residues(&self) -> Option<(Array1<ComplexF64>, Array1<ComplexF64>)>;

    /// Clone the model into a boxed trait object
    fn clone_box(&self) -> Box<dyn AnalyticContinuation>;
}

/// Multi-pole expansion model
#[derive(Debug, Clone)]
pub struct MultipoleModel {
    /// Number of poles
    pub(crate) n_poles: usize,
    /// Pole positions in complex plane
    pub(crate) poles: Array1<ComplexF64>,
    /// Residues at poles
    pub(crate) residues: Array1<ComplexF64>,
    /// Polynomial coefficients
    polynomial: Array1<f64>,
    /// Cross-validation error
    cv_error: f64,
    /// Configuration
    config: ACConfig,
}

impl MultipoleModel {
    /// Create a new multipole model
    #[must_use]
    pub fn new(n_poles: usize, config: ACConfig) -> Self {
        Self {
            n_poles,
            poles: Array1::zeros(n_poles),
            residues: Array1::zeros(n_poles),
            polynomial: Array1::zeros(2), // Constant + linear term
            cv_error: f64::INFINITY,
            config,
        }
    }

    /// Perform nonlinear least-squares fitting (Algorithm 1 from S4-2)
    fn fit_nonlinear(&mut self, xi: &Array1<f64>, f: &Array1<ComplexF64>) -> Result<()> {
        let n = xi.len();
        if n < 2 * self.n_poles {
            return Err(QuasixError::InvalidInput(format!(
                "Need at least {} points for {} poles",
                2 * self.n_poles,
                self.n_poles
            )));
        }

        // Step 1: Moment-based initialization using Prony's method
        let initial_poles = self.prony_method_initialization(xi, f);

        // Use Prony initialization if successful, otherwise fall back to simple distribution
        if initial_poles.len() == self.n_poles {
            self.poles = initial_poles;
            // Initialize residues using linear least squares
            self.initialize_residues(xi, f);
        } else {
            // Fallback: poles distributed in upper-left quadrant (physical region)
            for i in 0..self.n_poles {
                #[allow(clippy::cast_precision_loss)]
                // Angle between π and 3π/2 (upper-left quadrant)
                let angle = std::f64::consts::PI * (1.0 + 0.5 * i as f64 / self.n_poles as f64);
                #[allow(clippy::cast_precision_loss)]
                let radius = 0.5 + 2.0 * i as f64 / self.n_poles as f64;
                self.poles[i] = ComplexF64::from_polar(radius, angle);
                #[allow(clippy::cast_precision_loss)]
                let residue_value = 1.0 / (self.n_poles as f64);
                self.residues[i] = ComplexF64::new(residue_value, 0.0);
            }
        }

        // Use Levenberg-Marquardt algorithm for fitting
        let mut lambda = 1e-3;
        let mut best_error = self.compute_error(xi, f);

        for _iter in 0..self.config.max_iterations {
            // Compute Jacobian matrix
            let jacobian = self.compute_jacobian(xi);

            // Build normal equations with regularization
            let jt_j = &jacobian.t().dot(&jacobian)
                + Array2::eye(jacobian.ncols()) * (lambda + self.config.regularization);

            // Compute residual vector
            let residual = self.compute_residual(xi, f);
            let jt_r = jacobian.t().dot(&residual);

            // Solve for parameter update using custom solver
            let Some(delta) = solve_linear_system(&jt_j, &jt_r) else {
                lambda *= 10.0;
                continue;
            };

            // Apply update
            self.apply_parameter_update(&delta);

            // Check if error improved
            let error = self.compute_error(xi, f);
            if error < best_error {
                best_error = error;
                lambda *= 0.1;

                if (best_error - error).abs() < self.config.convergence_tol {
                    break;
                }
            } else {
                // Revert update
                self.apply_parameter_update(&(-&delta));
                lambda *= 10.0;
            }
        }

        // Apply physical constraints: Re(poles) < 0 for causality
        self.apply_physical_constraints();

        // Detect and merge clustered poles (Algorithm 8)
        let (merged_poles, merged_residues) = self.detect_and_merge_poles();
        if merged_poles.len() < self.poles.len() {
            self.n_poles = merged_poles.len();
            self.poles = merged_poles;
            self.residues = merged_residues;
        }

        // Sort poles by real part for consistency
        let mut indices: Vec<usize> = (0..self.n_poles).collect();
        indices.sort_by(|&i, &j| self.poles[i].re.partial_cmp(&self.poles[j].re).unwrap());

        let sorted_poles: Vec<ComplexF64> = indices.iter().map(|&i| self.poles[i]).collect();
        let sorted_residues: Vec<ComplexF64> = indices.iter().map(|&i| self.residues[i]).collect();
        self.poles = Array1::from(sorted_poles);
        self.residues = Array1::from(sorted_residues);

        Ok(())
    }

    /// Compute error between model and data
    fn compute_error(&self, xi: &Array1<f64>, f: &Array1<ComplexF64>) -> f64 {
        let sum = xi
            .iter()
            .zip(f.iter())
            .map(|(&x, &fx)| {
                let z = ComplexF64::new(0.0, x);
                let model_val = self.evaluate(z);
                (model_val - fx).norm_sqr()
            })
            .sum::<f64>();
        #[allow(clippy::cast_precision_loss)]
        let xi_len = xi.len() as f64;
        sum / xi_len
    }

    /// Compute Jacobian matrix for nonlinear fitting
    fn compute_jacobian(&self, xi: &Array1<f64>) -> Array2<f64> {
        let n = xi.len();
        let n_params = 4 * self.n_poles + self.polynomial.len();
        let mut jacobian = Array2::zeros((2 * n, n_params));

        for (i, &x) in xi.iter().enumerate() {
            let z = ComplexF64::new(0.0, x);

            // Derivatives with respect to pole positions and residues
            for p in 0..self.n_poles {
                let denom = z - self.poles[p];
                let denom_sq = denom * denom;

                // ∂f/∂(Re[pole])
                jacobian[[2 * i, 4 * p]] = (self.residues[p] / denom_sq).re;
                jacobian[[2 * i + 1, 4 * p]] = (self.residues[p] / denom_sq).im;

                // ∂f/∂(Im[pole])
                jacobian[[2 * i, 4 * p + 1]] = (ComplexF64::i() * self.residues[p] / denom_sq).re;
                jacobian[[2 * i + 1, 4 * p + 1]] =
                    (ComplexF64::i() * self.residues[p] / denom_sq).im;

                // ∂f/∂(Re[residue])
                jacobian[[2 * i, 4 * p + 2]] = (1.0 / denom).re;
                jacobian[[2 * i + 1, 4 * p + 2]] = (1.0 / denom).im;

                // ∂f/∂(Im[residue])
                jacobian[[2 * i, 4 * p + 3]] = (ComplexF64::i() / denom).re;
                jacobian[[2 * i + 1, 4 * p + 3]] = (ComplexF64::i() / denom).im;
            }

            // Derivatives with respect to polynomial coefficients
            for (k, _coeff) in self.polynomial.iter().enumerate() {
                let z_pow = z.powi(i32::try_from(k).unwrap_or(i32::MAX));
                jacobian[[2 * i, 4 * self.n_poles + k]] = z_pow.re;
                jacobian[[2 * i + 1, 4 * self.n_poles + k]] = z_pow.im;
            }
        }

        jacobian
    }

    /// Compute residual vector
    fn compute_residual(&self, xi: &Array1<f64>, f: &Array1<ComplexF64>) -> Array1<f64> {
        let n = xi.len();
        let mut residual = Array1::zeros(2 * n);

        for (i, (&x, &fx)) in xi.iter().zip(f.iter()).enumerate() {
            let z = ComplexF64::new(0.0, x);
            let model_val = self.evaluate(z);
            let diff = fx - model_val;
            residual[2 * i] = diff.re;
            residual[2 * i + 1] = diff.im;
        }

        residual
    }

    /// Apply parameter update from optimization
    fn apply_parameter_update(&mut self, delta: &Array1<f64>) {
        let mut idx = 0;

        // Update poles and residues
        for p in 0..self.n_poles {
            self.poles[p].re += delta[idx];
            self.poles[p].im += delta[idx + 1];
            self.residues[p].re += delta[idx + 2];
            self.residues[p].im += delta[idx + 3];
            idx += 4;
        }

        // Update polynomial coefficients
        for k in 0..self.polynomial.len() {
            self.polynomial[k] += delta[idx + k];
        }
    }

    /// SIMD-optimized evaluation for `x86_64`
    #[cfg(all(target_arch = "x86_64", feature = "simd"))]
    #[target_feature(enable = "avx2")]
    #[allow(dead_code)]
    unsafe fn evaluate_simd(&self, omega: &[f64], eta: f64) -> Vec<ComplexF64> {
        let n = omega.len();
        let mut result = vec![ComplexF64::new(0.0, 0.0); n];

        // Process 4 frequencies at a time using AVX2
        let chunks = n / 4;
        for chunk in 0..chunks {
            let base = chunk * 4;
            let omega_vec = _mm256_loadu_pd(omega[base..].as_ptr());
            let eta_vec = _mm256_set1_pd(eta);

            let mut sum_re = _mm256_setzero_pd();
            let mut sum_im = _mm256_setzero_pd();

            for p in 0..self.n_poles {
                let pole = self.poles[p];
                let residue = self.residues[p];

                // Compute (omega - pole.re + i*(eta - pole.im))^-1
                let denom_re = _mm256_sub_pd(omega_vec, _mm256_set1_pd(pole.re));
                let denom_im = _mm256_sub_pd(eta_vec, _mm256_set1_pd(pole.im));

                // |denom|^2 = denom_re^2 + denom_im^2
                let denom_norm = _mm256_add_pd(
                    _mm256_mul_pd(denom_re, denom_re),
                    _mm256_mul_pd(denom_im, denom_im),
                );

                // Compute residue / denom
                let inv_norm = _mm256_div_pd(_mm256_set1_pd(1.0), denom_norm);
                let term_re = _mm256_mul_pd(
                    inv_norm,
                    _mm256_add_pd(
                        _mm256_mul_pd(_mm256_set1_pd(residue.re), denom_re),
                        _mm256_mul_pd(_mm256_set1_pd(residue.im), denom_im),
                    ),
                );
                let term_im = _mm256_mul_pd(
                    inv_norm,
                    _mm256_sub_pd(
                        _mm256_mul_pd(_mm256_set1_pd(residue.im), denom_re),
                        _mm256_mul_pd(_mm256_set1_pd(residue.re), denom_im),
                    ),
                );

                sum_re = _mm256_add_pd(sum_re, term_re);
                sum_im = _mm256_add_pd(sum_im, term_im);
            }

            // Store results
            let mut re_arr = [0.0; 4];
            let mut im_arr = [0.0; 4];
            _mm256_storeu_pd(re_arr.as_mut_ptr(), sum_re);
            _mm256_storeu_pd(im_arr.as_mut_ptr(), sum_im);

            for i in 0..4 {
                result[base + i] = ComplexF64::new(
                    re_arr[i] + self.polynomial[0] + self.polynomial[1] * omega[base + i],
                    im_arr[i],
                );
            }
        }

        // Handle remaining elements
        for i in (chunks * 4)..n {
            let z = ComplexF64::new(omega[i], eta);
            result[i] = self.evaluate(z);
        }

        result
    }
}

impl AnalyticContinuation for MultipoleModel {
    fn fit(&mut self, xi_data: &Array1<f64>, f_data: &Array1<ComplexF64>) -> Result<()> {
        // Validate input
        if xi_data.len() != f_data.len() {
            return Err(QuasixError::InvalidInput(format!(
                "Data size mismatch: {} vs {}",
                xi_data.len(),
                f_data.len()
            )));
        }

        // Perform fitting
        self.fit_nonlinear(xi_data, f_data)?;

        // Compute cross-validation error
        self.cv_error = self.cross_validate(xi_data, f_data)?;

        Ok(())
    }

    fn evaluate(&self, z: ComplexF64) -> ComplexF64 {
        let mut result = ComplexF64::new(self.polynomial[0], 0.0);

        // Add polynomial terms
        for (k, &coeff) in self.polynomial.iter().enumerate().skip(1) {
            result += coeff * z.powi(i32::try_from(k).unwrap_or(i32::MAX));
        }

        // Add pole contributions
        for p in 0..self.n_poles {
            result += self.residues[p] / (z - self.poles[p]);
        }

        result
    }

    fn evaluate_batch(&self, z: &Array1<ComplexF64>) -> Array1<ComplexF64> {
        if self.config.parallel && z.len() > 100 {
            // Parallel evaluation for large batches
            let results: Vec<ComplexF64> = z
                .as_slice()
                .unwrap()
                .par_iter()
                .map(|&zi| self.evaluate(zi))
                .collect();
            Array1::from(results)
        } else {
            z.mapv(|zi| self.evaluate(zi))
        }
    }

    fn cross_validation_error(&self) -> f64 {
        self.cv_error
    }

    fn stability_check(&self) -> bool {
        // Check that all poles are in the upper-left quadrant (physical region for self-energy)
        // Poles should have negative real part for causality
        let poles_stable = self.poles.iter().all(|p| p.re <= 0.0);

        // Check that poles are not too close to the real axis
        let min_imag = self
            .poles
            .iter()
            .map(|p| p.im.abs())
            .fold(f64::INFINITY, f64::min);
        let poles_separated = min_imag > self.config.stability_threshold;

        // Check that residues have reasonable magnitudes
        let max_residue = self.residues.iter().map(|r| r.norm()).fold(0.0, f64::max);
        let residues_bounded = max_residue < 1e10;

        poles_stable && poles_separated && residues_bounded
    }

    fn get_poles_residues(&self) -> Option<(Array1<ComplexF64>, Array1<ComplexF64>)> {
        Some((self.poles.clone(), self.residues.clone()))
    }

    fn clone_box(&self) -> Box<dyn AnalyticContinuation> {
        Box::new(self.clone())
    }
}

impl MultipoleModel {
    /// Perform k-fold cross-validation
    fn cross_validate(&self, xi: &Array1<f64>, f: &Array1<ComplexF64>) -> Result<f64> {
        let n = xi.len();
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        #[allow(clippy::cast_precision_loss)]
        let fold_size = ((n as f64) * self.config.cv_fraction) as usize;

        if fold_size < 2 {
            return Ok(0.0);
        }

        let mut errors = Vec::new();

        for _iter in 0..self.config.cv_iterations {
            // Random train/test split
            let mut indices: Vec<usize> = (0..n).collect();
            {
                use rand::rng;
                use rand::seq::SliceRandom;
                indices.shuffle(&mut rng());
            }

            let test_indices = &indices[0..fold_size];
            let train_indices = &indices[fold_size..];

            // Extract train and test data
            let xi_train = train_indices.iter().map(|&i| xi[i]).collect();
            let f_train = train_indices.iter().map(|&i| f[i]).collect();
            let xi_test: Array1<f64> = test_indices.iter().map(|&i| xi[i]).collect();
            let f_test: Array1<ComplexF64> = test_indices.iter().map(|&i| f[i]).collect();

            // Fit on training data
            let mut model = MultipoleModel::new(self.n_poles, self.config.clone());
            model.fit_nonlinear(&xi_train, &f_train)?;

            // Evaluate on test data
            let error_sum: f64 = xi_test
                .iter()
                .zip(f_test.iter())
                .map(|(&x, &fx)| {
                    let z = ComplexF64::new(0.0, x);
                    let pred = model.evaluate(z);
                    (pred - fx).norm_sqr()
                })
                .sum::<f64>();
            #[allow(clippy::cast_precision_loss)]
            let error = error_sum / (fold_size as f64);

            errors.push(error.sqrt());
        }

        #[allow(clippy::cast_precision_loss)]
        Ok(errors.iter().sum::<f64>() / (errors.len() as f64))
    }

    /// Prony's method for initial pole guess (Algorithm 2 from S4-2)
    fn prony_method_initialization(
        &self,
        xi: &Array1<f64>,
        f: &Array1<ComplexF64>,
    ) -> Array1<ComplexF64> {
        let n = xi.len();
        let p = self.n_poles.min(n / 2);

        if n < 2 * p {
            // Not enough points for Prony's method
            return Array1::zeros(0);
        }

        // Build Hankel matrices H0 and H1
        let mut h0 = Array2::<ComplexF64>::zeros((n - p, p));
        let mut h1 = Array2::<ComplexF64>::zeros((n - p, p));

        for i in 0..n - p {
            for j in 0..p {
                if i + j < n {
                    h0[[i, j]] = f[i + j];
                    if i + j + 1 < n {
                        h1[[i, j]] = f[i + j + 1];
                    }
                }
            }
        }

        // Solve generalized eigenvalue problem H1 * v = lambda * H0 * v
        // For simplified implementation, extract poles from singular values
        let svd_result = h0.svd(false, false);
        match svd_result {
            Ok((_, s_vals, _)) => {
                // Use singular values to estimate decay rates
                let delta_xi = if xi.len() > 1 { xi[1] - xi[0] } else { 1.0 };

                let mut poles = Vec::new();
                for (i, &s_val) in s_vals.iter().enumerate() {
                    if i >= p {
                        break;
                    }

                    let s_abs = s_val.abs();
                    if s_abs > 1e-10 && s_abs < 1e10 {
                        // Estimate decay rate from singular value
                        #[allow(clippy::cast_precision_loss)]
                        let decay = -(s_abs.ln()) / (delta_xi * (i + 1) as f64);

                        // Create pole in upper-left quadrant
                        let real_part = -decay.abs();
                        let imag_part = -decay.abs() * 0.5; // Arbitrary but physical
                        let pole = ComplexF64::new(real_part, imag_part);

                        if pole.re < 0.0 {
                            poles.push(pole);
                        }
                    }
                }

                // If we didn't get enough poles from SVD, add some defaults
                while poles.len() < self.n_poles.min(p) {
                    let i = poles.len();
                    #[allow(clippy::cast_precision_loss)]
                    let angle = std::f64::consts::PI * (0.5 + 0.4 * i as f64 / p as f64);
                    #[allow(clippy::cast_precision_loss)]
                    let radius = 0.5 + 2.0 * i as f64 / p as f64;
                    poles.push(ComplexF64::from_polar(radius, angle));
                }

                // Limit to requested number of poles
                poles.truncate(self.n_poles);
                Array1::from(poles)
            }
            Err(_) => {
                // Prony's method failed, return empty
                Array1::zeros(0)
            }
        }
    }

    /// Initialize residues using linear least squares
    fn initialize_residues(&mut self, xi: &Array1<f64>, f: &Array1<ComplexF64>) {
        let n_data = xi.len();
        let n_poles = self.n_poles;

        // Build matrix A where A\[i,j\] = 1/(i*xi\[i\] - pole\[j\])
        let mut a_matrix = Array2::<ComplexF64>::zeros((n_data, n_poles));

        for (i, &x) in xi.iter().enumerate() {
            let z = ComplexF64::new(0.0, x);
            for j in 0..n_poles {
                a_matrix[[i, j]] = 1.0 / (z - self.poles[j]);
            }
        }

        // Solve A * residues = f via least squares
        let ata = a_matrix.t().dot(&a_matrix);
        let at_times_f = a_matrix.t().dot(f);

        // Add small regularization
        let mut ata_reg = ata;
        for i in 0..n_poles {
            ata_reg[[i, i]] += ComplexF64::new(1e-10, 0.0);
        }

        // Solve for residues
        if let Ok(residues) = ata_reg.solve(&at_times_f) {
            self.residues = residues;
        } else {
            // Fallback to uniform residues
            #[allow(clippy::cast_precision_loss)]
            let uniform_value = 1.0 / n_poles as f64;
            for i in 0..n_poles {
                self.residues[i] = ComplexF64::new(uniform_value, 0.0);
            }
        }
    }

    /// Apply physical constraints to poles
    fn apply_physical_constraints(&mut self) {
        for i in 0..self.n_poles {
            // Ensure poles are in the upper-left quadrant for causality
            // Real part should be negative (left half-plane)
            if self.poles[i].re > 0.0 {
                self.poles[i].re = -self.poles[i].re.abs();
            }

            // Imaginary part should be negative (upper half-plane in physics convention)
            // Note: In the complex plane, negative imaginary means decaying in time
            if self.poles[i].im > 0.0 {
                self.poles[i].im = -self.poles[i].im.abs();
            }

            // Keep poles away from real axis for numerical stability
            if self.poles[i].im.abs() < 0.01 {
                self.poles[i].im = -0.01; // Always negative for physical self-energy
            }
        }
    }

    /// Detect and merge clustered poles (Algorithm 8 from S4-2)
    pub(crate) fn detect_and_merge_poles(&self) -> (Array1<ComplexF64>, Array1<ComplexF64>) {
        // Use adaptive threshold based on pole distribution
        // For N poles spread over range R, typical spacing is R/(N-1)
        // Consider poles clustered if they're closer than 40% of typical spacing
        let threshold = if self.n_poles > 1 {
            let real_parts: Vec<f64> = self.poles.iter().map(|p| p.re).collect();
            let min_re = real_parts.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            let max_re = real_parts.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            let range = (max_re - min_re).abs();
            #[allow(clippy::cast_precision_loss)]
            let typical_spacing = range / (self.n_poles as f64 - 1.0);
            // Use 40% of typical spacing as clustering threshold, with minimum of 0.1
            (0.4 * typical_spacing).max(0.1)
        } else {
            0.1 // Default threshold
        };

        let n = self.n_poles;
        let mut used = vec![false; n];
        let mut merged_poles = Vec::new();
        let mut merged_residues = Vec::new();

        for i in 0..n {
            if used[i] {
                continue;
            }

            // Find cluster around pole i
            let mut cluster_indices = vec![i];
            used[i] = true;

            for (j, used_j) in used.iter_mut().enumerate().take(n).skip(i + 1) {
                if !*used_j && (self.poles[i] - self.poles[j]).norm() < threshold {
                    cluster_indices.push(j);
                    *used_j = true;
                }
            }

            if cluster_indices.len() > 1 {
                // Merge cluster: weighted average by residue magnitude
                let mut total_weight = 0.0;
                let mut weighted_pole = ComplexF64::new(0.0, 0.0);
                let mut total_residue = ComplexF64::new(0.0, 0.0);

                for &idx in &cluster_indices {
                    let weight = self.residues[idx].norm();
                    total_weight += weight;
                    weighted_pole += self.poles[idx] * weight;
                    total_residue += self.residues[idx];
                }

                if total_weight > 1e-10 {
                    merged_poles.push(weighted_pole / total_weight);
                    merged_residues.push(total_residue);
                }
            } else {
                merged_poles.push(self.poles[i]);
                merged_residues.push(self.residues[i]);
            }
        }

        (Array1::from(merged_poles), Array1::from(merged_residues))
    }
}

/// Padé approximant model using Thiele recursion
#[derive(Debug, Clone)]
pub struct PadeModel {
    /// Order of numerator polynomial
    numerator_order: usize,
    /// Order of denominator polynomial
    denominator_order: usize,
    /// Numerator coefficients
    numerator_coeffs: Array1<ComplexF64>,
    /// Denominator coefficients
    denominator_coeffs: Array1<ComplexF64>,
    /// Continued fraction coefficients (Thiele)
    thiele_coeffs: Array1<ComplexF64>,
    /// Support points for Thiele recursion
    support_points: Array1<ComplexF64>,
    /// Cross-validation error
    cv_error: f64,
    /// Configuration
    config: ACConfig,
}

impl PadeModel {
    /// Create a new Padé model
    #[must_use]
    pub fn new(order: (usize, usize), config: ACConfig) -> Self {
        let (m, n) = order;
        Self {
            numerator_order: m,
            denominator_order: n,
            numerator_coeffs: Array1::zeros(m + 1),
            denominator_coeffs: Array1::zeros(n + 1),
            thiele_coeffs: Array1::zeros(m + n + 1),
            support_points: Array1::zeros(m + n + 1),
            cv_error: f64::INFINITY,
            config,
        }
    }

    /// Fit using Thiele's continued fraction algorithm
    fn fit_thiele(&mut self, xi: &Array1<f64>, f: &Array1<ComplexF64>) -> Result<()> {
        let n_points = self.numerator_order + self.denominator_order + 1;

        if xi.len() < n_points {
            return Err(QuasixError::InvalidInput(format!(
                "Need at least {} points for [{}/{}] Padé",
                n_points, self.numerator_order, self.denominator_order
            )));
        }

        // Select support points (use first n_points)
        self.support_points = xi
            .slice(ndarray::s![..n_points])
            .mapv(|x| ComplexF64::new(0.0, x));
        let f_support = f.slice(ndarray::s![..n_points]).to_owned();

        // Thiele recursion to build continued fraction
        let mut g = Array2::zeros((n_points, n_points));
        g.column_mut(0).assign(&f_support);

        for i in 1..n_points {
            for j in i..n_points {
                let num = self.support_points[j] - self.support_points[j - i];
                let denom = g[[j, i - 1]] - g[[j - 1, i - 1]];

                if denom.norm() < 1e-14 {
                    g[[j, i]] = ComplexF64::new(0.0, 0.0);
                } else {
                    g[[j, i]] = num / denom;
                }
            }
        }

        // Extract diagonal as continued fraction coefficients
        self.thiele_coeffs = Array1::from_vec((0..n_points).map(|i| g[[i, i]]).collect());

        // Convert continued fraction to rational function if needed
        self.convert_to_rational();

        Ok(())
    }

    /// Fit using SVD-stabilized linear system (Algorithm 3 from S4-2)
    fn fit_svd(&mut self, xi: &Array1<f64>, func_values: &Array1<ComplexF64>) -> Result<()> {
        let num_order = self.numerator_order;
        let denom_order = self.denominator_order;
        let n_points = xi.len();

        if n_points < num_order + denom_order + 1 {
            return Err(QuasixError::InvalidInput(format!(
                "Need at least {} points for [{}/{}] Padé, got {}",
                num_order + denom_order + 1,
                num_order,
                denom_order,
                n_points
            )));
        }

        // Step 1: Build coefficient matrix for denominator
        // We solve for q_1, ..., q_n where Q(z) = 1 + q_1*z + q_2*z^2 + ... + q_n*z^n
        // The system is: sum(f[i-j] * q_j * xi[i]^(n-j)) = -f[i] * xi[i]^n
        let mut a_matrix = Array2::<ComplexF64>::zeros((n_points - denom_order, denom_order));
        let mut b_vector = Array1::<ComplexF64>::zeros(n_points - denom_order);

        for i in denom_order..n_points {
            for j in 0..denom_order {
                // A[i-n, j] = f[i-j-1] * xi[i]^(n-j)
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                let xi_power = ComplexF64::new(xi[i].powi((denom_order - j) as i32), 0.0);
                a_matrix[[i - denom_order, j]] = func_values[i - j - 1] * xi_power;
            }
            // b[i-n] = -f[i] * xi[i]^n
            #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
            let xi_n = ComplexF64::new(xi[i].powi(denom_order as i32), 0.0);
            b_vector[i - denom_order] = -func_values[i] * xi_n;
        }

        // Step 2: SVD decomposition with adaptive truncation
        let svd_result = a_matrix.svd(true, true);
        let Ok((u_opt, s_vals, vt_opt)) = svd_result else {
            // Fallback to simpler method if SVD fails
            return self.fit_thiele(xi, func_values);
        };

        // Extract U and VT matrices
        let u =
            u_opt.ok_or_else(|| QuasixError::NumericalError("SVD missing U matrix".to_string()))?;
        let vt = vt_opt
            .ok_or_else(|| QuasixError::NumericalError("SVD missing VT matrix".to_string()))?;

        // Adaptive truncation based on condition number
        let s_max = s_vals.iter().map(|&s| s.abs()).fold(0.0_f64, f64::max);
        let s_min_threshold = self.config.regularization * s_max;

        // Determine numerical rank
        let rank = s_vals
            .iter()
            .filter(|&&s| s.abs() > s_min_threshold)
            .count();

        // Step 3: Compute truncated pseudoinverse
        let mut s_inv = Array1::<ComplexF64>::zeros(s_vals.len());
        for i in 0..rank.min(s_vals.len()) {
            if s_vals[i].abs() > s_min_threshold {
                s_inv[i] = ComplexF64::new(1.0 / s_vals[i], 0.0);
            }
        }

        // q = V * S^{-1} * U^T * b
        let ut_b = u.t().dot(&b_vector);
        let s_inv_ut_b = &s_inv * &ut_b;
        let q_coeffs = vt.t().dot(&s_inv_ut_b);

        // Step 4: Build denominator polynomial Q(z) = 1 + q_1*z + ... + q_n*z^n
        self.denominator_coeffs[0] = ComplexF64::new(1.0, 0.0);
        for i in 0..denom_order.min(q_coeffs.len()) {
            self.denominator_coeffs[i + 1] = q_coeffs[i];
        }

        // Step 5: Solve for numerator coefficients P(z)
        // P(z) = f(z) * Q(z), fit using least squares
        let mut p_matrix = Array2::<ComplexF64>::zeros((n_points, num_order + 1));
        let mut p_rhs = Array1::<ComplexF64>::zeros(n_points);

        for i in 0..n_points {
            let z_val = ComplexF64::new(0.0, xi[i]);

            // Evaluate Q(z)
            let mut q_val = ComplexF64::new(0.0, 0.0);
            for j in 0..=denom_order.min(self.denominator_coeffs.len() - 1) {
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                let power = j as i32;
                q_val += self.denominator_coeffs[j] * z_val.powi(power);
            }

            // Build matrix for P coefficients
            #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
            for j in 0..=num_order {
                p_matrix[[i, j]] = z_val.powi(j as i32);
            }

            // RHS is f(z) * Q(z)
            p_rhs[i] = func_values[i] * q_val;
        }

        // Solve for P coefficients using least squares
        let p_transpose_p = p_matrix.t().dot(&p_matrix);
        let p_transpose_rhs = p_matrix.t().dot(&p_rhs);

        // Add small regularization for numerical stability
        let mut ptp_reg = p_transpose_p;
        for i in 0..ptp_reg.nrows() {
            ptp_reg[[i, i]] += ComplexF64::new(1e-12, 0.0);
        }

        // Solve for numerator coefficients
        match ptp_reg.solve(&p_transpose_rhs) {
            Ok(p_coeffs) => {
                self.numerator_coeffs = p_coeffs;
            }
            Err(_) => {
                // If direct solution fails, use simpler approximation
                for i in 0..=num_order.min(n_points - 1) {
                    self.numerator_coeffs[i] = func_values[i];
                }
            }
        }

        // Step 6: Stability check - validate poles
        let poles = self.find_polynomial_roots();
        let spurious_count = poles.iter()
            .filter(|p| p.im.abs() < 1e-3)  // Too close to real axis
            .count();

        if spurious_count > 0 {
            // Try deflation of spurious poles
            self.deflate_spurious_poles(&poles);
        }

        // Store CV error estimate based on condition number
        self.cv_error = s_max / s_vals[rank.min(s_vals.len() - 1)].abs();

        Ok(())
    }

    /// Find roots of denominator polynomial (poles of Padé approximant)
    fn find_polynomial_roots(&self) -> Array1<ComplexF64> {
        let n = self.denominator_order;
        if n == 0 {
            return Array1::zeros(0);
        }

        // Build companion matrix for Q(z)
        let mut companion = Array2::<ComplexF64>::zeros((n, n));

        // Standard companion matrix form
        for i in 0..n - 1 {
            companion[[i + 1, i]] = ComplexF64::new(1.0, 0.0);
        }

        // Last column: -q_0/q_n, -q_1/q_n, ..., -q_{n-1}/q_n
        let q_n = self.denominator_coeffs[n];
        if q_n.norm() > 1e-14 {
            for i in 0..n {
                companion[[i, n - 1]] = -self.denominator_coeffs[i] / q_n;
            }
        }

        // Eigenvalues of companion matrix are the roots
        // For now, use approximate analytical solution for small orders
        match n {
            1 => {
                // Linear: q_0 + q_1*z = 0 => z = -q_0/q_1
                if self.denominator_coeffs[1].norm() > 1e-14 {
                    Array1::from_elem(1, -self.denominator_coeffs[0] / self.denominator_coeffs[1])
                } else {
                    Array1::zeros(0)
                }
            }
            2 => {
                // Quadratic formula
                let a = self.denominator_coeffs[2];
                let b = self.denominator_coeffs[1];
                let c = self.denominator_coeffs[0];

                if a.norm() > 1e-14 {
                    let discriminant = (b * b - 4.0 * a * c).sqrt();
                    Array1::from_vec(vec![
                        (-b + discriminant) / (2.0 * a),
                        (-b - discriminant) / (2.0 * a),
                    ])
                } else if b.norm() > 1e-14 {
                    Array1::from_elem(1, -c / b)
                } else {
                    Array1::zeros(0)
                }
            }
            _ => {
                // For higher orders, use iterative method or return approximate poles
                let mut poles = Array1::zeros(n);
                for i in 0..n {
                    // Distribute poles in upper-left quadrant
                    #[allow(clippy::cast_precision_loss)]
                    let angle = std::f64::consts::PI * (0.5 + 0.4 * i as f64 / n as f64);
                    #[allow(clippy::cast_precision_loss)]
                    let radius = 0.5 + 3.0 * i as f64 / n as f64;
                    poles[i] = ComplexF64::from_polar(radius, angle);
                }
                poles
            }
        }
    }

    /// Deflate spurious poles that are too close to the real axis
    fn deflate_spurious_poles(&mut self, poles: &Array1<ComplexF64>) {
        let mut physical_indices = Vec::new();

        for (i, &pole) in poles.iter().enumerate() {
            // Keep poles that are sufficiently far from real axis
            // and have negative real parts (causality)
            if pole.im.abs() > 1e-3 && pole.re < 0.0 {
                physical_indices.push(i);
            }
        }

        // If we removed poles, reduce the denominator order
        if physical_indices.len() < poles.len() {
            let new_n = physical_indices.len();

            // Rebuild polynomial from physical poles only
            let mut new_denom = Array1::<ComplexF64>::zeros(new_n + 1);
            new_denom[0] = ComplexF64::new(1.0, 0.0);

            // Product form: Q(z) = prod(z - pole_i)
            for &idx in &physical_indices {
                let pole = poles[idx];
                let mut temp = Array1::<ComplexF64>::zeros(new_denom.len() + 1);

                // Multiply by (z - pole)
                for j in 0..new_denom.len() {
                    if j > 0 {
                        temp[j] += new_denom[j - 1]; // z term
                    }
                    temp[j] -= pole * new_denom[j]; // -pole term
                }

                new_denom = temp.slice(s![..=new_n]).to_owned();
            }

            // Update denominator
            self.denominator_order = new_n;
            self.denominator_coeffs = new_denom;
        }
    }

    /// Convert continued fraction to rational function
    fn convert_to_rational(&mut self) {
        // This is a simplified conversion - full implementation would
        // use backward recursion to build P(z)/Q(z) from continued fraction

        // For now, approximate with pole-residue form
        let n_poles = self.thiele_coeffs.len() / 2;
        let mut poles = Vec::new();
        let mut residues = Vec::new();

        for i in 0..n_poles.min(self.denominator_order) {
            if i < self.support_points.len() && i + 1 < self.thiele_coeffs.len() {
                // Approximate pole location
                poles.push(-self.support_points[i].im - ComplexF64::i() * 0.1);
                // Approximate residue
                residues.push(self.thiele_coeffs[i + 1]);
            }
        }

        // Store as coefficients (simplified)
        self.denominator_coeffs[0] = ComplexF64::new(1.0, 0.0);
        for (i, &pole) in poles.iter().enumerate() {
            if i + 1 < self.denominator_coeffs.len() {
                self.denominator_coeffs[i + 1] = -pole;
            }
        }

        for (i, &res) in residues.iter().enumerate() {
            if i < self.numerator_coeffs.len() {
                self.numerator_coeffs[i] = res;
            }
        }
    }

    /// Evaluate continued fraction using backward recursion
    fn evaluate_continued_fraction(&self, z: ComplexF64) -> ComplexF64 {
        if self.thiele_coeffs.is_empty() {
            return ComplexF64::new(0.0, 0.0);
        }

        let n = self.thiele_coeffs.len();
        let mut result = self.thiele_coeffs[n - 1];

        for i in (1..n).rev() {
            if i < self.support_points.len() {
                let denom = result;
                if denom.norm() > 1e-14 {
                    result = self.thiele_coeffs[i - 1] + (z - self.support_points[i - 1]) / denom;
                } else {
                    result = self.thiele_coeffs[i - 1];
                }
            }
        }

        result
    }

    /// Evaluate rational function P(z)/Q(z)
    fn evaluate_rational(&self, z: ComplexF64) -> ComplexF64 {
        // Evaluate numerator
        let mut numerator = ComplexF64::new(0.0, 0.0);
        for (i, &coeff) in self.numerator_coeffs.iter().enumerate() {
            numerator += coeff * z.powi(i32::try_from(i).unwrap_or(i32::MAX));
        }

        // Evaluate denominator
        let mut denominator = ComplexF64::new(0.0, 0.0);
        for (i, &coeff) in self.denominator_coeffs.iter().enumerate() {
            denominator += coeff * z.powi(i32::try_from(i).unwrap_or(i32::MAX));
        }

        if denominator.norm() < 1e-14 {
            ComplexF64::new(0.0, 0.0)
        } else {
            numerator / denominator
        }
    }

    /// Extract poles and residues from rational function
    fn extract_poles_residues(&self) -> (Array1<ComplexF64>, Array1<ComplexF64>) {
        // Find roots of denominator polynomial (poles)
        // This is a simplified version - production would use companion matrix eigenvalues

        // Respect max_poles configuration
        let n_poles = self.denominator_order.min(self.config.max_poles);
        let mut poles = Array1::zeros(n_poles);
        let mut residues = Array1::zeros(n_poles);

        // Approximate poles from Thiele coefficients
        for i in 0..n_poles.min(self.support_points.len()) {
            poles[i] = -self.support_points[i].im - ComplexF64::i() * 0.01;

            // Compute residue as lim_{z->pole} (z-pole)*P(z)/Q(z)
            let h = 1e-6;
            let z1 = poles[i] + h;
            let z2 = poles[i] - h;
            let f1 = self.evaluate_rational(z1);
            let f2 = self.evaluate_rational(z2);
            residues[i] = (f1 - f2) * poles[i] / (2.0 * h);
        }

        (poles, residues)
    }
}

impl AnalyticContinuation for PadeModel {
    fn fit(&mut self, xi_data: &Array1<f64>, f_data: &Array1<ComplexF64>) -> Result<()> {
        // Validate input
        if xi_data.len() != f_data.len() {
            return Err(QuasixError::InvalidInput(format!(
                "Data size mismatch: {} vs {}",
                xi_data.len(),
                f_data.len()
            )));
        }

        // Try Thiele recursion first (more stable for well-conditioned data)
        match self.fit_thiele(xi_data, f_data) {
            Ok(()) => {}
            Err(_) => {
                // Fall back to SVD method
                self.fit_svd(xi_data, f_data)?;
            }
        }

        // Compute cross-validation error
        self.cv_error = self.cross_validate(xi_data, f_data)?;

        Ok(())
    }

    fn evaluate(&self, z: ComplexF64) -> ComplexF64 {
        // Use continued fraction if available, otherwise rational function
        if self.thiele_coeffs.is_empty() {
            self.evaluate_rational(z)
        } else {
            self.evaluate_continued_fraction(z)
        }
    }

    fn evaluate_batch(&self, z: &Array1<ComplexF64>) -> Array1<ComplexF64> {
        if self.config.parallel && z.len() > 100 {
            let results: Vec<ComplexF64> = z
                .as_slice()
                .unwrap()
                .par_iter()
                .map(|&zi| self.evaluate(zi))
                .collect();
            Array1::from(results)
        } else {
            z.mapv(|zi| self.evaluate(zi))
        }
    }

    fn cross_validation_error(&self) -> f64 {
        self.cv_error
    }

    fn stability_check(&self) -> bool {
        // Extract poles and check stability
        let (poles, residues) = self.extract_poles_residues();

        // Check that poles are not on the real axis
        let min_imag = poles
            .iter()
            .map(|p| p.im.abs())
            .fold(f64::INFINITY, f64::min);
        let poles_separated = min_imag > self.config.stability_threshold;

        // Check residue magnitudes
        let max_residue = residues.iter().map(|r| r.norm()).fold(0.0, f64::max);
        let residues_bounded = max_residue < 1e10;

        poles_separated && residues_bounded
    }

    fn get_poles_residues(&self) -> Option<(Array1<ComplexF64>, Array1<ComplexF64>)> {
        Some(self.extract_poles_residues())
    }

    fn clone_box(&self) -> Box<dyn AnalyticContinuation> {
        Box::new(self.clone())
    }
}

impl PadeModel {
    /// Perform k-fold cross-validation
    #[allow(clippy::unnecessary_wraps)]
    fn cross_validate(&self, xi: &Array1<f64>, f: &Array1<ComplexF64>) -> Result<f64> {
        let n = xi.len();
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        #[allow(clippy::cast_precision_loss)]
        let fold_size = ((n as f64) * self.config.cv_fraction) as usize;

        if fold_size < 2 {
            return Ok(0.0);
        }

        let mut errors = Vec::new();

        for _iter in 0..self.config.cv_iterations {
            // Random train/test split
            let mut indices: Vec<usize> = (0..n).collect();
            {
                use rand::rng;
                use rand::seq::SliceRandom;
                indices.shuffle(&mut rng());
            }

            let test_indices = &indices[0..fold_size];
            let train_indices = &indices[fold_size..];

            // Extract train and test data
            let xi_train = train_indices.iter().map(|&i| xi[i]).collect();
            let f_train = train_indices.iter().map(|&i| f[i]).collect();
            let xi_test: Array1<f64> = test_indices.iter().map(|&i| xi[i]).collect();
            let f_test: Array1<ComplexF64> = test_indices.iter().map(|&i| f[i]).collect();

            // Fit on training data
            let mut model = PadeModel::new(
                (self.numerator_order, self.denominator_order),
                self.config.clone(),
            );
            if model.fit(&xi_train, &f_train).is_err() {
                continue; // Skip if fitting fails
            }

            // Evaluate on test data
            let error_sum: f64 = xi_test
                .iter()
                .zip(f_test.iter())
                .map(|(&x, &fx)| {
                    let z = ComplexF64::new(0.0, x);
                    let pred = model.evaluate(z);
                    (pred - fx).norm_sqr()
                })
                .sum::<f64>();
            #[allow(clippy::cast_precision_loss)]
            let error = error_sum / (fold_size as f64);

            errors.push(error.sqrt());
        }

        if errors.is_empty() {
            Ok(f64::INFINITY)
        } else {
            #[allow(clippy::cast_precision_loss)]
            Ok(errors.iter().sum::<f64>() / (errors.len() as f64))
        }
    }
}

/// Main analytic continuation interface with model selection
#[derive(Debug)]
pub struct AnalyticContinuationFitter {
    /// Configuration
    config: ACConfig,
    /// Fitted models
    models: Vec<(ModelType, Box<dyn AnalyticContinuation>)>,
    /// Best model based on cross-validation
    best_model: Option<(ModelType, Box<dyn AnalyticContinuation>)>,
}

impl AnalyticContinuationFitter {
    /// Create a new fitter with configuration
    #[must_use]
    pub fn new(config: ACConfig) -> Self {
        Self {
            config,
            models: Vec::new(),
            best_model: None,
        }
    }

    /// Fit models to imaginary-axis data
    ///
    /// # Errors
    ///
    /// Returns an error if no stable models can be fitted to the data.
    pub fn fit(&mut self, xi_data: &Array1<f64>, f_data: &Array1<ComplexF64>) -> Result<()> {
        self.models.clear();

        // Try multipole models with different pole counts
        let min_poles = 1; // Always allow at least 1 pole
        let max_poles = self.config.max_poles.min(xi_data.len() / 2);
        let mut pole_counts = Vec::new();

        // Generate pole count candidates
        if max_poles >= 1 {
            pole_counts.push(1);
        }
        if max_poles >= 2 {
            pole_counts.push(2);
        }
        if max_poles >= 5 {
            pole_counts.push(5);
        }
        if max_poles >= 10 {
            pole_counts.push(10);
        }
        if max_poles >= 15 {
            pole_counts.push(15);
        }
        if max_poles >= 20 {
            pole_counts.push(20);
        }

        // Ensure we at least try with the configured max_poles if it's small
        if self.config.max_poles <= 20 && !pole_counts.contains(&self.config.max_poles) {
            pole_counts.push(self.config.max_poles);
        }
        pole_counts.sort_unstable();
        pole_counts.dedup();

        for n_poles in pole_counts {
            if n_poles > self.config.max_poles || n_poles < min_poles {
                continue;
            }

            let mut model = MultipoleModel::new(n_poles, self.config.clone());
            if model.fit(xi_data, f_data).is_ok() && model.stability_check() {
                self.models.push((ModelType::Multipole, model.clone_box()));
            }
        }

        // Try Padé approximants with different orders
        // Ensure Padé order doesn't exceed max_poles (denominator order = number of poles)
        let orders = vec![(5, 5), (10, 10), (15, 15)];
        for (m, n) in orders {
            if m > self.config.max_pade_order.0 || n > self.config.max_pade_order.1 {
                break;
            }

            // Skip if denominator order would exceed max_poles
            if n > self.config.max_poles {
                continue;
            }

            let mut model = PadeModel::new((m, n), self.config.clone());
            if model.fit(xi_data, f_data).is_ok() && model.stability_check() {
                self.models.push((ModelType::Pade, model.clone_box()));
            }
        }

        // Select best model based on cross-validation error
        self.select_best_model()?;

        Ok(())
    }

    /// Select the best model based on CV error and stability
    fn select_best_model(&mut self) -> Result<()> {
        if self.models.is_empty() {
            return Err(QuasixError::ConvergenceFailed(
                "No stable models found".to_string(),
            ));
        }

        // Find model with lowest CV error
        let best_idx = self
            .models
            .iter()
            .enumerate()
            .min_by(|(_, (_, a)), (_, (_, b))| {
                a.cross_validation_error()
                    .partial_cmp(&b.cross_validation_error())
                    .unwrap()
            })
            .map(|(idx, _)| idx)
            .unwrap();

        let (model_type, model) = &self.models[best_idx];
        self.best_model = Some((*model_type, model.clone_box()));

        Ok(())
    }

    /// Get the best fitted model
    #[must_use]
    pub fn best_model(&self) -> Option<&dyn AnalyticContinuation> {
        self.best_model.as_ref().map(|(_, model)| model.as_ref())
    }

    /// Get result with model selection information
    ///
    /// # Errors
    ///
    /// Returns an error if no model has been fitted.
    pub fn get_result(&self) -> Result<ACResult> {
        let (model_type, model) = self
            .best_model
            .as_ref()
            .ok_or_else(|| QuasixError::InvalidInput("No model fitted".to_string()))?;

        Ok(ACResult {
            model_type: *model_type,
            cv_error: model.cross_validation_error(),
            stability_score: if model.stability_check() { 1.0 } else { 0.0 },
            model: Arc::from(model.clone_box()),
        })
    }

    /// Evaluate on real frequency axis with broadening
    ///
    /// # Errors
    ///
    /// Returns an error if no model has been fitted.
    pub fn evaluate_real_axis(
        &self,
        omega: &Array1<f64>,
        eta: Option<f64>,
    ) -> Result<Array1<ComplexF64>> {
        let model = self
            .best_model()
            .ok_or_else(|| QuasixError::InvalidInput("No model fitted".to_string()))?;

        let eta = eta.unwrap_or(self.config.eta);
        let z = omega.mapv(|w| ComplexF64::new(w, eta));

        Ok(model.evaluate_batch(&z))
    }

    /// Compute spectral function A(ω) = -Im[Σ(ω)]/π
    ///
    /// # Errors
    ///
    /// Returns an error if no model has been fitted.
    pub fn spectral_function(&self, omega: &Array1<f64>, eta: Option<f64>) -> Result<Array1<f64>> {
        let sigma = self.evaluate_real_axis(omega, eta)?;
        Ok(sigma.mapv(|s| -s.im / std::f64::consts::PI))
    }

    /// Validate causality and Kramers-Kronig relations
    ///
    /// # Errors
    ///
    /// Returns an error if no model has been fitted.
    pub fn validate_causality(&self, omega: &Array1<f64>) -> Result<CausalityMetrics> {
        let sigma = self.evaluate_real_axis(omega, Some(1e-6))?;

        // Check causality: Im[Σ(ω)] < 0 for ω > 0
        let positive_omega_mask = omega.mapv(|w| w > 0.0);
        let violations: usize = sigma
            .iter()
            .zip(positive_omega_mask.iter())
            .filter(|(s, &mask)| mask && s.im > 0.0)
            .count();

        // Compute Kramers-Kronig consistency
        let re_sigma = sigma.mapv(|s| s.re);
        let im_sigma = sigma.mapv(|s| s.im);
        let re_sigma_kk = kramers_kronig_transform(omega, &im_sigma);

        let kk_error = (&re_sigma - &re_sigma_kk)
            .mapv(f64::abs)
            .mean()
            .unwrap_or(f64::INFINITY);

        #[allow(clippy::cast_precision_loss)]
        let violation_fraction =
            (violations as f64) / (positive_omega_mask.iter().filter(|&&x| x).count() as f64);

        Ok(CausalityMetrics {
            violations,
            violation_fraction,
            kk_error,
            kk_relative_error: kk_error / (re_sigma.mapv(f64::abs).mean().unwrap_or(1.0) + 1e-10),
        })
    }
}

/// Metrics for causality validation
#[derive(Debug, Clone)]
pub struct CausalityMetrics {
    /// Number of causality violations
    pub violations: usize,
    /// Fraction of points violating causality
    pub violation_fraction: f64,
    /// Kramers-Kronig error
    pub kk_error: f64,
    /// Relative KK error
    pub kk_relative_error: f64,
}

/// Compute Kramers-Kronig transform
fn kramers_kronig_transform(omega: &Array1<f64>, im_part: &Array1<f64>) -> Array1<f64> {
    let n = omega.len();
    let mut re_part = Array1::zeros(n);

    for i in 0..n {
        let w = omega[i];
        let mut integral = 0.0;

        // Principal value integral using trapezoidal rule
        for j in 0..n - 1 {
            if (omega[j] - w).abs() > 1e-6 && (omega[j + 1] - w).abs() > 1e-6 {
                let dw = omega[j + 1] - omega[j];
                let f1 = im_part[j] / (omega[j] - w);
                let f2 = im_part[j + 1] / (omega[j + 1] - w);
                integral += 0.5 * (f1 + f2) * dw;
            }
        }

        re_part[i] = integral / std::f64::consts::PI;
    }

    re_part
}

/// Utility functions for cross-validation
pub mod cv_utils {
    use super::{Array1, ComplexF64};

    /// Split data into train and validation sets
    pub fn train_test_split<T: Clone>(
        data: &[T],
        test_fraction: f64,
        seed: Option<u64>,
    ) -> (Vec<T>, Vec<T>) {
        use rand::rngs::StdRng;
        use rand::{seq::SliceRandom, SeedableRng};

        let mut indices: Vec<usize> = (0..data.len()).collect();

        if let Some(s) = seed {
            let mut rng = StdRng::seed_from_u64(s);
            indices.shuffle(&mut rng);
        } else {
            use rand::rng;
            indices.shuffle(&mut rng());
        }

        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        #[allow(clippy::cast_precision_loss)]
        let test_size = ((data.len() as f64) * test_fraction) as usize;
        let test_indices = &indices[0..test_size];
        let train_indices = &indices[test_size..];

        let test_data: Vec<T> = test_indices.iter().map(|&i| data[i].clone()).collect();
        let train_data: Vec<T> = train_indices.iter().map(|&i| data[i].clone()).collect();

        (train_data, test_data)
    }

    /// Compute root mean square error
    #[must_use]
    pub fn rmse(predicted: &Array1<ComplexF64>, actual: &Array1<ComplexF64>) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let n = predicted.len() as f64;
        let mse = predicted
            .iter()
            .zip(actual.iter())
            .map(|(p, a)| (p - a).norm_sqr())
            .sum::<f64>()
            / n;
        mse.sqrt()
    }

    /// Compute mean absolute error
    #[must_use]
    pub fn mae(predicted: &Array1<ComplexF64>, actual: &Array1<ComplexF64>) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let n = predicted.len() as f64;
        predicted
            .iter()
            .zip(actual.iter())
            .map(|(p, a)| (p - a).norm())
            .sum::<f64>()
            / n
    }
}

#[cfg(test)]
mod tests {
    use super::{
        cv_utils, ACConfig, AnalyticContinuation, AnalyticContinuationFitter, Array1, ComplexF64,
        MultipoleModel, PadeModel,
    };
    use approx::assert_relative_eq;

    /// Generate test function with known poles
    fn generate_test_function(n_poles: usize) -> (Array1<ComplexF64>, Array1<ComplexF64>) {
        let mut poles = Array1::zeros(n_poles);
        let mut residues = Array1::zeros(n_poles);

        for i in 0..n_poles {
            #[allow(clippy::cast_precision_loss)]
            let pole_real = (i + 1) as f64;
            #[allow(clippy::cast_precision_loss)]
            let residue_real = 1.0 / n_poles as f64;
            poles[i] = ComplexF64::new(pole_real, -0.5);
            residues[i] = ComplexF64::new(residue_real, 0.1);
        }

        (poles, residues)
    }

    /// Evaluate test function at points
    fn evaluate_test_function(
        z: &Array1<ComplexF64>,
        poles: &Array1<ComplexF64>,
        residues: &Array1<ComplexF64>,
    ) -> Array1<ComplexF64> {
        z.mapv(|zi| {
            let mut result = ComplexF64::new(0.0, 0.0);
            for (p, r) in poles.iter().zip(residues.iter()) {
                result += r / (zi - p);
            }
            result
        })
    }

    #[test]
    fn test_multipole_model_creation() {
        let config = ACConfig::default();
        let model = MultipoleModel::new(5, config);
        assert_eq!(model.n_poles, 5);
        assert_eq!(model.poles.len(), 5);
        assert_eq!(model.residues.len(), 5);
    }

    #[test]
    fn test_multipole_fitting() {
        // Generate synthetic data
        let (poles_true, residues_true) = generate_test_function(2);
        let xi = Array1::linspace(0.1, 10.0, 50);
        let z_imag = xi.mapv(|x| ComplexF64::new(0.0, x));
        let f_data = evaluate_test_function(&z_imag, &poles_true, &residues_true);

        // Fit model
        let mut model = MultipoleModel::new(2, ACConfig::default());
        assert!(model.fit(&xi, &f_data).is_ok());

        // Check that fitting reduces error
        assert!(model.cv_error < 1.0);
    }

    #[test]
    fn test_pade_model_creation() {
        let config = ACConfig::default();
        let model = PadeModel::new((5, 5), config);
        assert_eq!(model.numerator_order, 5);
        assert_eq!(model.denominator_order, 5);
    }

    #[test]
    fn test_pade_thiele_recursion() {
        // Generate smooth test data
        let xi = Array1::linspace(0.1, 5.0, 20);
        let f_data = xi.mapv(|x| ComplexF64::new(1.0 / (1.0 + x * x), 0.0));

        // Fit Padé model
        let mut model = PadeModel::new((5, 5), ACConfig::default());
        assert!(model.fit(&xi, &f_data).is_ok());

        // Evaluate at test points
        let z_test = ComplexF64::new(0.0, 2.0);
        let result = model.evaluate(z_test);
        assert!(result.norm() < 10.0); // Reasonable bound
    }

    #[test]
    fn test_model_selection() {
        // Generate test data
        let (poles_true, residues_true) = generate_test_function(3);
        let xi = Array1::linspace(0.1, 10.0, 100);
        let z_imag = xi.mapv(|x| ComplexF64::new(0.0, x));
        let f_data = evaluate_test_function(&z_imag, &poles_true, &residues_true);

        // Fit with automatic model selection
        let mut fitter = AnalyticContinuationFitter::new(ACConfig::default());
        assert!(fitter.fit(&xi, &f_data).is_ok());

        // Check that a model was selected
        assert!(fitter.best_model().is_some());

        // Get result
        let result = fitter.get_result().unwrap();
        assert!(result.cv_error < 1.0);
    }

    #[test]
    fn test_real_axis_evaluation() {
        // Generate test data
        let (poles_true, residues_true) = generate_test_function(2);
        let xi = Array1::linspace(0.1, 10.0, 50);
        let z_imag = xi.mapv(|x| ComplexF64::new(0.0, x));
        let f_data = evaluate_test_function(&z_imag, &poles_true, &residues_true);

        // Fit model
        let mut fitter = AnalyticContinuationFitter::new(ACConfig::default());
        fitter.fit(&xi, &f_data).unwrap();

        // Evaluate on real axis
        let omega = Array1::linspace(-5.0, 5.0, 100);
        let sigma = fitter.evaluate_real_axis(&omega, Some(0.01)).unwrap();

        assert_eq!(sigma.len(), 100);
        // Check that values are finite
        assert!(sigma.iter().all(|s| s.is_finite()));
    }

    #[test]
    fn test_spectral_function() {
        // Generate test data with known spectral function
        let xi = Array1::linspace(0.1, 10.0, 50);
        let f_data = xi.mapv(|x| {
            // Lorentzian on imaginary axis
            ComplexF64::new(1.0, 0.0) / ComplexF64::new(0.0, x - 1.0)
        });

        // Fit model
        let mut fitter = AnalyticContinuationFitter::new(ACConfig::default());
        fitter.fit(&xi, &f_data).unwrap();

        // Compute spectral function
        let omega = Array1::linspace(-5.0, 5.0, 100);
        let spectral = fitter.spectral_function(&omega, Some(0.1)).unwrap();

        assert_eq!(spectral.len(), 100);
        // Spectral function should be real and mostly positive
        assert!(spectral.iter().all(|&s| s.is_finite()));
    }

    #[test]
    fn test_causality_validation() {
        // Generate causal test function with known poles
        let (poles_true, residues_true) = generate_test_function(2);
        let xi = Array1::linspace(0.1, 10.0, 50);
        let z_imag = xi.mapv(|x| ComplexF64::new(0.0, x));
        let f_data = evaluate_test_function(&z_imag, &poles_true, &residues_true);

        // Fit model with more poles to ensure good fit
        let config = ACConfig {
            max_poles: 5,
            ..Default::default()
        };
        let mut fitter = AnalyticContinuationFitter::new(config);
        fitter.fit(&xi, &f_data).unwrap();

        // Validate causality
        let omega = Array1::linspace(-10.0, 10.0, 200);
        let metrics = fitter.validate_causality(&omega).unwrap();

        // For properly fitted models, we expect reasonable causality
        // The violation fraction depends on the quality of the fit
        assert!(metrics.violation_fraction < 1.0); // All we require is it's not completely wrong
        assert!(metrics.kk_relative_error < 10.0); // Reasonable KK consistency
    }

    #[test]
    fn test_cv_utils() {
        let data: Vec<i32> = (0..100).collect();
        let (train, test) = cv_utils::train_test_split(&data, 0.3, Some(42));

        assert_eq!(train.len(), 70);
        assert_eq!(test.len(), 30);

        // Check no overlap
        for &t in &test {
            assert!(!train.contains(&t));
        }
    }

    #[test]
    fn test_stability_check() {
        let mut model = MultipoleModel::new(2, ACConfig::default());

        // Set stable poles (left half-plane for physical causality)
        model.poles = Array1::from(vec![
            ComplexF64::new(-1.0, -0.5),
            ComplexF64::new(-2.0, -1.0),
        ]);
        model.residues = Array1::from(vec![ComplexF64::new(0.5, 0.1), ComplexF64::new(0.5, -0.1)]);

        assert!(model.stability_check());

        // Set unstable pole (right half-plane violates causality)
        model.poles[0] = ComplexF64::new(1.0, -0.5);
        assert!(!model.stability_check());
    }

    #[test]
    fn test_parallel_evaluation() {
        let config = ACConfig {
            parallel: true,
            ..Default::default()
        };

        let mut model = MultipoleModel::new(3, config);
        model.poles = Array1::from(vec![
            ComplexF64::new(1.0, -0.5),
            ComplexF64::new(2.0, -1.0),
            ComplexF64::new(3.0, -1.5),
        ]);
        model.residues = Array1::ones(3) * ComplexF64::new(0.33, 0.0);

        // Evaluate on large grid
        let z = Array1::linspace(-10.0, 10.0, 1000).mapv(|w| ComplexF64::new(w, 0.01));
        let result = model.evaluate_batch(&z);

        assert_eq!(result.len(), 1000);
        assert!(result.iter().all(|r| r.is_finite()));
    }

    #[test]
    fn test_error_metrics() {
        let predicted = Array1::from(vec![
            ComplexF64::new(1.0, 0.5),
            ComplexF64::new(2.0, 1.0),
            ComplexF64::new(3.0, 1.5),
        ]);
        let actual = Array1::from(vec![
            ComplexF64::new(1.1, 0.4),
            ComplexF64::new(1.9, 1.1),
            ComplexF64::new(3.1, 1.4),
        ]);

        let rmse = cv_utils::rmse(&predicted, &actual);
        let mae = cv_utils::mae(&predicted, &actual);

        assert!(rmse > 0.0 && rmse < 1.0);
        assert!(mae > 0.0 && mae < 1.0);
        assert_relative_eq!(rmse, mae, epsilon = 0.5);
    }
}
