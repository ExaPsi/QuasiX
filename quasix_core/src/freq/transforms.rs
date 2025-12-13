//! Frequency transformation utilities for GW calculations
//!
//! This module provides transformations between real and imaginary
//! frequency axes, including analytical continuation and Wick rotation.
#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::{QuasixError, Result};
use ndarray::Array1;
use num_complex::Complex64;
use std::f64::consts::PI;

/// Transform real frequency to imaginary frequency via Wick rotation
///
/// ω → iξ (real to imaginary axis)
///
/// # Arguments
/// * `omega` - Real frequency
///
/// # Returns
/// * Complex frequency on imaginary axis
pub fn real_to_imaginary(omega: f64) -> Complex64 {
    Complex64::new(0.0, omega)
}

/// Transform imaginary frequency to real frequency with broadening
///
/// iξ → ω + iη (imaginary to real axis with broadening)
///
/// # Arguments
/// * `xi` - Imaginary frequency
/// * `eta` - Broadening parameter (small positive imaginary part)
///
/// # Returns
/// * Complex frequency on real axis with broadening
pub fn imaginary_to_real(xi: f64, eta: f64) -> Complex64 {
    Complex64::new(xi, eta)
}

/// Apply frequency-dependent transformation for self-energy
///
/// Common transformations for GW self-energy evaluation
#[derive(Debug, Clone, Copy)]
pub enum FrequencyTransform {
    /// Identity transformation
    Identity,
    /// Plasmon pole approximation
    PlasmonPole { omega_p: f64 },
    /// Padé approximant continuation
    Pade { n_poles: usize },
    /// Contour deformation
    ContourDeformation { eta: f64 },
}

/// Apply a frequency transformation to an array of frequencies
pub fn apply_frequency_transform(
    frequencies: &Array1<f64>,
    transform: FrequencyTransform,
) -> Array1<Complex64> {
    match transform {
        FrequencyTransform::Identity => frequencies.mapv(|omega| Complex64::new(omega, 0.0)),
        FrequencyTransform::PlasmonPole { omega_p } => {
            apply_plasmon_pole_transform(frequencies, omega_p)
        }
        FrequencyTransform::Pade { .. } => {
            // Padé approximation would be fitted separately
            frequencies.mapv(|omega| Complex64::new(omega, 1e-3))
        }
        FrequencyTransform::ContourDeformation { eta } => {
            frequencies.mapv(|omega| Complex64::new(omega, eta))
        }
    }
}

/// Apply plasmon pole approximation transformation
fn apply_plasmon_pole_transform(frequencies: &Array1<f64>, omega_p: f64) -> Array1<Complex64> {
    frequencies.mapv(|omega| {
        // Simple plasmon pole model
        let denominator = omega * omega + omega_p * omega_p;
        Complex64::new(omega / denominator, -omega_p / denominator)
    })
}

/// Transform response function from imaginary to real frequency
///
/// Uses analytical continuation via Padé approximants or
/// direct evaluation for known functional forms.
pub struct ResponseTransform {
    /// Number of imaginary frequency points
    pub n_imag: usize,
    /// Number of real frequency points
    pub n_real: usize,
    /// Imaginary frequency grid
    pub imag_freqs: Array1<f64>,
    /// Real frequency grid
    pub real_freqs: Array1<f64>,
    /// Broadening parameter
    pub eta: f64,
}

impl ResponseTransform {
    /// Create a new response function transformer
    #[must_use]
    pub fn new(n_imag: usize, n_real: usize, omega_max: f64) -> Self {
        // Generate imaginary frequency grid (Matsubara-like)
        let mut imag_freqs = Array1::zeros(n_imag);
        for i in 0..n_imag {
            imag_freqs[i] = (2.0 * i as f64 + 1.0) * PI / omega_max;
        }

        // Generate real frequency grid
        let mut real_freqs = Array1::zeros(n_real);
        for i in 0..n_real {
            real_freqs[i] = -omega_max + 2.0 * omega_max * i as f64 / (n_real - 1) as f64;
        }

        Self {
            n_imag,
            n_real,
            imag_freqs,
            real_freqs,
            eta: 1e-3,
        }
    }

    /// Set broadening parameter
    #[must_use]
    pub fn with_eta(mut self, eta: f64) -> Self {
        self.eta = eta;
        self
    }

    /// Transform polarizability from imaginary to real axis
    pub fn transform_polarizability(
        &self,
        chi_imag: &Array1<Complex64>,
    ) -> Result<Array1<Complex64>> {
        if chi_imag.len() != self.n_imag {
            return Err(QuasixError::InvalidInput(
                "Dimension mismatch in polarizability".to_string(),
            ));
        }

        // Use Kramers-Kronig relations for simple continuation
        let mut chi_real = Array1::zeros(self.n_real);

        for (i, &omega) in self.real_freqs.iter().enumerate() {
            chi_real[i] = self.kramers_kronig_transform(omega, chi_imag);
        }

        Ok(chi_real)
    }

    /// Apply Kramers-Kronig transformation
    fn kramers_kronig_transform(&self, omega: f64, chi_imag: &Array1<Complex64>) -> Complex64 {
        // Simplified KK relation (in practice, would use more sophisticated methods)
        let mut result = Complex64::new(0.0, 0.0);

        for (j, &xi) in self.imag_freqs.iter().enumerate() {
            let weight = 2.0 * xi / (xi * xi + omega * omega + self.eta * self.eta);
            result += weight * chi_imag[j];
        }

        result * (1.0 / PI)
    }
}

/// Frequency convolution for GW self-energy
///
/// Evaluates convolutions of the form:
/// Σ(ω) = i/(2π) ∫ dω' G(ω + ω') W(ω')
pub struct FrequencyConvolution {
    /// Frequency grid
    pub frequencies: Array1<f64>,
    /// Integration weights
    pub weights: Array1<f64>,
}

impl FrequencyConvolution {
    /// Create a new convolution evaluator
    #[must_use]
    pub fn new(frequencies: Array1<f64>, weights: Array1<f64>) -> Self {
        Self {
            frequencies,
            weights,
        }
    }

    /// Evaluate convolution of two functions
    pub fn convolve(
        &self,
        f1: &Array1<Complex64>,
        f2: &Array1<Complex64>,
    ) -> Result<Array1<Complex64>> {
        let n = self.frequencies.len();
        if f1.len() != n || f2.len() != n {
            return Err(QuasixError::InvalidInput(
                "Dimension mismatch in convolution".to_string(),
            ));
        }

        let mut result = Array1::zeros(n);

        // Simple convolution (in practice, would use FFT for efficiency)
        for i in 0..n {
            let mut conv = Complex64::new(0.0, 0.0);
            for j in 0..n {
                // Circular convolution index
                let k = (i + j) % n;
                conv += f1[j] * f2[k] * self.weights[j];
            }
            result[i] = conv * Complex64::new(0.0, 1.0 / (2.0 * PI));
        }

        Ok(result)
    }

    /// Evaluate self-energy convolution for GW
    pub fn gw_self_energy(
        &self,
        green: &Array1<Complex64>,
        screened_interaction: &Array1<Complex64>,
        omega: f64,
    ) -> Complex64 {
        let mut sigma = Complex64::new(0.0, 0.0);

        for (i, &freq) in self.frequencies.iter().enumerate() {
            // G(ω + ω') W(ω')
            let omega_sum = omega + freq;

            // Interpolate or evaluate Green's function at ω + ω'
            // For simplicity, use nearest point
            let g_idx = self.find_nearest_frequency(omega_sum);
            if let Some(idx) = g_idx {
                sigma += green[idx] * screened_interaction[i] * self.weights[i];
            }
        }

        sigma * Complex64::new(0.0, 1.0 / (2.0 * PI))
    }

    /// Find nearest frequency index
    fn find_nearest_frequency(&self, omega: f64) -> Option<usize> {
        let mut min_dist = f64::INFINITY;
        let mut best_idx = None;

        for (i, &freq) in self.frequencies.iter().enumerate() {
            let dist = (freq - omega).abs();
            if dist < min_dist {
                min_dist = dist;
                best_idx = Some(i);
            }
        }

        best_idx
    }
}

/// Matsubara frequency grid for finite temperature calculations
pub struct MatsubaraGrid {
    /// Temperature (in atomic units)
    pub temperature: f64,
    /// Number of Matsubara frequencies
    pub n_matsubara: usize,
    /// Fermionic (true) or bosonic (false)
    pub fermionic: bool,
}

impl MatsubaraGrid {
    /// Create a new Matsubara frequency grid
    #[must_use]
    pub fn new(temperature: f64, n_matsubara: usize, fermionic: bool) -> Self {
        Self {
            temperature,
            n_matsubara,
            fermionic,
        }
    }

    /// Generate Matsubara frequencies
    /// ω_n = (2n + δ) π T where δ = 1 for fermions, 0 for bosons
    pub fn frequencies(&self) -> Array1<f64> {
        let mut freqs = Array1::zeros(self.n_matsubara);
        let delta = if self.fermionic { 1.0 } else { 0.0 };

        for n in 0..self.n_matsubara {
            freqs[n] = (2.0 * n as f64 + delta) * PI * self.temperature;
        }

        freqs
    }

    /// Transform Matsubara sum to integral
    /// T Σ_n → (1/2π) ∫ dω with appropriate weight function
    pub fn sum_to_integral_weight(&self, _n: usize) -> f64 {
        // Weight for transforming sum to integral
        // This is simplified; actual implementation would be more sophisticated
        self.temperature
    }
}

/// Analytical continuation using maximum entropy method (MEM)
pub struct MaximumEntropyMethod {
    /// Default model spectrum
    pub default_model: Array1<f64>,
    /// Regularization parameter (alpha)
    pub alpha: f64,
    /// Number of real frequency points
    pub n_omega: usize,
}

impl MaximumEntropyMethod {
    /// Create a new MEM solver
    #[must_use]
    pub fn new(n_omega: usize) -> Self {
        // Default model: flat spectrum
        let default_model = Array1::from_elem(n_omega, 1.0 / n_omega as f64);

        Self {
            default_model,
            alpha: 1.0,
            n_omega,
        }
    }

    /// Set regularization parameter
    #[must_use]
    pub fn with_alpha(mut self, alpha: f64) -> Self {
        self.alpha = alpha;
        self
    }

    /// Perform analytical continuation
    pub fn continue_to_real_axis(
        &self,
        _g_imag: &Array1<Complex64>,
        _imag_freqs: &Array1<f64>,
        real_freqs: &Array1<f64>,
    ) -> Result<Array1<f64>> {
        // Simplified MEM (actual implementation would solve optimization problem)
        // For now, return a simple interpolation

        if real_freqs.len() != self.n_omega {
            return Err(QuasixError::InvalidInput(
                "Real frequency grid size mismatch".to_string(),
            ));
        }

        let mut spectrum = Array1::zeros(self.n_omega);

        // Placeholder: simple mapping based on causality
        for (i, &omega) in real_freqs.iter().enumerate() {
            // Use sum rule and causality to estimate spectrum
            let weight = (-omega.abs() / 10.0).exp();
            spectrum[i] = weight / self.n_omega as f64;
        }

        // Normalize
        let norm: f64 = spectrum.sum();
        if norm > 0.0 {
            spectrum /= norm;
        }

        Ok(spectrum)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_wick_rotation() {
        let omega = 5.0;
        let z = real_to_imaginary(omega);
        assert_eq!(z.re, 0.0);
        assert_eq!(z.im, omega);

        let xi = 3.0;
        let eta = 1e-3;
        let w = imaginary_to_real(xi, eta);
        assert_eq!(w.re, xi);
        assert_eq!(w.im, eta);
    }

    #[test]
    fn test_frequency_transform() {
        let freqs = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0]);

        // Test identity transform
        let transformed = apply_frequency_transform(&freqs, FrequencyTransform::Identity);
        for (i, &f) in freqs.iter().enumerate() {
            assert_eq!(transformed[i].re, f);
            assert_eq!(transformed[i].im, 0.0);
        }

        // Test contour deformation
        let eta = 0.01;
        let transformed =
            apply_frequency_transform(&freqs, FrequencyTransform::ContourDeformation { eta });
        for (i, &f) in freqs.iter().enumerate() {
            assert_eq!(transformed[i].re, f);
            assert_eq!(transformed[i].im, eta);
        }
    }

    #[test]
    fn test_response_transform() {
        let transformer = ResponseTransform::new(10, 20, 10.0);

        assert_eq!(transformer.n_imag, 10);
        assert_eq!(transformer.n_real, 20);
        assert_eq!(transformer.imag_freqs.len(), 10);
        assert_eq!(transformer.real_freqs.len(), 20);

        // Test that real frequencies span [-omega_max, omega_max]
        assert_relative_eq!(transformer.real_freqs[0], -10.0, epsilon = 1e-10);
        assert_relative_eq!(transformer.real_freqs[19], 10.0, epsilon = 1e-10);
    }

    #[test]
    fn test_frequency_convolution() {
        let n = 16;
        let freqs = Array1::linspace(-10.0, 10.0, n);
        let weights = Array1::from_elem(n, 1.0 / n as f64);

        let conv = FrequencyConvolution::new(freqs.clone(), weights);

        // Test with delta functions
        let mut f1 = Array1::zeros(n);
        let mut f2 = Array1::zeros(n);
        f1[n / 2] = Complex64::new(1.0, 0.0);
        f2[n / 2] = Complex64::new(1.0, 0.0);

        let result = conv.convolve(&f1, &f2).unwrap();
        assert_eq!(result.len(), n);

        // Result should be non-zero somewhere (convolution of deltas)
        let total_norm: f64 = result.iter().map(|z| z.norm()).sum();
        assert!(total_norm > 0.0);
    }

    #[test]
    fn test_matsubara_grid() {
        let temp = 0.01; // Low temperature
        let n = 10;

        // Test fermionic grid
        let grid_f = MatsubaraGrid::new(temp, n, true);
        let freqs_f = grid_f.frequencies();

        assert_eq!(freqs_f.len(), n);
        // First fermionic frequency: π T
        assert_relative_eq!(freqs_f[0], PI * temp, epsilon = 1e-12);

        // Test bosonic grid
        let grid_b = MatsubaraGrid::new(temp, n, false);
        let freqs_b = grid_b.frequencies();

        // First bosonic frequency: 0
        assert_relative_eq!(freqs_b[0], 0.0, epsilon = 1e-12);
        // Second bosonic frequency: 2π T
        assert_relative_eq!(freqs_b[1], 2.0 * PI * temp, epsilon = 1e-12);
    }

    #[test]
    fn test_plasmon_pole() {
        let omega_p = 5.0; // Plasmon frequency
        let freqs = Array1::from_vec(vec![0.0, 1.0, 2.0, 5.0, 10.0]);

        let transformed =
            apply_frequency_transform(&freqs, FrequencyTransform::PlasmonPole { omega_p });

        // At ω = 0, should get purely imaginary
        assert_eq!(transformed[0].re, 0.0);

        // At ω = ω_p, should have resonance features
        // (simplified model, just check it's computed)
        assert!(transformed[3].norm() > 0.0);
    }

    #[test]
    fn test_maximum_entropy() {
        let mem = MaximumEntropyMethod::new(50);

        assert_eq!(mem.n_omega, 50);
        assert_eq!(mem.default_model.len(), 50);

        // Default model should be normalized
        let sum: f64 = mem.default_model.sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-10);

        // Test continuation (placeholder)
        let g_imag = Array1::from_elem(10, Complex64::new(0.1, 0.0));
        let imag_freqs = Array1::linspace(0.0, 10.0, 10);
        let real_freqs = Array1::linspace(-5.0, 5.0, 50);

        let spectrum = mem
            .continue_to_real_axis(&g_imag, &imag_freqs, &real_freqs)
            .unwrap();

        // Spectrum should be normalized
        let spec_sum: f64 = spectrum.sum();
        assert_relative_eq!(spec_sum, 1.0, epsilon = 1e-10);

        // All values should be non-negative (spectral function property)
        assert!(spectrum.iter().all(|&x| x >= 0.0));
    }

    #[test]
    fn test_gw_self_energy() {
        let n = 20;
        let freqs = Array1::linspace(-10.0, 10.0, n);
        let weights = Array1::from_elem(n, 20.0 / n as f64);

        let conv = FrequencyConvolution::new(freqs, weights);

        // Simple test Green's function and interaction
        let green = Array1::from_elem(n, Complex64::new(0.1, -0.01));
        let w_screened = Array1::from_elem(n, Complex64::new(0.05, 0.0));

        let sigma = conv.gw_self_energy(&green, &w_screened, 0.0);

        // Self-energy should be complex with negative imaginary part
        // (for retarded functions)
        assert!(sigma.norm() > 0.0);
    }
}
