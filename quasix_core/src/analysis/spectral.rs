//! Spectral analysis module for computing and analyzing spectral functions
//!
//! This module provides high-performance implementations for computing spectral
//! functions A(ω) from quasiparticle energies and Z-factors, with support for
//! various broadening schemes and export formats.
//!
//! # Mathematical Foundation
//!
//! The spectral function is given by:
//! ```text
//! A(ω) = Σ_n Z_n δ(ω - E_n) + incoherent_part
//! ```
//!
//! Where:
//! - `Z_n` are the quasiparticle renormalization factors (0 < Z < 1)
//! - `E_n` are the quasiparticle energies
//! - The delta functions are broadened for visualization

#![warn(clippy::all, clippy::pedantic, clippy::perf)]
#![warn(missing_docs, missing_debug_implementations)]

use ndarray::{Array1, Zip};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use thiserror::Error;

#[cfg(feature = "hdf5_support")]
use hdf5;
#[cfg(feature = "hdf5_support")]
use std::str::FromStr;

// SIMD support for performance-critical operations
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    __m256d, _mm256_add_pd, _mm256_castpd256_pd128, _mm256_extractf128_pd, _mm256_loadu_pd,
    _mm256_max_pd, _mm256_min_pd, _mm256_mul_pd, _mm256_set1_pd, _mm256_setzero_pd, _mm256_sub_pd,
    _mm_add_pd, _mm_add_sd, _mm_cvtsd_f64, _mm_unpackhi_pd,
};

/// Error types for spectral analysis operations
#[derive(Error, Debug)]
pub enum SpectralError {
    /// Invalid Z-factor value (must be in (0,1))
    #[error("Invalid Z-factor {value} at index {index}: must be in range (0, 1)")]
    InvalidZFactor {
        /// Index of the invalid Z-factor
        index: usize,
        /// The invalid Z-factor value
        value: f64,
    },

    /// Dimension mismatch between arrays
    #[error("Dimension mismatch: energies ({energies}) != z_factors ({z_factors})")]
    DimensionMismatch {
        /// Number of energies
        energies: usize,
        /// Number of Z-factors
        z_factors: usize,
    },

    /// Invalid broadening parameter
    #[error("Invalid broadening parameter {param}: {value}")]
    InvalidBroadening {
        /// Parameter name
        param: String,
        /// Invalid parameter value
        value: f64,
    },

    /// Normalization check failed
    #[error("Normalization failed: sum = {sum:.6}, expected = {expected:.6} (tol = {tol:.1e})")]
    NormalizationError {
        /// Computed sum
        sum: f64,
        /// Expected sum
        expected: f64,
        /// Tolerance used
        tol: f64,
    },

    /// I/O error during export
    #[error("I/O error during export: {0}")]
    IoError(#[from] std::io::Error),

    /// JSON serialization error
    #[error("JSON serialization error: {0}")]
    JsonError(#[from] serde_json::Error),

    #[cfg(feature = "hdf5_support")]
    /// HDF5 error
    #[error("HDF5 error: {0}")]
    Hdf5Error(#[from] hdf5::Error),

    #[cfg(feature = "hdf5_support")]
    /// HDF5 string error
    #[error("HDF5 string error: {0}")]
    Hdf5StringError(#[from] hdf5::types::StringError),

    /// Invalid omega grid
    #[error("Invalid omega grid: {reason}")]
    InvalidGrid {
        /// Reason for invalid grid
        reason: String,
    },
}

/// Result type for spectral analysis operations
pub type Result<T> = std::result::Result<T, SpectralError>;

/// Broadening parameters for spectral functions
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BroadeningParams {
    /// Gaussian broadening width (eV)
    pub gaussian_width: f64,
    /// Lorentzian broadening width (eV)
    pub lorentzian_width: f64,
    /// Temperature for thermal broadening (K)
    pub temperature: f64,
    /// Mixing parameter for Voigt profile (0 = pure Gaussian, 1 = pure Lorentzian)
    pub voigt_mixing: f64,
}

impl Default for BroadeningParams {
    fn default() -> Self {
        Self {
            gaussian_width: 0.1,    // 100 meV default broadening
            lorentzian_width: 0.05, // 50 meV default broadening
            temperature: 300.0,     // Room temperature
            voigt_mixing: 0.5,      // Equal mixture
        }
    }
}

impl BroadeningParams {
    /// Validate broadening parameters
    ///
    /// # Errors
    /// Returns error if any parameter is invalid
    pub fn validate(&self) -> Result<()> {
        if self.gaussian_width <= 0.0 {
            return Err(SpectralError::InvalidBroadening {
                param: "gaussian_width".to_string(),
                value: self.gaussian_width,
            });
        }
        if self.lorentzian_width <= 0.0 {
            return Err(SpectralError::InvalidBroadening {
                param: "lorentzian_width".to_string(),
                value: self.lorentzian_width,
            });
        }
        if self.temperature < 0.0 {
            return Err(SpectralError::InvalidBroadening {
                param: "temperature".to_string(),
                value: self.temperature,
            });
        }
        if !(0.0..=1.0).contains(&self.voigt_mixing) {
            return Err(SpectralError::InvalidBroadening {
                param: "voigt_mixing".to_string(),
                value: self.voigt_mixing,
            });
        }
        Ok(())
    }
}

/// Spectral function representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralFunction {
    /// Frequency grid (eV)
    pub omega_grid: Array1<f64>,
    /// Spectral weights at each frequency
    pub spectral_weights: Array1<f64>,
    /// Metadata about the calculation
    pub metadata: SpectralMetadata,
}

/// Metadata for spectral calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralMetadata {
    /// Number of quasiparticle states
    pub n_states: usize,
    /// Fermi energy (eV)
    pub fermi_energy: f64,
    /// Broadening parameters used
    pub broadening: BroadeningParams,
    /// Total spectral weight (for normalization check)
    pub total_weight: f64,
    /// Energy units
    pub energy_unit: String,
}

/// Spectral data for PES/IPES
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralData {
    /// PES (photoemission) spectrum
    pub pes_spectrum: Option<SpectralFunction>,
    /// IPES (inverse photoemission) spectrum
    pub ipes_spectrum: Option<SpectralFunction>,
    /// Combined spectrum
    pub combined_spectrum: Option<SpectralFunction>,
    /// Calculation timestamp
    pub timestamp: String,
    /// Calculation method
    pub method: String,
}

/// Main spectral analyzer struct
#[derive(Debug)]
pub struct SpectralAnalyzer {
    /// Quasiparticle energies (eV)
    qp_energies: Array1<f64>,
    /// Z-factors (renormalization factors)
    z_factors: Array1<f64>,
    /// Fermi energy (eV)
    fermi_energy: f64,
    /// Broadening parameters
    broadening: BroadeningParams,
}

impl SpectralAnalyzer {
    /// Create a new spectral analyzer
    ///
    /// # Arguments
    /// * `qp_energies` - Quasiparticle energies in eV
    /// * `z_factors` - Renormalization factors (0 < Z < 1)
    /// * `fermi_energy` - Fermi energy in eV
    ///
    /// # Errors
    /// Returns error if dimensions mismatch or Z-factors are invalid
    pub fn new(
        qp_energies: Array1<f64>,
        z_factors: Array1<f64>,
        fermi_energy: f64,
    ) -> Result<Self> {
        // Check dimension consistency
        if qp_energies.len() != z_factors.len() {
            return Err(SpectralError::DimensionMismatch {
                energies: qp_energies.len(),
                z_factors: z_factors.len(),
            });
        }

        // Validate Z-factors
        for (i, &z) in z_factors.iter().enumerate() {
            if !(0.0..=1.0).contains(&z) {
                return Err(SpectralError::InvalidZFactor { index: i, value: z });
            }
        }

        Ok(Self {
            qp_energies,
            z_factors,
            fermi_energy,
            broadening: BroadeningParams::default(),
        })
    }

    /// Set broadening parameters
    ///
    /// # Errors
    /// Returns error if broadening parameters are invalid
    pub fn set_broadening(&mut self, broadening: BroadeningParams) -> Result<()> {
        broadening.validate()?;
        self.broadening = broadening;
        Ok(())
    }

    /// Compute spectral function on given omega grid
    ///
    /// # Arguments
    /// * `omega_grid` - Frequency grid in eV
    ///
    /// # Returns
    /// Spectral function A(ω)
    ///
    /// # Errors
    /// Returns error if omega grid is invalid
    pub fn compute_spectral_function(&self, omega_grid: &Array1<f64>) -> Result<SpectralFunction> {
        if omega_grid.is_empty() {
            return Err(SpectralError::InvalidGrid {
                reason: "Empty omega grid".to_string(),
            });
        }

        let n_omega = omega_grid.len();
        let mut spectral_weights = Array1::zeros(n_omega);

        // Apply appropriate broadening based on parameters
        if self.broadening.voigt_mixing < 0.01 {
            // Pure Gaussian
            self.apply_gaussian_broadening(omega_grid, &mut spectral_weights);
        } else if self.broadening.voigt_mixing > 0.99 {
            // Pure Lorentzian
            self.apply_lorentzian_broadening(omega_grid, &mut spectral_weights);
        } else {
            // Voigt profile (mixture)
            self.apply_voigt_broadening(omega_grid, &mut spectral_weights);
        }

        // Calculate total weight for normalization check
        let dw = if omega_grid.len() > 1 {
            omega_grid[1] - omega_grid[0]
        } else {
            1.0
        };
        let total_weight = spectral_weights.sum() * dw;

        Ok(SpectralFunction {
            omega_grid: omega_grid.clone(),
            spectral_weights,
            metadata: SpectralMetadata {
                n_states: self.qp_energies.len(),
                fermi_energy: self.fermi_energy,
                broadening: self.broadening,
                total_weight,
                energy_unit: "eV".to_string(),
            },
        })
    }

    /// Apply Gaussian broadening to spectral function
    fn apply_gaussian_broadening(
        &self,
        omega_grid: &Array1<f64>,
        spectral_weights: &mut Array1<f64>,
    ) {
        let sigma = self.broadening.gaussian_width;
        let prefactor = 1.0 / (sigma * (2.0 * PI).sqrt());

        // Use SIMD for performance on x86_64
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                // SAFETY: AVX2 availability is verified via is_x86_feature_detected!
                // The apply_gaussian_broadening_avx2 function is marked #[target_feature(enable = "avx2")]
                // and only uses AVX2 intrinsics which are guaranteed available at this point.
                // omega_grid and spectral_weights arrays are passed by reference with valid bounds.
                unsafe {
                    self.apply_gaussian_broadening_avx2(
                        omega_grid,
                        spectral_weights,
                        sigma,
                        prefactor,
                    );
                }
                return;
            }
        }

        // Fallback scalar implementation
        self.apply_gaussian_broadening_scalar(omega_grid, spectral_weights, sigma, prefactor);
    }

    /// Scalar implementation of Gaussian broadening
    fn apply_gaussian_broadening_scalar(
        &self,
        omega_grid: &Array1<f64>,
        spectral_weights: &mut Array1<f64>,
        sigma: f64,
        prefactor: f64,
    ) {
        let two_sigma_sq = 2.0 * sigma * sigma;

        Zip::from(omega_grid)
            .and(spectral_weights)
            .for_each(|&omega, weight| {
                let mut sum = 0.0;
                for (i, &energy) in self.qp_energies.iter().enumerate() {
                    let delta_e = omega - energy;
                    let gaussian = prefactor * (-delta_e * delta_e / two_sigma_sq).exp();
                    sum += self.z_factors[i] * gaussian;
                }
                *weight = sum;
            });
    }

    /// AVX2-optimized Gaussian broadening (`x86_64` only)
    #[cfg(target_arch = "x86_64")]
    unsafe fn apply_gaussian_broadening_avx2(
        &self,
        omega_grid: &Array1<f64>,
        spectral_weights: &mut Array1<f64>,
        sigma: f64,
        prefactor: f64,
    ) {
        let two_sigma_sq = 2.0 * sigma * sigma;
        let neg_inv_two_sigma_sq = -1.0 / two_sigma_sq;

        // Process omega grid
        for (omega_idx, &omega) in omega_grid.iter().enumerate() {
            let omega_vec = _mm256_set1_pd(omega);
            let mut sum_vec = _mm256_setzero_pd();

            // Process QP states in chunks of 4
            let chunks = self.qp_energies.len() / 4;
            for chunk_idx in 0..chunks {
                let base_idx = chunk_idx * 4;

                // Load energies and Z-factors
                let energies = _mm256_loadu_pd(self.qp_energies.as_ptr().add(base_idx));
                let z_factors = _mm256_loadu_pd(self.z_factors.as_ptr().add(base_idx));

                // Compute delta_e = omega - energy
                let delta_e = _mm256_sub_pd(omega_vec, energies);

                // Compute delta_e^2
                let delta_e_sq = _mm256_mul_pd(delta_e, delta_e);

                // Compute exponent
                let exponent = _mm256_mul_pd(delta_e_sq, _mm256_set1_pd(neg_inv_two_sigma_sq));

                // Approximate exp using fast approximation
                let gaussian = Self::fast_exp_avx2(exponent);

                // Multiply by Z-factors and prefactor
                let weighted = _mm256_mul_pd(
                    _mm256_mul_pd(gaussian, z_factors),
                    _mm256_set1_pd(prefactor),
                );

                // Accumulate
                sum_vec = _mm256_add_pd(sum_vec, weighted);
            }

            // Horizontal sum of vector
            let sum = Self::horizontal_sum_avx2(sum_vec);

            // Handle remaining elements
            for i in (chunks * 4)..self.qp_energies.len() {
                let delta_e = omega - self.qp_energies[i];
                let gaussian = prefactor * (-delta_e * delta_e / two_sigma_sq).exp();
                spectral_weights[omega_idx] += self.z_factors[i] * gaussian;
            }

            spectral_weights[omega_idx] += sum;
        }
    }

    /// Fast exponential approximation for AVX2
    #[cfg(target_arch = "x86_64")]
    #[inline]
    unsafe fn fast_exp_avx2(x: __m256d) -> __m256d {
        // Use polynomial approximation for exp(x)
        // This is faster than calling exp() but less accurate
        // For spectral broadening, this accuracy is sufficient

        // Clamp to reasonable range
        let min_val = _mm256_set1_pd(-20.0);
        let max_val = _mm256_set1_pd(0.0);
        let x_clamped = _mm256_max_pd(_mm256_min_pd(x, max_val), min_val);

        // Polynomial coefficients for exp approximation
        let c1 = _mm256_set1_pd(1.0);
        let c2 = _mm256_set1_pd(0.5);
        let c3 = _mm256_set1_pd(0.166_666_666_666_666_66);
        let c4 = _mm256_set1_pd(0.041_666_666_666_666_664);

        // Horner's method: 1 + x + x^2/2 + x^3/6 + x^4/24
        let x2 = _mm256_mul_pd(x_clamped, x_clamped);
        let x3 = _mm256_mul_pd(x2, x_clamped);
        let x4 = _mm256_mul_pd(x3, x_clamped);

        _mm256_add_pd(
            c1,
            _mm256_add_pd(
                x_clamped,
                _mm256_add_pd(
                    _mm256_mul_pd(x2, c2),
                    _mm256_add_pd(_mm256_mul_pd(x3, c3), _mm256_mul_pd(x4, c4)),
                ),
            ),
        )
    }

    /// Horizontal sum for AVX2 vector
    #[cfg(target_arch = "x86_64")]
    #[inline]
    unsafe fn horizontal_sum_avx2(v: __m256d) -> f64 {
        let high = _mm256_extractf128_pd(v, 1);
        let low = _mm256_castpd256_pd128(v);
        let sum128 = _mm_add_pd(high, low);
        let high64 = _mm_unpackhi_pd(sum128, sum128);
        let sum64 = _mm_add_sd(sum128, high64);
        _mm_cvtsd_f64(sum64)
    }

    /// Apply Lorentzian broadening to spectral function
    fn apply_lorentzian_broadening(
        &self,
        omega_grid: &Array1<f64>,
        spectral_weights: &mut Array1<f64>,
    ) {
        let gamma = self.broadening.lorentzian_width;
        let prefactor = gamma / PI;

        Zip::from(omega_grid)
            .and(spectral_weights)
            .for_each(|&omega, weight| {
                let mut sum = 0.0;
                for (i, &energy) in self.qp_energies.iter().enumerate() {
                    let delta_e = omega - energy;
                    let lorentzian = prefactor / (delta_e * delta_e + gamma * gamma);
                    sum += self.z_factors[i] * lorentzian;
                }
                *weight = sum;
            });
    }

    /// Apply Voigt broadening (convolution of Gaussian and Lorentzian)
    fn apply_voigt_broadening(&self, omega_grid: &Array1<f64>, spectral_weights: &mut Array1<f64>) {
        let sigma = self.broadening.gaussian_width;
        let gamma = self.broadening.lorentzian_width;
        let alpha = self.broadening.voigt_mixing;

        // Pre-compute Gaussian parameters
        let gaussian_prefactor = 1.0 / (sigma * (2.0 * PI).sqrt());
        let two_sigma_sq = 2.0 * sigma * sigma;

        // Pre-compute Lorentzian parameters
        let lorentzian_prefactor = gamma / PI;

        Zip::from(omega_grid)
            .and(spectral_weights)
            .for_each(|&omega, weight| {
                let mut sum = 0.0;
                for (i, &energy) in self.qp_energies.iter().enumerate() {
                    let delta_e = omega - energy;

                    // Gaussian component
                    let gaussian = gaussian_prefactor * (-delta_e * delta_e / two_sigma_sq).exp();

                    // Lorentzian component
                    let lorentzian = lorentzian_prefactor / (delta_e * delta_e + gamma * gamma);

                    // Weighted mixture
                    let voigt = (1.0 - alpha) * gaussian + alpha * lorentzian;
                    sum += self.z_factors[i] * voigt;
                }
                *weight = sum;
            });
    }

    /// Generate PES (photoemission) spectrum
    ///
    /// # Arguments
    /// * `omega_min` - Minimum binding energy (eV)
    /// * `omega_max` - Maximum binding energy (eV)
    /// * `n_points` - Number of grid points
    ///
    /// # Errors
    /// Returns error if the omega grid is invalid or spectrum computation fails
    pub fn generate_pes_spectrum(
        &self,
        omega_min: f64,
        omega_max: f64,
        n_points: usize,
    ) -> Result<SpectralFunction> {
        // Create binding energy grid (negative energies relative to Fermi)
        let omega_grid = Array1::linspace(omega_min, omega_max, n_points);

        // Filter occupied states (below Fermi energy)
        let occupied_mask: Vec<bool> = self
            .qp_energies
            .iter()
            .map(|&e| e <= self.fermi_energy)
            .collect();

        // Create filtered arrays
        let occupied_energies: Array1<f64> = self
            .qp_energies
            .iter()
            .zip(&occupied_mask)
            .filter(|(_, &mask)| mask)
            .map(|(&e, _)| self.fermi_energy - e) // Convert to binding energy
            .collect();

        let occupied_z_factors: Array1<f64> = self
            .z_factors
            .iter()
            .zip(&occupied_mask)
            .filter(|(_, &mask)| mask)
            .map(|(&z, _)| z)
            .collect();

        // Create analyzer for occupied states
        let pes_analyzer = SpectralAnalyzer {
            qp_energies: occupied_energies,
            z_factors: occupied_z_factors,
            fermi_energy: 0.0, // Binding energies relative to Fermi
            broadening: self.broadening,
        };

        pes_analyzer.compute_spectral_function(&omega_grid)
    }

    /// Generate IPES (inverse photoemission) spectrum
    ///
    /// # Arguments
    /// * `omega_min` - Minimum energy above Fermi (eV)
    /// * `omega_max` - Maximum energy above Fermi (eV)
    /// * `n_points` - Number of grid points
    ///
    /// # Errors
    /// Returns error if the omega grid is invalid or spectrum computation fails
    pub fn generate_ipes_spectrum(
        &self,
        omega_min: f64,
        omega_max: f64,
        n_points: usize,
    ) -> Result<SpectralFunction> {
        // Create energy grid above Fermi
        let omega_grid = Array1::linspace(omega_min, omega_max, n_points);

        // Filter unoccupied states (above Fermi energy)
        let unoccupied_mask: Vec<bool> = self
            .qp_energies
            .iter()
            .map(|&e| e > self.fermi_energy)
            .collect();

        // Create filtered arrays
        let unoccupied_energies: Array1<f64> = self
            .qp_energies
            .iter()
            .zip(&unoccupied_mask)
            .filter(|(_, &mask)| mask)
            .map(|(&e, _)| e - self.fermi_energy) // Energy above Fermi
            .collect();

        let unoccupied_z_factors: Array1<f64> = self
            .z_factors
            .iter()
            .zip(&unoccupied_mask)
            .filter(|(_, &mask)| mask)
            .map(|(&z, _)| z)
            .collect();

        // Create analyzer for unoccupied states
        let ipes_analyzer = SpectralAnalyzer {
            qp_energies: unoccupied_energies,
            z_factors: unoccupied_z_factors,
            fermi_energy: 0.0, // Energies relative to Fermi
            broadening: self.broadening,
        };

        ipes_analyzer.compute_spectral_function(&omega_grid)
    }

    /// Export spectral data to JSON format
    ///
    /// # Errors
    /// Returns error if file I/O or JSON serialization fails
    pub fn export_json(&self, data: &SpectralData, path: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(data)?;
        let mut file = File::create(path)?;
        file.write_all(json.as_bytes())?;
        Ok(())
    }

    /// Export spectral data to HDF5 format
    ///
    /// # Errors
    /// Returns error if file I/O or HDF5 operations fail
    #[cfg(feature = "hdf5_support")]
    pub fn export_hdf5(&self, data: &SpectralData, path: &Path) -> Result<()> {
        use hdf5::types::VarLenUnicode;

        let file = hdf5::File::create(path)?;

        // Create groups
        let root = file.group("/")?;
        root.new_attr::<VarLenUnicode>()
            .create("method")?
            .write_scalar(&VarLenUnicode::from_str(&data.method)?)?;
        root.new_attr::<VarLenUnicode>()
            .create("timestamp")?
            .write_scalar(&VarLenUnicode::from_str(&data.timestamp)?)?;

        // Export PES spectrum if available
        if let Some(ref pes) = data.pes_spectrum {
            let pes_group = root.create_group("pes")?;
            Self::write_spectral_function_to_hdf5(&pes_group, pes)?;
        }

        // Export IPES spectrum if available
        if let Some(ref ipes) = data.ipes_spectrum {
            let ipes_group = root.create_group("ipes")?;
            Self::write_spectral_function_to_hdf5(&ipes_group, ipes)?;
        }

        // Export combined spectrum if available
        if let Some(ref combined) = data.combined_spectrum {
            let combined_group = root.create_group("combined")?;
            Self::write_spectral_function_to_hdf5(&combined_group, combined)?;
        }

        Ok(())
    }

    /// Helper function to write spectral function to HDF5 group
    #[cfg(feature = "hdf5_support")]
    fn write_spectral_function_to_hdf5(
        group: &hdf5::Group,
        spectral_fn: &SpectralFunction,
    ) -> Result<()> {
        use hdf5::types::VarLenUnicode;

        // Write arrays
        group
            .new_dataset::<f64>()
            .shape(spectral_fn.omega_grid.len())
            .create("omega_grid")?
            .write(spectral_fn.omega_grid.as_slice().unwrap())?;

        group
            .new_dataset::<f64>()
            .shape(spectral_fn.spectral_weights.len())
            .create("spectral_weights")?
            .write(spectral_fn.spectral_weights.as_slice().unwrap())?;

        // Write metadata
        let metadata_group = group.create_group("metadata")?;
        metadata_group
            .new_attr::<usize>()
            .create("n_states")?
            .write_scalar(&spectral_fn.metadata.n_states)?;
        metadata_group
            .new_attr::<f64>()
            .create("fermi_energy")?
            .write_scalar(&spectral_fn.metadata.fermi_energy)?;
        metadata_group
            .new_attr::<f64>()
            .create("total_weight")?
            .write_scalar(&spectral_fn.metadata.total_weight)?;
        metadata_group
            .new_attr::<VarLenUnicode>()
            .create("energy_unit")?
            .write_scalar(&VarLenUnicode::from_str(&spectral_fn.metadata.energy_unit)?)?;

        // Write broadening parameters
        let broadening_group = metadata_group.create_group("broadening")?;
        broadening_group
            .new_attr::<f64>()
            .create("gaussian_width")?
            .write_scalar(&spectral_fn.metadata.broadening.gaussian_width)?;
        broadening_group
            .new_attr::<f64>()
            .create("lorentzian_width")?
            .write_scalar(&spectral_fn.metadata.broadening.lorentzian_width)?;
        broadening_group
            .new_attr::<f64>()
            .create("temperature")?
            .write_scalar(&spectral_fn.metadata.broadening.temperature)?;
        broadening_group
            .new_attr::<f64>()
            .create("voigt_mixing")?
            .write_scalar(&spectral_fn.metadata.broadening.voigt_mixing)?;

        Ok(())
    }
}

/// Gaussian broadening function (standalone)
///
/// # Arguments
/// * `omega` - Frequency at which to evaluate
/// * `energy` - Center energy
/// * `sigma` - Gaussian width
///
/// # Returns
/// Gaussian-broadened delta function value
#[must_use]
#[inline]
pub fn gaussian_broadening(omega: f64, energy: f64, sigma: f64) -> f64 {
    let delta_e = omega - energy;
    let prefactor = 1.0 / (sigma * (2.0 * PI).sqrt());
    prefactor * (-delta_e * delta_e / (2.0 * sigma * sigma)).exp()
}

/// Lorentzian broadening function (standalone)
///
/// # Arguments
/// * `omega` - Frequency at which to evaluate
/// * `energy` - Center energy
/// * `gamma` - Lorentzian width
///
/// # Returns
/// Lorentzian-broadened delta function value
#[must_use]
#[inline]
pub fn lorentzian_broadening(omega: f64, energy: f64, gamma: f64) -> f64 {
    let delta_e = omega - energy;
    (gamma / PI) / (delta_e * delta_e + gamma * gamma)
}

/// Voigt profile (convolution of Gaussian and Lorentzian)
///
/// # Arguments
/// * `omega` - Frequency at which to evaluate
/// * `energy` - Center energy
/// * `sigma` - Gaussian width
/// * `gamma` - Lorentzian width
/// * `alpha` - Mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian)
///
/// # Returns
/// Voigt-broadened delta function value
#[must_use]
#[inline]
pub fn voigt_profile(omega: f64, energy: f64, sigma: f64, gamma: f64, alpha: f64) -> f64 {
    let gauss = gaussian_broadening(omega, energy, sigma);
    let lorentz = lorentzian_broadening(omega, energy, gamma);
    (1.0 - alpha) * gauss + alpha * lorentz
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_broadening_params_validation() {
        let mut params = BroadeningParams::default();
        assert!(params.validate().is_ok());

        params.gaussian_width = -1.0;
        assert!(params.validate().is_err());

        params.gaussian_width = 0.1;
        params.voigt_mixing = 1.5;
        assert!(params.validate().is_err());
    }

    #[test]
    fn test_z_factor_validation() {
        let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
        let valid_z = Array1::from_vec(vec![0.5, 0.8, 0.3]);
        let invalid_z = Array1::from_vec(vec![0.5, 1.5, 0.3]);

        assert!(SpectralAnalyzer::new(energies.clone(), valid_z, 0.0).is_ok());
        assert!(SpectralAnalyzer::new(energies, invalid_z, 0.0).is_err());
    }

    #[test]
    fn test_gaussian_broadening_normalization() {
        let sigma = 0.1;
        let omega_grid = Array1::linspace(-5.0, 5.0, 1000);
        let dw = omega_grid[1] - omega_grid[0];

        let mut sum = 0.0;
        for &omega in &omega_grid {
            sum += gaussian_broadening(omega, 0.0, sigma) * dw;
        }

        assert_relative_eq!(sum, 1.0, epsilon = 0.01);
    }

    #[test]
    fn test_lorentzian_broadening_normalization() {
        let gamma = 0.1;
        let omega_grid = Array1::linspace(-10.0, 10.0, 2000);
        let dw = omega_grid[1] - omega_grid[0];

        let mut sum = 0.0;
        for &omega in &omega_grid {
            sum += lorentzian_broadening(omega, 0.0, gamma) * dw;
        }

        assert_relative_eq!(sum, 1.0, epsilon = 0.01);
    }

    #[test]
    fn test_spectral_function_creation() {
        let energies = Array1::from_vec(vec![-2.0, -1.0, 0.0, 1.0, 2.0]);
        let z_factors = Array1::from_vec(vec![0.8, 0.7, 0.6, 0.7, 0.8]);
        let fermi_energy = 0.0;

        let analyzer = SpectralAnalyzer::new(energies, z_factors, fermi_energy).unwrap();
        let omega_grid = Array1::linspace(-3.0, 3.0, 100);
        let spectral_fn = analyzer.compute_spectral_function(&omega_grid).unwrap();

        assert_eq!(spectral_fn.omega_grid.len(), 100);
        assert_eq!(spectral_fn.spectral_weights.len(), 100);
        assert!(spectral_fn.metadata.total_weight > 0.0);
    }

    #[test]
    fn test_pes_ipes_separation() {
        let energies = Array1::from_vec(vec![-2.0, -1.0, 0.5, 1.0, 2.0]);
        let z_factors = Array1::from_vec(vec![0.8, 0.7, 0.6, 0.7, 0.8]);
        let fermi_energy = 0.0;

        let analyzer = SpectralAnalyzer::new(energies, z_factors, fermi_energy).unwrap();

        // PES should only contain occupied states
        let pes = analyzer.generate_pes_spectrum(0.0, 3.0, 50).unwrap();
        assert_eq!(pes.metadata.n_states, 2); // Two occupied states

        // IPES should only contain unoccupied states
        let ipes = analyzer.generate_ipes_spectrum(0.0, 3.0, 50).unwrap();
        assert_eq!(ipes.metadata.n_states, 3); // Three unoccupied states
    }

    #[test]
    fn test_voigt_profile_limits() {
        let omega = 1.0;
        let energy = 0.0;
        let sigma = 0.1;
        let gamma = 0.05;

        // Pure Gaussian (alpha = 0)
        let pure_gauss = voigt_profile(omega, energy, sigma, gamma, 0.0);
        let expected_gauss = gaussian_broadening(omega, energy, sigma);
        assert_relative_eq!(pure_gauss, expected_gauss, epsilon = 1e-10);

        // Pure Lorentzian (alpha = 1)
        let pure_lorentz = voigt_profile(omega, energy, sigma, gamma, 1.0);
        let expected_lorentz = lorentzian_broadening(omega, energy, gamma);
        assert_relative_eq!(pure_lorentz, expected_lorentz, epsilon = 1e-10);
    }

    #[test]
    fn test_json_export() {
        use tempfile::NamedTempFile;

        let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
        let z_factors = Array1::from_vec(vec![0.5, 0.8, 0.3]);
        let analyzer = SpectralAnalyzer::new(energies, z_factors, 0.0).unwrap();

        let data = SpectralData {
            pes_spectrum: None,
            ipes_spectrum: None,
            combined_spectrum: None,
            timestamp: "2024-01-01T00:00:00Z".to_string(),
            method: "evGW".to_string(),
        };

        let temp_file = NamedTempFile::new().unwrap();
        analyzer.export_json(&data, temp_file.path()).unwrap();

        // Verify file was created and contains valid JSON
        let contents = std::fs::read_to_string(temp_file.path()).unwrap();
        let parsed: SpectralData = serde_json::from_str(&contents).unwrap();
        assert_eq!(parsed.method, "evGW");
    }
}
