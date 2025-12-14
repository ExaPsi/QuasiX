//! Unified parameter mapping for CD and AC methods
//!
//! This module provides parameter structures and mapping functions to ensure
//! consistent configuration across contour deformation and analytical continuation
//! methods for validation comparisons.

use serde::{Deserialize, Serialize};

/// Parameters for contour deformation method
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CDParameters {
    /// Broadening parameter (η) in eV
    pub eta: f64,
    /// Number of frequency points
    pub n_freq: usize,
    /// Maximum frequency for integration (eV)
    pub omega_max: f64,
    /// Integration tolerance
    pub tolerance: f64,
    /// Use adaptive grid refinement
    pub adaptive_grid: bool,
    /// Number of contour segments
    pub n_segments: usize,
    /// Complex plane rotation angle (radians)
    pub rotation_angle: f64,
}

impl Default for CDParameters {
    fn default() -> Self {
        Self {
            eta: 0.01,           // 10 meV broadening
            n_freq: 64,          // Number of frequency points
            omega_max: 100.0,    // 100 eV max frequency
            tolerance: 1e-6,     // Integration tolerance
            adaptive_grid: true, // Enable adaptive refinement
            n_segments: 4,       // Contour segments
            rotation_angle: 0.0, // No rotation by default
        }
    }
}

/// Parameters for analytical continuation method
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ACParameters {
    /// Number of imaginary frequency points
    pub n_imag: usize,
    /// Maximum imaginary frequency (eV)
    pub xi_max: f64,
    /// Number of poles in rational approximation
    pub n_poles: usize,
    /// Regularization parameter for fitting
    pub regularization: f64,
    /// Model type: "multipole" or "pade"
    pub model_type: String,
    /// Convergence threshold
    pub convergence_tol: f64,
    /// Maximum iterations for fitting
    pub max_iter: usize,
}

impl Default for ACParameters {
    fn default() -> Self {
        Self {
            n_imag: 32,                // Imaginary frequency points
            xi_max: 50.0,              // 50 eV max imaginary frequency
            n_poles: 20,               // Number of poles
            regularization: 1e-8,      // Regularization
            model_type: "pade".into(), // Use Padé by default
            convergence_tol: 1e-6,     // Convergence tolerance
            max_iter: 100,             // Max iterations
        }
    }
}

/// Unified parameters ensuring consistency between methods
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnifiedParameters {
    /// Common frequency grid size
    pub n_freq_common: usize,
    /// Common energy range (eV)
    pub energy_range: (f64, f64),
    /// Common broadening for comparison
    pub eta_common: f64,
    /// Temperature (eV)
    pub temperature: f64,
    /// Chemical potential (eV)
    pub chemical_potential: f64,
    /// CD-specific parameters
    pub cd_params: CDParameters,
    /// AC-specific parameters
    pub ac_params: ACParameters,
}

impl Default for UnifiedParameters {
    fn default() -> Self {
        Self {
            n_freq_common: 64,
            energy_range: (-30.0, 30.0),
            eta_common: 0.01,
            temperature: 0.025, // ~300K
            chemical_potential: 0.0,
            cd_params: CDParameters::default(),
            ac_params: ACParameters::default(),
        }
    }
}

impl UnifiedParameters {
    /// Create new unified parameters with validation
    pub fn new(n_freq: usize, energy_range: (f64, f64)) -> Self {
        let mut params = Self {
            n_freq_common: n_freq,
            energy_range,
            ..Default::default()
        };

        // Synchronize frequency grids
        params.cd_params.n_freq = n_freq;
        params.ac_params.n_imag = n_freq / 2; // AC needs fewer imaginary points

        // Set consistent energy scales
        let energy_width = (energy_range.1 - energy_range.0).abs();
        params.cd_params.omega_max = energy_width * 2.0;
        params.ac_params.xi_max = energy_width;

        params
    }

    /// Validate parameter consistency
    pub fn validate(&self) -> Result<(), String> {
        // Check energy range
        if self.energy_range.0 >= self.energy_range.1 {
            return Err("Invalid energy range".to_string());
        }

        // Check frequency grids
        if self.n_freq_common < 16 {
            return Err("Too few frequency points".to_string());
        }

        // Check broadening parameters
        if self.eta_common <= 0.0 || self.cd_params.eta <= 0.0 {
            return Err("Broadening must be positive".to_string());
        }

        // Check AC parameters
        if self.ac_params.n_poles < 2 {
            return Err("Too few poles for AC".to_string());
        }

        // Check temperature
        if self.temperature < 0.0 {
            return Err("Temperature must be non-negative".to_string());
        }

        Ok(())
    }

    /// Get frequency grid for CD method
    pub fn cd_frequency_grid(&self) -> Vec<f64> {
        let (e_min, e_max) = self.energy_range;
        let n = self.cd_params.n_freq;
        (0..n)
            .map(|i| e_min + (e_max - e_min) * (i as f64) / ((n - 1) as f64))
            .collect()
    }

    /// Get imaginary frequency grid for AC method
    pub fn ac_imaginary_grid(&self) -> Vec<f64> {
        let n = self.ac_params.n_imag;
        let xi_max = self.ac_params.xi_max;

        // Use Gauss-Legendre-like spacing for better convergence
        (0..n)
            .map(|i| {
                let t = (i as f64 + 0.5) / (n as f64);
                xi_max * t * t // Quadratic spacing concentrates points near zero
            })
            .collect()
    }

    /// Map CD parameters to frequency grid configuration
    pub fn cd_grid_config(&self) -> FrequencyGridConfig {
        FrequencyGridConfig {
            n_points: self.cd_params.n_freq,
            min_freq: self.energy_range.0,
            max_freq: self.energy_range.1,
            broadening: self.cd_params.eta,
            grid_type: GridType::Linear,
        }
    }

    /// Map AC parameters to frequency grid configuration
    pub fn ac_grid_config(&self) -> FrequencyGridConfig {
        FrequencyGridConfig {
            n_points: self.ac_params.n_imag,
            min_freq: 0.0,
            max_freq: self.ac_params.xi_max,
            broadening: self.eta_common,
            grid_type: GridType::GaussLegendre,
        }
    }
}

/// Frequency grid configuration
#[derive(Debug, Clone)]
pub struct FrequencyGridConfig {
    /// Number of grid points
    pub n_points: usize,
    /// Minimum frequency
    pub min_freq: f64,
    /// Maximum frequency
    pub max_freq: f64,
    /// Broadening parameter
    pub broadening: f64,
    /// Grid type
    pub grid_type: GridType,
}

/// Grid type enumeration
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GridType {
    /// Linear spacing
    Linear,
    /// Gauss-Legendre quadrature
    GaussLegendre,
    /// Logarithmic spacing
    Logarithmic,
    /// Chebyshev nodes
    Chebyshev,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unified_parameters_default() {
        let params = UnifiedParameters::default();
        assert_eq!(params.n_freq_common, 64);
        assert_eq!(params.energy_range, (-30.0, 30.0));
        assert!(params.validate().is_ok());
    }

    #[test]
    fn test_parameter_validation() {
        let mut params = UnifiedParameters::default();

        // Test invalid energy range
        params.energy_range = (10.0, -10.0);
        assert!(params.validate().is_err());

        // Test invalid broadening
        params.energy_range = (-30.0, 30.0);
        params.eta_common = -0.01;
        assert!(params.validate().is_err());

        // Test too few frequency points
        params.eta_common = 0.01;
        params.n_freq_common = 8;
        assert!(params.validate().is_err());
    }

    #[test]
    fn test_grid_generation() {
        let params = UnifiedParameters::new(32, (-20.0, 20.0));

        let cd_grid = params.cd_frequency_grid();
        assert_eq!(cd_grid.len(), 32);
        assert!((cd_grid[0] - (-20.0)).abs() < 1e-10);
        assert!((cd_grid[31] - 20.0).abs() < 1e-10);

        let ac_grid = params.ac_imaginary_grid();
        assert_eq!(ac_grid.len(), 16); // n_freq / 2
        assert!(ac_grid[0] > 0.0);
        assert!(ac_grid[15] <= params.ac_params.xi_max);
    }

    #[test]
    fn test_parameter_synchronization() {
        let params = UnifiedParameters::new(128, (-50.0, 50.0));

        assert_eq!(params.cd_params.n_freq, 128);
        assert_eq!(params.ac_params.n_imag, 64);
        assert_eq!(params.cd_params.omega_max, 200.0);
        assert_eq!(params.ac_params.xi_max, 100.0);
    }
}
