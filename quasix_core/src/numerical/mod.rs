// Numerical methods module for QuasiX
// High-precision algorithms for GW/BSE calculations

pub mod richardson;

pub use richardson::{RichardsonDerivative, ZFactorDerivative};

/// Numerical parameters for GW calculations
pub struct NumericalParams {
    /// Z-factor derivative computation
    pub z_factor_delta: f64,
    pub z_factor_levels: usize,
    pub z_factor_tolerance: f64,

    /// Physical bounds
    pub z_factor_min: f64,
    pub z_factor_max: f64,

    /// Convergence criteria
    pub derivative_max: f64,
}

impl Default for NumericalParams {
    fn default() -> Self {
        Self {
            // Richardson parameters (optimal for GW)
            z_factor_delta: 5e-4,
            z_factor_levels: 3,
            z_factor_tolerance: 1e-12,

            // Physical bounds
            z_factor_min: 0.01,  // Strong decay
            z_factor_max: 0.99,  // Weak correlation

            // Sanity checks
            derivative_max: 10.0,  // Detect instabilities
        }
    }
}

impl NumericalParams {
    /// Conservative parameters for difficult cases
    pub fn conservative() -> Self {
        Self {
            z_factor_delta: 1e-3,
            z_factor_levels: 2,
            z_factor_tolerance: 1e-10,
            ..Default::default()
        }
    }

    /// High-precision parameters for benchmarks
    pub fn high_precision() -> Self {
        Self {
            z_factor_delta: 1e-4,
            z_factor_levels: 4,
            z_factor_tolerance: 1e-14,
            ..Default::default()
        }
    }
}