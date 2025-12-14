//! Common utilities, types, and error handling for QuasiX
//!
//! This module provides shared functionality used across all QuasiX modules,
//! including error types, common constants, validation guards, and utility functions.

pub mod validation;

pub use validation::{
    validate_array1, validate_array1_complex, validate_array2, validate_array2_complex,
    validate_array3, validate_df_metric, validate_gw_dimensions, validate_mo_energies,
    validate_mo_occupations, ValidationConfig,
};

use thiserror::Error;

/// Main error type for QuasiX operations
#[derive(Error, Debug, Clone)]
pub enum QuasixError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("Dimension mismatch: {0}")]
    DimensionMismatch(String),

    #[error("Convergence failed: {0}")]
    ConvergenceFailed(String),

    #[error("Numerical error: {0}")]
    NumericalError(String),

    #[error("IO error: {0}")]
    IoError(String),

    #[error("Serialization error: {0}")]
    SerializationError(String),

    #[error("Not implemented: {0}")]
    NotImplemented(String),

    #[error("Convergence error: {0}")]
    ConvergenceError(String),

    #[error("Physics error: {0}")]
    PhysicsError(String),

    #[error("Python interop error: {0}")]
    PythonError(String),

    #[error("Parallelization error: {0}")]
    ParallelizationError(String),

    #[error("Rank deficiency detected: {message} (rank={rank}/{full_rank})")]
    RankDeficiency {
        message: String,
        rank: usize,
        full_rank: usize,
    },

    #[error("Runtime error: {0}")]
    RuntimeError(String),

    #[error("Validation error: {0}")]
    ValidationError(String),

    #[error("IO error occurred: {0}")]
    IOError(String),

    #[error("Deprecated method: {0}")]
    DeprecatedMethod(String),
}

impl From<anyhow::Error> for QuasixError {
    fn from(err: anyhow::Error) -> Self {
        QuasixError::NumericalError(err.to_string())
    }
}

impl From<std::io::Error> for QuasixError {
    fn from(err: std::io::Error) -> Self {
        QuasixError::IoError(err.to_string())
    }
}

impl From<serde_json::Error> for QuasixError {
    fn from(err: serde_json::Error) -> Self {
        QuasixError::SerializationError(err.to_string())
    }
}

impl From<ndarray_linalg::error::LinalgError> for QuasixError {
    fn from(err: ndarray_linalg::error::LinalgError) -> Self {
        QuasixError::NumericalError(format!("Linear algebra error: {}", err))
    }
}

/// `Result` type alias for QuasiX operations
pub type Result<T> = std::result::Result<T, QuasixError>;

/// Physical constants used in calculations
///
/// Note: Kept as crate-internal to avoid CFFI issues with float #define macros
/// These constants are used in various modules for unit conversions.
#[allow(dead_code)] // Constants may not be used in all configurations
pub(crate) mod constants {
    /// Hartree to eV conversion factor
    pub const HARTREE_TO_EV: f64 = 27.211_386_245_988;

    /// Bohr to Angstrom conversion factor
    pub const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_67;

    /// Default convergence threshold
    pub const DEFAULT_CONVERGENCE: f64 = 1e-6;
}

/// Get the version of the QuasiX core library
pub fn version() -> String {
    use tracing::info;
    let version = env!("CARGO_PKG_VERSION").to_string();
    info!(version = %version, "QuasiX version check");
    version
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        let v = version();
        assert!(!v.is_empty());
        assert_eq!(v, "0.6.0");
    }

    #[test]
    fn test_error_creation() {
        let err = QuasixError::InvalidInput("test".to_string());
        assert_eq!(format!("{err}"), "Invalid input: test");
    }

    #[test]
    fn test_constants() {
        assert!(constants::HARTREE_TO_EV > 27.0);
        assert!(constants::BOHR_TO_ANGSTROM > 0.5);
    }
}
