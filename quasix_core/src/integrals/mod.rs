//! Integral computation module for density fitting
//!
//! This module provides functionality for computing 2-center and 3-center
//! integrals required for density fitting/resolution-of-identity calculations.
//! Uses real libcint library for production-quality integral evaluation.
#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::{QuasixError, Result};
use ndarray::{Array2, Array3};
use tracing::{debug, info, warn};

pub mod basis;
pub mod molecule;
pub mod three_center;
pub mod three_center_optimized;
pub mod two_center;

// Real libcint implementation using the actual crate API
pub mod libcint_real;

// Re-export main types
pub use basis::{BasisFunction, BasisSet, BasisType};
pub use molecule::{Atom, Molecule};
pub use three_center::compute_3center_integrals;
pub use three_center_optimized::{
    compute_3center_integrals_optimized, OptimizedThreeCenterIntegrals, PerformanceMetrics,
};
pub use two_center::{
    compute_2center_integrals, compute_2center_integrals_with_config, TwoCenterConfig,
    TwoCenterMetrics,
};

/// Error types specific to integral calculations
#[derive(Debug, thiserror::Error)]
pub enum IntegralError {
    #[error("Invalid basis set: {0}")]
    InvalidBasis(String),

    #[error("Invalid molecule: {0}")]
    InvalidMolecule(String),

    #[error("Integral computation failed: {0}")]
    ComputationFailed(String),

    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: String, actual: String },
}

impl From<IntegralError> for QuasixError {
    fn from(err: IntegralError) -> Self {
        QuasixError::NumericalError(err.to_string())
    }
}

/// Main entry point for computing DF integrals
pub struct IntegralEngine {
    #[allow(dead_code)]
    molecule: Molecule,
    basis: BasisSet,
    aux_basis: BasisSet,
    #[cfg(feature = "real_libcint")]
    engine: libcint_real::LibcintRealEngine,
}

impl IntegralEngine {
    /// Create a new integral engine
    #[must_use]
    pub fn new(molecule: Molecule, basis: BasisSet, aux_basis: BasisSet) -> Self {
        info!(
            "Creating IntegralEngine for {} with {} basis functions and {} auxiliary functions",
            molecule.name(),
            basis.size(),
            aux_basis.size()
        );

        #[cfg(feature = "real_libcint")]
        let engine = libcint_real::LibcintRealEngine::new(&molecule, &basis, &aux_basis)
            .expect("Failed to initialize libcint engine");

        Self {
            molecule,
            basis,
            aux_basis,
            #[cfg(feature = "real_libcint")]
            engine,
        }
    }

    /// Compute 3-center integrals (μν|P)
    pub fn compute_3center(&self) -> Result<Array3<f64>> {
        info!("Computing 3-center integrals (μν|P)");
        debug!(
            "Dimensions: ({}, {}, {})",
            self.basis.size(),
            self.basis.size(),
            self.aux_basis.size()
        );

        #[cfg(feature = "real_libcint")]
        let result = self.engine.compute_3center()?;

        #[cfg(not(feature = "real_libcint"))]
        let result = compute_3center_integrals(&self.molecule, &self.basis, &self.aux_basis)?;

        // Validate dimensions - allow for differences due to spherical vs cartesian
        // When using real libcint, the actual size may differ from mock basis size
        // because libcint uses spherical harmonics (d=5, f=7) vs cartesian (d=6, f=10)
        let actual_shape = result.dim();

        #[cfg(not(feature = "real_libcint"))]
        {
            let expected_nbasis = self.basis.size();
            // Check that first two dimensions match expected basis size
            if actual_shape.0 != expected_nbasis || actual_shape.1 != expected_nbasis {
                return Err(IntegralError::DimensionMismatch {
                    expected: format!("({}, {}, _)", expected_nbasis, expected_nbasis),
                    actual: format!("{:?}", actual_shape),
                }
                .into());
            }
        }

        // For real libcint, just ensure dimensions are consistent (square in first two indices)
        #[cfg(feature = "real_libcint")]
        {
            if actual_shape.0 != actual_shape.1 {
                return Err(IntegralError::DimensionMismatch {
                    expected: "Square matrix in first two indices".to_string(),
                    actual: format!("{:?}", actual_shape),
                }
                .into());
            }
        }

        // For aux basis, allow for spherical vs cartesian differences
        // Spherical harmonics have fewer components: d=5 vs 6, f=7 vs 10
        let expected_naux = self.aux_basis.size();
        let actual_naux = actual_shape.2;
        if actual_naux.abs_diff(expected_naux) > 6 {
            warn!(
                "Large difference in auxiliary basis size: expected {}, got {}",
                expected_naux, actual_naux
            );
        }

        info!("3-center integrals computed successfully");
        Ok(result)
    }

    /// Compute 2-center integrals (P|Q)
    pub fn compute_2center(&self) -> Result<Array2<f64>> {
        info!("Computing 2-center integrals (P|Q)");
        debug!(
            "Dimensions: ({}, {})",
            self.aux_basis.size(),
            self.aux_basis.size()
        );

        #[cfg(feature = "real_libcint")]
        let result = self.engine.compute_2center()?;

        #[cfg(not(feature = "real_libcint"))]
        let result = compute_2center_integrals(&self.aux_basis)?;

        // Validate dimensions - allow for differences due to spherical vs cartesian basis
        let actual_shape = result.dim();

        // Check that matrix is square
        if actual_shape.0 != actual_shape.1 {
            return Err(IntegralError::DimensionMismatch {
                expected: "Square matrix".to_string(),
                actual: format!("{:?}", actual_shape),
            }
            .into());
        }

        // For aux basis, allow for spherical vs cartesian differences
        let expected_naux = self.aux_basis.size();
        let actual_naux = actual_shape.0;
        if actual_naux.abs_diff(expected_naux) > 6 {
            warn!(
                "Large difference in auxiliary basis size: expected {}, got {}",
                expected_naux, actual_naux
            );
        }

        // Validate positive definiteness
        if !is_positive_definite(&result) {
            warn!("2-center integral matrix is not positive definite");
        }

        info!("2-center integrals computed successfully");
        Ok(result)
    }

    /// Get the number of basis functions
    pub fn nbasis(&self) -> usize {
        self.basis.size()
    }

    /// Get the number of auxiliary basis functions
    pub fn naux(&self) -> usize {
        self.aux_basis.size()
    }
}

/// Check if a matrix is positive definite (simplified check)
fn is_positive_definite(matrix: &Array2<f64>) -> bool {
    // Simple check: diagonal elements should be positive
    // and matrix should be symmetric
    let n = matrix.nrows();
    if n != matrix.ncols() {
        return false;
    }

    // Check diagonal elements
    for i in 0..n {
        if matrix[[i, i]] <= 0.0 {
            return false;
        }
    }

    // Check symmetry (with tolerance)
    let tol = 1e-10;
    for i in 0..n {
        for j in i + 1..n {
            if (matrix[[i, j]] - matrix[[j, i]]).abs() > tol {
                return false;
            }
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integral_engine_creation() {
        let mol = Molecule::water();
        let basis = BasisSet::mock_sto3g(8);
        let aux_basis = BasisSet::mock_ri(18);

        let engine = IntegralEngine::new(mol, basis, aux_basis);
        assert_eq!(engine.nbasis(), 8);
        assert_eq!(engine.naux(), 18);
    }

    #[test]
    fn test_positive_definite_check() {
        let mut matrix = Array2::eye(3);
        assert!(is_positive_definite(&matrix));

        matrix[[0, 0]] = -1.0;
        assert!(!is_positive_definite(&matrix));
    }
}
