//! Density fitting and resolution-of-identity module
//!
//! This module handles the construction and manipulation of density fitting
//! tensors, including RI basis management and symmetrization operations.
#![allow(clippy::many_single_char_names)]
// Mathematical notation
// Note: into_shape is deprecated but we use it for stable array reshaping
// TODO: Migrate to into_shape_with_order when API stabilizes
#![allow(deprecated)]

pub mod linear_dependency;
pub mod metric;
pub mod mo_transform;
pub mod parallel;
pub mod simd_ops;

#[cfg(test)]
mod simd_validation_test;

use crate::common::{QuasixError, Result};
use crate::integrals::{BasisSet, IntegralEngine, Molecule};
use ndarray::{Array2, Array3};
use ndarray_linalg::Solve;
use tracing::{debug, info};

// Re-export key types
pub use linear_dependency::{
    detect_linear_dependencies, regularize_metric, remove_dependencies_eigenvalue,
    transform_df_tensors, transform_mo_df_tensors, DependencyAnalysis, DependencyMethod,
    LinearDependencyConfig, TransformationMatrix,
};
pub use metric::{
    compute_cholesky_v, compute_cholesky_v_with_lindep, compute_cholesky_v_with_tolerance,
    incremental_cholesky_update, CholeskyMethod, CholeskyMetric,
};
pub use mo_transform::{generate_mock_mo_coefficients, transform_mo_3center};
// build_cholesky_metric is deprecated - use compute_cholesky_v instead

/// Density fitting tensor builder
///
/// Implements Cholesky-decomposed density fitting following PySCF conventions:
/// 1. Compute 3-center integrals (μν|P) via libcint
/// 2. Compute 2-center metric (P|Q) and Cholesky factor L
/// 3. Solve L @ cderi = (μν|P) to get Cholesky DF tensor
///
/// Reference: pyscf/df/incore.py::cholesky_eri
pub struct DFBuilder {
    /// Molecular structure
    molecule: Molecule,
    /// AO basis set
    basis: BasisSet,
    /// Auxiliary (RI) basis set
    aux_basis: BasisSet,
    /// Linear dependency tolerance
    lindep_tol: f64,
}

impl DFBuilder {
    /// Create a new DF builder
    ///
    /// # Arguments
    /// * `molecule` - Molecular structure
    /// * `basis` - AO basis set
    /// * `aux_basis` - Auxiliary (RI) basis set
    ///
    /// # Returns
    /// A configured DFBuilder ready to compute DF tensors
    pub fn new(molecule: Molecule, basis: BasisSet, aux_basis: BasisSet) -> Self {
        Self {
            molecule,
            basis,
            aux_basis,
            lindep_tol: 1e-7, // Match PySCF default
        }
    }

    /// Set linear dependency tolerance
    ///
    /// # Arguments
    /// * `tol` - Linear dependency threshold (default: 1e-7, matching PySCF)
    pub fn with_lindep_tolerance(mut self, tol: f64) -> Self {
        self.lindep_tol = tol;
        self
    }

    /// Build the 3-index DF tensor in AO basis
    ///
    /// Computes Cholesky-decomposed DF tensor following PySCF algorithm:
    /// 1. j3c = (μν|P) via compute_3center_integrals
    /// 2. j2c = (P|Q) via compute_2center_integrals
    /// 3. L = cholesky(j2c)
    /// 4. cderi = L^(-1) @ j3c (solve triangular system)
    ///
    /// # Returns
    /// Array3 of shape (nbasis, nbasis, naux) containing Cholesky DF tensor
    ///
    /// # Errors
    /// Returns error if:
    /// - Integral computation fails
    /// - Cholesky decomposition fails (singular metric)
    /// - Triangular solve fails
    ///
    /// # Reference
    /// PySCF: `pyscf/df/incore.py::cholesky_eri()`
    pub fn build_df_tensor(&self) -> Result<Array3<f64>> {
        let nbasis = self.basis.size();
        let naux = self.aux_basis.size();

        info!(
            "Building DF tensor for {} with {} AO and {} aux functions",
            self.molecule.name(),
            nbasis,
            naux
        );

        // Step 1: Compute 3-center integrals (μν|P)
        debug!("Computing 3-center integrals (μν|P)");
        let engine = IntegralEngine::new(
            self.molecule.clone(),
            self.basis.clone(),
            self.aux_basis.clone(),
        );
        let j3c = engine.compute_3center()?; // Shape: (nbasis, nbasis, naux)

        // Step 2: Compute 2-center metric (P|Q)
        debug!("Computing 2-center metric (P|Q)");
        let j2c = engine.compute_2center()?; // Shape: (naux, naux)

        // Step 3: Cholesky factorization: j2c = L @ L^T
        debug!("Computing Cholesky factorization of metric");
        let lindep_config = Some(LinearDependencyConfig {
            eigenvalue_threshold: self.lindep_tol,
            ..Default::default()
        });
        let cholesky_result = compute_cholesky_v_with_lindep(&j2c, lindep_config)?;

        if cholesky_result.method != CholeskyMethod::Standard {
            info!(
                "Using {:?} for metric factorization (condition number: {:.2e})",
                cholesky_result.method, cholesky_result.condition_number
            );
        }

        let l_matrix = &cholesky_result.l_matrix; // Lower triangular

        // Step 4: Solve L @ cderi = j3c for cderi
        // Reshape j3c: (nbasis, nbasis, naux) → (nbasis*nbasis, naux)
        debug!("Reshaping 3-center integrals for triangular solve");
        let j3c_2d = j3c
            .into_shape((nbasis * nbasis, naux))
            .map_err(|e| QuasixError::InvalidInput(format!("Reshape failed: {}", e)))?;

        // Solve: L @ cderi.T = j3c_2d.T  =>  cderi.T = L^(-1) @ j3c_2d.T
        // We need to solve L @ X = B where B = j3c_2d.T
        // This means X = L^(-1) @ B
        debug!("Solving triangular system L @ X = j3c_2d.T");

        // Transpose j3c_2d to get (naux, nbasis*nbasis)
        let j3c_t = j3c_2d.t().to_owned();

        // Solve for each column separately
        // Result will be (naux, nbasis*nbasis)
        let mut cderi_t = Array2::zeros((naux, nbasis * nbasis));
        for (i, col) in j3c_t.axis_iter(ndarray::Axis(1)).enumerate() {
            let x = l_matrix.solve(&col.to_owned())?;
            cderi_t.column_mut(i).assign(&x);
        }

        // Transpose back to get (nbasis*nbasis, naux)
        let cderi = cderi_t.t().to_owned();

        // Step 5: Reshape back: (nbasis*nbasis, naux) → (nbasis, nbasis, naux)
        debug!("Reshaping result back to 3D tensor");
        let df_tensor = cderi
            .into_shape((nbasis, nbasis, naux))
            .map_err(|e| QuasixError::InvalidInput(format!("Final reshape failed: {}", e)))?;

        info!(
            "DF tensor computed: shape ({}, {}, {}), norm: {:.6e}",
            nbasis,
            naux,
            naux,
            df_tensor.iter().map(|x| x * x).sum::<f64>().sqrt()
        );

        Ok(df_tensor)
    }

    /// Build the DF metric matrix (P|Q)
    ///
    /// # Returns
    /// Array2 of shape (naux, naux) containing 2-center metric
    ///
    /// # Errors
    /// Returns error if integral computation fails
    pub fn build_metric(&self) -> Result<Array2<f64>> {
        let engine = IntegralEngine::new(
            self.molecule.clone(),
            self.basis.clone(),
            self.aux_basis.clone(),
        );
        engine.compute_2center()
    }

    /// Compute Cholesky factor of metric: V = L @ L^T
    ///
    /// # Arguments
    /// * `metric` - 2-center metric matrix (P|Q)
    ///
    /// # Returns
    /// CholeskyMetric containing L matrix and diagnostics
    ///
    /// # Errors
    /// Returns error if Cholesky decomposition fails
    pub fn compute_metric_cholesky(&self, metric: &Array2<f64>) -> Result<CholeskyMetric> {
        let lindep_config = Some(LinearDependencyConfig {
            eigenvalue_threshold: self.lindep_tol,
            ..Default::default()
        });
        compute_cholesky_v_with_lindep(metric, lindep_config)
    }
}

/// RI basis set information
pub struct RIBasis {
    /// Number of auxiliary functions
    pub size: usize,
    /// Basis type identifier
    pub basis_type: String,
}

impl RIBasis {
    /// Create a new RI basis
    #[must_use]
    pub fn new(size: usize, basis_type: String) -> Self {
        Self { size, basis_type }
    }

    /// Validate the RI basis
    pub fn validate(&self) -> Result<()> {
        todo!("Implement RI basis validation")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrals::{Atom, BasisSet, Molecule};

    #[test]
    fn test_ri_basis_creation() {
        let basis = RIBasis::new(100, "cc-pVDZ-RI".to_string());
        assert_eq!(basis.size, 100);
        assert_eq!(basis.basis_type, "cc-pVDZ-RI");
    }

    #[test]
    fn test_df_builder_creation() {
        // Create minimal H2 molecule for testing
        let atoms = vec![
            Atom::new("H", 1, [0.0, 0.0, 0.0]),
            Atom::new("H", 1, [0.0, 0.0, 1.4]),
        ];
        let mol = Molecule::new(atoms, 0, 1);

        // Create minimal basis sets using mock functions
        let basis = BasisSet::mock_sto3g(2);
        let aux_basis = BasisSet::mock_ri(10);

        let builder = DFBuilder::new(mol, basis, aux_basis);
        assert_eq!(builder.lindep_tol, 1e-7);
    }
}
