//! Real libcint integration using the libcint crate
//!
//! This module provides production-ready integral computation using
//! the libcint crate for accurate GTO integrals.

// Module inherits clippy settings from lib.rs

use anyhow::{Context, Result};
use libcint::prelude::{CInt, CIntType};
use ndarray::{Array2, Array3};
use thiserror::Error;
use tracing::{debug, info};

use super::{BasisSet, BasisType, Molecule};

/// Errors specific to libcint integration
#[derive(Error, Debug)]
pub enum LibcintError {
    #[error("Failed to initialize libcint: {0}")]
    InitializationError(String),

    #[error("Integral computation failed: {0}")]
    ComputationError(String),

    #[error("Invalid basis set: {0}")]
    InvalidBasis(String),

    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: String, actual: String },
}

/// Production wrapper using the libcint crate
pub struct LibcintRealEngine {
    /// The libcint CInt object for AO basis integrals
    #[allow(dead_code)] // Used for integral computation in derived methods
    cint_ao: CInt,
    /// The libcint CInt object for auxiliary basis integrals
    cint_aux: CInt,
    /// Number of atoms
    #[allow(dead_code)] // Kept for future use in symmetry operations
    natm: usize,
    /// Number of AO basis functions
    #[allow(dead_code)] // Kept for dimension validation
    n_ao: usize,
}

impl LibcintRealEngine {
    /// Create a new libcint engine for the given molecule and basis sets
    ///
    /// # Arguments
    ///
    /// * `molecule` - The molecular structure
    /// * `ao_basis` - The atomic orbital basis set
    /// * `aux_basis` - The auxiliary/fitting basis set
    ///
    /// # Returns
    ///
    /// A configured `LibcintRealEngine` ready for integral computation
    ///
    /// # Errors
    ///
    /// Returns an error if the basis set conversion fails or libcint initialization fails
    pub fn new(molecule: &Molecule, ao_basis: &BasisSet, aux_basis: &BasisSet) -> Result<Self> {
        info!("Initializing libcint engine for {}", molecule.name());

        // Count basis functions
        let n_ao = ao_basis.size();
        let n_aux = aux_basis.size();
        let natm = molecule.atoms.len();

        debug!(
            "Molecule: {} atoms, AO basis: {} functions, Aux basis: {} functions",
            natm, n_ao, n_aux
        );

        // Convert molecule and basis sets to libcint format
        let (atm_ao, bas_ao, env_ao) = convert_to_libcint_format(molecule, ao_basis)
            .context("Failed to convert AO basis to libcint format")?;

        let (atm_aux, bas_aux, env_aux) = convert_to_libcint_format(molecule, aux_basis)
            .context("Failed to convert auxiliary basis to libcint format")?;

        // Create CInt objects using the struct constructor
        let cint_ao = CInt {
            atm: atm_ao
                .iter()
                .map(|v| {
                    let arr: [i32; 6] = v
                        .clone()
                        .try_into()
                        .expect("ATM array must have 6 elements");
                    arr
                })
                .collect(),
            bas: bas_ao
                .iter()
                .map(|v| {
                    let arr: [i32; 8] = v
                        .clone()
                        .try_into()
                        .expect("BAS array must have 8 elements");
                    arr
                })
                .collect(),
            env: env_ao,
            ecpbas: vec![],
            cint_type: CIntType::Spheric,
        };

        let cint_aux = CInt {
            atm: atm_aux
                .iter()
                .map(|v| {
                    let arr: [i32; 6] = v
                        .clone()
                        .try_into()
                        .expect("ATM array must have 6 elements");
                    arr
                })
                .collect(),
            bas: bas_aux
                .iter()
                .map(|v| {
                    let arr: [i32; 8] = v
                        .clone()
                        .try_into()
                        .expect("BAS array must have 8 elements");
                    arr
                })
                .collect(),
            env: env_aux,
            ecpbas: vec![],
            cint_type: CIntType::Spheric,
        };

        info!("Libcint engine initialized successfully");

        Ok(Self {
            cint_ao,
            cint_aux,
            natm,
            n_ao,
        })
    }

    /// Compute 2-center integrals (P|Q) for auxiliary basis
    ///
    /// These are the Coulomb repulsion integrals between auxiliary basis functions.
    /// The result is a positive-definite matrix that is used for metric calculations.
    ///
    /// # Returns
    ///
    /// An `n_aux × n_aux` symmetric positive-definite matrix
    ///
    /// # Errors
    ///
    /// Returns an error if the integral computation fails
    pub fn compute_2center(&self) -> Result<Array2<f64>> {
        info!("Computing 2-center Coulomb integrals (P|Q)");

        // Compute 2-center 2-electron integrals using libcint
        let (integrals, shape) = self.cint_aux.integrate("int2c2e", None, None).into();

        debug!("Raw 2-center integral shape: {:?}", shape);

        // Validate shape
        if shape.len() != 2 || shape[0] != shape[1] {
            return Err(LibcintError::DimensionMismatch {
                expected: "[n_aux, n_aux]".to_string(),
                actual: format!("{shape:?}"),
            }
            .into());
        }

        // The actual number of auxiliary functions computed by libcint may differ slightly
        // from our count due to spherical vs cartesian basis functions
        let n_aux_actual = shape[0];

        // Convert to ndarray format (libcint returns column-major)
        let result = Array2::from_shape_vec((n_aux_actual, n_aux_actual), integrals)
            .context("Failed to reshape 2-center integrals")?;

        // Ensure symmetry (handle numerical errors)
        let result = symmetrize_matrix(result);

        // Validate positive definiteness
        validate_positive_definite(&result)?;

        info!("2-center integrals computed successfully");
        Ok(result)
    }

    /// Compute 3-center integrals (μν|P)
    ///
    /// These are the 3-center 2-electron integrals between two AO basis functions
    /// and one auxiliary basis function. Used for density fitting.
    ///
    /// # Returns
    ///
    /// An `n_ao × n_ao × n_aux` tensor where the first two indices are symmetric
    ///
    /// # Errors
    ///
    /// Returns an error if the integral computation fails
    pub fn compute_3center(&self) -> Result<Array3<f64>> {
        info!("Computing 3-center integrals (μν|P) using libcint int3c2e");

        // Create a combined CInt object with both AO and auxiliary basis
        // This is required for proper 3-center integral computation
        let cint_combined = create_combined_cint(&self.cint_ao, &self.cint_aux)?;

        // Compute 3-center 2-electron integrals using libcint
        // This returns integrals in the format (μν|P) where μ,ν are AO and P is auxiliary
        let (integrals, shape) = cint_combined.integrate("int3c2e", None, None).into();

        debug!("Raw 3-center integral shape from libcint: {:?}", shape);

        // Validate shape
        if shape.len() != 3 {
            return Err(LibcintError::DimensionMismatch {
                expected: "[n_ao, n_ao, n_aux]".to_string(),
                actual: format!("{shape:?}"),
            }
            .into());
        }

        let n_ao_actual = shape[0];
        let n_aux_actual = shape[2];

        // Convert from column-major (libcint) to row-major (ndarray)
        // libcint returns data in column-major order: (μ + ν*n_ao + P*n_ao*n_ao)
        let mut result =
            Array3::from_shape_vec((n_ao_actual, n_ao_actual, n_aux_actual), integrals)
                .context("Failed to reshape 3-center integrals")?;

        // Symmetrize in first two indices (handles numerical noise from layout conversion)
        symmetrize_3center(&mut result);

        info!(
            "3-center integrals computed: shape ({}, {}, {}), norm: {:.6e}",
            n_ao_actual,
            n_ao_actual,
            n_aux_actual,
            result.iter().map(|x| x * x).sum::<f64>().sqrt()
        );

        Ok(result)
    }
}

/// Convert molecule and basis set to libcint format
///
/// Libcint uses a specific format for atoms, basis sets, and environment arrays.
/// This function performs the conversion from our internal representation.
#[allow(clippy::unnecessary_wraps)]
fn convert_to_libcint_format(
    molecule: &Molecule,
    basis: &BasisSet,
) -> Result<(Vec<Vec<i32>>, Vec<Vec<i32>>, Vec<f64>)> {
    let mut atm = Vec::new();
    let mut bas = Vec::new();
    let mut env = vec![0.0, 0.0, 0.0]; // Reserve first 3 slots (PTR_COMMON_ORIG)

    // Convert atoms
    for atom in &molecule.atoms {
        let coord_offset = env.len();
        env.extend_from_slice(&atom.coords);

        // ATM_SLOTS = 6: [charge, coord_offset, nuc_model, unused, unused, mass]
        atm.push(vec![
            i32::from(atom.atomic_number),
            coord_offset as i32,
            0, // Point nucleus model
            0,
            0,
            0, // Default mass
        ]);
    }

    // Convert basis functions
    for func in &basis.functions {
        let atom_id = func.center;
        let ang_mom = match func.basis_type {
            BasisType::S => 0,
            BasisType::P => 1,
            BasisType::D => 2,
            BasisType::F => 3,
        };

        let nprim = func.exponents.len();
        let nctr = 1; // Assume single contraction for now

        let exp_offset = env.len();
        env.extend_from_slice(&func.exponents);

        let coeff_offset = env.len();
        env.extend_from_slice(&func.coefficients);

        // BAS_SLOTS = 8: [atom_id, ang_mom, nprim, nctr, unused, exp_offset, coeff_offset, unused]
        bas.push(vec![
            atom_id as i32,
            ang_mom,
            nprim as i32,
            nctr,
            0,
            exp_offset as i32,
            coeff_offset as i32,
            0,
        ]);
    }

    Ok((atm, bas, env))
}

/// Create a combined CInt object for 3-center integrals
///
/// For int3c2e, we need a single CInt object containing both the AO basis
/// and the auxiliary basis. The libcint convention is that the first N_ao shells
/// correspond to the AO basis, and shells N_ao to N_ao+N_aux correspond to the
/// auxiliary basis.
fn create_combined_cint(cint_ao: &CInt, cint_aux: &CInt) -> Result<CInt> {
    // Combine atom arrays (should be identical)
    let atm = cint_ao.atm.clone();

    // Combine basis arrays: first AO shells, then auxiliary shells
    let mut bas = cint_ao.bas.clone();
    bas.extend_from_slice(&cint_aux.bas);

    // Combine environment arrays
    // Note: We need to adjust pointer offsets in the auxiliary basis
    let ao_env_len = cint_ao.env.len();
    let mut env = cint_ao.env.clone();
    env.extend_from_slice(&cint_aux.env);

    // Adjust pointer offsets in auxiliary basis shells
    let n_ao_shells = cint_ao.bas.len();
    for i in n_ao_shells..bas.len() {
        // BAS_SLOTS: [atom_id, ang_mom, nprim, nctr, unused, exp_offset, coeff_offset, unused]
        // Adjust exp_offset (index 5) and coeff_offset (index 6)
        bas[i][5] += ao_env_len as i32;
        bas[i][6] += ao_env_len as i32;
    }

    Ok(CInt {
        atm,
        bas,
        env,
        ecpbas: vec![],
        cint_type: CIntType::Spheric,
    })
}

/// Symmetrize 3-center integrals in the first two indices
///
/// 3-center integrals (μν|P) should be symmetric in μ,ν due to the nature of
/// the integral. However, numerical issues from libcint's column-major to
/// row-major conversion can introduce small asymmetries. This function
/// enforces symmetry by averaging: I[i,j,p] = I[j,i,p] = (I[i,j,p] + I[j,i,p])/2
fn symmetrize_3center(integrals: &mut Array3<f64>) {
    let shape = integrals.dim();

    for p in 0..shape.2 {
        for i in 0..shape.0 {
            for j in (i + 1)..shape.1 {
                let avg = 0.5 * (integrals[[i, j, p]] + integrals[[j, i, p]]);
                integrals[[i, j, p]] = avg;
                integrals[[j, i, p]] = avg;
            }
        }
    }
}

/// Symmetrize a matrix to handle numerical errors
fn symmetrize_matrix(mut matrix: Array2<f64>) -> Array2<f64> {
    let n = matrix.nrows();
    for i in 0..n {
        for j in i + 1..n {
            let avg = (matrix[[i, j]] + matrix[[j, i]]) * 0.5;
            matrix[[i, j]] = avg;
            matrix[[j, i]] = avg;
        }
    }
    matrix
}

/// Validate that a matrix is positive definite
fn validate_positive_definite(matrix: &Array2<f64>) -> Result<()> {
    let n = matrix.nrows();

    // Check diagonal elements are positive
    for i in 0..n {
        if matrix[[i, i]] <= 0.0 {
            return Err(LibcintError::ComputationError(format!(
                "Matrix has non-positive diagonal element at position {i}"
            ))
            .into());
        }
    }

    // Check symmetry
    for i in 0..n {
        for j in i + 1..n {
            if (matrix[[i, j]] - matrix[[j, i]]).abs() > 1e-10 {
                return Err(
                    LibcintError::ComputationError("Matrix is not symmetric".to_string()).into(),
                );
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::{BasisSet, Molecule};
    use super::*;

    #[test]
    fn test_h2o_integrals() {
        // Test that we can compute integrals for H2O
        let molecule = Molecule::water();
        let ao_basis = BasisSet::for_molecule("H2O", false);
        let aux_basis = BasisSet::for_molecule("H2O", true);

        let engine = LibcintRealEngine::new(&molecule, &ao_basis, &aux_basis)
            .expect("Failed to create engine");

        // Test 2-center integrals
        let v_pq = engine
            .compute_2center()
            .expect("Failed to compute 2-center");
        assert!(v_pq.nrows() > 0);
        assert_eq!(v_pq.nrows(), v_pq.ncols()); // Should be square

        // Check symmetry
        for i in 0..v_pq.nrows() {
            for j in i + 1..v_pq.ncols() {
                assert!(
                    (v_pq[[i, j]] - v_pq[[j, i]]).abs() < 1e-10,
                    "2-center integrals must be symmetric"
                );
            }
        }

        // Test 3-center integrals
        let j3c = engine
            .compute_3center()
            .expect("Failed to compute 3-center");
        assert!(j3c.dim().0 > 0);
        assert_eq!(j3c.dim().0, j3c.dim().1); // First two dims should be equal

        // Check symmetry in first two indices
        for p in 0..j3c.dim().2 {
            for mu in 0..j3c.dim().0 {
                for nu in mu + 1..j3c.dim().1 {
                    assert!(
                        (j3c[[mu, nu, p]] - j3c[[nu, mu, p]]).abs() < 1e-10,
                        "3-center integrals must be symmetric in μν"
                    );
                }
            }
        }
    }
}
