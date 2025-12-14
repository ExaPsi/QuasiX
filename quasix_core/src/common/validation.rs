//! Comprehensive validation guards for tensor inputs and numerical sanity
//!
//! This module provides validation utilities to ensure:
//! 1. Input tensors are non-zero and physically meaningful
//! 2. Dimensions are consistent
//! 3. Numerical values are within acceptable ranges
//! 4. BLAS/LAPACK operations will not fail due to invalid inputs
//!
//! # Golden Rule Compliance
//! QuasiX is a HIGH QUALITY Quantum Chemistry package. ALL calculations MUST use REAL libraries
//! giving results comparable to PySCF. These validation guards ensure data quality before
//! expensive computations.

use crate::common::{QuasixError, Result};
use ndarray::{ArrayView1, ArrayView2, ArrayView3};
use num_complex::Complex64;

/// Validation configuration
#[derive(Debug, Clone)]
pub struct ValidationConfig {
    /// Minimum allowed non-zero value (for checking tensor sparsity)
    pub min_nonzero_threshold: f64,
    /// Maximum allowed absolute value (for detecting numerical overflow)
    pub max_abs_value: f64,
    /// Minimum fraction of non-zero elements
    pub min_nonzero_fraction: f64,
    /// Enable strict NaN/Inf checks
    pub check_finite: bool,
    /// Enable positive definiteness checks for metric matrices
    pub check_positive_definite: bool,
}

impl Default for ValidationConfig {
    fn default() -> Self {
        Self {
            min_nonzero_threshold: 1e-16,
            max_abs_value: 1e10,
            min_nonzero_fraction: 0.0001, // At least 0.01% non-zero
            check_finite: true,
            check_positive_definite: true,
        }
    }
}

/// Validate 1D array (energies, occupations, etc.)
pub fn validate_array1(arr: &ArrayView1<f64>, name: &str, config: &ValidationConfig) -> Result<()> {
    let len = arr.len();

    if len == 0 {
        return Err(QuasixError::InvalidInput(format!(
            "{} is empty (length 0)",
            name
        )));
    }

    // Check for NaN/Inf
    if config.check_finite {
        for (i, &val) in arr.iter().enumerate() {
            if !val.is_finite() {
                return Err(QuasixError::NumericalError(format!(
                    "{} contains non-finite value {} at index {}",
                    name, val, i
                )));
            }
        }
    }

    // Check for reasonable magnitude
    let max_abs = arr.iter().fold(0.0_f64, |max, &x| max.max(x.abs()));
    if max_abs > config.max_abs_value {
        return Err(QuasixError::NumericalError(format!(
            "{} has unreasonably large values (max_abs = {:.2e}, threshold = {:.2e})",
            name, max_abs, config.max_abs_value
        )));
    }

    // Check that array is not all zeros
    let nonzero_count = arr
        .iter()
        .filter(|&&x| x.abs() > config.min_nonzero_threshold)
        .count();
    if nonzero_count == 0 {
        return Err(QuasixError::InvalidInput(format!(
            "{} is all zeros (threshold = {:.2e})",
            name, config.min_nonzero_threshold
        )));
    }

    Ok(())
}

/// Validate complex 1D array
pub fn validate_array1_complex(
    arr: &ArrayView1<Complex64>,
    name: &str,
    config: &ValidationConfig,
) -> Result<()> {
    let len = arr.len();

    if len == 0 {
        return Err(QuasixError::InvalidInput(format!(
            "{} is empty (length 0)",
            name
        )));
    }

    // Check for NaN/Inf
    if config.check_finite {
        for (i, val) in arr.iter().enumerate() {
            if !val.re.is_finite() || !val.im.is_finite() {
                return Err(QuasixError::NumericalError(format!(
                    "{} contains non-finite value {:?} at index {}",
                    name, val, i
                )));
            }
        }
    }

    // Check for reasonable magnitude
    let max_abs = arr.iter().fold(0.0_f64, |max, x| max.max(x.norm()));
    if max_abs > config.max_abs_value {
        return Err(QuasixError::NumericalError(format!(
            "{} has unreasonably large values (max_abs = {:.2e}, threshold = {:.2e})",
            name, max_abs, config.max_abs_value
        )));
    }

    Ok(())
}

/// Validate 2D array (matrices)
pub fn validate_array2(arr: &ArrayView2<f64>, name: &str, config: &ValidationConfig) -> Result<()> {
    let shape = arr.shape();

    if shape[0] == 0 || shape[1] == 0 {
        return Err(QuasixError::InvalidInput(format!(
            "{} has zero dimension: shape = [{}, {}]",
            name, shape[0], shape[1]
        )));
    }

    // Check for NaN/Inf
    if config.check_finite {
        for (idx, &val) in arr.indexed_iter() {
            if !val.is_finite() {
                return Err(QuasixError::NumericalError(format!(
                    "{} contains non-finite value {} at index {:?}",
                    name, val, idx
                )));
            }
        }
    }

    // Check for reasonable magnitude
    let max_abs = arr.iter().fold(0.0_f64, |max, &x| max.max(x.abs()));
    if max_abs > config.max_abs_value {
        return Err(QuasixError::NumericalError(format!(
            "{} has unreasonably large values (max_abs = {:.2e}, threshold = {:.2e})",
            name, max_abs, config.max_abs_value
        )));
    }

    // Check that matrix is not all zeros
    let total_elements = shape[0] * shape[1];
    let nonzero_count = arr
        .iter()
        .filter(|&&x| x.abs() > config.min_nonzero_threshold)
        .count();
    let nonzero_fraction = nonzero_count as f64 / total_elements as f64;

    if nonzero_fraction < config.min_nonzero_fraction {
        return Err(QuasixError::InvalidInput(format!(
            "{} is too sparse: only {:.2}% non-zero elements (threshold = {:.2}%)",
            name,
            nonzero_fraction * 100.0,
            config.min_nonzero_fraction * 100.0
        )));
    }

    Ok(())
}

/// Validate complex 2D array
pub fn validate_array2_complex(
    arr: &ArrayView2<Complex64>,
    name: &str,
    config: &ValidationConfig,
) -> Result<()> {
    let shape = arr.shape();

    if shape[0] == 0 || shape[1] == 0 {
        return Err(QuasixError::InvalidInput(format!(
            "{} has zero dimension: shape = [{}, {}]",
            name, shape[0], shape[1]
        )));
    }

    // Check for NaN/Inf
    if config.check_finite {
        for (idx, val) in arr.indexed_iter() {
            if !val.re.is_finite() || !val.im.is_finite() {
                return Err(QuasixError::NumericalError(format!(
                    "{} contains non-finite value {:?} at index {:?}",
                    name, val, idx
                )));
            }
        }
    }

    // Check for reasonable magnitude
    let max_abs = arr.iter().fold(0.0_f64, |max, x| max.max(x.norm()));
    if max_abs > config.max_abs_value {
        return Err(QuasixError::NumericalError(format!(
            "{} has unreasonably large values (max_abs = {:.2e}, threshold = {:.2e})",
            name, max_abs, config.max_abs_value
        )));
    }

    Ok(())
}

/// Validate 3D array (tensors)
pub fn validate_array3(arr: &ArrayView3<f64>, name: &str, config: &ValidationConfig) -> Result<()> {
    let shape = arr.shape();

    if shape[0] == 0 || shape[1] == 0 || shape[2] == 0 {
        return Err(QuasixError::InvalidInput(format!(
            "{} has zero dimension: shape = [{}, {}, {}]",
            name, shape[0], shape[1], shape[2]
        )));
    }

    // Check for NaN/Inf
    if config.check_finite {
        for (idx, &val) in arr.indexed_iter() {
            if !val.is_finite() {
                return Err(QuasixError::NumericalError(format!(
                    "{} contains non-finite value {} at index {:?}",
                    name, val, idx
                )));
            }
        }
    }

    // Check for reasonable magnitude
    let max_abs = arr.iter().fold(0.0_f64, |max, &x| max.max(x.abs()));
    if max_abs > config.max_abs_value {
        return Err(QuasixError::NumericalError(format!(
            "{} has unreasonably large values (max_abs = {:.2e}, threshold = {:.2e})",
            name, max_abs, config.max_abs_value
        )));
    }

    // Check that tensor is not all zeros (sample-based for efficiency)
    let total_elements = shape[0] * shape[1] * shape[2];
    let sample_size = (total_elements / 100).max(1000).min(total_elements); // Sample 1% or 1000 elements
    let nonzero_count = arr
        .iter()
        .take(sample_size)
        .filter(|&&x| x.abs() > config.min_nonzero_threshold)
        .count();
    let nonzero_fraction = nonzero_count as f64 / sample_size as f64;

    if nonzero_fraction < config.min_nonzero_fraction {
        return Err(QuasixError::InvalidInput(format!(
            "{} is too sparse: only {:.2}% non-zero elements in sample (threshold = {:.2}%)",
            name,
            nonzero_fraction * 100.0,
            config.min_nonzero_fraction * 100.0
        )));
    }

    Ok(())
}

/// Validate MO energies (specific checks for eigenvalues)
pub fn validate_mo_energies(energies: &ArrayView1<f64>) -> Result<()> {
    let config = ValidationConfig::default();

    // Basic validation
    validate_array1(energies, "MO energies", &config)?;

    // Check for reasonable energy range (-100 to +100 Ha typical)
    let min_energy = energies.iter().fold(f64::INFINITY, |min, &x| min.min(x));
    let max_energy = energies
        .iter()
        .fold(f64::NEG_INFINITY, |max, &x| max.max(x));

    if min_energy < -200.0 || max_energy > 200.0 {
        return Err(QuasixError::NumericalError(format!(
            "MO energies outside reasonable range: [{:.2}, {:.2}] Ha (expected roughly [-100, +100] Ha)",
            min_energy, max_energy
        )));
    }

    Ok(())
}

/// Validate MO occupations (must be between 0 and 2 for RHF, 0 and 1 for UHF)
pub fn validate_mo_occupations(occupations: &ArrayView1<f64>, max_occ: f64) -> Result<()> {
    let config = ValidationConfig {
        min_nonzero_fraction: 0.0, // Occupations can have many zeros (virtual orbitals)
        ..Default::default()
    };

    // Basic validation
    validate_array1(occupations, "MO occupations", &config)?;

    // Check range
    for (i, &occ) in occupations.iter().enumerate() {
        if occ < 0.0 || occ > max_occ + 1e-10 {
            return Err(QuasixError::InvalidInput(format!(
                "MO occupation {} at index {} is out of valid range [0, {}]",
                occ, i, max_occ
            )));
        }
    }

    // Check that there are some occupied orbitals
    let n_occupied = occupations.iter().filter(|&&occ| occ > 0.1).count();
    if n_occupied == 0 {
        return Err(QuasixError::InvalidInput(
            "No occupied orbitals found (all occupations near zero)".to_string(),
        ));
    }

    Ok(())
}

/// Validate DF metric matrix (should be symmetric positive definite)
pub fn validate_df_metric(metric: &ArrayView2<f64>) -> Result<()> {
    let config = ValidationConfig::default();

    // Basic validation
    validate_array2(metric, "DF metric", &config)?;

    // Check square
    let shape = metric.shape();
    if shape[0] != shape[1] {
        return Err(QuasixError::InvalidInput(format!(
            "DF metric must be square, got shape [{}, {}]",
            shape[0], shape[1]
        )));
    }

    // Check symmetry (sample diagonal and off-diagonal for large matrices)
    let n = shape[0];
    let check_points = if n > 100 {
        // For large matrices, check diagonal + sample of off-diagonal
        n.min(100)
    } else {
        // For small matrices, check all
        n
    };

    for i in 0..check_points {
        for j in i..check_points {
            let diff = (metric[[i, j]] - metric[[j, i]]).abs();
            if diff > 1e-10 {
                return Err(QuasixError::NumericalError(format!(
                    "DF metric is not symmetric: |V[{},{}] - V[{},{}]| = {:.2e}",
                    i, j, j, i, diff
                )));
            }
        }
    }

    // Check diagonal is positive
    for i in 0..n {
        if metric[[i, i]] <= 0.0 {
            return Err(QuasixError::NumericalError(format!(
                "DF metric has non-positive diagonal element at index {}: {}",
                i,
                metric[[i, i]]
            )));
        }
    }

    Ok(())
}

/// Validate dimension consistency for GW calculation
pub fn validate_gw_dimensions(
    n_mo: usize,
    n_aux: usize,
    mo_energy: &ArrayView1<f64>,
    mo_occ: &ArrayView1<f64>,
    df_tensor: &ArrayView3<f64>,
) -> Result<()> {
    // Check MO energy dimension
    if mo_energy.len() != n_mo {
        return Err(QuasixError::InvalidInput(format!(
            "MO energy length {} does not match n_mo = {}",
            mo_energy.len(),
            n_mo
        )));
    }

    // Check MO occupation dimension
    if mo_occ.len() != n_mo {
        return Err(QuasixError::InvalidInput(format!(
            "MO occupation length {} does not match n_mo = {}",
            mo_occ.len(),
            n_mo
        )));
    }

    // Check DF tensor dimensions
    let df_shape = df_tensor.shape();
    if df_shape[2] != n_aux {
        return Err(QuasixError::InvalidInput(format!(
            "DF tensor auxiliary dimension {} does not match n_aux = {}",
            df_shape[2], n_aux
        )));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array1;

    #[test]
    fn test_validate_array1_good() {
        let arr = Array1::from_vec(vec![1.0, 2.0, 3.0]);
        let config = ValidationConfig::default();
        assert!(validate_array1(&arr.view(), "test", &config).is_ok());
    }

    #[test]
    fn test_validate_array1_all_zeros() {
        let arr = Array1::from_vec(vec![0.0, 0.0, 0.0]);
        let config = ValidationConfig::default();
        assert!(validate_array1(&arr.view(), "test", &config).is_err());
    }

    #[test]
    fn test_validate_mo_energies_good() {
        let energies = Array1::from_vec(vec![-10.0, -5.0, -1.0, 2.0, 5.0]);
        assert!(validate_mo_energies(&energies.view()).is_ok());
    }

    #[test]
    fn test_validate_mo_occupations_good() {
        let occs = Array1::from_vec(vec![2.0, 2.0, 2.0, 0.0, 0.0]);
        assert!(validate_mo_occupations(&occs.view(), 2.0).is_ok());
    }

    #[test]
    fn test_validate_mo_occupations_bad_range() {
        let occs = Array1::from_vec(vec![2.0, 3.0, 2.0, 0.0, 0.0]); // 3.0 is out of range
        assert!(validate_mo_occupations(&occs.view(), 2.0).is_err());
    }
}
