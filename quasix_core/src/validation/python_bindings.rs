//! PyO3 bindings for GW100 validation module
//!
//! Exposes high-performance validation functionality to Python with
//! efficient data transfer and error handling.

use numpy::PyReadonlyArray1;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::HashMap;

use super::gw100::{compute_validation_stats, BenchmarkResult, SimdStatistics};

/// Python-exposed benchmark result
#[pyclass(name = "BenchmarkResult")]
#[derive(Clone)]
pub struct PyBenchmarkResult {
    #[pyo3(get)]
    pub molecule: String,
    #[pyo3(get)]
    pub basis_set: String,
    #[pyo3(get)]
    pub aux_basis: String,
    #[pyo3(get)]
    pub qp_energies: Vec<f64>,
    #[pyo3(get)]
    pub z_factors: Vec<f64>,
    #[pyo3(get)]
    pub reference_energies: Vec<f64>,
    #[pyo3(get)]
    pub deviations: Vec<f64>,
    #[pyo3(get)]
    pub ip_error: f64,
    #[pyo3(get)]
    pub ea_error: f64,
    #[pyo3(get)]
    pub wall_time: f64,
    #[pyo3(get)]
    pub converged: bool,
}

#[pymethods]
impl PyBenchmarkResult {
    #[new]
    #[pyo3(signature = (molecule, basis_set, aux_basis, qp_energies, z_factors, reference_energies, homo_idx, lumo_idx, wall_time))]
    fn new(
        molecule: String,
        basis_set: String,
        aux_basis: String,
        qp_energies: Vec<f64>,
        z_factors: Vec<f64>,
        reference_energies: Vec<f64>,
        homo_idx: usize,
        lumo_idx: usize,
        wall_time: f64,
    ) -> Self {
        let result = BenchmarkResult::new(
            molecule.clone(),
            basis_set.clone(),
            aux_basis.clone(),
            qp_energies.clone(),
            z_factors.clone(),
            reference_energies.clone(),
            homo_idx,
            lumo_idx,
            wall_time,
        );

        Self {
            molecule,
            basis_set,
            aux_basis,
            qp_energies,
            z_factors,
            reference_energies,
            deviations: result.deviations,
            ip_error: result.ip_error,
            ea_error: result.ea_error,
            wall_time,
            converged: result.converged,
        }
    }

    /// Get mean absolute deviation
    fn mad(&self) -> f64 {
        if self.deviations.is_empty() {
            return 0.0;
        }
        self.deviations.iter().sum::<f64>() / self.deviations.len() as f64
    }

    /// Check if result passes threshold
    fn passes_threshold(&self, threshold_ev: f64) -> bool {
        self.ip_error <= threshold_ev && self.ea_error <= threshold_ev
    }

    /// String representation
    fn __repr__(&self) -> String {
        format!(
            "BenchmarkResult(molecule='{}', basis='{}', IP_error={:.3} eV, EA_error={:.3} eV, converged={})",
            self.molecule, self.basis_set, self.ip_error, self.ea_error, self.converged
        )
    }
}

/// Python-exposed validation statistics
#[pyclass(name = "ValidationStats")]
#[derive(Clone)]
pub struct PyValidationStats {
    #[pyo3(get)]
    pub mad: f64,
    #[pyo3(get)]
    pub rmse: f64,
    #[pyo3(get)]
    pub max_error: f64,
    #[pyo3(get)]
    pub correlation: f64,
    #[pyo3(get)]
    pub n_samples: usize,
    #[pyo3(get)]
    pub mad_ci_lower: f64,
    #[pyo3(get)]
    pub mad_ci_upper: f64,
    #[pyo3(get)]
    pub n_outliers: usize,
    #[pyo3(get)]
    pub outlier_indices: Vec<usize>,
    #[pyo3(get)]
    pub per_molecule_mad: HashMap<String, f64>,
}

#[pymethods]
impl PyValidationStats {
    /// Check if validation passes thresholds
    fn passes_validation(&self, mad_threshold: f64, correlation_threshold: f64) -> bool {
        self.mad <= mad_threshold && self.correlation >= correlation_threshold
    }

    /// Get summary string
    fn summary(&self) -> String {
        format!(
            "MAD: {:.3} eV, RMSE: {:.3} eV, R²: {:.3}, Max: {:.3} eV, N: {}, Outliers: {}",
            self.mad, self.rmse, self.correlation, self.max_error, self.n_samples, self.n_outliers
        )
    }

    /// Convert to dictionary
    fn to_dict(&self) -> HashMap<String, Py<PyAny>> {
        Python::attach(|py| {
            let mut dict = HashMap::new();
            dict.insert(
                "mad".to_string(),
                self.mad.into_pyobject(py).unwrap().unbind().into_any(),
            );
            dict.insert(
                "rmse".to_string(),
                self.rmse.into_pyobject(py).unwrap().unbind().into_any(),
            );
            dict.insert(
                "max_error".to_string(),
                self.max_error
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict.insert(
                "correlation".to_string(),
                self.correlation
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict.insert(
                "n_samples".to_string(),
                self.n_samples
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict.insert(
                "mad_ci_lower".to_string(),
                self.mad_ci_lower
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict.insert(
                "mad_ci_upper".to_string(),
                self.mad_ci_upper
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict.insert(
                "n_outliers".to_string(),
                self.n_outliers
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict.insert(
                "outlier_indices".to_string(),
                self.outlier_indices
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict.insert(
                "per_molecule_mad".to_string(),
                self.per_molecule_mad
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            dict
        })
    }

    fn __repr__(&self) -> String {
        self.summary()
    }
}

/// Compute MAD using SIMD acceleration
#[pyfunction]
#[pyo3(name = "compute_mad_simd")]
pub fn py_compute_mad_simd(
    _py: Python<'_>,
    calculated: PyReadonlyArray1<f64>,
    reference: PyReadonlyArray1<f64>,
) -> PyResult<f64> {
    let calc = calculated.as_slice()?;
    let ref_vals = reference.as_slice()?;

    SimdStatistics::compute_mad_simd(calc, ref_vals)
        .map_err(|e| PyRuntimeError::new_err(format!("MAD computation failed: {}", e)))
}

/// Compute RMSE
#[pyfunction]
#[pyo3(name = "compute_rmse")]
pub fn py_compute_rmse(
    _py: Python<'_>,
    calculated: PyReadonlyArray1<f64>,
    reference: PyReadonlyArray1<f64>,
) -> PyResult<f64> {
    let calc = calculated.as_slice()?;
    let ref_vals = reference.as_slice()?;

    SimdStatistics::compute_rmse(calc, ref_vals)
        .map_err(|e| PyRuntimeError::new_err(format!("RMSE computation failed: {}", e)))
}

/// Compute correlation coefficient R²
#[pyfunction]
#[pyo3(name = "compute_correlation")]
pub fn py_compute_correlation(
    _py: Python<'_>,
    x: PyReadonlyArray1<f64>,
    y: PyReadonlyArray1<f64>,
) -> PyResult<f64> {
    let x_slice = x.as_slice()?;
    let y_slice = y.as_slice()?;

    SimdStatistics::compute_correlation(x_slice, y_slice)
        .map_err(|e| PyRuntimeError::new_err(format!("Correlation computation failed: {}", e)))
}

/// Detect outliers in data
#[pyfunction]
#[pyo3(name = "detect_outliers")]
pub fn py_detect_outliers(
    _py: Python<'_>,
    values: PyReadonlyArray1<f64>,
    threshold: f64,
) -> PyResult<Vec<usize>> {
    let vals = values.as_slice()?;
    Ok(SimdStatistics::detect_outliers(vals, threshold))
}

/// Compute validation statistics from benchmark results
#[pyfunction]
#[pyo3(name = "compute_validation_stats")]
pub fn py_compute_validation_stats(results: Vec<PyBenchmarkResult>) -> PyResult<PyValidationStats> {
    // Convert Python results to Rust results
    let rust_results: Vec<BenchmarkResult> = results
        .iter()
        .map(|r| {
            let mut result = BenchmarkResult::new(
                r.molecule.clone(),
                r.basis_set.clone(),
                r.aux_basis.clone(),
                r.qp_energies.clone(),
                r.z_factors.clone(),
                r.reference_energies.clone(),
                0, // We don't need homo_idx for statistics
                0, // We don't need lumo_idx for statistics
                r.wall_time,
            );
            result.ip_error = r.ip_error;
            result.ea_error = r.ea_error;
            result.converged = r.converged;
            result.deviations = r.deviations.clone();
            result
        })
        .collect();

    // Compute statistics
    let stats = compute_validation_stats(&rust_results)
        .map_err(|e| PyRuntimeError::new_err(format!("Statistics computation failed: {}", e)))?;

    // Convert to Python stats
    Ok(PyValidationStats {
        mad: stats.mad,
        rmse: stats.rmse,
        max_error: stats.max_error,
        correlation: stats.correlation,
        n_samples: stats.n_samples,
        mad_ci_lower: stats.mad_ci_lower,
        mad_ci_upper: stats.mad_ci_upper,
        n_outliers: stats.n_outliers,
        outlier_indices: stats.outlier_indices,
        per_molecule_mad: stats.per_molecule_mad,
    })
}

/// Bootstrap confidence interval calculation
#[pyfunction]
#[pyo3(name = "bootstrap_confidence_interval")]
pub fn py_bootstrap_ci(
    _py: Python<'_>,
    calculated: PyReadonlyArray1<f64>,
    reference: PyReadonlyArray1<f64>,
    n_bootstrap: usize,
    confidence_level: f64,
) -> PyResult<(f64, f64)> {
    let calc = calculated.as_slice()?;
    let ref_vals = reference.as_slice()?;

    super::statistics::bootstrap_confidence_interval(calc, ref_vals, n_bootstrap, confidence_level)
        .map_err(|e| PyRuntimeError::new_err(format!("Bootstrap CI failed: {}", e)))
}

/// Register the validation module with Python
pub fn register_module(py: Python<'_>, parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let validation_module = PyModule::new(py, "validation")?;

    // Add classes
    validation_module.add_class::<PyBenchmarkResult>()?;
    validation_module.add_class::<PyValidationStats>()?;

    // Add functions
    validation_module.add_function(wrap_pyfunction!(py_compute_mad_simd, &validation_module)?)?;
    validation_module.add_function(wrap_pyfunction!(py_compute_rmse, &validation_module)?)?;
    validation_module.add_function(wrap_pyfunction!(
        py_compute_correlation,
        &validation_module
    )?)?;
    validation_module.add_function(wrap_pyfunction!(py_detect_outliers, &validation_module)?)?;
    validation_module.add_function(wrap_pyfunction!(
        py_compute_validation_stats,
        &validation_module
    )?)?;
    validation_module.add_function(wrap_pyfunction!(py_bootstrap_ci, &validation_module)?)?;

    parent_module.add_submodule(&validation_module)?;

    Ok(())
}
