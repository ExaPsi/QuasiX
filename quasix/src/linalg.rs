// Python bindings for linear algebra utilities

use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use quasix_core::linalg::MetricSqrt;

/// Compute the square root decomposition of a metric tensor
///
/// Parameters
/// ----------
/// metric : numpy.ndarray
///     Symmetric positive semi-definite matrix [n_aux, n_aux]
/// threshold : float, optional
///     Eigenvalue cutoff (default 1e-10)
///
/// Returns
/// -------
/// dict
///     Dictionary containing:
///     - 'sqrt_eigenvals': Square root of kept eigenvalues
///     - 'inv_sqrt_eigenvals': Inverse square root of kept eigenvalues
///     - 'eigenvecs_kept': Kept eigenvectors [n_aux, n_kept]
///     - 'condition_number': Condition number of the matrix
///     - 'n_kept': Number of kept eigenvalues
///     - 'n_total': Total dimension
///     - 'identity_error': Maximum error in v^{-1/2} v v^{-1/2} = I test
#[pyfunction]
#[pyo3(signature = (metric, threshold=None))]
pub fn compute_metric_sqrt<'py>(
    py: Python<'py>,
    metric: &Bound<'py, PyArray2<f64>>,
    threshold: Option<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    // Convert numpy array to ndarray
    let metric = metric.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to metric: {}",
            e
        ))
    })?;
    let metric_array = metric.as_array();
    let metric_owned = metric_array.to_owned();

    // Compute the decomposition
    let metric_sqrt = MetricSqrt::compute(&metric_owned, threshold).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Failed to compute metric sqrt: {}",
            e
        ))
    })?;

    // Verify identity
    let identity_error = metric_sqrt.verify_identity(&metric_owned).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to verify identity: {}", e))
    })?;

    // Create result dictionary
    let result = PyDict::new(py);
    result.set_item(
        "sqrt_eigenvals",
        PyArray1::from_owned_array(py, metric_sqrt.sqrt_eigenvals),
    )?;
    result.set_item(
        "inv_sqrt_eigenvals",
        PyArray1::from_owned_array(py, metric_sqrt.inv_sqrt_eigenvals),
    )?;
    result.set_item(
        "eigenvecs_kept",
        PyArray2::from_owned_array(py, metric_sqrt.eigenvecs_kept),
    )?;
    result.set_item(
        "eigenvals",
        PyArray1::from_owned_array(py, metric_sqrt.eigenvals),
    )?;
    result.set_item("condition_number", metric_sqrt.condition_number)?;
    result.set_item("n_kept", metric_sqrt.n_kept)?;
    result.set_item("n_total", metric_sqrt.n_total)?;
    result.set_item("n_negative", metric_sqrt.n_negative)?;
    result.set_item("n_zero", metric_sqrt.n_zero)?;
    result.set_item("threshold_used", metric_sqrt.threshold_used)?;
    result.set_item("reconstruction_error", metric_sqrt.reconstruction_error)?;
    result.set_item("identity_error", identity_error)?;

    Ok(result)
}

/// Apply v^{1/2} to a vector or matrix
///
/// Parameters
/// ----------
/// metric : numpy.ndarray
///     Symmetric positive semi-definite matrix [n_aux, n_aux]
/// x : numpy.ndarray
///     Vector [n_aux] or matrix [n_aux, n_cols] to transform
/// threshold : float, optional
///     Eigenvalue cutoff (default 1e-10)
///
/// Returns
/// -------
/// numpy.ndarray
///     Transformed vector/matrix v^{1/2} @ x
#[pyfunction]
#[pyo3(signature = (metric, x, threshold=None))]
pub fn apply_metric_sqrt<'py>(
    py: Python<'py>,
    metric: &Bound<'py, PyArray2<f64>>,
    x: &Bound<'py, PyArray2<f64>>,
    threshold: Option<f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    // Convert numpy arrays to ndarray
    let metric = metric.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to metric: {}",
            e
        ))
    })?;
    let metric_array = metric.as_array();
    let metric_owned = metric_array.to_owned();

    let x_readonly = x.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to x: {}",
            e
        ))
    })?;
    let x_array = x_readonly.as_array();

    // Compute the decomposition
    let metric_sqrt = MetricSqrt::compute(&metric_owned, threshold).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Failed to compute metric sqrt: {}",
            e
        ))
    })?;

    // Apply the transformation
    let result = metric_sqrt.apply_sqrt_matrix(&x_array).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to apply sqrt: {}", e))
    })?;

    Ok(PyArray2::from_owned_array(py, result))
}

/// Apply v^{-1/2} to a vector or matrix
///
/// Parameters
/// ----------
/// metric : numpy.ndarray
///     Symmetric positive semi-definite matrix [n_aux, n_aux]
/// x : numpy.ndarray
///     Vector [n_aux] or matrix [n_aux, n_cols] to transform
/// threshold : float, optional
///     Eigenvalue cutoff (default 1e-10)
///
/// Returns
/// -------
/// numpy.ndarray
///     Transformed vector/matrix v^{-1/2} @ x
#[pyfunction]
#[pyo3(signature = (metric, x, threshold=None))]
pub fn apply_metric_inv_sqrt<'py>(
    py: Python<'py>,
    metric: &Bound<'py, PyArray2<f64>>,
    x: &Bound<'py, PyArray2<f64>>,
    threshold: Option<f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    // Convert numpy arrays to ndarray
    let metric = metric.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to metric: {}",
            e
        ))
    })?;
    let metric_array = metric.as_array();
    let metric_owned = metric_array.to_owned();

    let x_readonly = x.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to x: {}",
            e
        ))
    })?;
    let x_array = x_readonly.as_array();

    // Compute the decomposition
    let metric_sqrt = MetricSqrt::compute(&metric_owned, threshold).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Failed to compute metric sqrt: {}",
            e
        ))
    })?;

    // Apply the transformation
    let result = metric_sqrt.apply_inv_sqrt_matrix(&x_array).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to apply inv_sqrt: {}", e))
    })?;

    Ok(PyArray2::from_owned_array(py, result))
}

/// Verify the identity v^{-1/2} v v^{-1/2} = I
///
/// Parameters
/// ----------
/// metric : numpy.ndarray
///     Symmetric positive semi-definite matrix [n_aux, n_aux]
/// threshold : float, optional
///     Eigenvalue cutoff (default 1e-10)
///
/// Returns
/// -------
/// float
///     Maximum absolute error in the identity test
#[pyfunction]
#[pyo3(signature = (metric, threshold=None))]
pub fn verify_metric_identity(
    metric: &Bound<'_, PyArray2<f64>>,
    threshold: Option<f64>,
) -> PyResult<f64> {
    // Convert numpy array to ndarray
    let metric = metric.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to metric: {}",
            e
        ))
    })?;
    let metric_array = metric.as_array();
    let metric_owned = metric_array.to_owned();

    // Compute the decomposition
    let metric_sqrt = MetricSqrt::compute(&metric_owned, threshold).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Failed to compute metric sqrt: {}",
            e
        ))
    })?;

    // Verify identity
    let identity_error = metric_sqrt.verify_identity(&metric_owned).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to verify identity: {}", e))
    })?;

    Ok(identity_error)
}

/// Register the linalg module with Python
#[allow(dead_code)] // Preserved for future module registration
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let linalg_module = PyModule::new(parent_module.py(), "linalg")?;
    linalg_module.add_function(wrap_pyfunction!(compute_metric_sqrt, &linalg_module)?)?;
    linalg_module.add_function(wrap_pyfunction!(apply_metric_sqrt, &linalg_module)?)?;
    linalg_module.add_function(wrap_pyfunction!(apply_metric_inv_sqrt, &linalg_module)?)?;
    linalg_module.add_function(wrap_pyfunction!(verify_metric_identity, &linalg_module)?)?;
    parent_module.add_submodule(&linalg_module)?;
    Ok(())
}
