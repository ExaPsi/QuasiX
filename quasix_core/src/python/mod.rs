//! Python bindings module for QuasiX core functionality
//!
//! This module provides PyO3 bindings to expose Rust functionality to Python,
//! enabling high-performance quantum chemistry calculations with a Python interface.
//!
//! # Clean Re-implementation
//!
//! After cleanup, this module exposes:
//! - `df_hdf5_bindings`: Validated S1-S2 DF tensor I/O (KEEP)
//! - `g0w0`: G₀W₀ Python bindings (to be re-implemented)

pub mod df_hdf5_bindings;
pub mod g0w0;

use numpy::PyArray1;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Python binding for Gauss-Legendre grid generation
///
/// Returns GL quadrature nodes and weights on [-1, 1]
///
/// # Arguments
/// * `nfreq` - Number of quadrature points
///
/// # Returns
/// * Dictionary with 'nodes' and 'weights' arrays
#[pyfunction]
pub fn gauss_legendre_grid<'py>(py: Python<'py>, nfreq: usize) -> PyResult<Bound<'py, PyDict>> {
    use crate::freq::imaginary_grid;
    use ndarray::Array1;
    use numpy::PyArray1;

    let (nodes, weights) = imaginary_grid::gauss_legendre_grid(nfreq).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("GL grid error: {}", e))
    })?;

    // Convert Vec to Array1, then to PyArray1
    let nodes_arr = Array1::from_vec(nodes);
    let weights_arr = Array1::from_vec(weights);

    let result = PyDict::new(py);
    result.set_item("nodes", PyArray1::from_array(py, &nodes_arr))?;
    result.set_item("weights", PyArray1::from_array(py, &weights_arr))?;

    Ok(result)
}

/// Python binding for imaginary axis transformation (linear - DEPRECATED)
///
/// **DEPRECATED**: Use `transform_to_imaginary_axis_pyscf()` for GW calculations!
///
/// Transforms GL nodes from [-1, 1] to [0, omega_max]
///
/// # Arguments
/// * `nodes` - GL nodes in [-1, 1]
/// * `omega_max` - Maximum imaginary frequency (Ha)
///
/// # Returns
/// * Transformed frequencies as numpy array
#[pyfunction]
pub fn transform_to_imaginary_axis<'py>(
    py: Python<'py>,
    nodes: numpy::PyReadonlyArray1<f64>,
    omega_max: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    use crate::freq::imaginary_grid;
    use ndarray::Array1;
    use numpy::PyArray1;

    let nodes_vec = nodes.as_slice()?.to_vec();
    #[allow(deprecated)]
    let freqs = imaginary_grid::transform_to_imaginary_axis(&nodes_vec, omega_max);

    let freqs_arr = Array1::from_vec(freqs);
    Ok(PyArray1::from_array(py, &freqs_arr))
}

/// Python binding for imaginary axis transformation (PySCF convention)
///
/// **CRITICAL**: This is the CORRECT transformation for GW calculations!
///
/// Transforms GL nodes from [-1, 1] to [0, ∞) using rational mapping
///
/// # Arguments
/// * `nodes` - GL nodes in [-1, 1]
/// * `x0` - Scaling parameter (default: 0.5 in PySCF)
///
/// # Returns
/// * Transformed frequencies as numpy array [0, ∞)
#[pyfunction]
pub fn transform_to_imaginary_axis_pyscf<'py>(
    py: Python<'py>,
    nodes: numpy::PyReadonlyArray1<f64>,
    x0: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    use crate::freq::imaginary_grid;
    use ndarray::Array1;
    use numpy::PyArray1;

    let nodes_vec = nodes.as_slice()?.to_vec();
    let freqs = imaginary_grid::transform_to_imaginary_axis_pyscf(&nodes_vec, x0);

    let freqs_arr = Array1::from_vec(freqs);
    Ok(PyArray1::from_array(py, &freqs_arr))
}

/// Python binding for weight transformation (linear - DEPRECATED)
///
/// **DEPRECATED**: Use `transform_weights_pyscf()` for GW calculations!
///
/// Applies Jacobian for linear transformation to [0, omega_max]
///
/// # Arguments
/// * `weights` - GL weights on [-1, 1]
/// * `omega_max` - Maximum imaginary frequency (Ha)
///
/// # Returns
/// * Transformed weights as numpy array
#[pyfunction]
pub fn transform_weights<'py>(
    py: Python<'py>,
    weights: numpy::PyReadonlyArray1<f64>,
    omega_max: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    use crate::freq::imaginary_grid;
    use ndarray::Array1;
    use numpy::PyArray1;

    let weights_vec = weights.as_slice()?.to_vec();
    #[allow(deprecated)]
    let transformed = imaginary_grid::transform_weights(&weights_vec, omega_max);

    let transformed_arr = Array1::from_vec(transformed);
    Ok(PyArray1::from_array(py, &transformed_arr))
}

/// Python binding for weight transformation (PySCF convention)
///
/// **CRITICAL**: This is the CORRECT weight transformation for GW calculations!
///
/// Applies Jacobian for rational transformation: dξ/dx = 2*x₀ / (1 - x)²
///
/// # Arguments
/// * `weights` - GL weights on [-1, 1]
/// * `nodes` - GL nodes on [-1, 1] (needed for Jacobian)
/// * `x0` - Scaling parameter (default: 0.5 in PySCF)
///
/// # Returns
/// * Transformed weights for [0, ∞) integration
#[pyfunction]
pub fn transform_weights_pyscf<'py>(
    py: Python<'py>,
    weights: numpy::PyReadonlyArray1<f64>,
    nodes: numpy::PyReadonlyArray1<f64>,
    x0: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    use crate::freq::imaginary_grid;
    use ndarray::Array1;
    use numpy::PyArray1;

    let weights_vec = weights.as_slice()?.to_vec();
    let nodes_vec = nodes.as_slice()?.to_vec();
    let transformed = imaginary_grid::transform_weights_pyscf(&weights_vec, &nodes_vec, x0);

    let transformed_arr = Array1::from_vec(transformed);
    Ok(PyArray1::from_array(py, &transformed_arr))
}

/// Initialize the Python module
pub fn init_python_module(py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    // Register G₀W₀ submodule (placeholder during re-implementation)
    g0w0::register_module(py, m)?;

    // Register DF HDF5 I/O bindings (validated S1-S2 infrastructure)
    df_hdf5_bindings::register_module(py, m)?;

    // Register frequency grid functions (S3-1)
    m.add_function(wrap_pyfunction!(gauss_legendre_grid, m)?)?;
    m.add_function(wrap_pyfunction!(transform_to_imaginary_axis, m)?)?; // DEPRECATED
    m.add_function(wrap_pyfunction!(transform_weights, m)?)?; // DEPRECATED
    m.add_function(wrap_pyfunction!(transform_to_imaginary_axis_pyscf, m)?)?; // CORRECT
    m.add_function(wrap_pyfunction!(transform_weights_pyscf, m)?)?; // CORRECT

    // Register validation submodule
    #[cfg(feature = "python")]
    crate::validation::python_bindings::register_module(py, m)?;

    // Add module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", "QuasiX Team")?;

    Ok(())
}
