//! PyO3 bindings for analytic continuation module
//!
//! This module exposes the high-performance Rust analytic continuation
//! implementation to Python, providing direct access to multipole and
//! Padé models with SIMD optimization and parallel processing.

use ndarray::Array1;
use num_complex::Complex64;
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use quasix_core::freq::{
    ACConfig, AnalyticContinuation, AnalyticContinuationFitter, ModelType, MultipoleModel,
    PadeModel,
};

/// Configuration for analytic continuation
///
/// This class configures the behavior of analytic continuation models.
///
/// Args:
///     model_selection (str): Model type selection ('multipole' or 'pade', default: 'multipole')
///     max_poles (int): Maximum number of poles for multipole model (default: 20)
///     max_pade_order (tuple): Maximum order for Padé approximant [M, N] (default: (15, 15))
///     cv_fraction (float): Fraction of data for cross-validation (default: 0.3)
///     cv_iterations (int): Number of CV iterations for error estimation (default: 5)
///     regularization (float): Regularization parameter for fitting (default: 1e-10)
///     stability_threshold (float): Stability threshold for pole analysis (default: 1e-3)
///     convergence_tol (float): Convergence tolerance for optimization (default: 1e-8)
///     max_iterations (int): Maximum iterations for fitting (default: 1000)
///     eta (float): Broadening parameter η for real-axis evaluation (default: 0.01)
///     parallel (bool): Enable parallel fitting over matrix elements (default: True)
#[pyclass(name = "ACConfig")]
#[derive(Clone)]
pub struct PyACConfig {
    inner: ACConfig,
}

#[pymethods]
impl PyACConfig {
    #[new]
    #[pyo3(signature = (
        model_selection="multipole",
        max_poles=20,
        max_pade_order=(15, 15),
        cv_fraction=0.3,
        cv_iterations=5,
        regularization=1e-10,
        stability_threshold=1e-3,
        convergence_tol=1e-8,
        max_iterations=1000,
        eta=0.01,
        parallel=true
    ))]
    fn new(
        model_selection: &str,
        max_poles: usize,
        max_pade_order: (usize, usize),
        cv_fraction: f64,
        cv_iterations: usize,
        regularization: f64,
        stability_threshold: f64,
        convergence_tol: f64,
        max_iterations: usize,
        eta: f64,
        parallel: bool,
    ) -> PyResult<Self> {
        let model_type = match model_selection.to_lowercase().as_str() {
            "multipole" => ModelType::Multipole,
            "pade" => ModelType::Pade,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Invalid model_selection: {}. Must be 'multipole' or 'pade'",
                    model_selection
                )))
            }
        };

        Ok(Self {
            inner: ACConfig {
                model_selection: model_type,
                max_poles,
                max_pade_order,
                cv_fraction,
                cv_iterations,
                regularization,
                stability_threshold,
                convergence_tol,
                max_iterations,
                eta,
                parallel,
            },
        })
    }

    /// Get configuration as dictionary
    fn as_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        let model_str = match self.inner.model_selection {
            ModelType::Multipole => "multipole",
            ModelType::Pade => "pade",
        };
        dict.set_item("model_selection", model_str)?;
        dict.set_item("max_poles", self.inner.max_poles)?;
        dict.set_item("max_pade_order", self.inner.max_pade_order)?;
        dict.set_item("cv_fraction", self.inner.cv_fraction)?;
        dict.set_item("cv_iterations", self.inner.cv_iterations)?;
        dict.set_item("regularization", self.inner.regularization)?;
        dict.set_item("stability_threshold", self.inner.stability_threshold)?;
        dict.set_item("convergence_tol", self.inner.convergence_tol)?;
        dict.set_item("max_iterations", self.inner.max_iterations)?;
        dict.set_item("eta", self.inner.eta)?;
        dict.set_item("parallel", self.inner.parallel)?;
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "ACConfig(max_poles={}, max_pade_order={:?}, cv_fraction={}, parallel={})",
            self.inner.max_poles,
            self.inner.max_pade_order,
            self.inner.cv_fraction,
            self.inner.parallel
        )
    }
}

/// Multi-pole expansion model for analytic continuation
///
/// This model fits data to a sum of poles: f(z) = Σ_p A_p/(z - ω_p) + polynomial
///
/// Args:
///     n_poles (int): Number of poles in the model
///     config (ACConfig): Configuration object (optional)
#[pyclass(name = "MultipoleModel")]
pub struct PyMultipoleModel {
    inner: Box<MultipoleModel>,
}

#[pymethods]
impl PyMultipoleModel {
    #[new]
    #[pyo3(signature = (n_poles, config=None))]
    fn new(n_poles: usize, config: Option<PyACConfig>) -> Self {
        let cfg = config.map(|c| c.inner).unwrap_or_default();
        Self {
            inner: Box::new(MultipoleModel::new(n_poles, cfg)),
        }
    }

    /// Fit the model to imaginary-axis data
    ///
    /// Args:
    ///     xi_data: Imaginary frequencies, shape (n_points,)
    ///     f_data_real: Real part of function values, shape (n_points,)
    ///     f_data_imag: Imaginary part of function values, shape (n_points,)
    ///
    /// Returns:
    ///     None (modifies model in-place)
    fn fit<'py>(
        &mut self,
        xi_data: &Bound<'py, PyArray1<f64>>,
        f_data_real: &Bound<'py, PyArray1<f64>>,
        f_data_imag: &Bound<'py, PyArray1<f64>>,
    ) -> PyResult<()> {
        // Convert numpy arrays
        let xi = xi_data.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read xi_data: {}",
                e
            ))
        })?;
        let xi_array = xi.as_array().to_owned();

        let f_real = f_data_real.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read f_data_real: {}",
                e
            ))
        })?;
        let f_real_array = f_real.as_array();

        let f_imag = f_data_imag.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read f_data_imag: {}",
                e
            ))
        })?;
        let f_imag_array = f_imag.as_array();

        // Combine into complex array
        let f_complex: Array1<Complex64> = f_real_array
            .iter()
            .zip(f_imag_array.iter())
            .map(|(&re, &im)| Complex64::new(re, im))
            .collect();

        // Perform fitting
        self.inner.fit(&xi_array, &f_complex).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Fitting failed: {:?}", e))
        })?;

        Ok(())
    }

    /// Evaluate the model at complex frequencies
    ///
    /// Args:
    ///     z_real: Real parts of evaluation points, shape (n_eval,)
    ///     z_imag: Imaginary parts of evaluation points, shape (n_eval,)
    ///
    /// Returns:
    ///     Tuple of (real, imag) arrays with function values
    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        z_real: &Bound<'py, PyArray1<f64>>,
        z_imag: &Bound<'py, PyArray1<f64>>,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
        // Convert input arrays
        let z_re = z_real.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read z_real: {}",
                e
            ))
        })?;
        let z_re_array = z_re.as_array();

        let z_im = z_imag.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read z_imag: {}",
                e
            ))
        })?;
        let z_im_array = z_im.as_array();

        // Create complex array
        let z_complex: Array1<Complex64> = z_re_array
            .iter()
            .zip(z_im_array.iter())
            .map(|(&re, &im)| Complex64::new(re, im))
            .collect();

        // Evaluate using batch method
        let result = self.inner.evaluate_batch(&z_complex);

        // Split into real and imaginary parts
        let result_real: Array1<f64> = result.mapv(|c| c.re);
        let result_imag: Array1<f64> = result.mapv(|c| c.im);

        Ok((
            PyArray1::from_owned_array(py, result_real),
            PyArray1::from_owned_array(py, result_imag),
        ))
    }

    /// Get cross-validation error
    fn cross_validation_error(&self) -> f64 {
        self.inner.cross_validation_error()
    }

    /// Check stability of the fitted model
    fn stability_check(&self) -> bool {
        self.inner.stability_check()
    }

    /// Get poles and residues
    ///
    /// Returns:
    ///     Dictionary with 'poles_real', 'poles_imag', 'residues_real', 'residues_imag'
    fn get_poles_residues<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        if let Some((poles, residues)) = self.inner.get_poles_residues() {
            let poles_real: Array1<f64> = poles.mapv(|c| c.re);
            let poles_imag: Array1<f64> = poles.mapv(|c| c.im);
            let residues_real: Array1<f64> = residues.mapv(|c| c.re);
            let residues_imag: Array1<f64> = residues.mapv(|c| c.im);

            dict.set_item("poles_real", PyArray1::from_owned_array(py, poles_real))?;
            dict.set_item("poles_imag", PyArray1::from_owned_array(py, poles_imag))?;
            dict.set_item(
                "residues_real",
                PyArray1::from_owned_array(py, residues_real),
            )?;
            dict.set_item(
                "residues_imag",
                PyArray1::from_owned_array(py, residues_imag),
            )?;
        }

        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "MultipoleModel(cv_error={:.3e})",
            self.inner.cross_validation_error()
        )
    }
}

/// Padé approximant model for analytic continuation
///
/// This model fits data to a rational function: f(z) ≈ P_M(z)/Q_N(z)
///
/// Args:
///     order (tuple): Order of approximant (M, N) for numerator and denominator
///     config (ACConfig): Configuration object (optional)
#[pyclass(name = "PadeModel")]
pub struct PyPadeModel {
    inner: Box<PadeModel>,
}

#[pymethods]
impl PyPadeModel {
    #[new]
    #[pyo3(signature = (order, config=None))]
    fn new(order: (usize, usize), config: Option<PyACConfig>) -> Self {
        let cfg = config.map(|c| c.inner).unwrap_or_default();
        Self {
            inner: Box::new(PadeModel::new(order, cfg)),
        }
    }

    /// Fit the model to imaginary-axis data
    ///
    /// Args:
    ///     xi_data: Imaginary frequencies, shape (n_points,)
    ///     f_data_real: Real part of function values, shape (n_points,)
    ///     f_data_imag: Imaginary part of function values, shape (n_points,)
    ///
    /// Returns:
    ///     None (modifies model in-place)
    fn fit<'py>(
        &mut self,
        xi_data: &Bound<'py, PyArray1<f64>>,
        f_data_real: &Bound<'py, PyArray1<f64>>,
        f_data_imag: &Bound<'py, PyArray1<f64>>,
    ) -> PyResult<()> {
        // Convert numpy arrays
        let xi = xi_data.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read xi_data: {}",
                e
            ))
        })?;
        let xi_array = xi.as_array().to_owned();

        let f_real = f_data_real.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read f_data_real: {}",
                e
            ))
        })?;
        let f_real_array = f_real.as_array();

        let f_imag = f_data_imag.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read f_data_imag: {}",
                e
            ))
        })?;
        let f_imag_array = f_imag.as_array();

        // Combine into complex array
        let f_complex: Array1<Complex64> = f_real_array
            .iter()
            .zip(f_imag_array.iter())
            .map(|(&re, &im)| Complex64::new(re, im))
            .collect();

        // Perform fitting
        self.inner.fit(&xi_array, &f_complex).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Fitting failed: {:?}", e))
        })?;

        Ok(())
    }

    /// Evaluate the model at complex frequencies
    ///
    /// Args:
    ///     z_real: Real parts of evaluation points, shape (n_eval,)
    ///     z_imag: Imaginary parts of evaluation points, shape (n_eval,)
    ///
    /// Returns:
    ///     Tuple of (real, imag) arrays with function values
    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        z_real: &Bound<'py, PyArray1<f64>>,
        z_imag: &Bound<'py, PyArray1<f64>>,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
        // Convert input arrays
        let z_re = z_real.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read z_real: {}",
                e
            ))
        })?;
        let z_re_array = z_re.as_array();

        let z_im = z_imag.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read z_imag: {}",
                e
            ))
        })?;
        let z_im_array = z_im.as_array();

        // Create complex array
        let z_complex: Array1<Complex64> = z_re_array
            .iter()
            .zip(z_im_array.iter())
            .map(|(&re, &im)| Complex64::new(re, im))
            .collect();

        // Evaluate using batch method
        let result = self.inner.evaluate_batch(&z_complex);

        // Split into real and imaginary parts
        let result_real: Array1<f64> = result.mapv(|c| c.re);
        let result_imag: Array1<f64> = result.mapv(|c| c.im);

        Ok((
            PyArray1::from_owned_array(py, result_real),
            PyArray1::from_owned_array(py, result_imag),
        ))
    }

    /// Get cross-validation error
    fn cross_validation_error(&self) -> f64 {
        self.inner.cross_validation_error()
    }

    /// Check stability of the fitted model
    fn stability_check(&self) -> bool {
        self.inner.stability_check()
    }

    /// Get poles and residues
    ///
    /// Returns:
    ///     Dictionary with 'poles_real', 'poles_imag', 'residues_real', 'residues_imag'
    fn get_poles_residues<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        if let Some((poles, residues)) = self.inner.get_poles_residues() {
            let poles_real: Array1<f64> = poles.mapv(|c| c.re);
            let poles_imag: Array1<f64> = poles.mapv(|c| c.im);
            let residues_real: Array1<f64> = residues.mapv(|c| c.re);
            let residues_imag: Array1<f64> = residues.mapv(|c| c.im);

            dict.set_item("poles_real", PyArray1::from_owned_array(py, poles_real))?;
            dict.set_item("poles_imag", PyArray1::from_owned_array(py, poles_imag))?;
            dict.set_item(
                "residues_real",
                PyArray1::from_owned_array(py, residues_real),
            )?;
            dict.set_item(
                "residues_imag",
                PyArray1::from_owned_array(py, residues_imag),
            )?;
        }

        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "PadeModel(cv_error={:.3e})",
            self.inner.cross_validation_error()
        )
    }
}

/// Main analytic continuation fitter with automatic model selection
///
/// This class performs analytic continuation with cross-validation to
/// automatically select the best model (multipole or Padé).
///
/// Args:
///     config (ACConfig): Configuration object (optional)
#[pyclass(name = "AnalyticContinuationFitter")]
pub struct PyAnalyticContinuationFitter {
    inner: AnalyticContinuationFitter,
}

#[pymethods]
impl PyAnalyticContinuationFitter {
    #[new]
    #[pyo3(signature = (config=None))]
    fn new(config: Option<PyACConfig>) -> Self {
        let cfg = config.map(|c| c.inner).unwrap_or_default();
        Self {
            inner: AnalyticContinuationFitter::new(cfg),
        }
    }

    /// Fit models to imaginary-axis data with automatic model selection
    ///
    /// Args:
    ///     xi_data: Imaginary frequencies, shape (n_points,)
    ///     f_data_real: Real part of function values, shape (n_points,) or (n_points, n_orbitals)
    ///     f_data_imag: Imaginary part of function values, shape (n_points,) or (n_points, n_orbitals)
    ///
    /// Returns:
    ///     Dictionary with fitting results
    fn fit<'py>(
        &mut self,
        py: Python<'py>,
        xi_data: &Bound<'py, PyArray1<f64>>,
        f_data_real: &Bound<'py, PyAny>,
        f_data_imag: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyDict>> {
        // Convert xi data
        let xi = xi_data.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read xi_data: {}",
                e
            ))
        })?;
        let xi_array = xi.as_array().to_owned();

        // Check if f_data is 1D or 2D
        let result_dict = PyDict::new(py);

        // Try to interpret as 1D array first
        if let (Ok(f_real_1d), Ok(f_imag_1d)) = (
            f_data_real.downcast::<PyArray1<f64>>(),
            f_data_imag.downcast::<PyArray1<f64>>(),
        ) {
            // Single orbital case
            let f_real = f_real_1d.try_readonly().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to read f_data_real: {}",
                    e
                ))
            })?;
            let f_real_array = f_real.as_array();

            let f_imag = f_imag_1d.try_readonly().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to read f_data_imag: {}",
                    e
                ))
            })?;
            let f_imag_array = f_imag.as_array();

            // Combine into complex array
            let f_complex: Array1<Complex64> = f_real_array
                .iter()
                .zip(f_imag_array.iter())
                .map(|(&re, &im)| Complex64::new(re, im))
                .collect();

            // Perform fitting
            self.inner.fit(&xi_array, &f_complex).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Fitting failed: {:?}",
                    e
                ))
            })?;

            // Get best model info
            if let Some(model) = self.inner.best_model() {
                result_dict.set_item("has_model", true)?;
                result_dict.set_item("cv_error", model.cross_validation_error())?;
                result_dict.set_item("is_stable", model.stability_check())?;
            } else {
                result_dict.set_item("has_model", false)?;
            }
        } else if let (Ok(f_real_2d), Ok(f_imag_2d)) = (
            f_data_real.downcast::<PyArray2<f64>>(),
            f_data_imag.downcast::<PyArray2<f64>>(),
        ) {
            // Multiple orbitals case
            let f_real = f_real_2d.try_readonly().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to read f_data_real: {}",
                    e
                ))
            })?;
            let f_real_array = f_real.as_array();

            let f_imag = f_imag_2d.try_readonly().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to read f_data_imag: {}",
                    e
                ))
            })?;
            let f_imag_array = f_imag.as_array();

            let n_orbitals = f_real_array.ncols();
            let results_list = PyList::empty(py);

            // Fit each orbital separately
            for orb in 0..n_orbitals {
                let f_real_orb = f_real_array.column(orb);
                let f_imag_orb = f_imag_array.column(orb);

                // Combine into complex array
                let f_complex: Array1<Complex64> = f_real_orb
                    .iter()
                    .zip(f_imag_orb.iter())
                    .map(|(&re, &im)| Complex64::new(re, im))
                    .collect();

                // Create new fitter for each orbital with same config
                let config = ACConfig::default(); // Use default config for each orbital
                let mut fitter = AnalyticContinuationFitter::new(config);

                if fitter.fit(&xi_array, &f_complex).is_ok() {
                    if let Some(model) = fitter.best_model() {
                        let orb_dict = PyDict::new(py);
                        orb_dict.set_item("orbital", orb)?;
                        orb_dict.set_item("cv_error", model.cross_validation_error())?;
                        orb_dict.set_item("is_stable", model.stability_check())?;
                        results_list.append(orb_dict)?;
                    }
                }
            }

            result_dict.set_item("orbital_results", results_list)?;
        } else {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "f_data must be 1D or 2D numpy array",
            ));
        }

        Ok(result_dict)
    }

    /// Evaluate the fitted model on real frequency axis
    ///
    /// Args:
    ///     omega: Real frequencies, shape (n_freq,)
    ///     eta: Broadening parameter (optional)
    ///
    /// Returns:
    ///     Tuple of (real, imag) arrays with Σ(ω)
    #[pyo3(signature = (omega, eta=None))]
    fn evaluate_real_axis<'py>(
        &self,
        py: Python<'py>,
        omega: &Bound<'py, PyArray1<f64>>,
        eta: Option<f64>,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
        let omega_arr = omega.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read omega: {}",
                e
            ))
        })?;
        let omega_array = omega_arr.as_array().to_owned();

        let result = self
            .inner
            .evaluate_real_axis(&omega_array, eta)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Evaluation failed: {:?}",
                    e
                ))
            })?;

        // Split into real and imaginary parts
        let result_real: Array1<f64> = result.mapv(|c| c.re);
        let result_imag: Array1<f64> = result.mapv(|c| c.im);

        Ok((
            PyArray1::from_owned_array(py, result_real),
            PyArray1::from_owned_array(py, result_imag),
        ))
    }

    /// Compute spectral function A(ω) = -Im[Σ(ω)]/π
    ///
    /// Args:
    ///     omega: Real frequencies, shape (n_freq,)
    ///     eta: Broadening parameter (optional)
    ///
    /// Returns:
    ///     Spectral function array
    #[pyo3(signature = (omega, eta=None))]
    fn spectral_function<'py>(
        &self,
        py: Python<'py>,
        omega: &Bound<'py, PyArray1<f64>>,
        eta: Option<f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let omega_arr = omega.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read omega: {}",
                e
            ))
        })?;
        let omega_array = omega_arr.as_array().to_owned();

        let spectral = self
            .inner
            .spectral_function(&omega_array, eta)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to compute spectral function: {:?}",
                    e
                ))
            })?;

        Ok(PyArray1::from_owned_array(py, spectral))
    }

    /// Validate causality and Kramers-Kronig relations
    ///
    /// Args:
    ///     omega: Real frequencies for validation, shape (n_freq,)
    ///
    /// Returns:
    ///     Dictionary with causality metrics
    fn validate_causality<'py>(
        &self,
        py: Python<'py>,
        omega: &Bound<'py, PyArray1<f64>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let omega_arr = omega.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to read omega: {}",
                e
            ))
        })?;
        let omega_array = omega_arr.as_array().to_owned();

        let metrics = self.inner.validate_causality(&omega_array).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Validation failed: {:?}", e))
        })?;

        let dict = PyDict::new(py);
        dict.set_item("violations", metrics.violations)?;
        dict.set_item("violation_fraction", metrics.violation_fraction)?;
        dict.set_item("kk_error", metrics.kk_error)?;
        dict.set_item("kk_relative_error", metrics.kk_relative_error)?;

        Ok(dict)
    }

    /// Get information about the best fitted model
    ///
    /// Returns:
    ///     Dictionary with model information
    fn get_best_model_info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);

        if let Some(model) = self.inner.best_model() {
            dict.set_item("has_model", true)?;
            dict.set_item("cv_error", model.cross_validation_error())?;
            dict.set_item("is_stable", model.stability_check())?;

            if let Some((poles, residues)) = model.get_poles_residues() {
                dict.set_item("n_poles", poles.len())?;
                dict.set_item("n_residues", residues.len())?;
            }
        } else {
            dict.set_item("has_model", false)?;
        }

        Ok(dict)
    }

    fn __repr__(&self) -> String {
        if let Some(model) = self.inner.best_model() {
            format!(
                "AnalyticContinuationFitter(fitted=true, cv_error={:.3e})",
                model.cross_validation_error()
            )
        } else {
            "AnalyticContinuationFitter(fitted=false)".to_string()
        }
    }
}

/// Register the analytic continuation module with PyO3
pub fn register_analytic_continuation(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Add classes
    m.add_class::<PyACConfig>()?;
    m.add_class::<PyMultipoleModel>()?;
    m.add_class::<PyPadeModel>()?;
    m.add_class::<PyAnalyticContinuationFitter>()?;

    Ok(())
}
