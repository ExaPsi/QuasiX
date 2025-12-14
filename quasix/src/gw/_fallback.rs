// PyO3 bindings for the Fallback Controller
// Provides Python interface to the high-performance Rust implementation

use pyo3::prelude::*;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use ndarray::{Array1, Array2, Axis};
use num_complex::Complex64;
use std::collections::HashMap;

// Import from quasix_core
use quasix_core::fallback::{
    FallbackController as RustController,
    FallbackConfig as RustConfig,
    FallbackDecision as RustDecision,
    QualityMetrics as RustMetrics,
    ACModel,
};

/// Python-exposed Fallback Configuration
#[pyclass(name = "FallbackConfig")]
#[derive(Clone)]
pub struct FallbackConfig {
    cv_error_threshold: f64,
    pole_stability_threshold: f64,
    causality_tolerance: f64,
    sum_rule_threshold: f64,
    enable_simd: bool,
    cache_size: usize,
    num_threads: usize,
}

#[pymethods]
impl FallbackConfig {
    #[new]
    #[pyo3(signature = (
        cv_error_threshold = 0.01,
        pole_stability_threshold = -1e-6,
        causality_tolerance = 1e-3,
        sum_rule_threshold = 0.05,
        enable_simd = true,
        cache_size = 1000,
        num_threads = None
    ))]
    fn new(
        cv_error_threshold: f64,
        pole_stability_threshold: f64,
        causality_tolerance: f64,
        sum_rule_threshold: f64,
        enable_simd: bool,
        cache_size: usize,
        num_threads: Option<usize>,
    ) -> PyResult<Self> {
        // Validate thresholds
        if cv_error_threshold <= 0.0 || cv_error_threshold >= 1.0 {
            return Err(PyValueError::new_err(
                "cv_error_threshold must be in (0, 1)"
            ));
        }

        if pole_stability_threshold >= 0.0 {
            return Err(PyValueError::new_err(
                "pole_stability_threshold must be negative"
            ));
        }

        Ok(Self {
            cv_error_threshold,
            pole_stability_threshold,
            causality_tolerance,
            sum_rule_threshold,
            enable_simd,
            cache_size,
            num_threads: num_threads.unwrap_or_else(num_cpus::get),
        })
    }

    fn to_rust(&self) -> RustConfig {
        RustConfig {
            cv_error_threshold: self.cv_error_threshold,
            pole_stability_threshold: self.pole_stability_threshold,
            causality_tolerance: self.causality_tolerance,
            sum_rule_threshold: self.sum_rule_threshold,
            enable_simd: self.enable_simd,
            cache_size: self.cache_size,
            num_threads: self.num_threads,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "FallbackConfig(cv_error={:.3}, pole_stability={:.1e}, threads={})",
            self.cv_error_threshold, self.pole_stability_threshold, self.num_threads
        )
    }
}

/// Python-exposed Quality Metrics
#[pyclass(name = "QualityMetrics")]
pub struct QualityMetrics {
    #[pyo3(get)]
    cv_error: f64,
    #[pyo3(get)]
    pole_stability: f64,
    #[pyo3(get)]
    causality_score: f64,
    #[pyo3(get)]
    sum_rule_error: f64,
    #[pyo3(get)]
    convergence_rate: f64,
    #[pyo3(get)]
    condition_number: f64,
    residue_magnitudes: Array1<f64>,
}

#[pymethods]
impl QualityMetrics {
    #[getter]
    fn residue_magnitudes<'py>(&self, py: Python<'py>) -> &'py PyArray1<f64> {
        PyArray1::from_slice(py, self.residue_magnitudes.as_slice().unwrap())
    }

    fn passes_thresholds(&self, config: &FallbackConfig) -> bool {
        self.cv_error <= config.cv_error_threshold &&
        self.pole_stability <= config.pole_stability_threshold &&
        self.causality_score <= config.causality_tolerance &&
        self.sum_rule_error <= config.sum_rule_threshold
    }

    fn to_dict(&self, py: Python) -> PyResult<PyObject> {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("cv_error", self.cv_error)?;
        dict.set_item("pole_stability", self.pole_stability)?;
        dict.set_item("causality_score", self.causality_score)?;
        dict.set_item("sum_rule_error", self.sum_rule_error)?;
        dict.set_item("convergence_rate", self.convergence_rate)?;
        dict.set_item("condition_number", self.condition_number)?;
        dict.set_item("residue_magnitudes",
            PyArray1::from_slice(py, self.residue_magnitudes.as_slice().unwrap()))?;
        Ok(dict.into())
    }

    fn __repr__(&self) -> String {
        format!(
            "QualityMetrics(cv_error={:.3}, pole_stability={:.1e}, causality={:.3})",
            self.cv_error, self.pole_stability, self.causality_score
        )
    }
}

/// Python-exposed Fallback Decision
#[pyclass(name = "FallbackDecision")]
pub struct FallbackDecision {
    #[pyo3(get)]
    use_fallback: bool,
    #[pyo3(get)]
    method_used: String,
    #[pyo3(get)]
    reason: Option<String>,
    #[pyo3(get)]
    confidence: f64,
    #[pyo3(get)]
    recommendation: String,
    metrics: Option<QualityMetrics>,
}

#[pymethods]
impl FallbackDecision {
    #[getter]
    fn metrics(&self) -> Option<QualityMetrics> {
        self.metrics.clone()
    }

    fn to_dict(&self, py: Python) -> PyResult<PyObject> {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("use_fallback", self.use_fallback)?;
        dict.set_item("method_used", &self.method_used)?;
        dict.set_item("reason", &self.reason)?;
        dict.set_item("confidence", self.confidence)?;
        dict.set_item("recommendation", &self.recommendation)?;
        if let Some(ref metrics) = self.metrics {
            dict.set_item("metrics", metrics.to_dict(py)?)?;
        }
        Ok(dict.into())
    }

    fn __repr__(&self) -> String {
        if self.use_fallback {
            format!("FallbackDecision(use_fallback=True, method={}, reason={:?})",
                self.method_used, self.reason)
        } else {
            format!("FallbackDecision(use_fallback=False, method={}, confidence={:.2})",
                self.method_used, self.confidence)
        }
    }
}

/// Main Python-exposed Fallback Controller
#[pyclass(name = "FallbackController")]
pub struct FallbackController {
    rust_controller: RustController,
    config: FallbackConfig,
}

#[pymethods]
impl FallbackController {
    #[new]
    fn new(config: FallbackConfig) -> PyResult<Self> {
        let rust_config = config.to_rust();
        let rust_controller = RustController::new(rust_config)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to create controller: {}", e)))?;

        Ok(Self {
            rust_controller,
            config,
        })
    }

    /// Assess AC model quality and make fallback decision
    #[pyo3(signature = (poles, residues, validation_data, frequencies))]
    fn assess_and_decide<'py>(
        &self,
        py: Python<'py>,
        poles: PyReadonlyArray1<Complex64>,
        residues: PyReadonlyArray1<Complex64>,
        validation_data: PyReadonlyArray2<Complex64>,
        frequencies: PyReadonlyArray1<f64>,
    ) -> PyResult<FallbackDecision> {
        // Convert inputs to Rust types
        let poles_vec: Vec<Complex64> = poles.as_array().to_vec();
        let residues_vec: Vec<Complex64> = residues.as_array().to_vec();
        let validation_array = validation_data.as_array().to_owned();
        let freq_vec: Vec<f64> = frequencies.as_array().to_vec();

        // Create AC model
        let ac_model = ACModel {
            poles: poles_vec,
            residues: residues_vec,
            n_states: poles.len(),
        };

        // Release GIL for computation
        let result = py.allow_threads(|| {
            self.rust_controller.assess_and_decide(
                &ac_model,
                &validation_array,
                &freq_vec,
            )
        });

        // Convert result to Python
        match result {
            Ok(rust_decision) => Ok(self.convert_decision(rust_decision)),
            Err(e) => Err(PyRuntimeError::new_err(format!("Assessment failed: {}", e)))
        }
    }

    /// Batch assessment for multiple states
    #[pyo3(signature = (models_poles, models_residues, validation_batch, frequencies))]
    fn assess_batch<'py>(
        &self,
        py: Python<'py>,
        models_poles: Vec<PyReadonlyArray1<Complex64>>,
        models_residues: Vec<PyReadonlyArray1<Complex64>>,
        validation_batch: Vec<PyReadonlyArray2<Complex64>>,
        frequencies: PyReadonlyArray1<f64>,
    ) -> PyResult<Vec<FallbackDecision>> {
        if models_poles.len() != models_residues.len() ||
           models_poles.len() != validation_batch.len() {
            return Err(PyValueError::new_err(
                "All input lists must have the same length"
            ));
        }

        let freq_vec: Vec<f64> = frequencies.as_array().to_vec();

        // Process each model
        let mut decisions = Vec::new();
        for (i, (poles, residues)) in models_poles.iter()
            .zip(models_residues.iter())
            .enumerate() {

            let validation = &validation_batch[i];

            // Create AC model
            let ac_model = ACModel {
                poles: poles.as_array().to_vec(),
                residues: residues.as_array().to_vec(),
                n_states: poles.len(),
            };

            let validation_array = validation.as_array().to_owned();

            // Assess with GIL released
            let result = py.allow_threads(|| {
                self.rust_controller.assess_and_decide(
                    &ac_model,
                    &validation_array,
                    &freq_vec,
                )
            });

            match result {
                Ok(rust_decision) => decisions.push(self.convert_decision(rust_decision)),
                Err(e) => return Err(PyRuntimeError::new_err(
                    format!("Assessment failed for model {}: {}", i, e)
                ))
            }
        }

        Ok(decisions)
    }

    /// Get statistics snapshot
    fn get_statistics(&self, py: Python) -> PyResult<PyObject> {
        let stats = self.rust_controller.get_statistics();
        let dict = pyo3::types::PyDict::new(py);

        dict.set_item("total_assessments", stats.total_assessments)?;
        dict.set_item("fallbacks_triggered", stats.fallbacks_triggered)?;
        dict.set_item("cache_hits", stats.cache_hits)?;
        dict.set_item("fallback_rate", stats.fallback_rate())?;

        // Add reason breakdown
        let reasons_dict = pyo3::types::PyDict::new(py);
        for (reason, count) in &stats.reasons_count {
            reasons_dict.set_item(reason, count)?;
        }
        dict.set_item("reasons_count", reasons_dict)?;

        Ok(dict.into())
    }

    /// Clear the decision cache
    fn clear_cache(&self) {
        self.rust_controller.clear_cache();
    }

    /// Update configuration (requires cache clear)
    fn update_config(&mut self, config: FallbackConfig) -> PyResult<()> {
        let rust_config = config.to_rust();
        self.rust_controller.update_config(rust_config)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to update config: {}", e)))?;
        self.config = config;
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("FallbackController(config={:?})", self.config)
    }
}

impl FallbackController {
    /// Convert Rust decision to Python decision
    fn convert_decision(&self, rust_decision: RustDecision) -> FallbackDecision {
        FallbackDecision {
            use_fallback: rust_decision.use_fallback,
            method_used: if rust_decision.use_fallback { "CD".to_string() } else { "AC".to_string() },
            reason: rust_decision.reason.map(|r| format!("{:?}", r)),
            confidence: rust_decision.confidence,
            recommendation: rust_decision.recommendation,
            metrics: rust_decision.metrics.map(|m| QualityMetrics {
                cv_error: m.cv_error,
                pole_stability: m.pole_stability,
                causality_score: m.causality_score,
                sum_rule_error: m.sum_rule_error,
                convergence_rate: m.convergence_rate,
                condition_number: m.condition_number,
                residue_magnitudes: m.residue_magnitudes,
            }),
        }
    }
}

/// Run AC GW calculation via Rust
///
/// # DEPRECATED - DO NOT USE IN PRODUCTION
///
/// This is a placeholder implementation that returns FAKE DATA.
/// For real analytic continuation GW calculations, use:
/// - `quasix.selfenergy.correlation.CorrelationSelfEnergyCD` with AC config
/// - Or use the full evGW driver which handles AC internally
///
/// This function exists only for API compatibility during refactoring.
#[pyfunction]
#[pyo3(signature = (
    mo_energy,
    mo_occ,
    iaP,
    chol_v,
    n_freq = 200,
    eta = 0.01,
    max_iter = 50,
    conv_tol = 1e-6
))]
#[deprecated(note = "Use proper AC implementation from selfenergy module")]
pub fn run_ac_gw<'py>(
    _py: Python<'py>,
    _mo_energy: PyReadonlyArray1<f64>,
    _mo_occ: PyReadonlyArray1<f64>,
    _iaP: PyReadonlyArray2<f64>,
    _chol_v: PyReadonlyArray2<f64>,
    _n_freq: usize,
    _eta: f64,
    _max_iter: usize,
    _conv_tol: f64,
) -> PyResult<HashMap<String, PyObject>> {
    Err(pyo3::exceptions::PyNotImplementedError::new_err(
        "run_ac_gw() is a deprecated placeholder. Use quasix_core selfenergy module with AC config instead."
    ))
}

/// Run CD GW calculation via Rust
///
/// # DEPRECATED - DO NOT USE IN PRODUCTION
///
/// This is a placeholder implementation that returns FAKE DATA.
/// For real contour deformation GW calculations, use:
/// - `quasix.selfenergy.correlation.CorrelationSelfEnergyCD`
/// - Or use the full evGW driver with freq_int="cd"
///
/// This function exists only for API compatibility during refactoring.
#[pyfunction]
#[pyo3(signature = (
    mo_energy,
    mo_occ,
    iaP,
    chol_v,
    n_freq = 60,
    eta = 0.01,
    max_iter = 50,
    conv_tol = 1e-6
))]
#[deprecated(note = "Use proper CD implementation from selfenergy module")]
pub fn run_cd_gw<'py>(
    _py: Python<'py>,
    _mo_energy: PyReadonlyArray1<f64>,
    _mo_occ: PyReadonlyArray1<f64>,
    _iaP: PyReadonlyArray2<f64>,
    _chol_v: PyReadonlyArray2<f64>,
    _n_freq: usize,
    _eta: f64,
    _max_iter: usize,
    _conv_tol: f64,
) -> PyResult<HashMap<String, PyObject>> {
    Err(pyo3::exceptions::PyNotImplementedError::new_err(
        "run_cd_gw() is a deprecated placeholder. Use quasix.gw.run_evgw() with freq_int='cd' instead."
    ))
}

/// Register the fallback module with Python
pub fn register_fallback_module(py: Python, m: &PyModule) -> PyResult<()> {
    let fallback_module = PyModule::new(py, "fallback")?;

    fallback_module.add_class::<FallbackConfig>()?;
    fallback_module.add_class::<QualityMetrics>()?;
    fallback_module.add_class::<FallbackDecision>()?;
    fallback_module.add_class::<FallbackController>()?;

    fallback_module.add_function(wrap_pyfunction!(run_ac_gw, fallback_module)?)?;
    fallback_module.add_function(wrap_pyfunction!(run_cd_gw, fallback_module)?)?;

    m.add_submodule(fallback_module)?;
    Ok(())
}