//! BSE-TDA Python bindings
//!
//! This module provides Python bindings for the BSE-TDA kernel and Davidson solver.
//!
//! # Usage from Python
//!
//! ```python
//! import quasix
//! import numpy as np
//!
//! # Create BSE kernel
//! kernel = quasix.bse.BSETDAKernel(
//!     nocc=2, nvirt=3, naux=10,
//!     delta_qp=np.array([0.3, 0.4, 0.5, 0.35, 0.45, 0.55]),
//!     df_ia=np.random.rand(6, 10),
//!     df_ij=np.random.rand(4, 10),
//!     df_ab=np.random.rand(9, 10),
//!     w0=np.eye(10),
//!     spin="singlet"
//! )
//!
//! # Solve for lowest 3 eigenvalues
//! result = kernel.solve_davidson(n_roots=3)
//! print("Excitation energies:", result["eigenvalues"])
//! ```

use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use quasix_core::bse::{
    analyze_exciton, compute_absorption_spectrum, compute_amplitude_heatmap,
    compute_electron_density, compute_hole_density, compute_nto_decomposition,
    compute_oscillator_strengths, compute_participation_ratio, find_dominant_transition,
    verify_svd_roundtrip, BSEKernelConfig, BSETDAKernel, BroadeningType, DavidsonConfig,
    DavidsonResult, OpticalConfig, SpinType,
};

/// Hartree to eV conversion factor
const HARTREE_TO_EV: f64 = 27.2114;

/// Python wrapper for BSE-TDA kernel configuration
#[pyclass(name = "DavidsonConfig")]
#[derive(Clone)]
pub struct PyDavidsonConfig {
    inner: DavidsonConfig,
}

#[pymethods]
impl PyDavidsonConfig {
    #[new]
    #[pyo3(signature = (n_roots=3, max_iter=100, max_space=None, tol_residual=1e-5, level_shift=0.001))]
    fn new(
        n_roots: usize,
        max_iter: usize,
        max_space: Option<usize>,
        tol_residual: f64,
        level_shift: f64,
    ) -> Self {
        Self {
            inner: DavidsonConfig {
                n_roots,
                max_iter,
                max_space: max_space.unwrap_or(20 * n_roots),
                tol_residual,
                lindep_threshold: 1e-12,
                level_shift,
            },
        }
    }

    /// Create configuration with default settings for n_roots
    #[staticmethod]
    fn with_n_roots(n_roots: usize) -> Self {
        Self {
            inner: DavidsonConfig::with_n_roots(n_roots),
        }
    }

    /// Create configuration for tight convergence
    #[staticmethod]
    fn tight() -> Self {
        Self {
            inner: DavidsonConfig::tight(),
        }
    }

    #[getter]
    fn n_roots(&self) -> usize {
        self.inner.n_roots
    }

    #[getter]
    fn max_iter(&self) -> usize {
        self.inner.max_iter
    }

    #[getter]
    fn max_space(&self) -> usize {
        self.inner.max_space
    }

    #[getter]
    fn tol_residual(&self) -> f64 {
        self.inner.tol_residual
    }

    fn __repr__(&self) -> String {
        format!(
            "DavidsonConfig(n_roots={}, max_iter={}, max_space={}, tol_residual={:.1e})",
            self.inner.n_roots, self.inner.max_iter, self.inner.max_space, self.inner.tol_residual
        )
    }
}

/// Python wrapper for BSE-TDA kernel
///
/// Matrix-free BSE-TDA kernel for computing optical excitations.
#[pyclass(name = "BSETDAKernel")]
pub struct PyBSETDAKernel {
    inner: BSETDAKernel,
}

#[pymethods]
impl PyBSETDAKernel {
    /// Create a new BSE-TDA kernel
    ///
    /// Args:
    ///     nocc: Number of occupied orbitals
    ///     nvirt: Number of virtual orbitals
    ///     naux: Number of auxiliary basis functions
    ///     delta_qp: QP energy differences [nocc * nvirt]
    ///     df_ia: DF tensor (ia|P) [nocc*nvirt, naux]
    ///     df_ij: DF tensor (ij|P) [nocc*nocc, naux]
    ///     df_ab: DF tensor (ab|P) [nvirt*nvirt, naux]
    ///     w0: Static screened interaction [naux, naux]
    ///     spin: "singlet" or "triplet"
    #[new]
    #[pyo3(signature = (nocc, nvirt, naux, delta_qp, df_ia, df_ij, df_ab, w0, spin="singlet"))]
    fn new(
        nocc: usize,
        nvirt: usize,
        naux: usize,
        delta_qp: PyReadonlyArray1<f64>,
        df_ia: PyReadonlyArray2<f64>,
        df_ij: PyReadonlyArray2<f64>,
        df_ab: PyReadonlyArray2<f64>,
        w0: PyReadonlyArray2<f64>,
        spin: &str,
    ) -> PyResult<Self> {
        // Convert spin type
        let spin_type = match spin.to_lowercase().as_str() {
            "singlet" => SpinType::Singlet,
            "triplet" => SpinType::Triplet,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Unknown spin type: {}. Use 'singlet' or 'triplet'",
                    spin
                )))
            }
        };

        // Convert numpy arrays to ndarray
        let delta_qp = delta_qp.as_array().to_owned();
        let df_ia = df_ia.as_array().to_owned();
        let df_ij = df_ij.as_array().to_owned();
        let df_ab = df_ab.as_array().to_owned();
        let w0 = w0.as_array().to_owned();

        // Create config
        let config = BSEKernelConfig {
            spin_type,
            ..BSEKernelConfig::default()
        };

        // Create kernel
        let kernel =
            BSETDAKernel::new(nocc, nvirt, naux, delta_qp, df_ia, df_ij, df_ab, w0, config)
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                        "Failed to create BSE-TDA kernel: {}",
                        e
                    ))
                })?;

        Ok(Self { inner: kernel })
    }

    /// Get the BSE Hamiltonian dimension (nocc * nvirt)
    #[getter]
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    /// Get number of occupied orbitals
    #[getter]
    fn nocc(&self) -> usize {
        self.inner.nocc
    }

    /// Get number of virtual orbitals
    #[getter]
    fn nvirt(&self) -> usize {
        self.inner.nvirt
    }

    /// Get number of auxiliary functions
    #[getter]
    fn naux(&self) -> usize {
        self.inner.naux
    }

    /// Apply BSE-TDA Hamiltonian to a trial vector
    ///
    /// Args:
    ///     x: Trial vector [nocc * nvirt]
    ///
    /// Returns:
    ///     y = H @ x [nocc * nvirt]
    fn apply_hamiltonian<'py>(
        &self,
        py: Python<'py>,
        x: PyReadonlyArray1<f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let x = x.as_array().to_owned();

        let y = self.inner.apply_tda_hamiltonian(&x).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to apply Hamiltonian: {}",
                e
            ))
        })?;

        Ok(PyArray1::from_vec(py, y.to_vec()))
    }

    /// Apply BSE-TDA Hamiltonian to multiple trial vectors
    ///
    /// Args:
    ///     x_batch: Trial vectors [nocc * nvirt, n_vectors]
    ///
    /// Returns:
    ///     Y = H @ X [nocc * nvirt, n_vectors]
    fn apply_hamiltonian_batch<'py>(
        &self,
        py: Python<'py>,
        x_batch: PyReadonlyArray2<f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let x_batch = x_batch.as_array().to_owned();

        let y = self
            .inner
            .apply_tda_hamiltonian_batch(&x_batch)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to apply Hamiltonian batch: {}",
                    e
                ))
            })?;

        // Convert to numpy array
        let (n_trans, n_vecs) = y.dim();
        Ok(PyArray2::from_vec2(
            py,
            &(0..n_trans)
                .map(|i| (0..n_vecs).map(|j| y[[i, j]]).collect())
                .collect::<Vec<Vec<f64>>>(),
        )
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to convert result to numpy: {}",
                e
            ))
        })?)
    }

    /// Solve BSE-TDA eigenvalue problem using Davidson iteration
    ///
    /// Args:
    ///     n_roots: Number of lowest eigenvalues to compute (default: 3)
    ///     max_iter: Maximum iterations (default: 100)
    ///     tol_residual: Convergence tolerance (default: 1e-5)
    ///
    /// Returns:
    ///     Dictionary with keys:
    ///         - eigenvalues: Excitation energies [n_roots]
    ///         - eigenvectors: Exciton wavefunctions [n_trans, n_roots]
    ///         - converged: List of convergence flags
    ///         - iterations: Number of iterations performed
    ///         - residual_norms: Final residual norms
    ///         - all_converged: True if all roots converged
    #[pyo3(signature = (n_roots=3, max_iter=100, tol_residual=1e-5))]
    fn solve_davidson<'py>(
        &self,
        py: Python<'py>,
        n_roots: usize,
        max_iter: usize,
        tol_residual: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let config = DavidsonConfig {
            n_roots,
            max_iter,
            max_space: 20 * n_roots,
            tol_residual,
            ..DavidsonConfig::default()
        };

        let result = self.inner.solve_davidson(&config).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Davidson solver failed: {}",
                e
            ))
        })?;

        davidson_result_to_pydict(py, result)
    }

    /// Solve with PyDavidsonConfig
    fn solve_davidson_with_config<'py>(
        &self,
        py: Python<'py>,
        config: &PyDavidsonConfig,
    ) -> PyResult<Bound<'py, PyDict>> {
        let result = self.inner.solve_davidson(&config.inner).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Davidson solver failed: {}",
                e
            ))
        })?;

        davidson_result_to_pydict(py, result)
    }

    /// Convenience method: solve with default configuration
    ///
    /// Args:
    ///     n_roots: Number of eigenvalues to compute
    ///
    /// Returns:
    ///     Dictionary with eigenvalues, eigenvectors, and convergence info
    fn solve<'py>(&self, py: Python<'py>, n_roots: usize) -> PyResult<Bound<'py, PyDict>> {
        let result = self.inner.solve(n_roots).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Davidson solver failed: {}",
                e
            ))
        })?;

        davidson_result_to_pydict(py, result)
    }

    fn __repr__(&self) -> String {
        format!(
            "BSETDAKernel(nocc={}, nvirt={}, naux={}, dimension={})",
            self.inner.nocc,
            self.inner.nvirt,
            self.inner.naux,
            self.inner.dimension()
        )
    }
}

/// Convert DavidsonResult to Python dictionary
fn davidson_result_to_pydict<'py>(
    py: Python<'py>,
    result: DavidsonResult,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);

    // Eigenvalues
    dict.set_item(
        "eigenvalues",
        PyArray1::from_vec(py, result.eigenvalues.to_vec()),
    )?;

    // Eigenvectors
    let (n_trans, n_roots) = result.eigenvectors.dim();
    let eigvecs_2d: Vec<Vec<f64>> = (0..n_trans)
        .map(|i| (0..n_roots).map(|j| result.eigenvectors[[i, j]]).collect())
        .collect();
    dict.set_item("eigenvectors", PyArray2::from_vec2(py, &eigvecs_2d)?)?;

    // Convergence info
    dict.set_item("converged", result.converged.clone())?;
    dict.set_item("iterations", result.iterations)?;
    dict.set_item("residual_norms", result.residual_norms.clone())?;
    dict.set_item("all_converged", result.all_converged())?;
    dict.set_item("n_converged", result.n_converged())?;
    dict.set_item("max_residual", result.max_residual())?;

    Ok(dict)
}

/// Compute oscillator strengths from BSE-TDA eigenvalues and eigenvectors
///
/// # Arguments
///
/// * `eigenvalues` - BSE-TDA excitation energies [n_roots] in Hartree
/// * `eigenvectors` - BSE-TDA eigenvectors [n_trans, n_roots]
/// * `transition_dipoles` - Transition dipole integrals <i|r|a> [n_trans, 3] for x, y, z
///
/// # Returns
///
/// Dictionary containing:
/// * `oscillator_strengths` - [n_roots] array
/// * `transition_dipoles` - [n_roots, 3] array (contracted with eigenvectors)
/// * `f_sum` - Float, sum of oscillator strengths (for sum rule check)
///
/// # Formula
///
/// f_n = (2/3) * E_n * sum_alpha |d_n^alpha|^2
/// d_n^alpha = sum_ia X_n(ia) * mu_{ia}^alpha
#[pyfunction]
#[pyo3(name = "compute_oscillator_strengths")]
fn py_compute_oscillator_strengths<'py>(
    py: Python<'py>,
    eigenvalues: PyReadonlyArray1<f64>,
    eigenvectors: PyReadonlyArray2<f64>,
    transition_dipoles: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    // Convert numpy arrays to ndarray views
    let eigenvalues = eigenvalues.as_array();
    let eigenvectors = eigenvectors.as_array();
    let transition_dipoles = transition_dipoles.as_array();

    // Call Rust implementation
    let result = compute_oscillator_strengths(eigenvalues, eigenvectors, transition_dipoles)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute oscillator strengths: {}",
                e
            ))
        })?;

    // Build result dictionary
    let dict = PyDict::new(py);

    // Oscillator strengths [n_roots]
    dict.set_item(
        "oscillator_strengths",
        PyArray1::from_vec(py, result.oscillator_strengths.to_vec()),
    )?;

    // Transition dipoles [n_roots, 3]
    let (n_roots, n_comp) = result.transition_dipoles.dim();
    let dipoles_2d: Vec<Vec<f64>> = (0..n_roots)
        .map(|i| {
            (0..n_comp)
                .map(|j| result.transition_dipoles[[i, j]])
                .collect()
        })
        .collect();
    dict.set_item("transition_dipoles", PyArray2::from_vec2(py, &dipoles_2d)?)?;

    // Sum of oscillator strengths
    dict.set_item("f_sum", result.f_sum)?;

    Ok(dict)
}

/// Compute absorption spectrum with broadening
///
/// # Arguments
///
/// * `eigenvalues` - BSE-TDA excitation energies [n_roots] in Hartree
/// * `oscillator_strengths` - Oscillator strengths [n_roots]
/// * `energy_range` - Optional (min, max) tuple in Hartree. If None, auto-detect.
/// * `n_points` - Number of grid points (default: 1000)
/// * `broadening` - FWHM in eV (default: 0.1, converted to Hartree internally)
/// * `broadening_type` - "lorentzian" or "gaussian" (default: "lorentzian")
///
/// # Returns
///
/// Dictionary containing:
/// * `energies` - [n_points] array in Hartree
/// * `energies_eV` - [n_points] array in eV (convenience)
/// * `intensities` - [n_points] array
/// * `broadening_eV` - Float, broadening used in eV
/// * `broadening_type` - String, broadening type used
#[pyfunction]
#[pyo3(name = "compute_absorption_spectrum")]
#[pyo3(signature = (eigenvalues, oscillator_strengths, energy_range=None, n_points=1000, broadening=0.1, broadening_type="lorentzian"))]
fn py_compute_absorption_spectrum<'py>(
    py: Python<'py>,
    eigenvalues: PyReadonlyArray1<f64>,
    oscillator_strengths: PyReadonlyArray1<f64>,
    energy_range: Option<(f64, f64)>,
    n_points: usize,
    broadening: f64,
    broadening_type: &str,
) -> PyResult<Bound<'py, PyDict>> {
    // Convert numpy arrays to ndarray views
    let eigenvalues = eigenvalues.as_array();
    let oscillator_strengths = oscillator_strengths.as_array();

    // Convert broadening from eV to Hartree
    let broadening_ha = broadening / HARTREE_TO_EV;

    // Parse broadening type
    let br_type = match broadening_type.to_lowercase().as_str() {
        "lorentzian" => BroadeningType::Lorentzian,
        "gaussian" => BroadeningType::Gaussian,
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Unknown broadening type: {}. Use 'lorentzian' or 'gaussian'",
                broadening_type
            )))
        }
    };

    // Build config
    let config = OpticalConfig {
        broadening_type: br_type,
        broadening: broadening_ha,
        energy_range,
        n_points,
    };

    // Call Rust implementation
    let result =
        compute_absorption_spectrum(eigenvalues, oscillator_strengths, &config).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute absorption spectrum: {}",
                e
            ))
        })?;

    // Build result dictionary
    let dict = PyDict::new(py);

    // Energies in Hartree [n_points]
    let energies_vec = result.energies.to_vec();
    dict.set_item("energies", PyArray1::from_vec(py, energies_vec.clone()))?;

    // Energies in eV [n_points] (convenience)
    let energies_ev: Vec<f64> = energies_vec.iter().map(|&e| e * HARTREE_TO_EV).collect();
    dict.set_item("energies_eV", PyArray1::from_vec(py, energies_ev))?;

    // Intensities [n_points]
    dict.set_item(
        "intensities",
        PyArray1::from_vec(py, result.intensities.to_vec()),
    )?;

    // Broadening in eV (as provided by user)
    dict.set_item("broadening_eV", broadening)?;

    // Broadening type as string
    let br_str = match result.broadening_type {
        BroadeningType::Lorentzian => "lorentzian",
        BroadeningType::Gaussian => "gaussian",
    };
    dict.set_item("broadening_type", br_str)?;

    Ok(dict)
}

// =============================================================================
// Exciton Analysis Functions (S6-4)
// =============================================================================

/// Analyze exciton character from BSE-TDA eigenvectors
///
/// Performs complete exciton analysis including:
/// - Amplitude heatmap |X^S_{ia}|²
/// - NTO decomposition via SVD
/// - Participation ratio
/// - Hole and electron density matrices
/// - Dominant transition identification
///
/// # Arguments
///
/// * `eigenvalues` - BSE-TDA excitation energies [n_roots] in Hartree
/// * `eigenvectors` - BSE-TDA eigenvectors [n_trans, n_roots]
/// * `nocc` - Number of occupied orbitals
/// * `nvirt` - Number of virtual orbitals
/// * `state_index` - Which state to analyze (0-indexed)
///
/// # Returns
///
/// Dictionary containing:
/// * `state_index` - Int, the analyzed state index
/// * `excitation_energy` - Float, energy in Hartree
/// * `amplitude_squared` - [nocc, nvirt] array, |X^S_{ia}|²
/// * `nto_occupations` - [k] array, NTO occupation numbers (sum = 1)
/// * `participation_ratio` - Float, PR = 1 / Σλ_n²
/// * `hole_density` - [nocc, nocc] array, hole density matrix
/// * `electron_density` - [nvirt, nvirt] array, electron density matrix
/// * `dominant_transition` - Tuple (i, a, weight), dominant transition
/// * `nto_u` - [nocc, k] array, hole NTO coefficients
/// * `nto_singular_values` - [k] array, singular values
/// * `nto_vt` - [k, nvirt] array, particle NTO coefficients
#[pyfunction]
#[pyo3(name = "analyze_exciton")]
fn py_analyze_exciton<'py>(
    py: Python<'py>,
    eigenvalues: PyReadonlyArray1<f64>,
    eigenvectors: PyReadonlyArray2<f64>,
    nocc: usize,
    nvirt: usize,
    state_index: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let eigenvalues = eigenvalues.as_array();
    let eigenvectors = eigenvectors.as_array();

    let result =
        analyze_exciton(eigenvalues, eigenvectors, nocc, nvirt, state_index).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to analyze exciton: {}",
                e
            ))
        })?;

    // Build result dictionary
    let dict = PyDict::new(py);

    // Basic info
    dict.set_item("state_index", result.state_index)?;
    dict.set_item("excitation_energy", result.excitation_energy)?;
    dict.set_item(
        "excitation_energy_eV",
        result.excitation_energy * HARTREE_TO_EV,
    )?;
    dict.set_item("participation_ratio", result.participation_ratio)?;

    // Amplitude squared [nocc, nvirt]
    let amp_sq = array2_to_vec2(&result.amplitude_squared);
    dict.set_item("amplitude_squared", PyArray2::from_vec2(py, &amp_sq)?)?;

    // NTO occupations [k]
    dict.set_item(
        "nto_occupations",
        PyArray1::from_vec(py, result.nto_occupations.to_vec()),
    )?;

    // NTO singular values [k]
    dict.set_item(
        "nto_singular_values",
        PyArray1::from_vec(py, result.nto_singular_values.to_vec()),
    )?;

    // Hole density [nocc, nocc]
    let hole_dm = array2_to_vec2(&result.hole_density);
    dict.set_item("hole_density", PyArray2::from_vec2(py, &hole_dm)?)?;

    // Electron density [nvirt, nvirt]
    let elec_dm = array2_to_vec2(&result.electron_density);
    dict.set_item("electron_density", PyArray2::from_vec2(py, &elec_dm)?)?;

    // Dominant transition (i, a, weight)
    dict.set_item("dominant_transition", result.dominant_transition)?;

    // NTO matrices
    let nto_u = array2_to_vec2(&result.nto_u);
    dict.set_item("nto_u", PyArray2::from_vec2(py, &nto_u)?)?;

    let nto_vt = array2_to_vec2(&result.nto_vt);
    dict.set_item("nto_vt", PyArray2::from_vec2(py, &nto_vt)?)?;

    Ok(dict)
}

/// Compute amplitude heatmap |X^S_{ia}|² for visualization
///
/// # Arguments
///
/// * `eigenvector` - BSE eigenvector X^S, shape [n_trans]
/// * `nocc` - Number of occupied orbitals
/// * `nvirt` - Number of virtual orbitals
///
/// # Returns
///
/// [nocc, nvirt] array containing |X^S_{ia}|²
#[pyfunction]
#[pyo3(name = "compute_amplitude_heatmap")]
fn py_compute_amplitude_heatmap<'py>(
    py: Python<'py>,
    eigenvector: PyReadonlyArray1<f64>,
    nocc: usize,
    nvirt: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let eigenvector = eigenvector.as_array();

    let heatmap = compute_amplitude_heatmap(eigenvector, nocc, nvirt).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute amplitude heatmap: {}",
            e
        ))
    })?;

    let heatmap_vec = array2_to_vec2(&heatmap);
    Ok(PyArray2::from_vec2(py, &heatmap_vec)?)
}

/// Compute NTO decomposition via SVD
///
/// # Arguments
///
/// * `amplitude_matrix` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Dictionary containing:
/// * `u` - [nocc, k] hole NTO coefficients
/// * `singular_values` - [k] singular values
/// * `vt` - [k, nvirt] particle NTO coefficients
/// * `occupations` - [k] NTO occupations (σ², normalized to sum=1)
#[pyfunction]
#[pyo3(name = "compute_nto_decomposition")]
fn py_compute_nto_decomposition<'py>(
    py: Python<'py>,
    amplitude_matrix: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let amplitude = amplitude_matrix.as_array();

    let (u, s, vt) = compute_nto_decomposition(amplitude).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute NTO decomposition: {}",
            e
        ))
    })?;

    let dict = PyDict::new(py);

    // U [nocc, k]
    let u_vec = array2_to_vec2(&u);
    dict.set_item("u", PyArray2::from_vec2(py, &u_vec)?)?;

    // Singular values [k]
    dict.set_item("singular_values", PyArray1::from_vec(py, s.to_vec()))?;

    // Vt [k, nvirt]
    let vt_vec = array2_to_vec2(&vt);
    dict.set_item("vt", PyArray2::from_vec2(py, &vt_vec)?)?;

    // Occupations = s² normalized
    let occ_raw: Vec<f64> = s.iter().map(|x| x * x).collect();
    let occ_sum: f64 = occ_raw.iter().sum();
    let occupations: Vec<f64> = if occ_sum > 1e-14 {
        occ_raw.iter().map(|x| x / occ_sum).collect()
    } else {
        occ_raw
    };
    dict.set_item("occupations", PyArray1::from_vec(py, occupations))?;

    Ok(dict)
}

/// Compute participation ratio from NTO occupations
///
/// PR = 1 / Σλ_n² where λ_n are NTO occupations summing to 1.
///
/// # Arguments
///
/// * `nto_occupations` - NTO occupation numbers [k]
///
/// # Returns
///
/// Participation ratio value (PR = 1 for single config, PR > 1 for multiconfigurational)
#[pyfunction]
#[pyo3(name = "compute_participation_ratio")]
fn py_compute_participation_ratio(nto_occupations: PyReadonlyArray1<f64>) -> f64 {
    let occupations = nto_occupations.as_array();
    compute_participation_ratio(occupations)
}

/// Compute hole density matrix in MO basis
///
/// ρ^h_{ij} = Σ_a X^S_{ia} (X^S_{ja})*
///
/// # Arguments
///
/// * `amplitude_matrix` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Hole density matrix, shape [nocc, nocc]
#[pyfunction]
#[pyo3(name = "compute_hole_density")]
fn py_compute_hole_density<'py>(
    py: Python<'py>,
    amplitude_matrix: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let amplitude = amplitude_matrix.as_array();
    let hole_dm = compute_hole_density(amplitude);

    let hole_vec = array2_to_vec2(&hole_dm);
    Ok(PyArray2::from_vec2(py, &hole_vec)?)
}

/// Compute electron density matrix in MO basis
///
/// ρ^e_{ab} = Σ_i (X^S_{ia})* X^S_{ib}
///
/// # Arguments
///
/// * `amplitude_matrix` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Electron density matrix, shape [nvirt, nvirt]
#[pyfunction]
#[pyo3(name = "compute_electron_density")]
fn py_compute_electron_density<'py>(
    py: Python<'py>,
    amplitude_matrix: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let amplitude = amplitude_matrix.as_array();
    let elec_dm = compute_electron_density(amplitude);

    let elec_vec = array2_to_vec2(&elec_dm);
    Ok(PyArray2::from_vec2(py, &elec_vec)?)
}

/// Find dominant i→a transition
///
/// # Arguments
///
/// * `amplitude_squared` - Squared amplitudes |X^S_{ia}|², shape [nocc, nvirt]
///
/// # Returns
///
/// Tuple (i, a, weight) for the dominant transition
#[pyfunction]
#[pyo3(name = "find_dominant_transition")]
fn py_find_dominant_transition(amplitude_squared: PyReadonlyArray2<f64>) -> (usize, usize, f64) {
    let amp_sq = amplitude_squared.as_array();
    find_dominant_transition(amp_sq)
}

/// Verify SVD roundtrip reconstruction
///
/// Checks that U @ diag(S) @ Vt reconstructs the original matrix.
///
/// # Arguments
///
/// * `u` - Left singular vectors, shape [m, k]
/// * `s` - Singular values, shape [k]
/// * `vt` - Right singular vectors (transposed), shape [k, n]
/// * `original` - Original matrix to compare against, shape [m, n]
/// * `tolerance` - Maximum allowed reconstruction error
///
/// # Returns
///
/// True if reconstruction error is within tolerance
#[pyfunction]
#[pyo3(name = "verify_svd_roundtrip")]
fn py_verify_svd_roundtrip(
    u: PyReadonlyArray2<f64>,
    s: PyReadonlyArray1<f64>,
    vt: PyReadonlyArray2<f64>,
    original: PyReadonlyArray2<f64>,
    tolerance: f64,
) -> bool {
    let u = u.as_array();
    let s = s.as_array();
    let vt = vt.as_array();
    let original = original.as_array();

    verify_svd_roundtrip(u, s, vt, original, tolerance)
}

/// Helper function to convert Array2 to Vec<Vec<f64>>
fn array2_to_vec2(arr: &ndarray::Array2<f64>) -> Vec<Vec<f64>> {
    let (nrows, ncols) = arr.dim();
    (0..nrows)
        .map(|i| (0..ncols).map(|j| arr[[i, j]]).collect())
        .collect()
}

// ==============================================================================
// Convenience functions matching test API expectations
// ==============================================================================

/// Compute NTOs directly from amplitude matrix (test-compatible API)
///
/// This is a convenience wrapper that returns an object-like dictionary
/// with the same API as expected by the tests.
///
/// # Arguments
///
/// * `amplitude` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Dictionary with attributes: occupations, participation_ratio, nto_occupations,
/// u, singular_values, vt
#[pyfunction]
#[pyo3(name = "compute_ntos")]
fn py_compute_ntos<'py>(
    py: Python<'py>,
    amplitude: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let amplitude_arr = amplitude.as_array();

    // Compute NTO decomposition
    let (u, s, vt) = compute_nto_decomposition(amplitude_arr).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute NTO decomposition: {}",
            e
        ))
    })?;

    // Compute occupations (s² normalized)
    let occ_raw: Vec<f64> = s.iter().map(|x| x * x).collect();
    let occ_sum: f64 = occ_raw.iter().sum();
    let occupations: Vec<f64> = if occ_sum > 1e-14 {
        occ_raw.iter().map(|x| x / occ_sum).collect()
    } else {
        occ_raw
    };

    // Compute participation ratio
    let occ_arr = ndarray::Array1::from(occupations.clone());
    let pr = compute_participation_ratio(occ_arr.view());

    // Build result dictionary
    let dict = PyDict::new(py);

    // Occupations (test expects both .occupations and .nto_occupations)
    dict.set_item("occupations", PyArray1::from_vec(py, occupations.clone()))?;
    dict.set_item("nto_occupations", PyArray1::from_vec(py, occupations))?;
    dict.set_item("participation_ratio", pr)?;

    // SVD components
    let u_vec = array2_to_vec2(&u);
    dict.set_item("u", PyArray2::from_vec2(py, &u_vec)?)?;
    dict.set_item("singular_values", PyArray1::from_vec(py, s.to_vec()))?;
    let vt_vec = array2_to_vec2(&vt);
    dict.set_item("vt", PyArray2::from_vec2(py, &vt_vec)?)?;

    Ok(dict)
}

/// Compute density matrices from amplitude matrix (test-compatible API)
///
/// # Arguments
///
/// * `amplitude` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Dictionary with attributes: hole_density, electron_density
#[pyfunction]
#[pyo3(name = "compute_density_matrices")]
fn py_compute_density_matrices<'py>(
    py: Python<'py>,
    amplitude: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let amplitude_arr = amplitude.as_array();

    let hole_dm = compute_hole_density(amplitude_arr);
    let elec_dm = compute_electron_density(amplitude_arr);

    let dict = PyDict::new(py);

    let hole_vec = array2_to_vec2(&hole_dm);
    dict.set_item("hole_density", PyArray2::from_vec2(py, &hole_vec)?)?;

    let elec_vec = array2_to_vec2(&elec_dm);
    dict.set_item("electron_density", PyArray2::from_vec2(py, &elec_vec)?)?;

    Ok(dict)
}

/// Analyze exciton from amplitude matrix with optional energy (test-compatible API)
///
/// This is a convenience wrapper that accepts a 2D amplitude matrix directly.
///
/// # Arguments
///
/// * `amplitude` - Amplitude matrix X^S, shape [nocc, nvirt]
/// * `energy` - Excitation energy in Hartree (optional, default 0.0)
///
/// # Returns
///
/// Dictionary with full exciton analysis
#[pyfunction]
#[pyo3(name = "analyze_exciton_simple", signature = (amplitude, energy=0.0))]
fn py_analyze_exciton_simple<'py>(
    py: Python<'py>,
    amplitude: PyReadonlyArray2<f64>,
    energy: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let amplitude_arr = amplitude.as_array();
    let (nocc, nvirt) = amplitude_arr.dim();

    // Make a contiguous copy and flatten for amplitude heatmap
    let amplitude_owned = amplitude_arr.to_owned();
    let flat_vec: Vec<f64> = amplitude_owned.iter().cloned().collect();
    let flat_arr = ndarray::Array1::from(flat_vec);

    // Compute amplitude squared
    let amp_squared = compute_amplitude_heatmap(flat_arr.view(), nocc, nvirt).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute amplitude heatmap: {}",
            e
        ))
    })?;

    // Compute NTOs
    let (u, s, vt) = compute_nto_decomposition(amplitude_arr).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute NTO decomposition: {}",
            e
        ))
    })?;

    // Compute occupations
    let occ_raw: Vec<f64> = s.iter().map(|x| x * x).collect();
    let occ_sum: f64 = occ_raw.iter().sum();
    let occupations: Vec<f64> = if occ_sum > 1e-14 {
        occ_raw.iter().map(|x| x / occ_sum).collect()
    } else {
        occ_raw
    };

    // Compute participation ratio
    let occ_arr = ndarray::Array1::from(occupations.clone());
    let pr = compute_participation_ratio(occ_arr.view());

    // Compute densities
    let hole_dm = compute_hole_density(amplitude_arr);
    let elec_dm = compute_electron_density(amplitude_arr);

    // Find dominant transition
    let (dom_i, dom_a, dom_weight) = find_dominant_transition(amp_squared.view());

    // Build result dictionary
    let dict = PyDict::new(py);

    dict.set_item("state_index", 0_usize)?;
    dict.set_item("excitation_energy", energy)?;
    dict.set_item("excitation_energy_eV", energy * HARTREE_TO_EV)?;
    dict.set_item("participation_ratio", pr)?;

    // For tests: dominant_i and dominant_a
    dict.set_item("dominant_i", dom_i)?;
    dict.set_item("dominant_a", dom_a)?;
    dict.set_item("dominant_weight", dom_weight)?;
    dict.set_item("dominant_transition", (dom_i, dom_a, dom_weight))?;

    // Amplitude squared
    let amp_sq_vec = array2_to_vec2(&amp_squared);
    dict.set_item("amplitude_squared", PyArray2::from_vec2(py, &amp_sq_vec)?)?;

    // NTO data
    dict.set_item("occupations", PyArray1::from_vec(py, occupations.clone()))?;
    dict.set_item("nto_occupations", PyArray1::from_vec(py, occupations))?;
    dict.set_item("nto_singular_values", PyArray1::from_vec(py, s.to_vec()))?;

    let u_vec = array2_to_vec2(&u);
    dict.set_item("nto_u", PyArray2::from_vec2(py, &u_vec)?)?;
    let vt_vec = array2_to_vec2(&vt);
    dict.set_item("nto_vt", PyArray2::from_vec2(py, &vt_vec)?)?;

    // Density matrices
    let hole_vec = array2_to_vec2(&hole_dm);
    dict.set_item("hole_density", PyArray2::from_vec2(py, &hole_vec)?)?;
    let elec_vec = array2_to_vec2(&elec_dm);
    dict.set_item("electron_density", PyArray2::from_vec2(py, &elec_vec)?)?;

    Ok(dict)
}

/// Register BSE submodule
pub fn register_bse_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let bse = PyModule::new(m.py(), "bse")?;

    // Add classes to submodule
    bse.add_class::<PyBSETDAKernel>()?;
    bse.add_class::<PyDavidsonConfig>()?;

    // Add optical property functions to submodule
    bse.add_function(wrap_pyfunction!(py_compute_oscillator_strengths, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_compute_absorption_spectrum, &bse)?)?;

    // Add exciton analysis functions to submodule (S6-4)
    bse.add_function(wrap_pyfunction!(py_analyze_exciton, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_compute_amplitude_heatmap, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_compute_nto_decomposition, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_compute_participation_ratio, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_compute_hole_density, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_compute_electron_density, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_find_dominant_transition, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_verify_svd_roundtrip, &bse)?)?;

    // Convenience functions (test-compatible API)
    bse.add_function(wrap_pyfunction!(py_compute_ntos, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_compute_density_matrices, &bse)?)?;
    bse.add_function(wrap_pyfunction!(py_analyze_exciton_simple, &bse)?)?;

    // Add submodule to parent module
    m.add_submodule(&bse)?;

    // Also register classes in sys.modules so they can be imported directly
    // This is needed for PyO3 submodules to work properly with Python's import system
    let sys = m.py().import("sys")?;
    let modules = sys.getattr("modules")?;
    modules.set_item("quasix.bse", &bse)?;

    Ok(())
}
