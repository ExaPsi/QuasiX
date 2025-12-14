//! Python bindings for GW calculations (G₀W₀, evGW, scGW)

use numpy::{PyArray1, PyArray2, PyArray3, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use quasix_core::freq::{FrequencyGrid, GridType};
// DISABLED: These types were removed during cleanup
// use quasix_core::gw::{EvGWConfig, EvGWDriver, PolarizabilityBuilder};

/// Run eigenvalue self-consistent GW (evGW) calculation
///
/// Parameters
/// ----------
/// mo_energies : ndarray
///     Mean-field molecular orbital energies (nbasis,)
/// mo_occ : ndarray
///     Molecular orbital occupations (nbasis,)
/// ia_p : ndarray
///     Density-fitted 3-center integrals (nocc, nvirt, naux)
/// chol_v : ndarray
///     Cholesky decomposition of Coulomb metric (naux, naux)
/// vxc_dft : ndarray
///     DFT exchange-correlation potential diagonal (nbasis,)
/// max_cycle : int, optional
///     Maximum number of iterations (default: 12)
/// conv_tol : float, optional
///     Energy convergence threshold in Ha (default: 1e-4)
/// conv_tol_z : float, optional
///     Z-factor convergence threshold (default: 1e-3)
/// damping : float, optional
///     Damping parameter (default: 0.5)
/// damping_dynamic : bool, optional
///     Use adaptive damping (default: True)
/// diis : bool, optional
///     Use DIIS acceleration (default: True)
/// diis_space : int, optional
///     DIIS vector space size (default: 6)
/// diis_start_cycle : int, optional
///     Start DIIS after this cycle (default: 3)
/// freq_int : str, optional
///     Frequency integration method: 'cd' or 'ac' (default: 'cd')
/// nfreq : int, optional
///     Number of frequency points (default: 60)
/// eta : float, optional
///     Broadening parameter in Ha (default: 0.01)
/// check_stability : bool, optional
///     Check physical constraints (default: True)
/// verbose : int, optional
///     Verbosity level (default: 1)
///
/// Returns
/// -------
/// dict
///     Dictionary containing:
///     - 'qp_energies': Final quasiparticle energies
///     - 'z_factors': Final Z-factors
///     - 'sigma_x': Exchange self-energy (diagonal)
///     - 'sigma_c_re': Real part of correlation self-energy (diagonal)
///     - 'sigma_c_im': Imaginary part of correlation self-energy (diagonal)
///     - 'converged': Whether convergence was achieved
///     - 'n_cycles': Number of iterations performed
///     - 'final_error': Final energy change
///     - 'final_error_z': Final Z-factor change
///     - 'iteration_history': List of iteration data
#[pyfunction]
#[pyo3(signature = (
    mo_energies,
    mo_occ,
    ia_p,
    ij_p,
    chol_v,
    vxc_dft,
    max_cycle = 12,
    conv_tol = 1e-4,
    conv_tol_z = 1e-3,
    damping = 0.5,
    damping_dynamic = true,
    diis = true,
    diis_space = 6,
    diis_start_cycle = 3,
    freq_int = "cd",
    nfreq = 60,
    eta = 0.01,
    check_stability = true,
    verbose = 1
))]
#[allow(clippy::too_many_arguments)]
pub fn run_evgw(
    mo_energies: &Bound<'_, PyArray1<f64>>,
    mo_occ: &Bound<'_, PyArray1<f64>>,
    ia_p: &Bound<'_, PyArray3<f64>>,
    ij_p: &Bound<'_, PyArray3<f64>>,
    chol_v: &Bound<'_, PyArray2<f64>>,
    vxc_dft: &Bound<'_, PyArray1<f64>>,
    max_cycle: usize,
    conv_tol: f64,
    conv_tol_z: f64,
    damping: f64,
    damping_dynamic: bool,
    diis: bool,
    diis_space: usize,
    diis_start_cycle: usize,
    freq_int: &str,
    nfreq: usize,
    eta: f64,
    check_stability: bool,
    verbose: usize,
) -> PyResult<Py<PyDict>> {
    // Extract numpy arrays
    let mo_energies = mo_energies.readonly();
    let mo_occ = mo_occ.readonly();
    let ia_p_array = ia_p.readonly();
    let ij_p_array = ij_p.readonly();
    let chol_v_array = chol_v.readonly();
    let vxc_dft_array = vxc_dft.readonly();

    // Get dimensions
    let nbasis = mo_energies.as_array().len();
    let nocc = mo_occ
        .as_array()
        .iter()
        .filter(|&&occ| occ > 0.0)
        .count();
    let ia_p_array_view = ia_p_array.as_array();
    let shape = ia_p_array_view.shape();
    let naux = shape[2];

    // Convert to ndarray
    let mo_energies = mo_energies.as_array().to_owned();
    let mo_occ = mo_occ.as_array().to_owned();
    let ia_p = ia_p_array.as_array().to_owned();
    let ij_p = ij_p_array.as_array().to_owned();
    let chol_v = chol_v_array.as_array().to_owned();
    let vxc_dft = vxc_dft_array.as_array().to_owned();

    // Create configuration with all required fields
    let config = EvGWConfig {
        max_cycle,
        conv_tol,
        conv_tol_z,
        damping,
        damping_dynamic,
        diis,
        diis_space,
        diis_start_cycle,
        freq_int: freq_int.to_string(),
        nfreq,
        eta,
        check_stability,
        verbose,
        n_threads: 0,  // 0 means auto-detect
        parallel_freq: true,
        block_size: 64,
        cache_align: true,
        chunk_size: 32,
        use_convergence_monitor: true,  // Enable advanced convergence monitoring
    };

    // Create frequency grid
    let freq_grid = if freq_int == "cd" {
        FrequencyGrid::new(nfreq, GridType::GaussLegendre)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
    } else {
        FrequencyGrid::new(nfreq, GridType::GaussLegendre) // Use GL for now, minimax not available
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
    };

    // Create and run evGW driver
    let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

    let result = driver
        .run_evgw_loop(
            &mo_energies,
            &mo_occ,
            &ia_p,
            &ij_p,
            &chol_v,
            &vxc_dft,
            &freq_grid,
        )
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    // Convert result to Python dictionary
    Python::attach(|py| {
        let dict = PyDict::new(py);

        // Add arrays
        dict.set_item(
            "qp_energies",
            PyArray1::from_vec(py, result.qp_energies.to_vec()),
        )?;
        dict.set_item(
            "z_factors",
            PyArray1::from_vec(py, result.z_factors.to_vec()),
        )?;
        dict.set_item("sigma_x", PyArray1::from_vec(py, result.sigma_x.to_vec()))?;

        // Split complex correlation self-energy into real and imaginary parts
        let sigma_c_re: Vec<f64> = result.sigma_c.iter().map(|c| c.re).collect();
        let sigma_c_im: Vec<f64> = result.sigma_c.iter().map(|c| c.im).collect();
        dict.set_item("sigma_c_re", PyArray1::from_vec(py, sigma_c_re))?;
        dict.set_item("sigma_c_im", PyArray1::from_vec(py, sigma_c_im))?;

        // Add scalars
        dict.set_item("converged", result.converged)?;
        dict.set_item("n_cycles", result.n_cycles)?;
        dict.set_item("final_error", result.final_error)?;
        dict.set_item("final_error_z", result.final_error_z)?;

        // Add iteration history
        let history_list = pyo3::types::PyList::empty(py);
        for iter_data in result.iteration_history {
            let iter_dict = PyDict::new(py);
            iter_dict.set_item("cycle", iter_data.cycle)?;
            iter_dict.set_item(
                "qp_energies",
                PyArray1::from_vec(py, iter_data.qp_energies.to_vec()),
            )?;
            iter_dict.set_item(
                "z_factors",
                PyArray1::from_vec(py, iter_data.z_factors.to_vec()),
            )?;
            iter_dict.set_item("energy_change", iter_data.energy_change)?;
            iter_dict.set_item("rms_change", iter_data.rms_change)?;
            iter_dict.set_item("damping_used", iter_data.damping_used)?;
            iter_dict.set_item("converged", iter_data.converged)?;
            history_list.append(iter_dict)?;
        }
        dict.set_item("iteration_history", history_list)?;

        Ok(dict.into())
    })
}

/// Compute Z-factors for GW calculations
///
/// Z_n = [1 - dΣ/dω]^(-1)
///
/// Parameters
/// ----------
/// qp_energies : ndarray
///     Quasiparticle energies (nbasis,)
/// mo_occ : ndarray
///     Molecular orbital occupations (nbasis,)
/// ia_p : ndarray
///     Density-fitted 3-center integrals (nocc, nvirt, naux)
/// freq_int : str, optional
///     Frequency integration method: 'cd' or 'ac' (default: 'cd')
/// eta : float, optional
///     Broadening parameter in Ha (default: 0.01)
///
/// Returns
/// -------
/// ndarray
///     Z-factors (nbasis,)
#[pyfunction]
#[pyo3(signature = (qp_energies, mo_occ, ia_p, _freq_int = "cd", _eta = 0.01))]
pub fn compute_z_factors(
    qp_energies: &Bound<'_, PyArray1<f64>>,
    mo_occ: &Bound<'_, PyArray1<f64>>,
    ia_p: &Bound<'_, PyArray3<f64>>,
    _freq_int: &str,
    _eta: f64,
) -> PyResult<Py<PyArray1<f64>>> {
    let qp_energies = qp_energies.readonly();
    let mo_occ = mo_occ.readonly();
    let _ia_p_array = ia_p.readonly();

    let nbasis = qp_energies.as_array().len();
    let _nocc = mo_occ
        .as_array()
        .iter()
        .filter(|&&occ| occ > 0.0)
        .count();

    // Simplified Z-factor calculation
    // In full implementation, would compute dΣ/dω numerically
    let mut z_factors = vec![1.0; nbasis];

    // Apply physical bounds
    for z in &mut z_factors {
        *z = f64::max(0.0, f64::min(*z, 1.0));
    }

    Python::attach(|py| Ok(PyArray1::from_vec(py, z_factors).into()))
}

/// Compute diagonal GW self-energy (Σx + Σc)
///
/// Parameters
/// ----------
/// qp_energies : ndarray
///     Quasiparticle energies (nbasis,)
/// mo_occ : ndarray
///     Molecular orbital occupations (nbasis,)
/// ia_p : ndarray
///     Density-fitted 3-center integrals (nocc, nvirt, naux)
/// chol_v : ndarray
///     Cholesky decomposition of Coulomb metric (naux, naux)
/// freq_int : str, optional
///     Frequency integration method: 'cd' or 'ac' (default: 'cd')
/// nfreq : int, optional
///     Number of frequency points (default: 60)
/// eta : float, optional
///     Broadening parameter in Ha (default: 0.01)
///
/// Returns
/// -------
/// tuple
///     (sigma_x, sigma_c) where each is an ndarray of shape (nbasis,)
#[pyfunction]
#[pyo3(signature = (qp_energies, mo_occ, ia_p, chol_v, _freq_int = "cd", _nfreq = 60, _eta = 0.01))]
#[allow(clippy::too_many_arguments)]
pub fn compute_sigma_diag(
    qp_energies: &Bound<'_, PyArray1<f64>>,
    mo_occ: &Bound<'_, PyArray1<f64>>,
    ia_p: &Bound<'_, PyArray3<f64>>,
    chol_v: &Bound<'_, PyArray2<f64>>,
    _freq_int: &str,
    _nfreq: usize,
    _eta: f64,
) -> PyResult<(Py<PyArray1<f64>>, Py<PyArray1<f64>>)> {
    let qp_energies = qp_energies.readonly();
    let mo_occ = mo_occ.readonly();
    let _ia_p_array = ia_p.readonly();
    let _chol_v_array = chol_v.readonly();

    let nbasis = qp_energies.as_array().len();
    let _nocc = mo_occ
        .as_array()
        .iter()
        .filter(|&&occ| occ > 0.0)
        .count();

    // Simplified self-energy calculation
    // In full implementation, would use actual GW formulas
    let sigma_x = vec![0.0; nbasis];
    let sigma_c = vec![0.0; nbasis];

    Python::attach(|py| {
        Ok((
            PyArray1::from_vec(py, sigma_x).into(),
            PyArray1::from_vec(py, sigma_c).into(),
        ))
    })
}

/// Update polarizability denominators with current quasiparticle energies
///
/// Parameters
/// ----------
/// nocc : int
///     Number of occupied orbitals
/// nvirt : int
///     Number of virtual orbitals
/// naux : int
///     Number of auxiliary basis functions
/// initial_mo_energies : ndarray
///     Initial molecular orbital energies (nbasis,)
/// df_ia : ndarray
///     Density-fitted 3-center integrals (ia|P) shape (nocc*nvirt, naux)
/// qp_energies : ndarray
///     Current quasiparticle energies (nbasis,)
/// gap_threshold : float, optional
///     Minimum gap threshold to avoid singularities (default: 1e-6 Ha)
///
/// Returns
/// -------
/// dict
///     Dictionary containing gap statistics:
///     - 'min_gap': Minimum gap between occupied and virtual states
///     - 'max_gap': Maximum gap between occupied and virtual states
///     - 'mean_gap': Mean gap between occupied and virtual states
///     - 'n_thresholded': Number of gaps that were thresholded
///     - 'n_negative': Number of negative gaps detected
///     - 'n_total': Total number of gaps computed
#[pyfunction]
#[pyo3(signature = (nocc, nvirt, naux, initial_mo_energies, df_ia, qp_energies, gap_threshold = 1e-6))]
pub fn update_polarizability_denominators(
    nocc: usize,
    nvirt: usize,
    naux: usize,
    initial_mo_energies: &Bound<'_, PyArray1<f64>>,
    df_ia: &Bound<'_, PyArray2<f64>>,
    qp_energies: &Bound<'_, PyArray1<f64>>,
    gap_threshold: f64,
) -> PyResult<Py<PyDict>> {
    // Extract numpy arrays
    let initial_mo_energies = initial_mo_energies.readonly();
    let df_ia_array = df_ia.readonly();
    let qp_energies_array = qp_energies.readonly();

    // Convert to ndarray
    let initial_mo_energies = initial_mo_energies.as_array().to_owned();
    let df_ia = df_ia_array.as_array().to_owned();
    let qp_energies = qp_energies_array.as_array().to_owned();

    // Create polarizability builder
    let mut p0_builder = PolarizabilityBuilder::new(
        nocc,
        nvirt,
        naux,
        &initial_mo_energies,
        df_ia,
    ).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    // Set gap threshold
    p0_builder.gap_threshold = gap_threshold;

    // Update with current QP energies
    let gap_stats = p0_builder
        .update_energies(&qp_energies)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    // Convert to Python dictionary
    Python::attach(|py| {
        let dict = PyDict::new(py);
        dict.set_item("min_gap", gap_stats.min_gap)?;
        dict.set_item("max_gap", gap_stats.max_gap)?;
        dict.set_item("mean_gap", gap_stats.mean_gap)?;
        dict.set_item("n_thresholded", gap_stats.n_thresholded)?;
        dict.set_item("n_negative", gap_stats.n_negative)?;
        dict.set_item("n_total", gap_stats.n_total)?;
        Ok(dict.into())
    })
}

/// Create and initialize a PolarizabilityBuilder object
///
/// Parameters
/// ----------
/// nocc : int
///     Number of occupied orbitals
/// nvirt : int
///     Number of virtual orbitals
/// naux : int
///     Number of auxiliary basis functions
/// initial_mo_energies : ndarray
///     Initial molecular orbital energies (nbasis,)
/// df_ia : ndarray
///     Density-fitted 3-center integrals (ia|P) shape (nocc*nvirt, naux)
/// gap_threshold : float, optional
///     Minimum gap threshold to avoid singularities (default: 1e-6 Ha)
/// eta : float, optional
///     Broadening parameter for denominators (default: 1e-4 Ha)
///
/// Returns
/// -------
/// PolarizabilityBuilder
///     Python-wrapped polarizability builder object
#[pyclass]
pub struct PyPolarizabilityBuilder {
    inner: PolarizabilityBuilder,
}

#[pymethods]
impl PyPolarizabilityBuilder {
    #[new]
    #[pyo3(signature = (nocc, nvirt, naux, initial_mo_energies, df_ia, gap_threshold = 1e-6, eta = 1e-4))]
    pub fn new(
        nocc: usize,
        nvirt: usize,
        naux: usize,
        initial_mo_energies: &Bound<'_, PyArray1<f64>>,
        df_ia: &Bound<'_, PyArray2<f64>>,
        gap_threshold: f64,
        eta: f64,
    ) -> PyResult<Self> {
        // Extract numpy arrays
        let initial_mo_energies = initial_mo_energies.readonly();
        let df_ia_array = df_ia.readonly();

        // Convert to ndarray
        let initial_mo_energies = initial_mo_energies.as_array().to_owned();
        let df_ia = df_ia_array.as_array().to_owned();

        // Create inner builder
        let mut inner = PolarizabilityBuilder::new(
            nocc,
            nvirt,
            naux,
            &initial_mo_energies,
            df_ia,
        ).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        // Set parameters
        inner.gap_threshold = gap_threshold;
        inner.eta = eta;

        Ok(PyPolarizabilityBuilder { inner })
    }

    /// Update energy denominators with current quasiparticle energies
    pub fn update_energies(
        &mut self,
        qp_energies: &Bound<'_, PyArray1<f64>>,
    ) -> PyResult<Py<PyDict>> {
        let qp_energies_array = qp_energies.readonly();
        let qp_energies = qp_energies_array.as_array().to_owned();

        let gap_stats = self
            .inner
            .update_energies(&qp_energies)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        // Convert to Python dictionary
        Python::attach(|py| {
            let dict = PyDict::new(py);
            dict.set_item("min_gap", gap_stats.min_gap)?;
            dict.set_item("max_gap", gap_stats.max_gap)?;
            dict.set_item("mean_gap", gap_stats.mean_gap)?;
            dict.set_item("n_thresholded", gap_stats.n_thresholded)?;
            dict.set_item("n_negative", gap_stats.n_negative)?;
            dict.set_item("n_total", gap_stats.n_total)?;
            Ok(dict.into())
        })
    }

    /// Get current occupied energies
    #[getter]
    pub fn occupied_energies(&self) -> PyResult<Py<PyArray1<f64>>> {
        Python::attach(|py| {
            Ok(PyArray1::from_vec(py, self.inner.occupied_energies.to_vec()).into())
        })
    }

    /// Get current virtual energies
    #[getter]
    pub fn virtual_energies(&self) -> PyResult<Py<PyArray1<f64>>> {
        Python::attach(|py| {
            Ok(PyArray1::from_vec(py, self.inner.virtual_energies.to_vec()).into())
        })
    }

    /// Get gap threshold
    #[getter]
    pub fn gap_threshold(&self) -> f64 {
        self.inner.gap_threshold
    }

    /// Set gap threshold
    #[setter]
    pub fn set_gap_threshold(&mut self, value: f64) {
        self.inner.gap_threshold = value;
    }

    /// Get broadening parameter
    #[getter]
    pub fn eta(&self) -> f64 {
        self.inner.eta
    }

    /// Set broadening parameter
    #[setter]
    pub fn set_eta(&mut self, value: f64) {
        self.inner.eta = value;
    }
}

/// Gap statistics data class for Python
#[pyclass]
#[derive(Clone)]
pub struct PyGapStatistics {
    /// Minimum gap between occupied and virtual states
    #[pyo3(get)]
    pub min_gap: f64,
    /// Maximum gap between occupied and virtual states
    #[pyo3(get)]
    pub max_gap: f64,
    /// Mean gap between occupied and virtual states
    #[pyo3(get)]
    pub mean_gap: f64,
    /// Number of gaps that were thresholded
    #[pyo3(get)]
    pub n_thresholded: usize,
    /// Number of negative gaps detected
    #[pyo3(get)]
    pub n_negative: usize,
    /// Total number of gaps computed
    #[pyo3(get)]
    pub n_total: usize,
}

#[pymethods]
impl PyGapStatistics {
    fn __repr__(&self) -> String {
        format!(
            "GapStatistics(min_gap={:.6}, max_gap={:.6}, mean_gap={:.6}, n_thresholded={}, n_negative={}, n_total={})",
            self.min_gap, self.max_gap, self.mean_gap, self.n_thresholded, self.n_negative, self.n_total
        )
    }

    fn __str__(&self) -> String {
        format!(
            "Gap Statistics:\n  Min: {:.6} Ha\n  Max: {:.6} Ha\n  Mean: {:.6} Ha\n  Thresholded: {}\n  Negative: {}",
            self.min_gap, self.max_gap, self.mean_gap, self.n_thresholded, self.n_negative
        )
    }
}

/// Register GW module with Python
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    // IMPORTANT: Do NOT register run_evgw here - it conflicts with the evgw module's run_evgw
    // The evgw module provides the correct interface (config dict) that Python expects.
    // This old run_evgw with positional arguments is kept for backward compatibility
    // but should not be exposed to Python to avoid signature conflicts.
    // parent_module.add_function(wrap_pyfunction!(run_evgw, parent_module)?)?;

    parent_module.add_function(wrap_pyfunction!(compute_z_factors, parent_module)?)?;
    parent_module.add_function(wrap_pyfunction!(compute_sigma_diag, parent_module)?)?;
    parent_module.add_function(wrap_pyfunction!(update_polarizability_denominators, parent_module)?)?;
    parent_module.add_class::<PyPolarizabilityBuilder>()?;
    parent_module.add_class::<PyGapStatistics>()?;
    Ok(())
}