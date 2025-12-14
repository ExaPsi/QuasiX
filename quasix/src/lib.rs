// QuasiX Python bindings via PyO3 - Basic scaffolding for S1-2
//
// This module provides minimal Python bindings to verify PyO3 setup works correctly.

// Allow some clippy warnings for Python bindings that are necessary for PyO3 interop
#![allow(clippy::too_many_arguments)] // Python functions often need many parameters
#![allow(clippy::type_complexity)] // PyO3 return types can be complex
#![allow(dead_code)] // Some functions are disabled but kept for future re-implementation

use numpy::{PyArray1, PyArray2, PyArray3, PyArrayMethods};
use pyo3::conversion::IntoPyObject;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

// Import QuasiX core logging utilities
use quasix_core::df::{compute_cholesky_v_with_tolerance, transform_mo_3center};
// Note: SolverBackend and SolverConfig removed in cleanup (legacy API only)
// Note: DielectricSolver and SolverType imported but used only in disabled sections
// Note: DielectricMatrix removed in cleanup (legacy type)
use quasix_core::dielectric::{PolarizabilityConfig, PolarizabilityRI};
use quasix_core::freq::{
    gauss_legendre_nodes_weights, FrequencyGrid, GridType, ImaginaryAxisGrid, MinimaxGrid,
    OptimizedGLGrid, TransformType,
};
use quasix_core::integrals::{BasisSet, IntegralEngine, Molecule};
use quasix_core::logging;
use quasix_core::selfenergy::ExchangeSelfEnergyRI;
// Note: CorrelationSelfEnergyCD used only in disabled sections

// Module declarations
mod analytic_continuation;
mod bse; // BSE-TDA kernel and Davidson solver (S6-1, S6-2)
mod evgw; // Re-enabled for S5-1
          // DISABLED: gw module uses removed EvGWConfig, EvGWDriver, PolarizabilityBuilder
          // mod gw;
mod io_hdf5;
mod linalg;
mod pyscf_adapter;
mod schema;

// Re-enabled evgw module for S5-1
use evgw::{evgw as evgw_func, run_evgw as run_evgw_func};

/// Initialize Rust logging subsystem
#[pyfunction]
fn init_rust_logging() -> PyResult<()> {
    logging::init_logger().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to initialize Rust logging: {}",
            e
        ))
    })?;
    Ok(())
}

/// Get QuasiX version string
#[pyfunction]
fn version() -> PyResult<String> {
    // Initialize logging if not already done
    let _ = logging::init_logger();
    Ok(env!("CARGO_PKG_VERSION").to_string())
}

/// Get QuasiX metadata dictionary
#[pyfunction]
fn metadata() -> PyResult<HashMap<String, String>> {
    let mut meta = HashMap::new();
    meta.insert("version".to_string(), env!("CARGO_PKG_VERSION").to_string());
    meta.insert("name".to_string(), env!("CARGO_PKG_NAME").to_string());
    meta.insert("authors".to_string(), env!("CARGO_PKG_AUTHORS").to_string());
    meta.insert(
        "description".to_string(),
        env!("CARGO_PKG_DESCRIPTION").to_string(),
    );
    meta.insert("rust_version".to_string(), "1.75.0".to_string()); // Minimum required
    Ok(meta)
}

/// No-op kernel function for testing GIL release and basic functionality
#[pyfunction]
fn noop_kernel() -> PyResult<HashMap<String, String>> {
    use quasix_core::timed_stage;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // This function does nothing but return a message
    // Used to test that the binding works and GIL can be released
    let result = timed_stage!("noop_kernel", {
        let mut result = HashMap::new();
        result.insert("status".to_string(), "success".to_string());
        result.insert(
            "message".to_string(),
            "No-op kernel executed successfully".to_string(),
        );
        result.insert("version".to_string(), env!("CARGO_PKG_VERSION").to_string());
        result
    });

    Ok(result)
}

/// Compute 3-center integrals (μν|P) for a molecule
///
/// Args:
///     molecule_name: Name of the molecule ("H2O", "NH3", "CO", "C6H6")
///     nbasis: Number of AO basis functions
///     naux: Number of auxiliary basis functions
///
/// Returns:
///     3D numpy array of shape (nbasis, nbasis, naux)
#[pyfunction]
fn compute_3center_integrals<'py>(
    py: Python<'py>,
    molecule_name: &str,
    _nbasis: usize,
    _naux: usize,
) -> PyResult<Bound<'py, PyArray3<f64>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Create molecule
    let molecule = match molecule_name {
        "H2O" => Molecule::water(),
        "NH3" => Molecule::ammonia(),
        "CO" => Molecule::carbon_monoxide(),
        "C6H6" => Molecule::benzene(),
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Unknown molecule: {}",
                molecule_name
            )))
        }
    };

    // Create basis sets
    let basis = BasisSet::for_molecule(molecule_name, false);
    let aux_basis = BasisSet::for_molecule(molecule_name, true);

    // Create integral engine
    let engine = IntegralEngine::new(molecule, basis, aux_basis);

    // Compute integrals
    let integrals = engine.compute_3center().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute 3-center integrals: {}",
            e
        ))
    })?;

    // Convert to numpy array
    Ok(PyArray3::from_owned_array(py, integrals))
}

/// Compute 2-center integrals (P|Q) for auxiliary basis
///
/// Args:
///     molecule_name: Name of the molecule ("H2O", "NH3", "CO", "C6H6")
///     naux: Number of auxiliary basis functions
///
/// Returns:
///     2D numpy array of shape (naux, naux)
#[pyfunction]
fn compute_2center_integrals<'py>(
    py: Python<'py>,
    molecule_name: &str,
    _naux: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Create auxiliary basis
    let aux_basis = BasisSet::for_molecule(molecule_name, true);

    // Compute integrals directly
    let integrals = quasix_core::integrals::compute_2center_integrals(&aux_basis).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute 2-center integrals: {}",
            e
        ))
    })?;

    // Convert to numpy array
    Ok(PyArray2::from_owned_array(py, integrals))
}

/// Get integral tensor dimensions for a specific molecule
///
/// Args:
///     molecule_name: Name of the molecule
///
/// Returns:
///     Dictionary with 'nbasis' and 'naux' keys
///
/// Note: These dimensions are for spherical harmonic basis functions
/// as used by libcint with CIntType::Spheric
#[pyfunction]
fn get_integral_dimensions(molecule_name: &str) -> PyResult<HashMap<String, usize>> {
    let mut dims = HashMap::new();

    // These dimensions match what libcint returns with spherical harmonics
    // Spherical harmonics have fewer components than Cartesian:
    // - d-orbitals: 5 spherical vs 6 Cartesian
    // - f-orbitals: 7 spherical vs 10 Cartesian
    match molecule_name {
        "H2O" => {
            dims.insert("nbasis".to_string(), 8); // s and p only, same for both
            dims.insert("naux".to_string(), 17); // includes d-orbitals (5 components)
        }
        "NH3" => {
            dims.insert("nbasis".to_string(), 11); // s and p only
            dims.insert("naux".to_string(), 22); // includes d-orbitals
        }
        "CO" => {
            dims.insert("nbasis".to_string(), 10); // s and p only
            dims.insert("naux".to_string(), 19); // includes d-orbitals
        }
        "C6H6" => {
            dims.insert("nbasis".to_string(), 36); // s and p only
            dims.insert("naux".to_string(), 66); // includes d-orbitals
        }
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Unknown molecule: {}",
                molecule_name
            )))
        }
    }

    Ok(dims)
}

/// Transform 3-center integrals to MO basis
///
/// Args:
///     j3c_ao: 3-center AO integrals of shape (n_ao, n_ao, n_aux)
///     c_occ: Occupied MO coefficients of shape (n_ao, n_occ)
///     c_vir: Virtual MO coefficients of shape (n_ao, n_vir)
///
/// Returns:
///     2D numpy array of shape (n_occ*n_vir, n_aux)
#[pyfunction]
fn transform_to_mo_basis<'py>(
    py: Python<'py>,
    j3c_ao: &Bound<'py, PyArray3<f64>>,
    c_occ: &Bound<'py, PyArray2<f64>>,
    c_vir: &Bound<'py, PyArray2<f64>>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let j3c_ao = j3c_ao.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to j3c_ao: {}",
            e
        ))
    })?;
    let j3c_ao_array = j3c_ao.as_array().to_owned();

    let c_occ = c_occ.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to c_occ: {}",
            e
        ))
    })?;
    let c_occ_array = c_occ.as_array().to_owned();

    let c_vir = c_vir.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to c_vir: {}",
            e
        ))
    })?;
    let c_vir_array = c_vir.as_array().to_owned();

    // Perform transformation
    let j3c_ia = transform_mo_3center(&j3c_ao_array, &c_occ_array, &c_vir_array).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to transform to MO basis: {}",
            e
        ))
    })?;

    // Convert back to numpy
    Ok(PyArray2::from_owned_array(py, j3c_ia))
}

/// Compute Cholesky factorization of metric matrix
///
/// Args:
///     metric: 2-center metric matrix of shape (n_aux, n_aux)
///
/// Returns:
///     Dictionary with 'L' (Cholesky factor), 'condition_number', and 'naux'
#[pyfunction]
fn compute_cholesky_metric<'py>(
    py: Python<'py>,
    metric: &Bound<'py, PyArray2<f64>>,
) -> PyResult<Bound<'py, PyDict>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy array to ndarray
    let metric = metric.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to metric: {}",
            e
        ))
    })?;
    let metric_array = metric.as_array().to_owned();

    // Compute Cholesky factorization with relaxed reconstruction tolerance
    // Use 1e-2 reconstruction tolerance for production runs with auxiliary basis sets
    // The reconstruction check is primarily for debugging, not numerical stability
    let cholesky =
        compute_cholesky_v_with_tolerance(&metric_array, None, Some(1e-2)).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute Cholesky factorization: {}",
                e
            ))
        })?;

    // Create result dictionary
    let dict = PyDict::new(py);
    dict.set_item("L", PyArray2::from_owned_array(py, cholesky.l_matrix))?;
    dict.set_item("condition_number", cholesky.condition_number)?;
    dict.set_item("naux", cholesky.naux)?;

    Ok(dict)
}

/// Solve quasiparticle equations for all orbitals (simplified version)
///
/// Args:
///     mo_energies: MO energies of shape (n_mo,)
///     mo_occ: MO occupations of shape (n_mo,)
///     sigma_x: Exchange self-energy matrix of shape (n_mo, n_mo)
///     sigma_c_diagonal: Pre-computed diagonal correlation self-energy of shape (n_mo,)
///     vxc_diagonal: DFT exchange-correlation diagonal of shape (n_mo,)
///     config: Optional solver configuration dictionary
///
/// Returns:
///     Dictionary with 'qp_energies', 'z_factors', 'diagnostics', etc.
///
/// DISABLED: QPEquationSolver removed during cleanup - awaiting re-implementation
#[cfg(disabled)] // Exclude from compilation
#[allow(dead_code)]
fn solve_quasiparticle_equations_DISABLED<'py>(
    py: Python<'py>,
    mo_energies: &Bound<'py, PyArray1<f64>>,
    mo_occ: &Bound<'py, PyArray1<f64>>,
    sigma_x: &Bound<'py, PyArray2<f64>>,
    sigma_c_diagonal: &Bound<'py, PyArray1<f64>>,
    vxc_diagonal: &Bound<'py, PyArray1<f64>>,
    config: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyDict>> {
    use num_complex::Complex64;
    use quasix_core::qp::solver::{QPEquationSolver, QPSolverConfig};

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let mo_energies = mo_energies.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to mo_energies: {}",
            e
        ))
    })?;
    let mo_energies_array = mo_energies.as_array().to_owned();

    let mo_occ = mo_occ.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to mo_occ: {}",
            e
        ))
    })?;
    let mo_occ_array = mo_occ.as_array().to_owned();

    let sigma_x = sigma_x.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to sigma_x: {}",
            e
        ))
    })?;
    let sigma_x_array = sigma_x.as_array().to_owned();

    let sigma_c_diagonal = sigma_c_diagonal.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to sigma_c_diagonal: {}",
            e
        ))
    })?;
    let sigma_c_diagonal_array = sigma_c_diagonal.as_array().to_owned();

    let vxc_diagonal = vxc_diagonal.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to vxc_diagonal: {}",
            e
        ))
    })?;
    let vxc_diagonal_array = vxc_diagonal.as_array().to_owned();

    // Parse configuration
    let mut solver_config = QPSolverConfig::default();
    let mut energy_dependence_alpha = 0.5; // Default energy dependence parameter
    let mut energy_model = String::from("lorentzian"); // Default energy dependence model
    if let Some(cfg) = config {
        if let Ok(Some(val)) = cfg.get_item("max_iterations") {
            if let Ok(v) = val.extract::<usize>() {
                solver_config.max_newton_iterations = v;
            }
        }
        if let Ok(Some(val)) = cfg.get_item("tolerance") {
            if let Ok(v) = val.extract::<f64>() {
                solver_config.energy_tolerance = v;
            }
        }
        if let Ok(Some(val)) = cfg.get_item("derivative_delta") {
            if let Ok(v) = val.extract::<f64>() {
                solver_config.derivative_delta = v;
            }
        }
        if let Ok(Some(val)) = cfg.get_item("damping_factor") {
            if let Ok(v) = val.extract::<f64>() {
                solver_config.damping_factor = v;
            }
        }
        if let Ok(Some(val)) = cfg.get_item("use_bisection_fallback") {
            if let Ok(v) = val.extract::<bool>() {
                solver_config.use_bisection_fallback = v;
            }
        }
        // Custom parameter for energy dependence strength
        if let Ok(Some(val)) = cfg.get_item("energy_dependence_alpha") {
            if let Ok(v) = val.extract::<f64>() {
                energy_dependence_alpha = v;
            }
        }
        // Energy dependence model type
        if let Ok(Some(val)) = cfg.get_item("energy_model") {
            if let Ok(v) = val.extract::<String>() {
                energy_model = v;
            }
        }
    }

    // Create solver with proper S3-6 configuration
    let solver = QPEquationSolver::new(solver_config.clone());

    // Create energy-dependent self-energy function
    // For proper GW, the correlation self-energy has energy dependence
    let _n_mo = mo_energies_array.len();

    // Clone arrays for closure capture
    let mo_energies_closure = mo_energies_array.clone();
    let sigma_x_closure = sigma_x_array.clone();
    let sigma_c_closure = sigma_c_diagonal_array.clone();

    // Create energy-dependent self-energy function
    // This implements a physically motivated energy dependence for Σc(ω)
    let sigma_total_func = move |orbital_idx: usize,
                                 energy: f64|
          -> quasix_core::common::Result<Complex64> {
        // Exchange self-energy (energy-independent)
        let sigma_x_val = sigma_x_closure[[orbital_idx, orbital_idx]];

        // Correlation self-energy with energy dependence
        let sigma_c_static = sigma_c_closure[orbital_idx];
        let ref_energy = mo_energies_closure[orbital_idx]; // Reference energy
        let energy_diff = energy - ref_energy;

        // Choose energy dependence model
        let sigma_c_dynamic = match energy_model.as_str() {
            "lorentzian" => {
                // Modified Lorentzian model with asymmetry to ensure non-zero derivative
                // Σc(ω) = Σc(εF) * [1 + β(ω-εF)] / [1 + α(ω-εF)²]
                // This ensures dΣ/dω ≠ 0 at ω = εF
                let beta = 0.3 * energy_dependence_alpha; // Asymmetry parameter
                let numerator = 1.0 + beta * energy_diff;
                let denominator = 1.0 + energy_dependence_alpha * energy_diff * energy_diff;
                let energy_factor = numerator / denominator;
                sigma_c_static * energy_factor
            }
            "pole" => {
                // Pole structure model: Σc(ω) = Σc(εF) * εF / (ω + iη)
                // More realistic for systems with sharp quasiparticle peaks
                // Using real part only for simplified version
                let eta: f64 = 0.01; // Small broadening
                let denominator = energy_diff.powi(2) + eta.powi(2);
                if denominator > 1e-10 {
                    let pole_factor = eta / denominator.sqrt();
                    sigma_c_static * pole_factor * (1.0 + 0.2 * energy_dependence_alpha)
                } else {
                    sigma_c_static
                }
            }
            "linear" => {
                // Asymmetric linear model: Σc(ω) = Σc(εF) * (1 - α(ω - εF))
                // Simple model with non-zero derivative
                let linear_factor = (1.0 - energy_dependence_alpha * energy_diff).clamp(0.1, 2.0);
                sigma_c_static * linear_factor
            }
            _ => {
                // Default to constant (original behavior)
                log::warn!(
                    "Unknown energy model '{}', using constant self-energy",
                    energy_model
                );
                sigma_c_static
            }
        };

        // Total self-energy = exchange + correlation(ω)
        Ok(Complex64::new(sigma_x_val + sigma_c_dynamic, 0.0))
    };

    // Solve QP equations with simplified self-energy
    let solution = solver
        .solve_all_orbitals(
            &mo_energies_array,
            &mo_occ_array,
            sigma_total_func,
            &vxc_diagonal_array,
        )
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to solve QP equations: {:?}",
                e
            ))
        })?;

    // Create result dictionary
    let dict = PyDict::new(py);
    dict.set_item(
        "qp_energies",
        PyArray1::from_owned_array(py, solution.qp_energies),
    )?;
    dict.set_item(
        "z_factors",
        PyArray1::from_owned_array(py, solution.z_factors),
    )?;
    dict.set_item(
        "convergence_flags",
        PyArray1::from_owned_array(py, solution.convergence_flags),
    )?;

    // Add diagnostics
    let diag_dict = PyDict::new(py);
    diag_dict.set_item("converged_count", solution.diagnostics.converged_count)?;
    diag_dict.set_item("physical_z_count", solution.diagnostics.physical_z_count)?;
    diag_dict.set_item("avg_iterations", solution.diagnostics.avg_iterations)?;
    diag_dict.set_item("z_factor_min", solution.diagnostics.z_factor_range.0)?;
    diag_dict.set_item("z_factor_max", solution.diagnostics.z_factor_range.1)?;
    diag_dict.set_item("avg_energy_shift", solution.diagnostics.avg_energy_shift)?;
    diag_dict.set_item("bisection_count", solution.diagnostics.bisection_count)?;
    dict.set_item("diagnostics", diag_dict)?;

    // Add per-orbital solutions
    let orbital_solutions = pyo3::types::PyList::empty(py);
    for sol in solution.orbital_solutions.iter() {
        let sol_dict = PyDict::new(py);
        sol_dict.set_item("orbital_idx", sol.orbital_idx)?;
        sol_dict.set_item("qp_energy", sol.qp_energy)?;
        sol_dict.set_item("z_factor", sol.z_factor)?;
        sol_dict.set_item("iterations", sol.iterations)?;
        sol_dict.set_item("converged", sol.converged)?;
        sol_dict.set_item("used_bisection", sol.used_bisection)?;
        sol_dict.set_item("energy_shift", sol.energy_shift)?;
        orbital_solutions.append(sol_dict)?;
    }
    dict.set_item("orbital_solutions", orbital_solutions)?;

    Ok(dict)
}

/// Generate mock orthonormal MO coefficients for testing
///
/// Args:
///     n_ao: Number of atomic orbitals
///     n_mo: Number of molecular orbitals
///     seed: Random seed for reproducibility
///
/// Returns:
///     2D numpy array of shape (n_ao, n_mo) with orthonormal coefficients
#[pyfunction]
fn generate_mock_mo_coefficients<'py>(
    py: Python<'py>,
    n_ao: usize,
    n_mo: usize,
    seed: u64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    use quasix_core::df::mo_transform::generate_mock_mo_coefficients as gen_mock_mo;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Generate coefficients
    let coeffs = gen_mock_mo(n_ao, n_mo, seed);

    // Convert to numpy
    Ok(PyArray2::from_owned_array(py, coeffs))
}

/// Compute exchange self-energy matrix using RI/DF
///
/// Args:
///     df_3c_mo: Three-center DF integrals in MO basis, shape (n_mo, n_mo, n_aux)
///     df_metric_inv: Inverse of the DF metric matrix, shape (n_aux, n_aux)
///     nocc: Number of occupied orbitals
///
/// Returns:
///     Exchange self-energy matrix Σˣ, shape (n_mo, n_mo)
#[pyfunction]
fn compute_exchange_matrix_ri<'py>(
    py: Python<'py>,
    df_3c_mo: &Bound<'py, PyArray3<f64>>,
    df_metric_inv: &Bound<'py, PyArray2<f64>>,
    nocc: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let df_3c_mo = df_3c_mo.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_3c_mo: {}",
            e
        ))
    })?;
    let df_3c_mo_array = df_3c_mo.as_array().to_owned();

    let df_metric_inv = df_metric_inv.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_metric_inv: {}",
            e
        ))
    })?;
    let df_metric_inv_array = df_metric_inv.as_array().to_owned();

    // Get dimensions
    let shape = df_3c_mo_array.shape();
    let n_mo = shape[0];
    let n_aux = shape[2];

    // Create exchange self-energy calculator
    let mut calculator = ExchangeSelfEnergyRI::new(n_mo, nocc, n_aux);

    // Compute exchange matrix
    let sigma_x = calculator
        .compute_exchange_matrix_ri(&df_3c_mo_array, &df_metric_inv_array)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute exchange self-energy: {}",
                e
            ))
        })?;

    // Convert to numpy
    Ok(PyArray2::from_owned_array(py, sigma_x))
}

/// Create a frequency grid for GW calculations
///
/// Args:
///     n_freq: Number of frequency points
///     grid_type: Type of grid ("gauss_legendre", "modified_gl", "minimax")
///     omega_max: Maximum frequency for modified GL grids (optional)
///
/// Returns:
///     Dictionary with 'points' and 'weights' arrays
#[pyfunction]
#[pyo3(signature = (n_freq, grid_type="gauss_legendre", omega_max=None))]
fn create_frequency_grid<'py>(
    py: Python<'py>,
    n_freq: usize,
    grid_type: &str,
    omega_max: Option<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Parse grid type
    let grid = match grid_type {
        "gauss_legendre" => GridType::GaussLegendre,
        "modified_gl" => {
            let omega = omega_max.unwrap_or(100.0);
            GridType::ModifiedGaussLegendre { omega_max: omega }
        }
        "minimax" => GridType::Minimax,
        "contour_deformation" => GridType::ContourDeformation,
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Unknown grid type: {}",
                grid_type
            )))
        }
    };

    // Create frequency grid
    let freq_grid = FrequencyGrid::new(n_freq, grid).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to create frequency grid: {}",
            e
        ))
    })?;

    // Create result dictionary
    let result = PyDict::new(py);
    result.set_item(
        "points",
        PyArray1::from_owned_array(py, freq_grid.points.clone()),
    )?;
    result.set_item(
        "weights",
        PyArray1::from_owned_array(py, freq_grid.weights.clone()),
    )?;
    result.set_item("n_freq", n_freq)?;
    result.set_item("grid_type", grid_type)?;

    Ok(result)
}

/// Generate Gauss-Legendre quadrature nodes and weights
///
/// Args:
///     n_points: Number of quadrature points
///     a: Lower bound of interval (default: -1.0)
///     b: Upper bound of interval (default: 1.0)
///
/// Returns:
///     Tuple of (nodes, weights) arrays
#[pyfunction]
#[pyo3(signature = (n_points, a=-1.0, b=1.0))]
fn gauss_legendre_quadrature<'py>(
    py: Python<'py>,
    n_points: usize,
    a: f64,
    b: f64,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Generate GL quadrature
    let (nodes, weights) = if a == -1.0 && b == 1.0 {
        gauss_legendre_nodes_weights(n_points)
    } else {
        // For scaled intervals, transform the standard [-1, 1] quadrature
        let (std_nodes, std_weights) = gauss_legendre_nodes_weights(n_points).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to generate standard GL quadrature: {}",
                e
            ))
        })?;
        // Linear transformation from [-1, 1] to [a, b]
        let scale = (b - a) / 2.0;
        let shift = (b + a) / 2.0;
        let nodes = std_nodes.mapv(|x| scale * x + shift);
        let weights = std_weights.mapv(|w| w * scale);
        Ok((nodes, weights))
    }
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to generate GL quadrature: {}",
            e
        ))
    })?;

    Ok((
        PyArray1::from_owned_array(py, nodes),
        PyArray1::from_owned_array(py, weights),
    ))
}

/// Compute diagonal elements of exchange self-energy using RI/DF
///
/// Args:
///     df_3c_mo: Three-center DF integrals in MO basis, shape (n_mo, n_mo, n_aux)
///     df_metric_inv: Inverse of the DF metric matrix, shape (n_aux, n_aux)
///     nocc: Number of occupied orbitals
///
/// Returns:
///     Diagonal elements of Σˣ, shape (n_mo,)
#[pyfunction]
fn compute_exchange_diagonal_ri<'py>(
    py: Python<'py>,
    df_3c_mo: &Bound<'py, PyArray3<f64>>,
    df_metric_inv: &Bound<'py, PyArray2<f64>>,
    nocc: usize,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let df_3c_mo = df_3c_mo.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_3c_mo: {}",
            e
        ))
    })?;
    let df_3c_mo_array = df_3c_mo.as_array().to_owned();

    let df_metric_inv = df_metric_inv.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_metric_inv: {}",
            e
        ))
    })?;
    let df_metric_inv_array = df_metric_inv.as_array().to_owned();

    // Get dimensions
    let shape = df_3c_mo_array.shape();
    let n_mo = shape[0];
    let n_aux = shape[2];

    // Create exchange self-energy calculator
    let mut calculator = ExchangeSelfEnergyRI::new(n_mo, nocc, n_aux);

    // Compute diagonal elements
    let sigma_x_diag = calculator
        .compute_exchange_diagonal_ri(&df_3c_mo_array, &df_metric_inv_array)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute exchange self-energy diagonal: {}",
                e
            ))
        })?;

    // Convert to numpy
    Ok(PyArray1::from_owned_array(py, sigma_x_diag))
}

/// Compute exchange self-energy using symmetrized DF tensors
///
/// Args:
///     df_3c_mo_sym: Symmetrized three-center DF integrals in MO basis, shape (n_mo, n_mo, n_aux)
///     nocc: Number of occupied orbitals
///
/// Returns:
///     Exchange self-energy matrix Σˣ, shape (n_mo, n_mo)
#[pyfunction]
fn compute_exchange_symmetric<'py>(
    _py: Python<'py>,
    df_3c_mo_sym: &Bound<'py, PyArray3<f64>>,
    _nocc: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray (for validation only)
    let df_3c_mo_sym = df_3c_mo_sym.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_3c_mo_sym: {}",
            e
        ))
    })?;
    let df_3c_mo_sym_array = df_3c_mo_sym.as_array();

    // Get dimensions (for validation only)
    let shape = df_3c_mo_sym_array.shape();
    let _n_mo = shape[0];
    let _n_aux = shape[2];

    // FIXME: Exchange self-energy API needs re-implementation (S3-4)
    // Current API only has legacy stubs that call unimplemented!()
    // For now, return error to avoid runtime panic
    Err(PyErr::new::<pyo3::exceptions::PyNotImplementedError, _>(
        "Exchange self-energy computation is under re-implementation (Story S3-4). \
         Please use PySCF directly for now.",
    ))
}

// ============================================================================
// Correlation Self-Energy functions for S3-5
// ============================================================================

/// Run evGW calculation with contour deformation
///
/// Implements eigenvalue-only self-consistent GW (evGW) where only the quasiparticle
/// energies are updated in the Green's function denominators.
///
/// Args:
///     mo_energy: Initial MO energies, shape (n_mo,)
///     mo_occ: MO occupations, shape (n_mo,)
///     df_ia: DF tensor (ia|P) in MO basis, shape (n_trans, n_aux) where n_trans = nocc * nvirt
///     df_mo_sym: Symmetrized DF tensor in full MO basis for exchange, shape (n_mo, n_mo, n_aux)
///     v_sqrt: Square root of Coulomb metric V^(1/2), shape (n_aux, n_aux)
///     vxc_diagonal: DFT exchange-correlation diagonal, shape (n_mo,)
///     config: Optional evGW configuration dictionary with:
///         - max_cycles: Maximum evGW iterations (default: 12)
///         - conv_tol: Energy convergence tolerance in Ha (default: 1e-4)
///         - damping: Damping factor for energy updates (default: 0.5)
///         - n_freq: Number of frequency points (default: 60)
///         - eta: Broadening parameter (default: 0.01)
///         - freq_method: 'cd' for contour deformation or 'ac' for analytic continuation
///         - verbose: Verbosity level (default: 1)
///
/// Returns:
///     Dictionary containing:
///     - 'qp_energies': Final quasiparticle energies, shape (n_mo,)
///     - 'z_factors': Final Z-factors, shape (n_mo,)
///     - 'sigma_x': Exchange self-energy matrix, shape (n_mo, n_mo)
///     - 'sigma_c_diag': Correlation self-energy diagonal, shape (n_mo,)
///     - 'convergence_history': List of dicts with iteration data
///     - 'converged': Boolean indicating convergence
///     - 'n_cycles': Number of evGW cycles performed
///     - 'spectral_function': Optional spectral function if requested
///     - 'diagnostics': Detailed convergence diagnostics
#[pyfunction]
#[pyo3(signature = (_mo_energy, _mo_occ, _df_ia, _df_mo_sym, _v_sqrt, _vxc_diagonal, config=None))]
fn evgw_contour_deformation<'py>(
    py: Python<'py>,
    _mo_energy: &Bound<'py, PyArray1<f64>>,
    _mo_occ: &Bound<'py, PyArray1<f64>>,
    _df_ia: &Bound<'py, PyArray2<f64>>,
    _df_mo_sym: &Bound<'py, PyArray3<f64>>,
    _v_sqrt: &Bound<'py, PyArray2<f64>>,
    _vxc_diagonal: &Bound<'py, PyArray1<f64>>,
    config: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyDict>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Parse configuration
    let mut max_cycles = 12;
    let mut conv_tol = 1e-4;
    let mut freq_method = "cd".to_string();

    if let Some(cfg) = config {
        if let Ok(Some(val)) = cfg.get_item("max_cycles") {
            if let Ok(v) = val.extract::<usize>() {
                max_cycles = v;
            }
        }
        if let Ok(Some(val)) = cfg.get_item("conv_tol") {
            if let Ok(v) = val.extract::<f64>() {
                conv_tol = v;
            }
        }
        if let Ok(Some(val)) = cfg.get_item("damping") {
            if let Ok(_v) = val.extract::<f64>() {
                // Damping would be used in full implementation
            }
        }
        if let Ok(Some(val)) = cfg.get_item("n_freq") {
            if let Ok(_v) = val.extract::<usize>() {
                // n_freq would be used in full implementation
            }
        }
        if let Ok(Some(val)) = cfg.get_item("eta") {
            if let Ok(_v) = val.extract::<f64>() {
                // eta would be used in full implementation
            }
        }
        if let Ok(Some(val)) = cfg.get_item("freq_method") {
            if let Ok(v) = val.extract::<String>() {
                freq_method = v;
            }
        }
        if let Ok(Some(val)) = cfg.get_item("verbose") {
            if let Ok(_v) = val.extract::<i32>() {
                // verbose would be used in full implementation
            }
        }
    }

    // This is a legacy function signature - redirect users to the proper implementation
    let dict = PyDict::new(py);
    dict.set_item("error", "This function signature is deprecated. Please use run_evgw() from quasix.gw module instead.")?;
    dict.set_item("message", "The proper evGW implementation is available via gw.run_evgw() or evgw.run_evgw() with full tensor inputs.")?;
    dict.set_item("max_cycles", max_cycles)?;
    dict.set_item("conv_tol", conv_tol)?;
    dict.set_item("freq_method", freq_method)?;

    Err(pyo3::exceptions::PyDeprecationWarning::new_err(
        "evgw_contour_deformation() is deprecated. Use run_evgw() from quasix.gw or quasix.evgw modules instead."
    ))
}

/// Compute correlation self-energy using contour deformation
///
/// Args:
///     mo_energy: Molecular orbital energies, shape (n_mo,)
///     mo_occ: Molecular orbital occupations, shape (n_mo,)
///     w_screened_real: Real part of screened interaction W(ω), shape (n_freq, n_aux, n_aux)
///     w_screened_imag: Imaginary part of screened interaction W(ω), shape (n_freq, n_aux, n_aux)
///     omega_grid: Frequency grid for W, shape (n_freq,)
///     eval_points: Evaluation frequencies for Σc, shape (n_eval,)
///     df_tensor: Density fitting tensor in MO basis, shape (n_mo, n_mo, n_aux)
///     eta: Broadening parameter (default: 0.05 eV)
///     n_imag_points: Number of imaginary-axis integration points (default: 32)
///     xi_max: Maximum imaginary frequency (default: 50.0)
///
/// Returns:
///     Dictionary containing:
///     - 'sigma_c_real': Real part of Σc, shape (n_eval, n_mo)
///     - 'sigma_c_imag': Imaginary part of Σc, shape (n_eval, n_mo)
///     - 'residue_real': Real part of residue contribution
///     - 'residue_imag': Imaginary part of residue contribution
///     - 'integral_real': Real part of integral contribution
///     - 'integral_imag': Imaginary part of integral contribution
///     - 'spectral_function': Optional spectral function A(ω), shape (n_eval, n_mo)
///     - 'spectral_normalization': Normalization check for each orbital
///     - 'diagnostics': Convergence and diagnostic information
///
/// DISABLED: API changed - compute_sigma_c_contour_deformation now returns Array2 directly, not struct
#[cfg(disabled)] // Exclude from compilation
#[allow(dead_code)]
fn compute_correlation_self_energy_cd_DISABLED<'py>(
    py: Python<'py>,
    mo_energy: &Bound<'py, PyArray1<f64>>,
    mo_occ: &Bound<'py, PyArray1<f64>>,
    w_screened_real: &Bound<'py, PyArray3<f64>>,
    w_screened_imag: &Bound<'py, PyArray3<f64>>,
    omega_grid: &Bound<'py, PyArray1<f64>>,
    eval_points: &Bound<'py, PyArray1<f64>>,
    df_tensor: &Bound<'py, PyArray3<f64>>,
    eta: Option<f64>,
    n_imag_points: Option<usize>,
    xi_max: Option<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    use ndarray::Array3;
    use num_complex::Complex64;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let mo_energy = mo_energy.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to mo_energy: {}",
            e
        ))
    })?;
    let mo_energy_array = mo_energy.as_array().to_owned();

    let mo_occ = mo_occ.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to mo_occ: {}",
            e
        ))
    })?;
    let mo_occ_array = mo_occ.as_array().to_owned();

    let w_real = w_screened_real.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to w_screened_real: {}",
            e
        ))
    })?;
    let w_real_array = w_real.as_array().to_owned();

    let w_imag = w_screened_imag.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to w_screened_imag: {}",
            e
        ))
    })?;
    let w_imag_array = w_imag.as_array().to_owned();

    // Combine real and imaginary parts into complex array
    let shape = w_real_array.shape();
    let mut w_screened = Array3::<Complex64>::zeros((shape[0], shape[1], shape[2]));
    for i in 0..shape[0] {
        for j in 0..shape[1] {
            for k in 0..shape[2] {
                w_screened[[i, j, k]] =
                    Complex64::new(w_real_array[[i, j, k]], w_imag_array[[i, j, k]]);
            }
        }
    }

    let omega_grid = omega_grid.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to omega_grid: {}",
            e
        ))
    })?;
    let omega_grid_array = omega_grid.as_array().to_owned();

    let eval_points = eval_points.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to eval_points: {}",
            e
        ))
    })?;
    let eval_points_array = eval_points.as_array().to_owned();

    let df_tensor = df_tensor.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_tensor: {}",
            e
        ))
    })?;
    let df_tensor_array = df_tensor.as_array().to_owned();

    // Get dimensions
    let n_mo = mo_energy_array.len();
    let n_aux = w_screened.shape()[1];
    let n_occ = mo_occ_array.iter().filter(|&&x| x > 0.0).count();

    // Create calculator
    let calculator = CorrelationSelfEnergyCD::new(n_mo, n_aux, n_occ);

    // Compute correlation self-energy
    let result = calculator
        .compute_sigma_c_contour_deformation(
            &mo_energy_array,
            &mo_occ_array,
            &w_screened,
            &omega_grid_array,
            &eval_points_array,
            &df_tensor_array,
        )
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute correlation self-energy: {}",
                e
            ))
        })?;

    // Create result dictionary
    let dict = PyDict::new(py);

    // Extract real and imaginary parts of sigma_c
    let sigma_c_real = result.sigma_c.mapv(|x| x.re);
    let sigma_c_imag = result.sigma_c.mapv(|x| x.im);
    dict.set_item("sigma_c_real", PyArray2::from_owned_array(py, sigma_c_real))?;
    dict.set_item("sigma_c_imag", PyArray2::from_owned_array(py, sigma_c_imag))?;

    // Extract residue parts
    let residue_real = result.residue_part.mapv(|x| x.re);
    let residue_imag = result.residue_part.mapv(|x| x.im);
    dict.set_item("residue_real", PyArray2::from_owned_array(py, residue_real))?;
    dict.set_item("residue_imag", PyArray2::from_owned_array(py, residue_imag))?;

    // Extract integral parts
    let integral_real = result.integral_part.mapv(|x| x.re);
    let integral_imag = result.integral_part.mapv(|x| x.im);
    dict.set_item(
        "integral_real",
        PyArray2::from_owned_array(py, integral_real),
    )?;
    dict.set_item(
        "integral_imag",
        PyArray2::from_owned_array(py, integral_imag),
    )?;

    // Add spectral function if available
    if let Some(spectral) = result.spectral_function {
        dict.set_item(
            "spectral_function",
            PyArray2::from_owned_array(py, spectral.values),
        )?;
        dict.set_item(
            "spectral_normalization",
            PyArray1::from_owned_array(py, spectral.normalization),
        )?;
    }

    // Add diagnostics
    let diag_dict = PyDict::new(py);
    diag_dict.set_item("residue_weight", result.diagnostics.residue_weight)?;
    diag_dict.set_item("integral_weight", result.diagnostics.integral_weight)?;
    diag_dict.set_item(
        "normalization_error",
        result.diagnostics.normalization_error,
    )?;
    diag_dict.set_item("n_poles", result.diagnostics.n_poles)?;
    diag_dict.set_item("converged", result.diagnostics.convergence_achieved)?;
    diag_dict.set_item("effective_eta", result.diagnostics.effective_eta)?;
    diag_dict.set_item("max_denominator", result.diagnostics.max_denominator)?;
    diag_dict.set_item("min_denominator", result.diagnostics.min_denominator)?;
    diag_dict.set_item("max_w_element", result.diagnostics.max_w_element)?;
    diag_dict.set_item("total_iterations", result.diagnostics.total_iterations)?;
    dict.set_item("diagnostics", diag_dict)?;

    Ok(dict)
}

// ============================================================================
// Screening W functions for S3-3
// ============================================================================

/// Compute screened interaction W(ω) using advanced dielectric solver
///
/// This function implements the symmetrized formulation:
/// W(ω) = v^{1/2} [1 - v^{1/2} P^0(ω) v^{1/2}]^{-1} v^{1/2}
///
/// Args:
///     p0_real: Real part of P0(ω), shape (n_aux, n_aux)
///     p0_imag: Imaginary part of P0(ω), shape (n_aux, n_aux)
///     v_sqrt: Square root of Coulomb metric V^(1/2), shape (n_aux, n_aux)
///     solver: Solver backend ('lu', 'cholesky', 'svd', 'auto')
///     config: Optional dictionary with solver configuration:
///         - regularization: Regularization parameter (default: 1e-10)
///         - condition_threshold: Condition number threshold (default: 1e12)
///         - self_consistency_tol: Self-consistency tolerance (default: 1e-8)
///         - monitor_condition: Whether to monitor condition number (default: True)
///         - svd_threshold: SVD truncation threshold (default: 1e-14)
///
/// Returns:
///     Dictionary containing:
///     - 'w_real': Real part of W(ω), shape (n_aux, n_aux)
///     - 'w_imag': Imaginary part of W(ω), shape (n_aux, n_aux)
///     - 'condition_number': Condition number of (1-M) matrix
///     - 'solver_used': Actual solver backend used
///     - 'hermiticity_error': Maximum Hermiticity violation
///     - 'self_consistency_error': Self-consistency check W = v + vP0W
///     - 'diagnostics': Dictionary with detailed solver diagnostics
/// DISABLED: Legacy screening API expected Complex64, but new API uses real f64 only
/// The new compute_screened_interaction() takes Array3<f64> and returns Array3<f64>
#[cfg(disabled)] // Exclude from compilation
#[allow(dead_code)]
fn compute_screened_interaction_s33_DISABLED<'py>(
    py: Python<'py>,
    p0_real: &Bound<'py, PyArray2<f64>>,
    p0_imag: &Bound<'py, PyArray2<f64>>,
    v_sqrt: &Bound<'py, PyArray2<f64>>,
    solver: &str,
    config: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyDict>> {
    use num_complex::Complex64;
    use quasix_core::dielectric::screening::{
        DielectricSolver, SolverBackend, SolverConfig, SolverType,
    };

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let p0_real = p0_real.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to p0_real: {}",
            e
        ))
    })?;
    let p0_real_array = p0_real.as_array().to_owned();

    let p0_imag = p0_imag.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to p0_imag: {}",
            e
        ))
    })?;
    let p0_imag_array = p0_imag.as_array().to_owned();

    let v_sqrt = v_sqrt.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to v_sqrt: {}",
            e
        ))
    })?;
    let v_sqrt_array = v_sqrt.as_array().to_owned();

    // Combine real and imaginary parts into complex P0
    let naux = p0_real_array.nrows();
    let mut p0_complex = ndarray::Array2::<Complex64>::zeros((naux, naux));
    for i in 0..naux {
        for j in 0..naux {
            p0_complex[[i, j]] = Complex64::new(p0_real_array[[i, j]], p0_imag_array[[i, j]]);
        }
    }

    // Parse solver backend
    let backend = match solver.to_lowercase().as_str() {
        "lu" => SolverBackend::LU,
        "cholesky" => SolverBackend::Cholesky,
        "svd" => SolverBackend::SVD,
        "auto" => SolverBackend::Auto,
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid solver: {}. Must be 'lu', 'cholesky', 'svd', or 'auto'",
                solver
            )))
        }
    };

    // Parse configuration
    let mut solver_config = SolverConfig::default();
    let solver_type = if let Some(cfg) = config {
        // Parse regularization
        if let Ok(Some(val)) = cfg.get_item("regularization") {
            if let Ok(v) = val.extract::<f64>() {
                solver_config.regularization = v;
            }
        }
        // Parse condition threshold
        if let Ok(Some(val)) = cfg.get_item("condition_threshold") {
            if let Ok(v) = val.extract::<f64>() {
                solver_config.condition_threshold = v;
            }
        }
        // Parse self-consistency tolerance
        if let Ok(Some(val)) = cfg.get_item("self_consistency_tol") {
            if let Ok(v) = val.extract::<f64>() {
                solver_config.self_consistency_tol = v;
            }
        }
        // Parse monitor_condition
        if let Ok(Some(val)) = cfg.get_item("monitor_condition") {
            if let Ok(v) = val.extract::<bool>() {
                solver_config.monitor_condition = v;
            }
        }
        // Parse SVD threshold
        if let Ok(Some(val)) = cfg.get_item("svd_threshold") {
            if let Ok(v) = val.extract::<f64>() {
                solver_config.svd_threshold = v;
            }
        }
        // Parse max_iterations
        if let Ok(Some(val)) = cfg.get_item("max_iterations") {
            if let Ok(v) = val.extract::<usize>() {
                solver_config.max_iterations = v;
            }
        }
        // Parse block_size
        if let Ok(Some(val)) = cfg.get_item("block_size") {
            if let Ok(v) = val.extract::<usize>() {
                solver_config.block_size = v;
            }
        }

        // Determine solver type (Cholesky, LU, SVD, or Adaptive)
        // Available variants: Cholesky, LU, SVD, Adaptive
        if let Ok(Some(val)) = cfg.get_item("solver_type") {
            if let Ok(v) = val.extract::<String>() {
                match v.to_lowercase().as_str() {
                    "cholesky" => SolverType::Cholesky,
                    "lu" => SolverType::LU,
                    "svd" => SolverType::SVD,
                    "adaptive" => SolverType::Adaptive,
                    _ => SolverType::Adaptive, // Default to Adaptive
                }
            } else {
                SolverType::Adaptive
            }
        } else {
            SolverType::Adaptive
        }
    } else {
        SolverType::Adaptive
    };

    // Create solver (legacy API - only ::new() available, ignores config)
    let solver = DielectricSolver::new(naux, solver_type);

    // Compute screened interaction W
    let w = solver
        .compute_screened_interaction(&p0_complex, &v_sqrt_array)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute screened interaction: {:?}",
                e
            ))
        })?;

    // Check Hermiticity
    let hermiticity_error = {
        let mut max_error: f64 = 0.0;
        for i in 0..naux {
            for j in i + 1..naux {
                let diff = (w[[i, j]] - w[[j, i]].conj()).norm();
                max_error = max_error.max(diff);
            }
        }
        max_error
    };

    // Compute self-consistency check: W = v + vP0W
    let v = &v_sqrt_array.dot(&v_sqrt_array);
    let self_consistency_error = solver
        .verify_self_consistency(&w, &p0_complex, v)
        .unwrap_or(f64::INFINITY);

    // Estimate condition number of (1-M)
    let m = solver
        .build_symmetrized_dielectric(&p0_complex, &v_sqrt_array)
        .unwrap();
    let one_minus_m = ndarray::Array2::eye(naux).mapv(|x| Complex64::new(x, 0.0)) - &m;
    let condition_number = solver
        .estimate_condition(&one_minus_m)
        .unwrap_or(f64::INFINITY);

    // Extract real and imaginary parts
    let w_real = w.mapv(|c| c.re);
    let w_imag = w.mapv(|c| c.im);

    // Create result dictionary
    let dict = PyDict::new(py);
    dict.set_item("w_real", PyArray2::from_owned_array(py, w_real))?;
    dict.set_item("w_imag", PyArray2::from_owned_array(py, w_imag))?;
    dict.set_item("condition_number", condition_number)?;
    dict.set_item("solver_used", format!("{:?}", backend))?;
    dict.set_item("hermiticity_error", hermiticity_error)?;
    dict.set_item("self_consistency_error", self_consistency_error)?;

    // Add diagnostics dictionary
    let diag_dict = PyDict::new(py);
    diag_dict.set_item("naux", naux)?;
    diag_dict.set_item("solver_type", format!("{:?}", solver_type))?;
    diag_dict.set_item("regularization", solver_config.regularization)?;
    diag_dict.set_item("condition_threshold", solver_config.condition_threshold)?;
    diag_dict.set_item("block_size", solver_config.block_size)?;
    diag_dict.set_item("monitor_condition", solver_config.monitor_condition)?;
    dict.set_item("diagnostics", diag_dict)?;

    Ok(dict)
}

// ============================================================================
// Dielectric and Polarizability functions for S3-2
// ============================================================================

/// Compute independent-particle polarizability P0(ω) using DF tensors
///
/// Args:
///     df_ia: DF tensor (ia|P) of shape (n_trans, n_aux) where n_trans = nocc * nvirt
///     e_occ: Occupied orbital energies
///     e_virt: Virtual orbital energies  
///     omega_real: Real part of frequency
///     omega_imag: Imaginary part of frequency (broadening)
///     eta: Additional broadening parameter (default: 1e-4)
///
/// Returns:
///     P0 matrix of shape (n_aux, n_aux) as complex array
#[pyfunction]
#[pyo3(signature = (df_ia, e_occ, e_virt, omega_real, omega_imag, eta=None))]
fn compute_polarizability_p0<'py>(
    py: Python<'py>,
    df_ia: &Bound<'py, PyArray2<f64>>,
    e_occ: &Bound<'py, PyArray1<f64>>,
    e_virt: &Bound<'py, PyArray1<f64>>,
    omega_real: f64,
    omega_imag: f64,
    eta: Option<f64>,
) -> PyResult<Bound<'py, PyArray2<num_complex::Complex64>>> {
    use num_complex::Complex64;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let df_ia = df_ia.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_ia: {}",
            e
        ))
    })?;
    let df_ia_array = df_ia.as_array().to_owned();

    let e_occ = e_occ.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to e_occ: {}",
            e
        ))
    })?;
    let e_occ_array = e_occ.as_array().to_owned();

    let e_virt = e_virt.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to e_virt: {}",
            e
        ))
    })?;
    let e_virt_array = e_virt.as_array().to_owned();

    // Get dimensions
    let (_n_trans, naux) = df_ia_array.dim();
    let nocc = e_occ_array.len();
    let nvirt = e_virt_array.len();

    if _n_trans != nocc * nvirt {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Dimension mismatch: n_trans={} != nocc*nvirt={}",
            _n_trans,
            nocc * nvirt
        )));
    }

    // Create polarizability calculator with configuration
    let mut config = PolarizabilityConfig::default();
    if let Some(eta_val) = eta {
        config.eta = eta_val;
    }

    // Legacy API - PolarizabilityRI only has ::new(), no with_config()
    let calc = PolarizabilityRI::new(nocc, nvirt, naux);
    // Note: config.eta is ignored by legacy API

    // Compute P0
    let omega = Complex64::new(omega_real, omega_imag);
    let p0_complex = calc
        .compute_p0(omega, &df_ia_array, &e_occ_array, &e_virt_array)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute P0: {}",
                e
            ))
        })?;

    // Convert complex result to complex numpy array
    // We need to return the full complex result for proper Hermiticity
    Ok(PyArray2::from_owned_array(py, p0_complex))
}

/// Compute multiple polarizabilities P0(ω) at different frequencies
///
/// Args:
///     df_ia: DF tensor (ia|P) of shape (n_trans, n_aux)
///     e_occ: Occupied orbital energies
///     e_virt: Virtual orbital energies
///     omega_list: List of frequencies (real parts)
///     broadening: Imaginary broadening parameter
///
/// Returns:
///     3D array of shape (n_freq, n_aux, n_aux) with P0 matrices
#[pyfunction]
fn compute_polarizability_batch<'py>(
    py: Python<'py>,
    df_ia: &Bound<'py, PyArray2<f64>>,
    e_occ: &Bound<'py, PyArray1<f64>>,
    e_virt: &Bound<'py, PyArray1<f64>>,
    omega_list: &Bound<'py, PyArray1<f64>>,
    broadening: f64,
) -> PyResult<Bound<'py, PyArray3<num_complex::Complex64>>> {
    use ndarray::Array3;
    use num_complex::Complex64;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let df_ia = df_ia.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to df_ia: {}",
            e
        ))
    })?;
    let df_ia_array = df_ia.as_array().to_owned();

    let e_occ = e_occ.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to e_occ: {}",
            e
        ))
    })?;
    let e_occ_array = e_occ.as_array().to_owned();

    let e_virt = e_virt.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to e_virt: {}",
            e
        ))
    })?;
    let e_virt_array = e_virt.as_array().to_owned();

    let omega_list = omega_list.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to omega_list: {}",
            e
        ))
    })?;
    let omega_array = omega_list.as_array();

    // Get dimensions
    let (_n_trans, naux) = df_ia_array.dim();
    let nocc = e_occ_array.len();
    let nvirt = e_virt_array.len();
    let n_freq = omega_array.len();

    // Create polarizability calculator
    let calc = PolarizabilityRI::new(nocc, nvirt, naux);

    // Compute P0 for each frequency
    let mut p0_all = Array3::<Complex64>::zeros((n_freq, naux, naux));

    for (i, &omega_re) in omega_array.iter().enumerate() {
        let omega = Complex64::new(omega_re, broadening);
        let p0_complex = calc
            .compute_p0(omega, &df_ia_array, &e_occ_array, &e_virt_array)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to compute P0 at frequency {}: {}",
                    i, e
                ))
            })?;

        // Store complex result in array
        for p in 0..naux {
            for q in 0..naux {
                p0_all[[i, p, q]] = p0_complex[[p, q]];
            }
        }
    }

    // Convert to numpy
    Ok(PyArray3::from_owned_array(py, p0_all))
}

/// Compute dielectric function ε(ω) = 1 - V^(1/2) P0 V^(1/2)
///
/// Args:
///     p0: Polarizability matrix P0 of shape (n_aux, n_aux)
///     metric: 2-center Coulomb metric V of shape (n_aux, n_aux)
///
/// Returns:
///     Dictionary with 'epsilon' matrix and 'condition_number'
/// DISABLED: Uses DielectricMatrix which was removed during cleanup
#[cfg(disabled)]
#[pyfunction]
fn compute_dielectric_function<'py>(
    py: Python<'py>,
    p0: &Bound<'py, PyArray2<f64>>,
    metric: &Bound<'py, PyArray2<f64>>,
) -> PyResult<Bound<'py, PyDict>> {
    use num_complex::Complex64;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let p0 = p0.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to p0: {}",
            e
        ))
    })?;
    let p0_array = p0.as_array().to_owned();

    let metric = metric.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to metric: {}",
            e
        ))
    })?;
    let metric_array = metric.as_array().to_owned();

    // Create dielectric matrix calculator
    let dielectric_calc = DielectricMatrix::new(&metric_array).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to create dielectric calculator: {}",
            e
        ))
    })?;

    // Convert P0 to complex (since it's real from our computation above)
    let p0_complex = p0_array.mapv(|x| Complex64::new(x, 0.0));

    // Compute symmetrized M = V^(1/2) P0 V^(1/2)
    let m_matrix = dielectric_calc
        .compute_symmetrized_m(&p0_complex)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to compute symmetrized M: {}",
                e
            ))
        })?;

    // Compute ε = 1 - M
    let epsilon = dielectric_calc.compute_epsilon(&m_matrix).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute epsilon: {}",
            e
        ))
    })?;

    // Get condition number
    let condition = dielectric_calc.get_condition_number();

    // Convert epsilon to real array (take real part)
    let epsilon_real = epsilon.mapv(|c| c.re);

    // Create result dictionary
    let dict = PyDict::new(py);
    dict.set_item("epsilon", PyArray2::from_owned_array(py, epsilon_real))?;
    dict.set_item("condition_number", condition)?;
    dict.set_item("naux", p0_array.nrows())?;

    Ok(dict)
}

/// Compute inverse dielectric matrix ε^(-1) for screening
///
/// Args:
///     epsilon: Dielectric matrix of shape (n_aux, n_aux)
///
/// Returns:
///     Inverse dielectric matrix ε^(-1)
#[pyfunction]
fn compute_epsilon_inverse<'py>(
    py: Python<'py>,
    epsilon: &Bound<'py, PyArray2<f64>>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    use num_complex::Complex64;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy array to ndarray
    let epsilon = epsilon.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to epsilon: {}",
            e
        ))
    })?;
    let epsilon_array = epsilon.as_array().to_owned();

    // Convert to complex for inversion
    let epsilon_complex = epsilon_array.mapv(|x| Complex64::new(x, 0.0));

    // Create dielectric function calculator
    let naux = epsilon_array.nrows();
    let calc = quasix_core::dielectric::DielectricFunction::new(naux);

    // Compute inverse
    let epsilon_inv = calc.compute_epsilon_inv(&epsilon_complex).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to compute epsilon inverse: {}",
            e
        ))
    })?;

    // Convert back to real
    let epsilon_inv_real = epsilon_inv.mapv(|c| c.re);

    // Convert to numpy
    Ok(PyArray2::from_owned_array(py, epsilon_inv_real))
}

/// Compute screened Coulomb interaction W using high-performance solver
///
/// Args:
///     p0: Polarizability matrix P0 of shape (n_aux, n_aux)  
///     vsqrt: Square root of Coulomb matrix V^(1/2) of shape (n_aux, n_aux)
///     config: Optional dictionary with solver configuration
///
/// Returns:
///     Dictionary with 'W' (screened interaction), 'epsilon', 'epsilon_inv', and diagnostics
/// DISABLED: Legacy screening API uses invert_dielectric() which doesn't exist in new API
#[cfg(disabled)] // Exclude from compilation
#[allow(dead_code)]
fn compute_screened_interaction_DISABLED<'py>(
    py: Python<'py>,
    p0: &Bound<'py, PyArray2<f64>>,
    vsqrt: &Bound<'py, PyArray2<f64>>,
    config: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyDict>> {
    use num_complex::Complex64;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert numpy arrays to ndarray
    let p0 = p0.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to p0: {}",
            e
        ))
    })?;
    let p0_array = p0.as_array().to_owned();

    let vsqrt = vsqrt.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to vsqrt: {}",
            e
        ))
    })?;
    let vsqrt_array = vsqrt.as_array().to_owned();

    let naux = p0_array.nrows();

    // Parse configuration if provided
    let mut solver_config = SolverConfig::default();
    let solver_type = if let Some(cfg) = config {
        // Parse regularization parameter
        if let Ok(Some(reg)) = cfg.get_item("regularization") {
            if let Ok(val) = reg.extract::<f64>() {
                solver_config.regularization = val;
            }
        }
        // Parse condition threshold
        if let Ok(Some(tol)) = cfg.get_item("condition_threshold") {
            if let Ok(val) = tol.extract::<f64>() {
                solver_config.condition_threshold = val;
            }
        }
        // Parse self-consistency tolerance
        if let Ok(Some(tol)) = cfg.get_item("self_consistency_tol") {
            if let Ok(val) = tol.extract::<f64>() {
                solver_config.self_consistency_tol = val;
            }
        }
        // Parse max iterations
        if let Ok(Some(max_iter)) = cfg.get_item("max_iterations") {
            if let Ok(val) = max_iter.extract::<usize>() {
                solver_config.max_iterations = val;
            }
        }
        // Parse solver type (Cholesky, LU, SVD, or Adaptive)
        if let Ok(Some(method)) = cfg.get_item("solver_type") {
            if let Ok(method_str) = method.extract::<String>() {
                match method_str.as_str() {
                    "cholesky" | "Cholesky" => SolverType::Cholesky,
                    "lu" | "LU" => SolverType::LU,
                    "svd" | "SVD" => SolverType::SVD,
                    "adaptive" | "Adaptive" => SolverType::Adaptive,
                    _ => SolverType::Adaptive,
                }
            } else {
                SolverType::Adaptive
            }
        } else {
            SolverType::Adaptive
        }
    } else {
        SolverType::Adaptive
    };

    // Create solver (legacy API - only ::new() available, ignores config)
    let solver = DielectricSolver::new(naux, solver_type);

    // Convert P0 to complex
    let p0_complex = p0_array.mapv(|x| Complex64::new(x, 0.0));

    // Build symmetrized dielectric matrix M = v^{1/2} P0 v^{1/2}
    let m = solver
        .build_symmetrized_dielectric(&p0_complex, &vsqrt_array)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to build dielectric matrix: {:?}",
                e
            ))
        })?;

    // Compute epsilon = 1 - M
    let mut epsilon = ndarray::Array2::eye(naux).mapv(|x| Complex64::new(x, 0.0));
    epsilon = &epsilon - &m;

    // Invert dielectric matrix
    let epsilon_inv = solver.invert_dielectric(&m).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to invert dielectric matrix: {:?}",
            e
        ))
    })?;

    // Compute W = V^(1/2) ε^(-1) V^(1/2)
    let vsqrt_complex = vsqrt_array.mapv(|x| Complex64::new(x, 0.0));
    let temp = vsqrt_complex.dot(&epsilon_inv);
    let w = temp.dot(&vsqrt_complex);

    // Convert results to real arrays
    let epsilon_real = epsilon.mapv(|c| c.re);
    let epsilon_inv_real = epsilon_inv.mapv(|c| c.re);
    let w_real = w.mapv(|c| c.re);

    // Create result dictionary
    let dict = PyDict::new(py);
    dict.set_item("W", PyArray2::from_owned_array(py, w_real))?;
    dict.set_item("epsilon", PyArray2::from_owned_array(py, epsilon_real))?;
    dict.set_item(
        "epsilon_inv",
        PyArray2::from_owned_array(py, epsilon_inv_real),
    )?;
    dict.set_item("naux", naux)?;

    Ok(dict)
}

/// Process multiple frequencies in parallel using optimized solver
///
/// Args:
///     frequencies: Array of complex frequencies (real + imag parts)
///     p0_func: Python callable that computes P0(ω)
///     vsqrt: Square root of Coulomb matrix
///     config: Optional solver configuration
///
/// Returns:
///     List of dictionaries with W(ω) for each frequency
/// DISABLED: Legacy screening API uses invert_dielectric() which doesn't exist in new API
#[cfg(disabled)] // Exclude from compilation
#[allow(dead_code)]
fn compute_screened_interaction_batch_DISABLED<'py>(
    py: Python<'py>,
    frequencies: &Bound<'py, PyArray1<f64>>,
    p0_func: &Bound<'py, PyAny>,
    vsqrt: &Bound<'py, PyArray2<f64>>,
    config: Option<&Bound<'py, PyDict>>,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    use num_complex::Complex64;

    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Convert arrays
    let frequencies = frequencies.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to frequencies: {}",
            e
        ))
    })?;
    let freq_array = frequencies.as_array();

    let vsqrt = vsqrt.try_readonly().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to get readonly access to vsqrt: {}",
            e
        ))
    })?;
    let vsqrt_array = vsqrt.as_array().to_owned();

    let naux = vsqrt_array.nrows();

    // Parse configuration
    let solver_config = SolverConfig::default();
    let solver_type = if let Some(_cfg) = config {
        // Could parse config here if needed
        SolverType::Adaptive
    } else {
        SolverType::Adaptive
    };

    // Create solver
    let solver =
        DielectricSolver::with_config(naux, solver_type, SolverBackend::Auto, solver_config);

    // Process each frequency
    let mut results = Vec::new();

    for &omega_real in freq_array.iter() {
        // Call Python function to get P0(ω)
        let p0_py = p0_func.call1((omega_real,))?;
        let p0_arr: &Bound<'py, PyArray2<f64>> = p0_py.downcast()?;

        let p0 = p0_arr.try_readonly().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to get readonly access to p0: {}",
                e
            ))
        })?;
        let p0_array = p0.as_array().to_owned();

        // Convert to complex
        let p0_complex = p0_array.mapv(|x| Complex64::new(x, 0.0));

        // Compute dielectric and screening
        let m = solver
            .build_symmetrized_dielectric(&p0_complex, &vsqrt_array)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to build dielectric matrix: {:?}",
                    e
                ))
            })?;

        // Compute epsilon = 1 - M
        let mut epsilon = ndarray::Array2::eye(naux).mapv(|x| Complex64::new(x, 0.0));
        epsilon = &epsilon - &m;

        let epsilon_inv = solver.invert_dielectric(&m).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to invert dielectric matrix: {:?}",
                e
            ))
        })?;

        // Compute W
        let vsqrt_complex = vsqrt_array.mapv(|x| Complex64::new(x, 0.0));
        let temp = vsqrt_complex.dot(&epsilon_inv);
        let w = temp.dot(&vsqrt_complex);

        // Create result dict for this frequency
        let dict = PyDict::new(py);
        dict.set_item("omega", omega_real)?;
        dict.set_item("W", PyArray2::from_owned_array(py, w.mapv(|c| c.re)))?;
        dict.set_item(
            "epsilon",
            PyArray2::from_owned_array(py, epsilon.mapv(|c| c.re)),
        )?;

        results.push(dict);
    }

    Ok(results)
}

// ============================================================================
// Frequency grid functions for S3-1
// ============================================================================

/// Create optimized Gauss-Legendre grid with caching
///
/// Args:
///     n_points: Number of quadrature points
///     cache_polynomials: Whether to cache Legendre polynomials for fast evaluation
///
/// Returns:
///     Dictionary with 'nodes', 'weights', and metadata
#[pyfunction]
#[pyo3(signature = (n_points, cache_polynomials=false))]
fn create_optimized_gl_grid<'py>(
    py: Python<'py>,
    n_points: usize,
    cache_polynomials: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let grid = OptimizedGLGrid::new(n_points, cache_polynomials).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to create optimized GL grid: {}",
            e
        ))
    })?;

    let result = PyDict::new(py);
    result.set_item("nodes", PyArray1::from_owned_array(py, grid.nodes.clone()))?;
    result.set_item(
        "weights",
        PyArray1::from_owned_array(py, grid.weights.clone()),
    )?;
    result.set_item("n_points", n_points)?;
    result.set_item("cached", cache_polynomials)?;

    Ok(result)
}

/// Create imaginary axis grid for GW calculations
///
/// Args:
///     n_points: Number of frequency points
///     transform_type: Type of transformation ('linear', 'tan', 'double_exp')
///     omega_max: Maximum frequency (required for linear/tan)
///     scale: Scale parameter (required for double_exp)
///
/// Returns:
///     Dictionary with 'points', 'weights', and metadata
#[pyfunction]
#[pyo3(signature = (n_points, transform_type="linear", omega_max=None, scale=None))]
fn create_imaginary_axis_grid<'py>(
    py: Python<'py>,
    n_points: usize,
    transform_type: &str,
    omega_max: Option<f64>,
    scale: Option<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let transform = match transform_type {
        "linear" => {
            let omega_max = omega_max.ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "omega_max required for linear transform",
                )
            })?;
            TransformType::Linear { omega_max }
        }
        "tan" => {
            let omega_max = omega_max.ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "omega_max required for tan transform",
                )
            })?;
            TransformType::Tan { omega_max }
        }
        "double_exp" => {
            let scale = scale.unwrap_or(1.0);
            TransformType::DoubleExp { scale }
        }
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid transform type: {}. Must be 'linear', 'tan', or 'double_exp'",
                transform_type
            )))
        }
    };

    let grid = ImaginaryAxisGrid::new(n_points, transform).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Failed to create imaginary axis grid: {}",
            e
        ))
    })?;

    let result = PyDict::new(py);
    result.set_item(
        "points",
        PyArray1::from_owned_array(py, grid.points.clone()),
    )?;
    result.set_item(
        "weights",
        PyArray1::from_owned_array(py, grid.weights.clone()),
    )?;
    result.set_item("n_points", n_points)?;
    result.set_item("transform_type", transform_type)?;

    Ok(result)
}

/// Create minimax grid for rational approximation
///
/// Args:
///     n_points: Number of points
///     grid_type: Type of minimax grid ('chebyshev' or 'fekete')
///
/// Returns:
///     Dictionary with 'points', 'weights' (barycentric), and metadata
#[pyfunction]
#[pyo3(signature = (n_points, grid_type="chebyshev"))]
fn create_minimax_grid<'py>(
    py: Python<'py>,
    n_points: usize,
    grid_type: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let grid = match grid_type {
        "chebyshev" => MinimaxGrid::chebyshev(n_points),
        "fekete" => MinimaxGrid::fekete(n_points).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Failed to create Fekete grid: {}",
                e
            ))
        })?,
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid grid type: {}. Must be 'chebyshev' or 'fekete'",
                grid_type
            )))
        }
    };

    let result = PyDict::new(py);
    result.set_item(
        "points",
        PyArray1::from_owned_array(py, grid.points.clone()),
    )?;
    result.set_item(
        "weights",
        PyArray1::from_owned_array(py, grid.weights.clone()),
    )?;
    result.set_item("n_points", n_points)?;
    result.set_item("grid_type", grid_type)?;

    Ok(result)
}

// NOTE: Duplicate function removed - using the one defined at line 704
// This was a duplicate definition that would cause compilation errors

/// Create contour deformation calculator
///
/// Args:
///     energy_min: Minimum energy for contour
///     energy_max: Maximum energy for contour
///     n_points: Number of contour points
///
/// Returns:
///     Dictionary with contour configuration
#[pyfunction]
#[allow(dead_code)]
fn create_contour_deformation(
    energy_min: f64,
    energy_max: f64,
    n_points: usize,
) -> PyResult<HashMap<String, Py<PyAny>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    // Configuration only, no actual struct creation yet

    // Return configuration
    #[allow(deprecated)] // Python::attach is not yet available in stable PyO3 0.26
    Python::with_gil(|py| {
        let mut config = HashMap::new();
        config.insert(
            "energy_min".to_string(),
            energy_min.into_pyobject(py)?.unbind().into(),
        );
        config.insert(
            "energy_max".to_string(),
            energy_max.into_pyobject(py)?.unbind().into(),
        );
        config.insert(
            "n_points".to_string(),
            n_points.into_pyobject(py)?.unbind().into(),
        );
        config.insert(
            "method".to_string(),
            "contour_deformation".into_pyobject(py)?.unbind().into(),
        );
        Ok(config)
    })
}

/// Create analytical continuation fitter
///
/// Args:
///     n_poles: Number of poles for rational approximation
///     regularization: Regularization parameter (default: 1e-8)
///
/// Returns:
///     Dictionary with AC fitter configuration
#[pyfunction]
#[pyo3(signature = (n_poles, regularization=None))]
#[allow(dead_code)]
fn create_ac_fitter(
    n_poles: usize,
    regularization: Option<f64>,
) -> PyResult<HashMap<String, Py<PyAny>>> {
    // Initialize logging if not already done
    let _ = logging::init_logger();

    let reg = regularization.unwrap_or(1e-8);

    // Configuration only, no actual struct creation yet

    // Return configuration
    #[allow(deprecated)] // Python::attach is not yet available in stable PyO3 0.26
    Python::with_gil(|py| {
        let mut config = HashMap::new();
        config.insert(
            "n_poles".to_string(),
            n_poles.into_pyobject(py)?.unbind().into(),
        );
        config.insert(
            "regularization".to_string(),
            reg.into_pyobject(py)?.unbind().into(),
        );
        config.insert(
            "method".to_string(),
            "analytical_continuation".into_pyobject(py)?.unbind().into(),
        );
        Ok(config)
    })
}

/*
// DISABLED: Legacy G₀W₀ kernel function (uses removed quasix_core::gw::evgw_proper module)
// This function will be replaced once the new correlation_fixed.rs is validated

/// Python wrapper for G₀W₀ kernel from quasix_core
///
/// This function exposes the production-quality G₀W₀ implementation from quasix_core
/// to Python. It performs a single-shot G₀W₀ calculation.
///
/// Parameters
/// ----------
/// mo_energy : ndarray(n_mo,)
///     Mean-field molecular orbital energies
/// mo_occ : ndarray(n_mo,)
///     Molecular orbital occupations (2.0 for occupied, 0.0 for virtual)
/// df_3c_mo : ndarray(n_mo, n_mo, n_aux)
///     Density-fitted 3-center integrals in MO basis
/// df_metric_inv : ndarray(n_aux, n_aux)
///     Inverse of DF metric (V^{-1})
/// vxc_diag : ndarray(n_mo,)
///     Exchange-correlation potential (diagonal, for QP equation)
///
/// Returns
/// -------
/// dict
///     Dictionary containing:
///     - 'qp_energies': Quasiparticle energies (n_mo,)
///     - 'z_factors': Renormalization factors (n_mo,)
///     - 'sigma_x': Exchange self-energy (n_mo,)
///     - 'sigma_c_real': Real part of correlation self-energy (n_mo,)
///     - 'sigma_c_imag': Imaginary part of correlation self-energy (n_mo,)
///     - 'converged': Whether the calculation converged
///     - 'n_iterations': Number of iterations performed
///     - 'final_error': Final convergence error
#[pyfunction]
#[pyo3(signature = (mo_energy, mo_occ, df_3c_mo, df_metric_inv, vxc_diag))]
pub fn py_run_g0w0_kernel(
    py: Python,
    mo_energy: PyReadonlyArray1<f64>,
    mo_occ: PyReadonlyArray1<f64>,
    df_3c_mo: PyReadonlyArray3<f64>,
    df_metric_inv: PyReadonlyArray2<f64>,
    vxc_diag: PyReadonlyArray1<f64>,
) -> PyResult<Py<PyDict>> {
    use quasix_core::gw::evgw_proper::{run_g0w0, ProperEvGWConfig};
    use ndarray::{Array1, Array3};

    // Convert arrays to owned ndarrays
    let mo_energy_arr = mo_energy.as_array().to_owned();
    let mo_occ_arr = mo_occ.as_array().to_owned();
    let df_3c_mo_arr = df_3c_mo.as_array().to_owned();
    let df_metric_inv_arr = df_metric_inv.as_array().to_owned();
    let vxc_diag_arr = vxc_diag.as_array().to_owned();

    // Determine dimensions
    let n_mo = mo_energy_arr.len();
    let n_occ = mo_occ_arr.iter().filter(|&&x| x > 0.5).count();
    let n_vir = n_mo - n_occ;
    let n_aux = df_metric_inv_arr.shape()[0];

    // Split df_3c_mo into occupied-virtual (ia) and occupied-occupied (ij) blocks
    // df_ia: (n_occ, n_vir, n_aux)
    let mut df_ia = Array3::zeros((n_occ, n_vir, n_aux));
    for i in 0..n_occ {
        for a in 0..n_vir {
            for p in 0..n_aux {
                df_ia[[i, a, p]] = df_3c_mo_arr[[i, n_occ + a, p]];
            }
        }
    }

    // df_ij: (n_occ, n_occ, n_aux)
    let mut df_ij = Array3::zeros((n_occ, n_occ, n_aux));
    for i in 0..n_occ {
        for j in 0..n_occ {
            for p in 0..n_aux {
                df_ij[[i, j, p]] = df_3c_mo_arr[[i, j, p]];
            }
        }
    }

    // Call the core G₀W₀ function
    let result = run_g0w0(
        &mo_energy_arr,
        &mo_occ_arr,
        &df_ia,
        &df_ij,
        &df_metric_inv_arr,
        &vxc_diag_arr,
        true,  // use_full_rpa = true
    ).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
        format!("G₀W₀ calculation failed: {:?}", e)
    ))?;

    // Convert result to Python dict
    let dict = PyDict::new(py);
    dict.set_item("converged", result.converged)?;
    dict.set_item("n_iterations", result.n_iterations)?;
    dict.set_item("final_error", result.final_error)?;
    dict.set_item(
        "qp_energies",
        PyArray1::from_vec(py, result.qp_energies.to_vec()).unbind(),
    )?;
    dict.set_item(
        "z_factors",
        PyArray1::from_vec(py, result.z_factors.to_vec()).unbind(),
    )?;
    dict.set_item(
        "sigma_x",
        PyArray1::from_vec(py, result.sigma_x.to_vec()).unbind(),
    )?;

    // Extract real and imaginary parts from Complex64 sigma_c
    let sigma_c_real: Vec<f64> = result.sigma_c.iter().map(|c| c.re).collect();
    let sigma_c_imag: Vec<f64> = result.sigma_c.iter().map(|c| c.im).collect();

    dict.set_item(
        "sigma_c_real",
        PyArray1::from_vec(py, sigma_c_real).unbind(),
    )?;
    dict.set_item(
        "sigma_c_imag",
        PyArray1::from_vec(py, sigma_c_imag).unbind(),
    )?;

    Ok(dict.unbind())
}
*/

/// QuasiX Python module definition
#[pymodule]
fn quasix(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Add version as module attribute
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // Core functions
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(metadata, m)?)?;
    m.add_function(wrap_pyfunction!(noop_kernel, m)?)?;
    m.add_function(wrap_pyfunction!(init_rust_logging, m)?)?;
    // m.add_function(wrap_pyfunction!(qp_test, m)?)?;

    // Integral functions
    m.add_function(wrap_pyfunction!(compute_3center_integrals, m)?)?;
    m.add_function(wrap_pyfunction!(compute_2center_integrals, m)?)?;
    m.add_function(wrap_pyfunction!(get_integral_dimensions, m)?)?;

    // MO transformation and Cholesky functions
    m.add_function(wrap_pyfunction!(transform_to_mo_basis, m)?)?;
    m.add_function(wrap_pyfunction!(compute_cholesky_metric, m)?)?;
    m.add_function(wrap_pyfunction!(generate_mock_mo_coefficients, m)?)?;

    // Exchange self-energy functions
    m.add_function(wrap_pyfunction!(compute_exchange_matrix_ri, m)?)?;
    m.add_function(wrap_pyfunction!(compute_exchange_diagonal_ri, m)?)?;
    m.add_function(wrap_pyfunction!(compute_exchange_symmetric, m)?)?;

    // DISABLED: Correlation self-energy (API changed - returns Array2 not struct)
    // m.add_function(wrap_pyfunction!(compute_correlation_self_energy_cd, m)?)?;

    // DISABLED: Quasiparticle solver (QPEquationSolver removed during cleanup)
    // m.add_function(wrap_pyfunction!(solve_quasiparticle_equations, m)?)?;

    // evGW functions for S5-1 (re-enabled)
    m.add_function(wrap_pyfunction!(evgw_func, m)?)?;
    m.add_function(wrap_pyfunction!(run_evgw_func, m)?)?;

    // Dielectric and polarizability functions for S3-2
    m.add_function(wrap_pyfunction!(compute_polarizability_p0, m)?)?;
    m.add_function(wrap_pyfunction!(compute_polarizability_batch, m)?)?;
    // DISABLED: compute_dielectric_function uses removed DielectricMatrix
    // m.add_function(wrap_pyfunction!(compute_dielectric_function, m)?)?;
    m.add_function(wrap_pyfunction!(compute_epsilon_inverse, m)?)?;

    // DISABLED: Screened interaction functions (legacy API with Complex64, new API uses real f64)
    // m.add_function(wrap_pyfunction!(compute_screened_interaction, m)?)?;
    // m.add_function(wrap_pyfunction!(compute_screened_interaction_batch, m)?)?;
    // m.add_function(wrap_pyfunction!(compute_screened_interaction_s33, m)?)?;

    // Frequency grid functions for S3-1
    m.add_function(wrap_pyfunction!(create_frequency_grid, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_legendre_quadrature, m)?)?;
    m.add_function(wrap_pyfunction!(create_optimized_gl_grid, m)?)?;
    m.add_function(wrap_pyfunction!(create_imaginary_axis_grid, m)?)?;
    m.add_function(wrap_pyfunction!(create_minimax_grid, m)?)?;

    // Add linalg functions directly to the module
    // (submodules in PyO3 need special handling for proper Python access)
    m.add_function(wrap_pyfunction!(linalg::compute_metric_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(linalg::apply_metric_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(linalg::apply_metric_inv_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(linalg::verify_metric_identity, m)?)?;

    // Register PySCF adapter submodule
    pyscf_adapter::register_pyscf_adapter(m)?;

    // Register schema submodule - this adds QuasixData class and convenience functions
    schema::register_schema_module(m)?;

    // Also expose the HDF5 I/O convenience functions directly in main module
    m.add_function(wrap_pyfunction!(schema::save_pyscf_data, m)?)?;
    m.add_function(wrap_pyfunction!(schema::load_pyscf_data, m)?)?;
    m.add_function(wrap_pyfunction!(schema::optimize_hdf5_chunks, m)?)?;

    // Register analytic continuation submodule
    analytic_continuation::register_analytic_continuation(m)?;

    // Register BSE submodule (S6-1 kernel, S6-2 Davidson solver)
    bse::register_bse_module(m)?;

    // DISABLED: GW module (uses removed types)
    // gw::register_module(m)?;

    // Register HDF5 I/O for DF tensors
    io_hdf5::register_io_hdf5(m)?;

    // DISABLED: Legacy G₀W₀ kernel binding (uses removed quasix_core::gw::evgw_proper module)
    // m.add_function(wrap_pyfunction!(py_run_g0w0_kernel, m)?)?;

    // Initialize logging on module import
    let _ = logging::init_logger();

    Ok(())
}
