//! evGW Python bindings for the quasix package
//!
//! This module provides the Python interface for evGW calculations,
//! bridging between Python/PySCF and the Rust implementation.

// Allow non-snake_case for tensor variables (iaP, ijP, abP, mnP follow standard QC notation)
#![allow(non_snake_case)]

use numpy::{
    PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3, PyUntypedArrayMethods,
};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt};
use std::collections::HashMap;

// Import from quasix_core
use quasix_core::qp::{ContourDeformationConfig, EvGWConfig, EvGWDriver, QPSolverConfig};

/// Main evGW calculation function exposed to Python
///
/// This function accepts MO data and DF tensors from Python and runs the evGW calculation.
#[pyfunction(name = "evgw")]
#[pyo3(signature = (
    mo_energy,
    mo_occ,
    iaP,
    ijP,
    abP,
    chol_v,
    vxc_dft,
    max_cycle=12,
    conv_tol=0.0001,
    conv_tol_z=0.001,
    damping=0.5,
    damping_dynamic=false,
    diis=false,
    diis_space=6,
    _diis_start_cycle=3,
    _freq_int="cd",
    nfreq=48,
    eta=0.01,
    _check_stability=true,
    verbose=1,
    energy_window=2.0,
    frontier_convergence=true,
    frontier_window=1.0,
    n_occ_active=None,
    n_vir_active=None,
    qp_solver="linearized",
    newton_max_iterations=50,
    newton_energy_tolerance=1e-8
))]
#[allow(clippy::too_many_arguments)]
pub fn evgw<'py>(
    py: Python<'py>,
    mo_energy: PyReadonlyArray1<f64>,
    mo_occ: PyReadonlyArray1<f64>,
    iaP: PyReadonlyArray2<f64>,
    ijP: PyReadonlyArray2<f64>,
    abP: PyReadonlyArray2<f64>,
    chol_v: PyReadonlyArray2<f64>,
    vxc_dft: PyReadonlyArray1<f64>,
    max_cycle: usize,
    conv_tol: f64,
    conv_tol_z: f64,
    damping: f64,
    damping_dynamic: bool,
    diis: bool,
    diis_space: usize,
    _diis_start_cycle: usize, // Unused for now
    _freq_int: &str,          // Unused for now
    nfreq: usize,
    eta: f64,
    _check_stability: bool, // Unused for now
    verbose: u32,
    // Production-quality evGW settings (VASP/FHI-aims style)
    energy_window: f64,          // Energy window for self-consistency (Ha)
    frontier_convergence: bool,  // Check convergence on frontier orbitals only
    frontier_window: f64,        // Window for frontier orbital convergence check (Ha)
    n_occ_active: Option<usize>, // Number of occupied orbitals to include in convergence
    n_vir_active: Option<usize>, // Number of virtual orbitals to include in convergence
    // QP solver selection
    qp_solver: &str,              // "linearized" (default) or "newton"
    newton_max_iterations: usize, // Max Newton iterations per orbital (default: 50)
    newton_energy_tolerance: f64, // Newton energy tolerance in Ha (default: 1e-8)
) -> PyResult<HashMap<String, Py<PyAny>>> {
    // Convert NumPy arrays to ndarray views
    let mo_energy_arr = mo_energy.as_array();
    let mo_occ_arr = mo_occ.as_array();
    let vxc_dft_arr = vxc_dft.as_array();

    // Get dimensions
    let n_mo = mo_energy_arr.len();
    let n_occ = mo_occ_arr.iter().filter(|&&x| x > 0.0).count();
    let n_vir = n_mo - n_occ;
    let n_aux = chol_v.shape()[0];

    // Convert DF tensors
    let iaP_arr = iaP.as_array();
    let ijP_arr = ijP.as_array();
    let chol_v_arr = chol_v.as_array();

    // Reshape iaP and ijP to 3D tensors
    // iaP: occupied-virtual block (n_occ, n_vir, n_aux) for polarizability
    let iaP_3d = ndarray::Array3::from_shape_vec(
        (n_occ, n_vir, n_aux),
        iaP_arr.as_slice().unwrap().to_vec(),
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Failed to reshape iaP tensor: {}",
            e
        ))
    })?;

    // miP is the full (m|iP) tensor: (n_mo, n_occ, n_aux)
    // Required for exchange self-energy which needs ALL m orbitals
    // Python flattens it to (n_mo*n_occ, n_aux)
    let miP_3d =
        ndarray::Array3::from_shape_vec((n_mo, n_occ, n_aux), ijP_arr.as_slice().unwrap().to_vec())
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Failed to reshape miP tensor: {}",
                    e
                ))
            })?;

    // CRITICAL FIX (2025-11-18): Reshape abP to 3D: (nvir, nvir, naux)
    // Virtual-virtual block needed for full DF tensor reconstruction
    let abP_arr = abP.as_array();
    let abP_3d = ndarray::Array3::from_shape_vec(
        (n_vir, n_vir, n_aux),
        abP_arr.as_slice().unwrap().to_vec(),
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Failed to reshape abP tensor: {}",
            e
        ))
    })?;

    // For exchange self-energy, we only need the occupied blocks
    // Pass iaP_3d directly to evGW driver, which will handle the block structure

    // Create mock MO coefficients (identity matrix for now)
    let mo_coeff = ndarray::Array2::eye(n_mo);

    // ===== BUILD FULL DF TENSOR (SAFE IMPLEMENTATION) =====
    // The full DF tensor (n_mo, n_mo, n_aux) is required for correct Σc computation
    // This needs to be constructed from iaP and miP blocks
    // CRITICAL FIX (2025-11-25): Replace unsafe pointer arithmetic with safe ndarray operations
    // This fixes segfault for larger auxiliary basis sets (e.g., CO with def2-tzvpp-ri)
    println!("Building full DF tensor from blocks (optimized)...");

    // Reconstruct full (n_mo, n_mo, n_aux) tensor from blocks
    // iaP: (n_occ, n_vir, n_aux) - occupied-virtual transitions
    // miP: (n_mo, n_occ, n_aux) - all-occupied blocks
    // abP: (n_vir, n_vir, n_aux) - virtual-virtual block
    // We need: (n_mo, n_mo, n_aux) - all combinations

    let mut mnP_3d = ndarray::Array3::<f64>::zeros((n_mo, n_mo, n_aux));

    // Fill (m|iP) block from miP_3d using safe ndarray slicing
    // miP_3d has shape (n_mo, n_occ, n_aux)
    // We fill mnP_3d[0:n_mo, 0:n_occ, :]
    for m in 0..n_mo {
        for i in 0..n_occ {
            for p in 0..n_aux {
                mnP_3d[[m, i, p]] = miP_3d[[m, i, p]];
            }
        }
    }

    // Fill occupied-virtual block from iaP
    // iaP has shape (n_occ, n_vir, n_aux)
    // We fill mnP_3d[0:n_occ, n_occ:n_mo, :]
    for i in 0..n_occ {
        for a in 0..n_vir {
            for p in 0..n_aux {
                mnP_3d[[i, n_occ + a, p]] = iaP_3d[[i, a, p]];
            }
        }
    }

    // Fill virtual-virtual block from abP
    // abP has shape (n_vir, n_vir, n_aux)
    // We fill mnP_3d[n_occ:n_mo, n_occ:n_mo, :]
    for a in 0..n_vir {
        for b in 0..n_vir {
            for p in 0..n_aux {
                mnP_3d[[n_occ + a, n_occ + b, p]] = abP_3d[[a, b, p]];
            }
        }
    }

    // Symmetrize: mnP[m, n, P] should equal mnP[n, m, P]
    for m in 0..n_mo {
        for n in (m + 1)..n_mo {
            for p in 0..n_aux {
                let avg = (mnP_3d[[m, n, p]] + mnP_3d[[n, m, p]]) * 0.5;
                mnP_3d[[m, n, p]] = avg;
                mnP_3d[[n, m, p]] = avg;
            }
        }
    }

    let mnP_norm = mnP_3d.iter().map(|&x| x * x).sum::<f64>().sqrt();
    println!(
        "✓ Full DF tensor built from blocks (optimized): shape {:?}, norm {:.2e}",
        mnP_3d.shape(),
        mnP_norm
    );

    if mnP_norm < 1e-10 {
        log::warn!("Full DF tensor has very small norm - may indicate missing data");
    }

    // Parse QP solver type from string parameter
    let use_solved_qp = match qp_solver.to_lowercase().as_str() {
        "linearized" | "linear" => false,
        "newton" | "newton-raphson" | "solved" => true,
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid qp_solver '{}'. Valid options: 'linearized' (default) or 'newton'",
                qp_solver
            )));
        }
    };

    // Create configuration
    let config = EvGWConfig {
        max_iterations: max_cycle,
        energy_tolerance: conv_tol,
        density_tolerance: conv_tol_z,
        damping_factor: damping,
        z_min: 0.1,
        z_max: 0.999,
        use_diis: diis,
        diis_space_dim: diis_space,
        print_level: verbose,
        adaptive_damping: damping_dynamic,
        core_threshold: 10.0, // Use frozen-core approximation for |ε| > 10 Ha
        // Production-quality evGW settings (VASP/FHI-aims style)
        energy_window,        // Energy window for self-consistency (Ha)
        frontier_convergence, // Check convergence on frontier orbitals only
        frontier_window,      // Window for frontier orbital convergence check (Ha)
        n_occ_active,         // Number of occupied orbitals to include
        n_vir_active,         // Number of virtual orbitals to include
        cd_config: ContourDeformationConfig {
            eta,
            n_imag_points: nfreq,
            n_imag_freq: nfreq, // Backward compatibility alias
            xi_max: 5.0,        // CRITICAL FIX: Match PySCF's iw_cutoff=5.0 default
            omega_max: 5.0,     // CRITICAL FIX: Match PySCF default
            use_gl_quadrature: true,
            convergence_tol: 1e-10,
            regularization: 1e-10,
            use_simd: true,
            n_threads: None,
            verbose: verbose as usize,
            pole_threshold: 0.01,
            compute_spectral: false,
        },
        qp_solver_config: QPSolverConfig {
            energy_tolerance: newton_energy_tolerance,
            residual_tolerance: 1e-6,
            max_newton_iterations: newton_max_iterations,
            max_bisection_iterations: 20,
            initial_damping: damping,
            min_damping: 0.1,
            max_damping: 0.9,
            max_energy_step: 1.0,
            derivative_delta: 1e-6,
            use_line_search: false,
            core_threshold: Some(10.0), // Freeze core orbitals with |ε| > 10 Ha
            use_richardson: false,
            use_bisection_fallback: true,
            z_bounds: (0.1, 0.999),
            n_threads: None,
            damping_factor: damping,
        },
        compute_final_selfenergy: false, // Skip expensive final recalculation to avoid hanging
        // Newton-Raphson solved QP equation settings (from Python parameters)
        use_solved_qp,           // Parsed from qp_solver string parameter
        newton_energy_tolerance, // From Python parameter (default: 1e-8)
        newton_max_iterations,   // From Python parameter (default: 50)
        // CRITICAL FIX (2025-12-09): Frequency caching MUST be disabled for evGW!
        // When caching is enabled, P0/W are computed once with HF energies and reused,
        // which means evGW gives identical results to G0W0. For true evGW self-consistency,
        // P0 must be recomputed with QP energies at each iteration.
        use_frequency_caching: false, // DISABLED for correct evGW (was: true)
    };

    // Create evGW driver
    let mut driver = EvGWDriver::new(n_mo, n_aux, n_occ, config);

    // Convert arrays to owned versions
    let mo_energy_owned = mo_energy_arr.to_owned();
    let mo_occ_owned = mo_occ_arr.to_owned();
    let vxc_dft_owned = vxc_dft_arr.to_owned();
    let chol_v_2d = chol_v_arr.to_owned();

    // Extract df_mi (m,i,P) for ALL m and occupied i from the full tensor
    // This is needed for exchange self-energy computation
    let df_mi = mnP_3d.slice(ndarray::s![.., 0..n_occ, ..]).to_owned();

    // Run evGW calculation with block-structured DF tensors
    // Note: df_mi has shape (nmo, nocc, naux) - ALL m with occupied i
    // mnP_3d is the full DF tensor (nmo, nmo, naux) built from blocks
    let result = driver
        .run_evgw_blocks(
            &mo_energy_owned,
            &mo_occ_owned,
            &mo_coeff,
            &iaP_3d,
            &df_mi, // (nmo, nocc, naux) - ALL m with occupied i for exchange
            &chol_v_2d,
            &vxc_dft_owned,
            &mnP_3d, // Full DF tensor for correct Σc calculation (nmo, nmo, naux)
        )
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "evGW calculation failed: {}",
                e
            ))
        })?;

    // Exchange self-energy is now computed inside the evGW driver
    // Extract it from the result
    let sigma_x = result.sigma_x.diag().to_vec();

    // Extract correlation self-energy (diagonal elements)
    // result.sigma_c is already an Array1<Complex64> (diagonal elements)
    // Split into real and imaginary parts as expected by Python
    let sigma_c_re: Vec<f64> = result.sigma_c.iter().map(|c| c.re).collect();
    let sigma_c_im: Vec<f64> = result.sigma_c.iter().map(|c| c.im).collect();

    // Prepare output dictionary
    let mut output = HashMap::new();

    output.insert(
        "qp_energies".to_string(),
        PyArray1::from_vec(py, result.qp_energies.to_vec()).into(),
    );
    output.insert(
        "z_factors".to_string(),
        PyArray1::from_vec(py, result.z_factors.to_vec()).into(),
    );
    output.insert(
        "converged".to_string(),
        PyBool::new(py, result.converged).to_owned().unbind().into(),
    );
    output.insert(
        "n_cycles".to_string(),
        PyInt::new(py, result.n_iterations as i64)
            .to_owned()
            .unbind()
            .into(),
    );
    // Get final error from convergence history
    let final_error = result
        .convergence_history
        .max_energy_changes
        .last()
        .copied()
        .unwrap_or(0.0);
    output.insert(
        "final_error".to_string(),
        PyFloat::new(py, final_error).to_owned().unbind().into(),
    );
    output.insert(
        "sigma_x".to_string(),
        PyArray1::from_vec(py, sigma_x).into(),
    );
    // Add correlation self-energy components
    output.insert(
        "sigma_c_re".to_string(),
        PyArray1::from_vec(py, sigma_c_re).into(),
    );
    output.insert(
        "sigma_c_im".to_string(),
        PyArray1::from_vec(py, sigma_c_im).into(),
    );

    // Add iteration history for monitoring
    use pyo3::types::{PyDict, PyList};
    let history_list = PyList::empty(py);
    for (cycle, energies) in result.convergence_history.energies.iter().enumerate() {
        let iter_dict = PyDict::new(py);
        iter_dict.set_item("cycle", cycle + 1)?;
        iter_dict.set_item("qp_energies", PyArray1::from_vec(py, energies.to_vec()))?;
        // Use final z_factors for all iterations (simplified)
        iter_dict.set_item(
            "z_factors",
            PyArray1::from_vec(py, result.z_factors.to_vec()),
        )?;
        let energy_change = result
            .convergence_history
            .max_energy_changes
            .get(cycle)
            .copied()
            .unwrap_or(0.0);
        let rms_change = result
            .convergence_history
            .rms_energy_changes
            .get(cycle)
            .copied()
            .unwrap_or(0.0);
        iter_dict.set_item("energy_change", energy_change)?;
        iter_dict.set_item("rms_change", rms_change)?;
        iter_dict.set_item("damping_used", 0.5)?; // Simplified - use config damping
        iter_dict.set_item(
            "converged",
            cycle + 1 == result.n_iterations && result.converged,
        )?;
        history_list.append(iter_dict)?;
    }
    output.insert("iteration_history".to_string(), history_list.into());

    // Add additional output for compatibility
    output.insert(
        "ip".to_string(),
        PyFloat::new(py, -result.qp_energies[n_occ - 1])
            .to_owned()
            .unbind()
            .into(),
    );

    if n_vir > 0 {
        output.insert(
            "ea".to_string(),
            PyFloat::new(py, -result.qp_energies[n_occ])
                .to_owned()
                .unbind()
                .into(),
        );
        output.insert(
            "gap".to_string(),
            PyFloat::new(
                py,
                result.qp_energies[n_occ] - result.qp_energies[n_occ - 1],
            )
            .to_owned()
            .unbind()
            .into(),
        );
    } else {
        output.insert("ea".to_string(), py.None());
        output.insert("gap".to_string(), py.None());
    }

    Ok(output)
}

/// Alternative evGW interface using a configuration dictionary
#[pyfunction(name = "run_evgw")]
pub fn run_evgw(py: Python, config: &Bound<'_, PyDict>) -> PyResult<HashMap<String, Py<PyAny>>> {
    // Extract required arrays
    let mo_energy: PyReadonlyArray1<f64> = config
        .get_item("mo_energy")?
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>("mo_energy not found"))?
        .extract()?;

    let mo_occ: PyReadonlyArray1<f64> = config
        .get_item("mo_occ")?
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>("mo_occ not found"))?
        .extract()?;

    // CRITICAL FIX: Python passes iaP as 3D array (nocc, nvir, naux)
    // but the main evgw function expects 2D (nocc*nvir, naux)
    // We need to accept 3D and reshape to 2D
    let iaP_3d: PyReadonlyArray3<f64> = config
        .get_item("iaP")?
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>("iaP not found"))?
        .extract()?;

    // Get dimensions from the 3D array
    let iaP_shape = iaP_3d.shape();
    let (nocc_check, nvir_check, naux_check) = (iaP_shape[0], iaP_shape[1], iaP_shape[2]);

    // Reshape to 2D for the main evgw function
    // We'll reshape within the ndarray and pass it directly
    let iaP_arr = iaP_3d.as_array();
    let iaP_2d_owned = iaP_arr
        .to_shape((nocc_check * nvir_check, naux_check))
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Failed to reshape iaP from 3D to 2D: {}",
                e
            ))
        })?
        .to_owned();

    // Create a PyArray from the reshaped ndarray
    use numpy::ToPyArray;
    let iaP_pyarray = iaP_2d_owned.to_pyarray(py);
    let iaP_readonly: PyReadonlyArray2<f64> = iaP_pyarray.extract()?;

    let ijP: PyReadonlyArray2<f64> = config
        .get_item("ijP")?
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>("ijP not found"))?
        .extract()?;

    // CRITICAL FIX (2025-11-18): Extract abP for virtual-virtual block
    // Make it optional for backward compatibility with older code
    let abP_opt: Option<PyReadonlyArray2<f64>> = config
        .get_item("abP")?
        .map(|x| x.extract::<PyReadonlyArray2<f64>>())
        .transpose()?;

    // If abP not provided, create zeros array (backward compatibility)
    let abP_readonly = if let Some(abP_arr) = abP_opt {
        abP_arr
    } else {
        // Create dummy abP with zeros (nvir*nvir, naux)
        // This maintains backward compatibility but won't fix the symmetry issue
        let zeros = ndarray::Array2::<f64>::zeros((nvir_check * nvir_check, naux_check));
        use numpy::ToPyArray;
        let zeros_pyarray = zeros.to_pyarray(py);
        zeros_pyarray.extract()?
    };

    let chol_v: PyReadonlyArray2<f64> = config
        .get_item("chol_v")?
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>("chol_v not found"))?
        .extract()?;

    let vxc_dft: PyReadonlyArray1<f64> = config
        .get_item("vxc_dft")?
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyKeyError, _>("vxc_dft not found"))?
        .extract()?;

    // Extract full DF tensor (optional - may not be present in older code)
    // Currently unused but reserved for future full-matrix evGW implementation
    let _mnP_opt: Option<PyReadonlyArray2<f64>> = config
        .get_item("mnP")?
        .map(|x| x.extract::<PyReadonlyArray2<f64>>())
        .transpose()?;

    // Extract optional parameters with defaults
    let max_cycle = config
        .get_item("max_cycle")?
        .map(|x| x.extract::<usize>())
        .transpose()?
        .unwrap_or(12);

    let conv_tol = config
        .get_item("conv_tol")?
        .map(|x| x.extract::<f64>())
        .transpose()?
        .unwrap_or(0.0001);

    let damping = config
        .get_item("damping")?
        .map(|x| x.extract::<f64>())
        .transpose()?
        .unwrap_or(0.5);

    let freq_int = config
        .get_item("freq_int")?
        .map(|x| x.extract::<String>())
        .transpose()?
        .unwrap_or_else(|| "cd".to_string());

    let nfreq = config
        .get_item("nfreq")?
        .map(|x| x.extract::<usize>())
        .transpose()?
        .unwrap_or(48);

    let verbose = config
        .get_item("verbose")?
        .map(|x| x.extract::<u32>())
        .transpose()?
        .unwrap_or(1);

    // DIIS parameters - read from config
    let diis = config
        .get_item("diis")?
        .map(|x| x.extract::<bool>())
        .transpose()?
        .unwrap_or(true); // DIIS enabled by default for faster convergence

    let diis_space = config
        .get_item("diis_space")?
        .map(|x| x.extract::<usize>())
        .transpose()?
        .unwrap_or(6);

    // Production-quality evGW settings (VASP/FHI-aims style)
    let energy_window = config
        .get_item("energy_window")?
        .map(|x| x.extract::<f64>())
        .transpose()?
        .unwrap_or(2.0); // ~54 eV window for self-consistency

    let frontier_convergence = config
        .get_item("frontier_convergence")?
        .map(|x| x.extract::<bool>())
        .transpose()?
        .unwrap_or(true); // Check convergence on frontier orbitals only by default

    let frontier_window = config
        .get_item("frontier_window")?
        .map(|x| x.extract::<f64>())
        .transpose()?
        .unwrap_or(0.1); // ~2.7 eV window for convergence check (HOMO+LUMO)

    // Default: check convergence on HOMO + 2 nearby occupied and LUMO + 2 nearby virtual
    // This ensures robust convergence for most systems
    let n_occ_active: Option<usize> = config
        .get_item("n_occ_active")?
        .map(|x| x.extract::<usize>())
        .transpose()?
        .or(Some(3)); // Default: 3 occupied orbitals (HOMO-2, HOMO-1, HOMO)

    let n_vir_active: Option<usize> = config
        .get_item("n_vir_active")?
        .map(|x| x.extract::<usize>())
        .transpose()?
        .or(Some(3)); // Default: 3 virtual orbitals (LUMO, LUMO+1, LUMO+2)

    // QP solver selection: "linearized" (default) or "newton"
    let qp_solver = config
        .get_item("qp_solver")?
        .map(|x| x.extract::<String>())
        .transpose()?
        .unwrap_or_else(|| "linearized".to_string());

    // Newton solver parameters (only used when qp_solver="newton")
    let newton_max_iterations = config
        .get_item("newton_max_iterations")?
        .map(|x| x.extract::<usize>())
        .transpose()?
        .unwrap_or(50); // Default: 50 iterations per orbital

    let newton_energy_tolerance = config
        .get_item("newton_energy_tolerance")?
        .map(|x| x.extract::<f64>())
        .transpose()?
        .unwrap_or(1e-8); // Default: sub-uev precision

    // Call the main evgw function with the reshaped iaP and abP
    evgw(
        py,
        mo_energy,
        mo_occ,
        iaP_readonly,
        ijP,
        abP_readonly, // CRITICAL FIX (2025-11-18): Add abP parameter
        chol_v,
        vxc_dft,
        max_cycle,
        conv_tol,
        0.001, // conv_tol_z
        damping,
        false,      // damping_dynamic
        diis,       // DIIS acceleration from config
        diis_space, // DIIS space dimension from config
        3,          // diis_start_cycle
        &freq_int,
        nfreq,
        0.01, // eta
        true, // check_stability
        verbose,
        // Production-quality evGW settings (VASP/FHI-aims style)
        energy_window,
        frontier_convergence,
        frontier_window,
        n_occ_active,
        n_vir_active,
        // QP solver selection
        &qp_solver,
        newton_max_iterations,
        newton_energy_tolerance,
    )
}

/// Calculate exchange self-energy (diagonal elements) - simplified version
fn calculate_exchange_selfenergy(
    _ijP: &ndarray::ArrayView2<f64>,
    _chol_v: &ndarray::ArrayView2<f64>,
    n_occ: usize,
    n_mo: usize,
) -> Vec<f64> {
    let mut sigma_x = vec![0.0; n_mo];

    // Simplified calculation - just return zeros for now
    // A proper implementation would compute the exchange from ijP and chol_v
    // For occupied orbitals, compute exchange from ijP
    for i in 0..n_occ {
        // Placeholder: proper implementation would use ijP and chol_v
        sigma_x[i] = -0.1 * (i as f64 + 1.0); // Mock values
    }

    sigma_x
}
