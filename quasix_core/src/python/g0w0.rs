//! Python bindings for G₀W₀ calculation
//!
//! **Status**: Placeholder for re-implementation
//!
//! This module will expose the clean G₀W₀ implementation to Python following
//! the type signatures in `docs/G0W0/04_TYPE_SIGNATURES.md`.
//!
//! # Target API
//!
//! ```python
//! result = quasix_core.py_run_g0w0(
//!     mo_energy, mo_occ, df_3c_mo, df_metric_inv, vxc_diag,
//!     nfreq=32, eta=0.01, verbose=1
//! )
//! ```
//!
//! Returns dictionary with:
//! - qp_energies: np.ndarray[float64]
//! - z_factors: np.ndarray[float64]
//! - sigma_x: np.ndarray[float64]
//! - sigma_c_real: np.ndarray[float64]
//! - sigma_c_imag: np.ndarray[float64] (should be ~0)
//! - ip, ea, gap: float (in eV)

use numpy::{PyArray1, PyArray2, PyArray3, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3};
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Compute RPA polarizability P₀(iω) on imaginary frequency axis (CD method)
///
/// **Formula**: P⁰_PQ(iξ) = 4 Σᵢₐ (ia|P)(ia|Q) * εᵢₐ / (ω² + εᵢₐ²)
///
/// # Arguments
/// * `ia_P` - DF 3-center integrals [n_occ, n_virt, n_aux]
/// * `mo_energy` - MO energies [n_mo]
/// * `mo_occ` - MO occupations [n_mo]
/// * `freqs` - Imaginary frequencies [nfreq]
///
/// # Returns
/// * P₀(iω) - [nfreq, n_aux, n_aux] (REAL-valued from CD method!)
///
/// # Validation
/// * Element-wise match PySCF CD < 1e-10
/// * Symmetric: |P₀ - P₀.T| < 1e-12
/// * Real-valued (no imaginary part)
///
/// # PySCF Reference
/// `pyscf/gw/gw_cd.py::get_rho_response()` lines 135-145
#[pyfunction]
pub fn compute_polarizability_p0<'py>(
    py: Python<'py>,
    ia_p: PyReadonlyArray3<f64>,
    mo_energy: PyReadonlyArray1<f64>,
    mo_occ: PyReadonlyArray1<f64>,
    freqs: PyReadonlyArray1<f64>,
) -> PyResult<Bound<'py, PyArray3<f64>>> {
    // CHANGED: Returns REAL f64!
    use numpy::PyArray3;

    // Convert to ndarray
    let ia_p_arr = ia_p.as_array();
    let mo_energy_arr = mo_energy.as_array();
    let mo_occ_arr = mo_occ.as_array();
    let freqs_arr = freqs.as_array();

    // Call Rust implementation (returns Array3<f64> for CD method)
    let p0 = crate::dielectric::polarizability::compute_polarizability_p0(
        &ia_p_arr.to_owned(),
        &mo_energy_arr.to_owned(),
        &mo_occ_arr.to_owned(),
        &freqs_arr.to_owned(),
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("P₀ computation error: {}", e))
    })?;

    // Convert to Python (now returns REAL f64 array)
    Ok(PyArray3::from_owned_array(py, p0))
}

/// Compute screened interaction W(iω) via RPA (CD method)
///
/// **Formula**: W(iξ) = v^{1/2} [I - v^{1/2} P₀ v^{1/2}]^{-1} v^{1/2}
///
/// # Arguments
/// * `p0` - Polarizability [nfreq, n_aux, n_aux] (REAL-valued from CD!)
/// * `coulomb_matrix` - Bare Coulomb matrix [n_aux, n_aux] (real-valued)
/// * `freqs` - Imaginary frequency grid [nfreq]
///
/// # Returns
/// * W(iω) - Screened interaction [nfreq, n_aux, n_aux] (REAL-valued!)
///
/// # Validation
/// * Element-wise match PySCF CD < 1e-8
/// * Symmetric: |W - W.T| < 1e-12
/// * Real-valued (no imaginary part)
/// * Positive definite eigenvalues
///
/// # PySCF Reference
/// `pyscf/gw/gw_cd.py` lines 180-195
#[pyfunction]
#[pyo3(signature = (p0, coulomb_matrix, freqs, v_sqrt=None))]
pub fn compute_screening_w<'py>(
    py: Python<'py>,
    p0: PyReadonlyArray3<f64>, // CHANGED: REAL f64!
    coulomb_matrix: PyReadonlyArray2<f64>,
    freqs: PyReadonlyArray1<f64>,
    v_sqrt: Option<PyReadonlyArray2<f64>>,
) -> PyResult<Bound<'py, PyArray3<f64>>> {
    // CHANGED: Returns REAL f64!
    use numpy::PyArray3;

    // Convert to ndarray
    let p0_arr = p0.as_array();
    let v_aux_arr = coulomb_matrix.as_array();
    let freqs_arr = freqs.as_array();

    // Convert optional v_sqrt
    let v_sqrt_opt = v_sqrt.map(|arr| arr.as_array().to_owned());

    // Call Rust implementation (now handles REAL f64 P0 from CD method)
    let w_result = crate::dielectric::screening::compute_screened_interaction_with_vsqrt(
        &p0_arr.to_owned(),
        &v_aux_arr.to_owned(),
        &freqs_arr.to_owned(),
        v_sqrt_opt.as_ref(),
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("W(iω) computation error: {}", e))
    })?;

    // Convert to Python (now returns REAL f64 array)
    Ok(PyArray3::from_owned_array(py, w_result))
}

/// Export intermediate values for debugging W(iω) calculation
///
/// This function provides detailed diagnostics for the screening calculation
/// to identify where QuasiX diverges from PySCF.
///
/// # Arguments
/// * `coulomb_matrix` - DF Coulomb metric v [n_aux, n_aux]
///
/// # Returns
/// Dictionary with keys:
/// * cholesky_success - bool: Whether Cholesky succeeded
/// * eigvals - np.ndarray[float64]: Raw eigenvalues
/// * eigvecs - np.ndarray[float64]: Raw eigenvectors
/// * eigvals_thresholded - np.ndarray[float64]: Thresholded eigenvalues
/// * chol_v - np.ndarray[float64]: Computed v^{1/2} matrix
/// * chol_reconstruction - np.ndarray[float64]: Reconstructed v = chol_v @ chol_v.T
/// * max_reconstruction_error - float: Maximum reconstruction error
#[pyfunction]
pub fn debug_compute_screening_intermediates<'py>(
    py: Python<'py>,
    coulomb_matrix: PyReadonlyArray2<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let v_aux = coulomb_matrix.as_array().to_owned();

    // Call Rust debug function
    let (cholesky_success, eigvals, eigvecs, eigvals_thresh, chol_v, reconstructed, max_error) =
        crate::dielectric::screening::debug_compute_screening_intermediates(&v_aux).map_err(
            |e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Debug screening failed: {}",
                    e
                ))
            },
        )?;

    // Create Python dict with results
    let result = PyDict::new(py);
    result.set_item("cholesky_success", cholesky_success)?;
    result.set_item("eigvals", PyArray1::from_owned_array(py, eigvals))?;
    result.set_item("eigvecs", PyArray2::from_owned_array(py, eigvecs))?;
    result.set_item(
        "eigvals_thresholded",
        PyArray1::from_owned_array(py, eigvals_thresh),
    )?;
    result.set_item("chol_v", PyArray2::from_owned_array(py, chol_v))?;
    result.set_item(
        "chol_reconstruction",
        PyArray2::from_owned_array(py, reconstructed),
    )?;
    result.set_item("max_reconstruction_error", max_error)?;

    Ok(result)
}

/// Compute exchange self-energy Σˣ diagonal
///
/// **Formula**: Σˣ_nn' = -Σᵢ Σ_PQ (ni|P) v⁻¹_PQ (Q|n'i)
///
/// For diagonal elements (n=n'):
/// Σˣ_nn = -Σᵢ Σ_PQ (ni|P) v⁻¹_PQ (Q|ni)
///
/// # Arguments
/// * `ia_P` - DF 3-center integrals [n_occ, n_virt, n_aux]
/// * `coulomb_matrix_inv` - Inverse DF metric v⁻¹ [n_aux, n_aux]
/// * `n_mo` - Total number of molecular orbitals
///
/// # Returns
/// * Σˣ diagonal - [n_mo] (occupied orbitals have values, virtual orbitals are zero)
///
/// # Validation
/// * Element-wise match PySCF < 1e-10
/// * Real-valued (enforced by type)
/// * Occupied orbitals: negative values
/// * Virtual orbitals: zero
///
/// # PySCF Reference
/// `pyscf/gw/gw_cd.py` line 77
#[pyfunction]
pub fn compute_exchange_selfenergy<'py>(
    py: Python<'py>,
    ia_p: PyReadonlyArray3<f64>,
    coulomb_matrix_inv: PyReadonlyArray2<f64>,
    n_mo: usize,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    // Convert to ndarray
    let ia_p_arr = ia_p.as_array();
    let v_aux_inv_arr = coulomb_matrix_inv.as_array();

    // Call Rust implementation
    let sigma_x = crate::selfenergy::exchange::compute_exchange_diagonal(
        &ia_p_arr.to_owned(),
        &v_aux_inv_arr.to_owned(),
        n_mo,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Exchange self-energy computation error: {}",
            e
        ))
    })?;

    // Convert to Python
    Ok(PyArray1::from_owned_array(py, sigma_x))
}

/// Compute Green's function G(ω) at complex frequencies
///
/// **Formula**: G_n(ω) = 1 / (ω + μ - ε_n)
///
/// # Arguments
/// * `omega` - Complex frequencies [nfreq] (can be ±iω for occupied/virtual)
/// * `mo_energy` - HF orbital energies [nmo] (Hartree)
/// * `chemical_potential` - Fermi level μ (Hartree)
///
/// # Returns
/// * G(ω) - Green's function [nmo, nfreq] (complex-valued, transposed for PySCF)
///
/// # Validation
/// * All frequencies must be finite (no NaN/Inf)
/// * Result must be finite (no NaN/Inf)
///
/// # PySCF Convention
/// * Occupied orbitals: ω = -iω (negative imaginary)
/// * Virtual orbitals: ω = +iω (positive imaginary)
/// * Shape: [nmo, nfreq] (Rust internally uses [nfreq, nmo], transposed here)
///
/// # PySCF Reference
/// `pyscf/gw/gw_ac.py::get_sigma_element()` lines 156-212
#[pyfunction]
pub fn compute_green_function<'py>(
    py: Python<'py>,
    omega: PyReadonlyArray1<num_complex::Complex<f64>>,
    mo_energy: PyReadonlyArray1<f64>,
    chemical_potential: f64,
) -> PyResult<Bound<'py, PyArray2<num_complex::Complex<f64>>>> {
    let omega_arr = omega.as_array();
    let mo_energy_arr = mo_energy.as_array();

    let green = crate::selfenergy::green_function::compute_green_function(
        &omega_arr.to_owned(),
        &mo_energy_arr.to_owned(),
        chemical_potential,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Green's function error: {}", e))
    })?;

    // Transpose from Rust shape [nfreq, nmo] to PySCF shape [nmo, nfreq]
    let green_transposed = green.t().to_owned();
    Ok(PyArray2::from_owned_array(py, green_transposed))
}

/// Compute correlation self-energy Σᶜ diagonal (FIXED VERSION)
///
/// **CRITICAL**: This function now uses the FIXED implementation that matches PySCF!
/// The old version had a bug where HOMO correlation was always zero.
///
/// **Formula**: Σᶜ_n = -1/π ∫ dω Σᶜ_n(iω)
///
/// where Σᶜ_n(iω) = Σ_m (n|W(iω)|m) G_m(iω)
/// and W(iω) is computed internally from P₀(iω)
///
/// # Arguments
/// * `green_function` - Green's function G(iω) [nfreq, nmo] (complex)
/// * `p0_iw` - RPA polarizability P₀(iω) [nfreq, naux, naux] (REAL!) - CHANGED from w_iw!
/// * `lpq` - DF tensor (pq|P) [naux, nmo, nmo] (real)
/// * `freqs` - Imaginary frequency grid [nfreq] (positive real values)
/// * `weights` - Quadrature weights [nfreq] (Gauss-Legendre)
/// * `nocc` - Number of occupied orbitals
///
/// # Returns
/// * Σᶜ diagonal - [nmo] (complex, should be nearly real)
///
/// # Validation
/// * Element-wise match PySCF < 1e-10
/// * |Im(Σᶜ)| < 1e-5 (warns if violated)
///
/// # PySCF Reference
/// `pyscf/gw/gw_ac.py::get_sigma_diag()` lines 203-213
#[pyfunction]
pub fn compute_sigma_c_diagonal<'py>(
    py: Python<'py>,
    green_function: PyReadonlyArray2<num_complex::Complex<f64>>,
    p0_iw: PyReadonlyArray3<f64>, // CHANGED: Now takes P0, not W!
    lpq: PyReadonlyArray3<f64>,
    freqs: PyReadonlyArray1<f64>,
    weights: PyReadonlyArray1<f64>,
    nocc: usize,
) -> PyResult<Bound<'py, PyArray1<num_complex::Complex<f64>>>> {
    // Convert arrays
    let green_arr = green_function.as_array().to_owned();
    let p0_arr = p0_iw.as_array().to_owned(); // CHANGED: p0, not w!
    let lpq_arr = lpq.as_array().to_owned();
    let freqs_arr = freqs.as_array().to_owned();
    let weights_arr = weights.as_array().to_owned();

    // CRITICAL: Use the FIXED version that matches PySCF!
    // The old correlation::compute_sigma_c_diagonal had a bug
    // where HOMO correlation was zero.
    let sigma_c = crate::selfenergy::compute_sigma_c_diagonal_fixed(
        &green_arr,
        &p0_arr, // Pass P0 - W is computed internally!
        &lpq_arr,
        &freqs_arr,
        &weights_arr,
        nocc,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Correlation self-energy error: {}",
            e
        ))
    })?;

    Ok(PyArray1::from_owned_array(py, sigma_c))
}

/// Extract correlation self-energy at ω=0 from imaginary frequency data
///
/// **CRITICAL**: PySCF baseline sigma_c_iw ALREADY includes -1/π normalization!
/// This function just extracts the ω=0 value (first frequency point).
///
/// # Formula
/// Σᶜ_n(ω=0) = sigma_c_iw[n, 0]  (extract first frequency point)
///
/// # Arguments
/// * `sigma_c_iw` - Σᶜ(iω) with -1/π already applied [nmo, nfreq] (complex, PySCF shape)
/// * `freqs` - Imaginary frequencies [nfreq] (unused, kept for API compatibility)
/// * `weights` - Quadrature weights [nfreq] (unused, kept for API compatibility)
///
/// # Returns
/// * Extracted Σᶜ [nmo] (complex, should be nearly real)
///
/// # Physical Checks
/// * |Im(Σᶜ)| should be < 1e-5 (warns if violated)
/// * Result must be finite
///
/// # Note
/// Input is transposed from PySCF shape [nmo, nfreq] to Rust shape [nfreq, nmo]
#[pyfunction]
pub fn integrate_sigma_c<'py>(
    py: Python<'py>,
    sigma_c_iw: PyReadonlyArray2<num_complex::Complex<f64>>,
    freqs: PyReadonlyArray1<f64>,
    weights: PyReadonlyArray1<f64>,
) -> PyResult<Bound<'py, PyArray1<num_complex::Complex<f64>>>> {
    // Transpose from PySCF shape [nmo, nfreq] to Rust shape [nfreq, nmo]
    let sigma_transposed = sigma_c_iw.as_array().t().to_owned();
    let freqs_arr = freqs.as_array().to_owned();
    let weights_arr = weights.as_array().to_owned();

    let sigma_c = crate::selfenergy::correlation::integrate_sigma_c(
        &sigma_transposed,
        &freqs_arr,
        &weights_arr,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Σᶜ extraction error: {}", e))
    })?;

    Ok(PyArray1::from_owned_array(py, sigma_c))
}

/// Solve quasiparticle equation (linearized G₀W₀)
///
/// **Formula**: E_QP = ε_HF + Σˣ + Σᶜ - V_xc
///
/// # Arguments
/// * `mo_energy` - HF eigenvalues [n_mo] (Hartree)
/// * `sigma_x` - Exchange self-energy diagonal [n_mo] (Hartree)
/// * `sigma_c` - Correlation self-energy diagonal [n_mo] (Hartree, real part)
/// * `vxc` - Exchange-correlation potential diagonal [n_mo] (Hartree)
///
/// # Returns
/// Dictionary with keys:
/// * qp_energies - np.ndarray[float64] [n_mo] (Hartree)
/// * z_factors - np.ndarray[float64] [n_mo] (all 1.0 for G₀W₀)
/// * converged - bool (always True for linearized G₀W₀)
///
/// # Validation
/// * Element-wise match PySCF < 1e-8 Ha (< 0.01 eV)
/// * All energies finite
/// * HOMO < LUMO (positive gap)
///
/// # PySCF Reference
/// `pyscf/gw/gw_ac.py::kernel()` lines 112-125 (linearized branch)
#[pyfunction]
pub fn solve_qp_equation<'py>(
    py: Python<'py>,
    mo_energy: PyReadonlyArray1<f64>,
    sigma_x: PyReadonlyArray1<f64>,
    sigma_c: PyReadonlyArray1<f64>,
    vxc: PyReadonlyArray1<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    // Convert to ndarray
    let mo_energy_arr = mo_energy.as_array();
    let sigma_x_arr = sigma_x.as_array();
    let sigma_c_arr = sigma_c.as_array();
    let vxc_arr = vxc.as_array();

    // Call Rust implementation
    let (qp_energies, z_factors) = crate::qp::solver::solve_quasiparticle_linearized(
        &mo_energy_arr.to_owned(),
        &sigma_x_arr.to_owned(),
        &sigma_c_arr.to_owned(),
        &vxc_arr.to_owned(),
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("QP solver error: {}", e))
    })?;

    // Create result dictionary
    let result = PyDict::new(py);
    result.set_item("qp_energies", PyArray1::from_owned_array(py, qp_energies))?;
    result.set_item("z_factors", PyArray1::from_owned_array(py, z_factors))?;
    result.set_item("converged", true)?; // Linearized G₀W₀ always converges (one-shot)

    Ok(result)
}

/// Compute spectroscopic properties from QP energies
///
/// # Arguments
/// * `qp_energies` - Quasiparticle energies [n_mo] (Hartree)
/// * `n_occ` - Number of occupied orbitals
///
/// # Returns
/// Dictionary with keys:
/// * ip - Ionization potential (eV)
/// * ea - Electron affinity (eV)
/// * gap - Fundamental gap (eV)
/// * homo_index - HOMO orbital index (int)
/// * lumo_index - LUMO orbital index (int)
///
/// # PySCF Compatibility
/// Results match PySCF definition:
/// * IP = -E_HOMO (positive value)
/// * EA = -E_LUMO (negative value for stable molecules)
/// * gap = IP - EA = E_LUMO - E_HOMO (positive)
#[pyfunction]
pub fn compute_spectroscopic_properties<'py>(
    py: Python<'py>,
    qp_energies: PyReadonlyArray1<f64>,
    n_occ: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let qp_arr = qp_energies.as_array();

    let (ip, ea, gap) =
        crate::qp::solver::compute_spectroscopic_properties(&qp_arr.to_owned(), n_occ).map_err(
            |e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Spectroscopic properties error: {}",
                    e
                ))
            },
        )?;

    let result = PyDict::new(py);
    result.set_item("ip", ip)?;
    result.set_item("ea", ea)?;
    result.set_item("gap", gap)?;
    result.set_item("homo_index", n_occ - 1)?;
    result.set_item("lumo_index", n_occ)?;

    Ok(result)
}

/// Compute correlation self-energy Σᶜ diagonal (FIXED version matching PySCF AC)
///
/// **CRITICAL BUG FIX**: This is the corrected implementation that matches PySCF AC exactly.
/// Use this instead of `compute_sigma_c_diagonal` until the original is fixed.
///
/// # Formula
///
/// Σᶜ_n = -1/π Σ_w weight_w * Σ_m (n|W(iω_w)|m) * G_m(iω_w)
///
/// where (n|W|m) = Σ_PQ (n|P) [ε^{-1}_PQ(iω) - δ_PQ] (Q|m)
///
/// # Arguments
/// * `green_function` - Green's function G(iω) [nfreq, nmo] (complex)
/// * `p0_iw` - Polarizability P₀(iω) [nfreq, naux, naux] (REAL!)
/// * `lpq` - DF tensor (P|pq) [naux, nmo, nmo] (real)
/// * `freqs` - Imaginary frequency grid [nfreq] (positive real values)
/// * `weights` - Quadrature weights [nfreq] (includes Jacobian)
/// * `nocc` - Number of occupied orbitals (unused, for API compatibility)
///
/// # Returns
/// * Σᶜ diagonal - [nmo] (complex, should be nearly real)
///
/// # Validation
/// * Element-wise match PySCF AC < 1e-10
/// * |Im(Σᶜ)| < 1e-5
///
/// # PySCF Reference
/// `pyscf/gw/gw_ac.py::get_sigma_diag()` lines 203-213
#[pyfunction]
pub fn compute_sigma_c_diagonal_fixed<'py>(
    py: Python<'py>,
    green_function: PyReadonlyArray2<num_complex::Complex<f64>>,
    p0_iw: PyReadonlyArray3<f64>,
    lpq: PyReadonlyArray3<f64>,
    freqs: PyReadonlyArray1<f64>,
    weights: PyReadonlyArray1<f64>,
    nocc: usize,
) -> PyResult<Bound<'py, PyArray1<num_complex::Complex<f64>>>> {
    // Convert arrays
    let green_arr = green_function.as_array().to_owned();
    let p0_arr = p0_iw.as_array().to_owned();
    let lpq_arr = lpq.as_array().to_owned();
    let freqs_arr = freqs.as_array().to_owned();
    let weights_arr = weights.as_array().to_owned();

    let sigma_c = crate::selfenergy::correlation_fixed::compute_sigma_c_diagonal_fixed(
        &green_arr,
        &p0_arr,
        &lpq_arr,
        &freqs_arr,
        &weights_arr,
        nocc,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Correlation self-energy (FIXED AC) error: {}",
            e
        ))
    })?;

    Ok(PyArray1::from_owned_array(py, sigma_c))
}

/// Compute correlation self-energy Σᶜ using PySCF CD formula
///
/// **MISSION**: Implement exact PySCF contour deformation formula
/// with target-dependent Green's function.
///
/// # Formula (PySCF CD)
///
/// For each target orbital n:
///   Σᶜ_n(ω=0) = -1/π Σ_freq Σ_m (n|W(iω)|m) * G_m^(n)(iω)
///
/// where Green's function depends on TARGET:
///   emo = -1j*η*sign(ε_n - ε_m) - ε_m
///   G_m^(n)(iω) = weight * emo / (emo² + ω²)
///
/// # Arguments
/// * `mo_energy` - Orbital energies [nmo] (Hartree)
/// * `eta` - Broadening parameter (Ha, typically 0.01)
/// * `p0_iw` - Polarizability P₀(iω) [nfreq, naux, naux]
/// * `lpq` - DF tensor (P|pq) [naux, nmo, nmo]
/// * `freqs` - Imaginary frequency grid [nfreq] (Ha)
/// * `weights` - Quadrature weights [nfreq]
///
/// # Returns
/// * Σᶜ diagonal - [nmo] (complex, should be nearly real)
///
/// # Target
/// * Must match PySCF CD within < 0.01 eV for H₂/STO-3G
///
/// # PySCF Reference
/// `pyscf/gw/gw_cd.py::get_sigmaI_diag()` lines 169-178
#[pyfunction]
pub fn compute_sigma_c_pyscf_cd<'py>(
    py: Python<'py>,
    mo_energy: PyReadonlyArray1<f64>,
    eta: f64,
    p0_iw: PyReadonlyArray3<f64>,
    lpq: PyReadonlyArray3<f64>,
    freqs: PyReadonlyArray1<f64>,
    weights: PyReadonlyArray1<f64>,
) -> PyResult<Bound<'py, PyArray1<num_complex::Complex<f64>>>> {
    // Convert arrays
    let mo_energy_arr = mo_energy.as_array().to_owned();
    let p0_arr = p0_iw.as_array().to_owned();
    let lpq_arr = lpq.as_array().to_owned();
    let freqs_arr = freqs.as_array().to_owned();
    let weights_arr = weights.as_array().to_owned();

    let sigma_c = crate::selfenergy::compute_sigma_c_pyscf_cd(
        &mo_energy_arr,
        eta,
        &p0_arr,
        &lpq_arr,
        &freqs_arr,
        &weights_arr,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "Correlation self-energy (PySCF CD) error: {}",
            e
        ))
    })?;

    Ok(PyArray1::from_owned_array(py, sigma_c))
}

/// Compute self-energy on imaginary axis using AC method (PySCF-validated)
///
/// This is the validated implementation matching PySCF's gw_ac.py exactly.
///
/// # Arguments
/// * `mo_energy` - MO energies [n_mo]
/// * `lpq` - DF 3-center integrals [n_aux, n_mo, n_mo]
/// * `nocc` - Number of occupied orbitals (CRITICAL: pass from Python, not inferred)
/// * `orbs` - Which orbitals to compute (default: all)
/// * `nw` - Number of frequency points (default 100)
/// * `x0` - Frequency scaling parameter (default 0.5)
/// * `iw_cutoff` - Cutoff for imaginary frequency (default 5.0 Ha)
///
/// # Returns
/// Dictionary with:
/// * `sigma_c` - Correlation self-energy [n_orbs, nw_sigma] (complex)
/// * `omega` - Frequencies [n_orbs, nw_sigma] (complex)
/// * `freqs` - Quadrature frequencies [nw]
/// * `wts` - Quadrature weights [nw]
///
/// # Note
/// The `nocc` parameter MUST be passed explicitly from Python (based on electron count).
/// Previously we inferred nocc from eigenvalue signs, but this fails for DFT calculations
/// where LUMO can have negative energy (e.g., H2O/PBE, CO/PBE).
#[pyfunction]
#[pyo3(signature = (mo_energy, lpq, nocc, orbs=None, nw=100, x0=0.5, iw_cutoff=5.0))]
pub fn compute_sigma_c_ac<'py>(
    py: Python<'py>,
    mo_energy: PyReadonlyArray1<f64>,
    lpq: PyReadonlyArray3<f64>,
    nocc: usize,
    orbs: Option<Vec<usize>>,
    nw: usize,
    x0: f64,
    iw_cutoff: f64,
) -> PyResult<Bound<'py, PyDict>> {
    use crate::selfenergy::correlation_ac::{get_scaled_legendre_roots, get_sigma_diag, ACConfig};
    use numpy::{PyArray1, PyArray2};

    let mo_energy_arr = mo_energy.as_array().to_owned();
    let lpq_arr = lpq.as_array().to_owned();
    let nmo = mo_energy_arr.len();

    // Validate nocc parameter
    if nocc == 0 || nocc >= nmo {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Invalid nocc={}: must be in range [1, {})",
            nocc, nmo
        )));
    }

    // Default orbs: all orbitals
    let orbs_vec = orbs.unwrap_or_else(|| (0..nmo).collect());

    // Create config
    let config = ACConfig {
        nw,
        iw_cutoff,
        x0,
        eta: 0.0,
    };

    // Compute self-energy (with explicit nocc from Python)
    let (sigma, omega) = get_sigma_diag(
        &mo_energy_arr,
        &lpq_arr,
        &orbs_vec,
        nocc, // CRITICAL: Use explicit nocc from Python, not inferred from eigenvalue signs
        &config,
    )
    .map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("AC self-energy error: {}", e))
    })?;

    // Get frequency grid for output
    let (freqs, wts) = get_scaled_legendre_roots(nw, x0);

    // Create result dictionary
    let dict = PyDict::new(py);
    dict.set_item("sigma_c", PyArray2::from_owned_array(py, sigma))?;
    dict.set_item("omega", PyArray2::from_owned_array(py, omega))?;
    dict.set_item("freqs", PyArray1::from_owned_array(py, freqs))?;
    dict.set_item("wts", PyArray1::from_owned_array(py, wts))?;
    dict.set_item("orbs", orbs_vec)?;
    dict.set_item("nocc", nocc)?;

    Ok(dict)
}

/// Register G₀W₀ Python bindings
pub fn register_module(py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    let _ = py; // Suppress unused warning
                // S3-2: RPA polarizability P₀(iω)
    m.add_function(wrap_pyfunction!(compute_polarizability_p0, m)?)?;

    // S3-3: Screened interaction W(iω)
    m.add_function(wrap_pyfunction!(compute_screening_w, m)?)?;

    // S3-4: Exchange self-energy Σˣ
    m.add_function(wrap_pyfunction!(compute_exchange_selfenergy, m)?)?;

    // S3-5: Green's function and correlation self-energy Σᶜ
    m.add_function(wrap_pyfunction!(compute_green_function, m)?)?;
    m.add_function(wrap_pyfunction!(compute_sigma_c_diagonal, m)?)?; // Original (buggy)
    m.add_function(wrap_pyfunction!(compute_sigma_c_diagonal_fixed, m)?)?; // FIXED version (AC)
    m.add_function(wrap_pyfunction!(compute_sigma_c_pyscf_cd, m)?)?; // PySCF CD formula
    m.add_function(wrap_pyfunction!(integrate_sigma_c, m)?)?;

    // AC method (validated, matches PySCF exactly)
    m.add_function(wrap_pyfunction!(compute_sigma_c_ac, m)?)?;

    // S3-6: Quasiparticle solver
    m.add_function(wrap_pyfunction!(solve_qp_equation, m)?)?;
    m.add_function(wrap_pyfunction!(compute_spectroscopic_properties, m)?)?;

    // Debug function for S3-3 validation
    m.add_function(wrap_pyfunction!(debug_compute_screening_intermediates, m)?)?;

    // Placeholder - will add run_g0w0() function during S3-1 to S3-6 implementation
    // See: docs/G0W0/04_TYPE_SIGNATURES.md lines 96-108 for exact signature
    Ok(())
}
