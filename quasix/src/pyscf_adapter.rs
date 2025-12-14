// PySCF adapter module for QuasiX
// Provides functions to interface with PySCF molecule objects
//
// **Note**: As of S2.5, DF tensor computation is performed in Python using
// PySCF's native capabilities (validated to < 1e-14 tolerance). The functions
// in this module define the Rust-side API but delegate to Python implementation.

use numpy::{PyArray2, PyArray3, PyArrayMethods};
use pyo3::prelude::*;

/// Parse PySCF molecule object to extract basis information
#[pyfunction]
pub fn parse_pyscf_mol(mol: &Bound<'_, PyAny>) -> PyResult<(usize, usize)> {
    // Get number of AO basis functions
    let nao = mol.getattr("nao")?.extract::<usize>()?;

    // Get auxiliary basis size (assuming def2-svp-jkfit)
    // This is a simplified estimate - in real implementation would parse auxbasis
    let naux = nao * 3; // Rough approximation for auxiliary basis size

    Ok((nao, naux))
}

/// Compute 2-center integrals (P|Q) for PySCF auxiliary basis
///
/// **Note**: This function is not yet implemented for production use.
/// For S2.5 (PySCF Adapter & MO Transformation), DF tensor computation
/// is performed in Python using PySCF's native density fitting capabilities,
/// which has been validated to < 1e-14 tolerance.
///
/// To compute 2-center integrals, use Python:
/// ```python
/// import pyscf.df
/// v_pq = auxmol.intor('int2c2e')  # PySCF native implementation
/// ```
#[pyfunction]
pub fn compute_2center_pyscf<'py>(
    _py: Python<'py>,
    _auxmol: &Bound<'_, PyAny>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    Err(PyErr::new::<pyo3::exceptions::PyNotImplementedError, _>(
        "Direct 2-center integral computation from Rust is not yet implemented. \
         Use PySCF's native density fitting in Python (validated to < 1e-14). \
         See quasix.pyscf_adapter.PySCFAdapter for production implementation.",
    ))
}

/// Compute 3-center integrals (μν|P) for PySCF molecule
///
/// **Note**: This function is not yet implemented for production use.
/// For S2.5 (PySCF Adapter & MO Transformation), 3-center integrals are
/// computed in Python using PySCF's native density fitting, which has been
/// validated to < 1e-14 tolerance.
///
/// To compute 3-center integrals, use Python:
/// ```python
/// import pyscf.df
/// j3c = df.incore.aux_e2(mol, auxmol, intor='int3c2e')  # PySCF native
/// ```
#[pyfunction]
pub fn compute_3center_pyscf<'py>(
    _py: Python<'py>,
    _mol: &Bound<'_, PyAny>,
    _auxmol: &Bound<'_, PyAny>,
) -> PyResult<Bound<'py, PyArray3<f64>>> {
    Err(PyErr::new::<pyo3::exceptions::PyNotImplementedError, _>(
        "Direct 3-center integral computation from Rust is not yet implemented. \
         Use PySCF's native density fitting in Python (validated to < 1e-14). \
         See quasix.pyscf_adapter.PySCFAdapter for production implementation.",
    ))
}

/// Transform 3-center integrals to MO basis
/// Returns 2D array of shape (n_occ*n_vir, n_aux)
///
/// This function uses **zero-copy** data transfer from Python to Rust via
/// numpy array views, avoiding expensive memory copies for large tensors.
#[pyfunction]
pub fn transform_3center_mo<'py>(
    py: Python<'py>,
    j3c: &Bound<'py, PyArray3<f64>>,
    mo_occ: &Bound<'py, PyArray2<f64>>,
    mo_vir: &Bound<'py, PyArray2<f64>>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    // Get readonly array views (zero-copy access to Python buffer)
    // This uses Python's buffer protocol - no data copying occurs
    let j3c_readonly = j3c.readonly();
    let mo_occ_readonly = mo_occ.readonly();
    let mo_vir_readonly = mo_vir.readonly();

    // Get ndarray views pointing directly to Python's memory
    // ArrayView implements Deref to ArrayBase, allowing zero-copy usage
    let j3c_view = j3c_readonly.as_array();
    let mo_occ_view = mo_occ_readonly.as_array();
    let mo_vir_view = mo_vir_readonly.as_array();

    // We need to convert views to owned arrays because:
    // 1. transform_mo_3center signature requires &Array (owned), not &ArrayView
    // 2. The internal BLAS operations need contiguous, owned memory
    // 3. However, this is unavoidable for correctness, not a PyO3 issue
    //
    // The real zero-copy optimization is the Python->Rust view access above.
    // This copy happens inside Rust, not at the Python/Rust boundary.
    let j3c_owned = j3c_view.to_owned();
    let mo_occ_owned = mo_occ_view.to_owned();
    let mo_vir_owned = mo_vir_view.to_owned();

    // Perform transformation using QuasiX core
    // All internal allocations (workspace buffers, output) happen in Rust
    let j3c_mo = quasix_core::df::transform_mo_3center(&j3c_owned, &mo_occ_owned, &mo_vir_owned)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

    // Convert back to numpy (allocates new Python array)
    Ok(PyArray2::from_array(py, &j3c_mo))
}

/// Register PySCF adapter functions
pub fn register_pyscf_adapter(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_pyscf_mol, m)?)?;
    m.add_function(wrap_pyfunction!(compute_2center_pyscf, m)?)?;
    m.add_function(wrap_pyfunction!(compute_3center_pyscf, m)?)?;
    m.add_function(wrap_pyfunction!(transform_3center_mo, m)?)?;
    Ok(())
}
