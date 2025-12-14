//! HDF5 I/O bindings for the quasix Python module
//!
//! This module re-exports and wraps the S2.4 HDF5 functionality from quasix_core

use numpy::{PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use quasix_core::io::df_tensors::{DFMetadata, DFTensorsS24};
use std::path::PathBuf;

/// Python wrapper for DF tensor metadata
#[pyclass(name = "DFMetadata")]
#[derive(Clone)]
pub struct PyDFMetadata {
    inner: DFMetadata,
}

#[pymethods]
impl PyDFMetadata {
    #[new]
    fn new(
        naux: usize,
        nao: usize,
        nocc: usize,
        nvir: usize,
        auxbasis: String,
        basis: String,
    ) -> Self {
        Self {
            inner: DFMetadata {
                naux,
                nao,
                nocc,
                nvir,
                auxbasis,
                basis,
                version: "2.4".to_string(),
                timestamp: None,
                git_commit: None,
            },
        }
    }

    #[getter]
    fn naux(&self) -> usize {
        self.inner.naux
    }

    #[getter]
    fn nao(&self) -> usize {
        self.inner.nao
    }

    #[getter]
    fn nocc(&self) -> usize {
        self.inner.nocc
    }

    #[getter]
    fn nvir(&self) -> usize {
        self.inner.nvir
    }

    #[getter]
    fn auxbasis(&self) -> String {
        self.inner.auxbasis.clone()
    }

    #[getter]
    fn basis(&self) -> String {
        self.inner.basis.clone()
    }

    #[getter]
    fn version(&self) -> String {
        self.inner.version.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "DFMetadata(naux={}, nao={}, nocc={}, nvir={}, auxbasis='{}', basis='{}')",
            self.inner.naux,
            self.inner.nao,
            self.inner.nocc,
            self.inner.nvir,
            self.inner.auxbasis,
            self.inner.basis
        )
    }
}

/// Python wrapper for S2.4 DF tensor HDF5 I/O
#[pyclass(name = "DFTensorsHDF5")]
pub struct PyDFTensorsHDF5 {
    inner: DFTensorsS24,
}

#[pymethods]
impl PyDFTensorsHDF5 {
    #[new]
    #[pyo3(signature = (metric_2c, cderi_3c, ia_p, ij_p, naux, nao, nocc, nvir, auxbasis, basis))]
    fn new(
        metric_2c: PyReadonlyArray2<f64>,
        cderi_3c: PyReadonlyArray2<f64>,
        ia_p: PyReadonlyArray2<f64>,
        ij_p: PyReadonlyArray2<f64>,
        naux: usize,
        nao: usize,
        nocc: usize,
        nvir: usize,
        auxbasis: String,
        basis: String,
    ) -> PyResult<Self> {
        // Convert numpy arrays to ndarray::Array2
        let metric_2c_array = metric_2c.as_array().to_owned();
        let cderi_3c_array = cderi_3c.as_array().to_owned();
        let ia_p_array = ia_p.as_array().to_owned();
        let ij_p_array = ij_p.as_array().to_owned();

        // Create metadata
        let metadata = DFMetadata {
            naux,
            nao,
            nocc,
            nvir,
            auxbasis,
            basis,
            version: "1.0".to_string(), // S2.4 schema version 1.0
            timestamp: None,
            git_commit: None,
        };

        // Create DFTensorsS24 directly
        let df_tensors = DFTensorsS24 {
            metric_2c: metric_2c_array,
            cderi_3c: cderi_3c_array,
            ia_p: ia_p_array,
            ij_p: ij_p_array,
            metadata,
        };

        Ok(PyDFTensorsHDF5 { inner: df_tensors })
    }

    /// Save DF tensors to HDF5 file
    fn save_hdf5(&self, filename: String) -> PyResult<()> {
        use pyo3::exceptions::PyIOError;

        let path = PathBuf::from(filename);
        self.inner
            .write_hdf5(&path)
            .map_err(|e| PyIOError::new_err(format!("HDF5 write failed: {}", e)))
    }

    /// Load DF tensors from HDF5 file
    #[staticmethod]
    fn load_hdf5(filename: String) -> PyResult<Self> {
        use pyo3::exceptions::PyIOError;

        let path = PathBuf::from(filename);
        let df = DFTensorsS24::read_hdf5(&path)
            .map_err(|e| PyIOError::new_err(format!("HDF5 read failed: {}", e)))?;
        Ok(PyDFTensorsHDF5 { inner: df })
    }

    /// Get metric_2c as numpy array
    fn get_metric_2c<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_owned_array(py, self.inner.metric_2c.clone())
    }

    /// Get cderi_3c as numpy array
    fn get_cderi_3c<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_owned_array(py, self.inner.cderi_3c.clone())
    }

    /// Get ia_p as numpy array
    fn get_ia_p<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_owned_array(py, self.inner.ia_p.clone())
    }

    /// Get ij_p as numpy array
    fn get_ij_p<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_owned_array(py, self.inner.ij_p.clone())
    }

    /// Get metadata as Python dict
    fn get_metadata(&self, py: Python) -> PyResult<Py<PyAny>> {
        let dict = PyDict::new(py);
        dict.set_item("naux", self.inner.metadata.naux)?;
        dict.set_item("nao", self.inner.metadata.nao)?;
        dict.set_item("nocc", self.inner.metadata.nocc)?;
        dict.set_item("nvir", self.inner.metadata.nvir)?;
        dict.set_item("auxbasis", &self.inner.metadata.auxbasis)?;
        dict.set_item("basis", &self.inner.metadata.basis)?;
        dict.set_item("version", &self.inner.metadata.version)?;
        Ok(dict.into())
    }

    /// Get metadata object
    fn get_metadata_obj(&self) -> PyDFMetadata {
        PyDFMetadata {
            inner: self.inner.metadata.clone(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "DFTensorsHDF5(naux={}, nao={}, nocc={}, nvir={}, auxbasis='{}', basis='{}')",
            self.inner.metadata.naux,
            self.inner.metadata.nao,
            self.inner.metadata.nocc,
            self.inner.metadata.nvir,
            self.inner.metadata.auxbasis,
            self.inner.metadata.basis
        )
    }
}

/// Register HDF5 I/O classes with Python module
pub fn register_io_hdf5(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDFTensorsHDF5>()?;
    m.add_class::<PyDFMetadata>()?;
    Ok(())
}
