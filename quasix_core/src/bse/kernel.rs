//! Matrix-free BSE-TDA kernel implementation
//!
//! This module provides the core matrix-free BSE-TDA (Bethe-Salpeter Equation with
//! Tamm-Dancoff Approximation) kernel for computing optical excitations.
//!
//! # BSE-TDA Hamiltonian
//!
//! The BSE-TDA Hamiltonian in the singlet channel is:
//!
//! ```text
//! H^{BSE-TDA}_{ia,jb} = delta_ij * delta_ab * Delta_qp_{ia} + K^x_{ia,jb} + K^d_{ia,jb}
//! ```
//!
//! Where:
//! - `delta_ij`, `delta_ab` are Kronecker deltas
//! - `Delta_qp_{ia} = e_a^QP - e_i^QP` is the QP energy difference
//! - `K^x` is the exchange kernel (bare Coulomb)
//! - `K^d` is the direct kernel (screened Coulomb with W(omega=0))
//!
//! # Matrix-Free Philosophy
//!
//! The full BSE matrix has dimension (nocc * nvirt)^2 which can exceed 10^8 elements
//! for medium-sized molecules. Instead of storing the matrix explicitly, we implement
//! matrix-vector products y = H * x, enabling iterative eigensolvers (Davidson).
//!
//! # Memory Scaling
//!
//! - Full matrix: O(n_occ^2 * n_virt^2) - prohibitive for large systems
//! - Matrix-free: O(n_aux * n_occ * n_virt) - feasible with density fitting

#![allow(clippy::many_single_char_names)] // Mathematical notation (i,j,a,b,P,Q)

use crate::common::{QuasixError, Result};
use ndarray::{Array1, Array2, Axis};
use rayon::prelude::*;
use std::sync::Arc;

/// Module version for compatibility tracking
pub const MODULE_VERSION: &str = "0.6.1";

/// Numerical threshold for symmetry validation
pub const SYMMETRY_THRESHOLD: f64 = 1e-10;

/// Spin type for BSE calculations
///
/// Determines the prefactor for the exchange kernel:
/// - Singlet: 2 * K^x + K^d (spin-allowed optical transitions)
/// - Triplet: K^d only (spin-forbidden, no exchange contribution)
///
/// # Physical Interpretation
///
/// In singlet excitations (S=0), the electron and hole have opposite spins,
/// allowing the exchange interaction. In triplet excitations (S=1), they have
/// parallel spins, making the direct Coulomb term the only contribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SpinType {
    /// Singlet excitations (S=0, spin-allowed)
    /// Exchange kernel prefactor: 2.0
    #[default]
    Singlet,

    /// Triplet excitations (S=1, spin-forbidden)
    /// Exchange kernel prefactor: 0.0
    Triplet,
}

impl SpinType {
    /// Get the exchange kernel prefactor
    ///
    /// Returns 2.0 for singlet (accounts for spin summation) and 0.0 for triplet.
    #[inline]
    #[must_use]
    pub const fn exchange_prefactor(self) -> f64 {
        match self {
            Self::Singlet => 2.0,
            Self::Triplet => 0.0,
        }
    }
}

/// Configuration for BSE kernel calculations
///
/// Provides sensible defaults optimized for the target hardware
/// (2x Intel Xeon Silver 4314 with 24MB L3 cache per socket).
#[derive(Debug, Clone)]
pub struct BSEKernelConfig {
    /// Spin type (singlet or triplet)
    pub spin_type: SpinType,

    /// Use Tamm-Dancoff approximation (always true for this implementation)
    pub use_tda: bool,

    /// Number of parallel threads (defaults to available cores)
    pub n_threads: usize,

    /// Cache intermediate tensors for multiple applications
    pub cache_intermediates: bool,

    /// Block size for cache-optimal operations (L3 cache optimization)
    /// Default: 256 (optimized for 24MB L3 on Xeon Silver 4314)
    pub block_size: usize,

    /// Disable BLAS threading (use Rayon instead for outer parallelism)
    pub disable_blas_threads: bool,
}

impl Default for BSEKernelConfig {
    fn default() -> Self {
        Self {
            spin_type: SpinType::Singlet,
            use_tda: true,
            n_threads: rayon::current_num_threads(),
            cache_intermediates: true,
            block_size: 256,
            disable_blas_threads: true,
        }
    }
}

impl BSEKernelConfig {
    /// Configuration for small systems (< 100 transitions)
    ///
    /// Reduces block size and disables caching for minimal memory overhead.
    #[must_use]
    pub fn small_system() -> Self {
        Self {
            block_size: 64,
            cache_intermediates: false,
            ..Self::default()
        }
    }

    /// Configuration for large systems (> 10000 transitions)
    ///
    /// Increases block size and enables caching for better performance.
    #[must_use]
    pub fn large_system() -> Self {
        Self {
            block_size: 512,
            cache_intermediates: true,
            ..Self::default()
        }
    }

    /// Create configuration for triplet calculations
    #[must_use]
    pub fn triplet() -> Self {
        Self {
            spin_type: SpinType::Triplet,
            ..Self::default()
        }
    }
}

/// BSE-TDA kernel with matrix-free application
///
/// Implements the BSE Hamiltonian in TDA:
///
/// ```text
/// H^{BSE-TDA}_{ia,jb} = delta_ij * delta_ab * Delta_qp_{ia} + K^x_{ia,jb} + K^d_{ia,jb}
/// ```
///
/// Where:
/// - `Delta_qp_{ia} = e^{QP}_a - e^{QP}_i` (QP energy difference from evGW)
/// - `K^x_{ia,jb} = (ia|jb)` (exchange kernel, bare Coulomb)
/// - `K^d_{ia,jb} = -(ij|W|ab)` (direct kernel, screened Coulomb at omega=0)
///
/// # Matrix-Free Philosophy
///
/// Instead of forming the full (nocc*nvirt)^2 matrix, we compute matrix-vector
/// products y = H * x directly. This enables:
/// 1. Memory efficiency: O(n_aux * n_trans) instead of O(n_trans^2)
/// 2. Iterative eigensolvers (Davidson) for lowest eigenvalues
/// 3. Parallelization via Rayon for large systems
///
/// # References
///
/// - PySCF: pyscf/tdscf/bse.py
/// - Theory: docs/derivations/s6-1/theory.md
#[derive(Debug)]
pub struct BSETDAKernel {
    /// Number of occupied orbitals
    pub nocc: usize,

    /// Number of virtual orbitals
    pub nvirt: usize,

    /// Number of auxiliary basis functions
    pub naux: usize,

    /// Configuration
    pub config: BSEKernelConfig,

    /// QP energy differences: delta_qp[i*nvirt + a] = e_a^QP - e_i^QP
    delta_qp: Array1<f64>,

    /// DF tensor for occ-virt transitions: (ia|P), shape [nocc*nvirt, naux]
    /// Used for exchange kernel: K^x = B_ia @ B_jb.T
    df_ia: Arc<Array2<f64>>,

    /// DF tensor for occ-occ pairs: (ij|P), shape [nocc*nocc, naux]
    /// Used for direct kernel half-transform
    df_ij: Arc<Array2<f64>>,

    /// DF tensor for virt-virt pairs: (ab|P), shape [nvirt*nvirt, naux]
    /// Used for direct kernel half-transform
    df_ab: Arc<Array2<f64>>,

    /// Static screened interaction W(omega=0), shape [naux, naux]
    /// From S3-3: W = v [I - v P_0(0)]^{-1}
    w0: Arc<Array2<f64>>,

    /// Cached intermediate: B_ij @ W0, shape [nocc*nocc, naux]
    /// Pre-computed if config.cache_intermediates is true
    df_ij_w0: Option<Array2<f64>>,
}

impl BSETDAKernel {
    /// Create a new BSE-TDA kernel
    ///
    /// # Arguments
    ///
    /// * `nocc` - Number of occupied orbitals
    /// * `nvirt` - Number of virtual orbitals
    /// * `naux` - Number of auxiliary basis functions
    /// * `delta_qp` - QP energy differences, shape [nocc * nvirt]
    /// * `df_ia` - DF tensor (ia|P), shape [nocc*nvirt, naux]
    /// * `df_ij` - DF tensor (ij|P), shape [nocc*nocc, naux]
    /// * `df_ab` - DF tensor (ab|P), shape [nvirt*nvirt, naux]
    /// * `w0` - Static screened interaction W(0), shape [naux, naux]
    /// * `config` - Kernel configuration
    ///
    /// # Errors
    ///
    /// Returns error if dimensions are inconsistent or inputs contain non-finite values.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use quasix_core::bse::{BSETDAKernel, BSEKernelConfig};
    /// use ndarray::{Array1, Array2};
    ///
    /// let nocc = 2;
    /// let nvirt = 3;
    /// let naux = 10;
    ///
    /// let delta_qp = Array1::from_vec(vec![0.5; nocc * nvirt]);
    /// let df_ia = Array2::zeros((nocc * nvirt, naux));
    /// let df_ij = Array2::zeros((nocc * nocc, naux));
    /// let df_ab = Array2::zeros((nvirt * nvirt, naux));
    /// let w0 = Array2::eye(naux);
    ///
    /// let kernel = BSETDAKernel::new(
    ///     nocc, nvirt, naux,
    ///     delta_qp, df_ia, df_ij, df_ab, w0,
    ///     BSEKernelConfig::default()
    /// )?;
    /// ```
    pub fn new(
        nocc: usize,
        nvirt: usize,
        naux: usize,
        delta_qp: Array1<f64>,
        df_ia: Array2<f64>,
        df_ij: Array2<f64>,
        df_ab: Array2<f64>,
        w0: Array2<f64>,
        config: BSEKernelConfig,
    ) -> Result<Self> {
        // Validate dimensions
        validate_bse_inputs(nocc, nvirt, naux, &delta_qp, &df_ia, &df_ij, &df_ab, &w0)?;

        // Disable BLAS threading if configured (Rayon handles parallelism)
        if config.disable_blas_threads {
            // Note: This affects the process globally; consider using scoped settings
            std::env::set_var("OPENBLAS_NUM_THREADS", "1");
            std::env::set_var("MKL_NUM_THREADS", "1");
        }

        let df_ia = Arc::new(df_ia);
        let df_ij = Arc::new(df_ij);
        let df_ab = Arc::new(df_ab);
        let w0 = Arc::new(w0);

        // Pre-compute cached intermediate if enabled
        // B_ij @ W0: [nocc*nocc, naux] @ [naux, naux] -> [nocc*nocc, naux]
        let df_ij_w0 = if config.cache_intermediates {
            Some(df_ij.dot(&*w0))
        } else {
            None
        };

        Ok(Self {
            nocc,
            nvirt,
            naux,
            config,
            delta_qp,
            df_ia,
            df_ij,
            df_ab,
            w0,
            df_ij_w0,
        })
    }

    /// Get the BSE Hamiltonian dimension (nocc * nvirt)
    #[inline]
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.nocc * self.nvirt
    }

    /// Apply exchange kernel to a trial vector
    ///
    /// Computes: y^{Kx}_{ia} = sum_{jb} (ia|jb) * x_{jb}
    ///
    /// Using density fitting:
    /// ```text
    /// y^{Kx}_{ia} = sum_{jb} sum_P (ia|P) (jb|P) * x_{jb}
    ///            = sum_P (ia|P) * [sum_{jb} (jb|P) * x_{jb}]
    ///            = sum_P (ia|P) * z_P
    /// ```
    ///
    /// Where `z_P = sum_{jb} (jb|P) * x_{jb} = B_jb^T @ x`
    ///
    /// # Algorithm (Two-Step BLAS)
    ///
    /// 1. `z = B.T @ x`  => [naux] = [naux, n_trans].T @ [n_trans]
    /// 2. `y = B @ z`    => [n_trans] = [n_trans, naux] @ [naux]
    ///
    /// # Complexity
    ///
    /// - Time: O(2 * n_trans * n_aux)
    /// - Memory: O(n_aux) for intermediate z
    ///
    /// # Arguments
    ///
    /// * `x` - Trial vector, shape [nocc * nvirt]
    ///
    /// # Returns
    ///
    /// Exchange contribution y_kx, shape [nocc * nvirt]
    pub fn apply_exchange_kernel(&self, x: &Array1<f64>) -> Result<Array1<f64>> {
        if x.len() != self.dimension() {
            return Err(QuasixError::DimensionMismatch(format!(
                "Input vector length {} != BSE dimension {}",
                x.len(),
                self.dimension()
            )));
        }

        // Step 1: z_P = sum_{jb} B_{jb,P} * x_{jb}
        // z = B.T @ x  =>  [naux] = [naux, n_trans] @ [n_trans]
        let z = self.df_ia.t().dot(x);

        // Step 2: y_{ia} = sum_P B_{ia,P} * z_P
        // y = B @ z  =>  [n_trans] = [n_trans, naux] @ [naux]
        let y_kx = self.df_ia.dot(&z);

        Ok(y_kx)
    }

    /// Apply direct (screened) kernel to a trial vector
    ///
    /// Computes: y^{Kd}_{ia} = -sum_{jb} (ij|W|ab) * x_{jb}
    ///
    /// The direct kernel involves the screened Coulomb interaction W(omega=0):
    /// ```text
    /// K^d_{ia,jb} = -(ij|W(0)|ab) = -sum_{PQ} (ij|P) W_{PQ} (ab|Q)
    /// ```
    ///
    /// # Algorithm
    ///
    /// For each occupied pair (i,j) and virtual pair (a,b):
    /// 1. Compute `(ij|W|ab) = sum_{PQ} B_{ij,P} W_{PQ} B_{ab,Q}`
    /// 2. Contract with `x_{jb}`
    ///
    /// # Complexity
    ///
    /// - Time: O(nocc^2 * nvirt^2 * naux) - parallelized over i
    /// - Memory: O(naux) per thread
    ///
    /// # Arguments
    ///
    /// * `x` - Trial vector, shape [nocc * nvirt]
    ///
    /// # Returns
    ///
    /// Direct kernel contribution y_kd, shape [nocc * nvirt]
    pub fn apply_direct_kernel(&self, x: &Array1<f64>) -> Result<Array1<f64>> {
        if x.len() != self.dimension() {
            return Err(QuasixError::DimensionMismatch(format!(
                "Input vector length {} != BSE dimension {}",
                x.len(),
                self.dimension()
            )));
        }

        let nocc = self.nocc;
        let nvirt = self.nvirt;

        // Reshape x from [nocc*nvirt] to view as (j, b) indexing
        // x[j*nvirt + b] corresponds to transition j -> b

        // Get the B_ij @ W intermediate (cached or computed)
        let df_ij_w0 = match &self.df_ij_w0 {
            Some(cached) => cached.clone(),
            None => self.df_ij.dot(&*self.w0),
        };

        // Process in parallel over occupied index i
        let y_slices: Vec<Array1<f64>> = (0..nocc)
            .into_par_iter()
            .map(|i| {
                let mut y_i = Array1::<f64>::zeros(nvirt);

                for j in 0..nocc {
                    // Index for (i,j) pair: i*nocc + j
                    let ij_idx = i * nocc + j;

                    // Get [B_ij @ W]_Q for this (i,j) pair
                    let b_ij_w = df_ij_w0.row(ij_idx);

                    for a in 0..nvirt {
                        for b in 0..nvirt {
                            // Index for (a,b) pair: a*nvirt + b
                            let ab_idx = a * nvirt + b;

                            // Get B_{ab,Q}
                            let b_ab = self.df_ab.row(ab_idx);

                            // Compute (ij|W|ab) = sum_Q [B_ij @ W]_Q * B_{ab,Q}
                            let ij_w_ab: f64 =
                                b_ij_w.iter().zip(b_ab.iter()).map(|(&w, &b)| w * b).sum();

                            // Get x_{jb} = x[j*nvirt + b]
                            let x_jb = x[j * nvirt + b];

                            // Accumulate: y[i,a] += (ij|W|ab) * x[j,b]
                            y_i[a] += ij_w_ab * x_jb;
                        }
                    }
                }

                y_i
            })
            .collect();

        // Combine results into flat array
        let mut y = Array1::<f64>::zeros(nocc * nvirt);
        for (i, y_i) in y_slices.into_iter().enumerate() {
            let start = i * nvirt;
            for (a, &val) in y_i.iter().enumerate() {
                y[start + a] = val;
            }
        }

        // Apply negative sign for direct kernel: K^d = -(ij|W|ab)
        y.mapv_inplace(|v| -v);

        Ok(y)
    }

    /// Apply BSE-TDA Hamiltonian to a trial vector (matrix-free)
    ///
    /// Computes: y = H^{BSE-TDA} * x
    ///
    /// Where:
    /// ```text
    /// H^{BSE-TDA}_{ia,jb} = delta_ij * delta_ab * Delta_qp_{ia}
    ///                    + alpha_x * K^x_{ia,jb}
    ///                    + K^d_{ia,jb}
    /// ```
    ///
    /// - `alpha_x = 2.0` for singlet, `0.0` for triplet
    /// - `K^d` already includes the negative sign
    ///
    /// # Algorithm
    ///
    /// 1. Diagonal contribution: `y = Delta_qp * x` (element-wise)
    /// 2. Exchange contribution: `y += alpha_x * K^x(x)`
    /// 3. Direct contribution: `y += K^d(x)` (K^d already has negative sign)
    ///
    /// # Complexity
    ///
    /// - Time: O(n_trans * n_aux) for exchange + O(nocc^2 * nvirt^2 * naux) for direct
    /// - Memory: O(n_aux) temporary workspace
    ///
    /// # Arguments
    ///
    /// * `x` - Trial vector, shape [nocc * nvirt]
    ///
    /// # Returns
    ///
    /// Result vector y = H * x, shape [nocc * nvirt]
    ///
    /// # Errors
    ///
    /// Returns error if input dimension is wrong or computation produces non-finite values.
    pub fn apply_tda_hamiltonian(&self, x: &Array1<f64>) -> Result<Array1<f64>> {
        // Validate input dimension
        if x.len() != self.dimension() {
            return Err(QuasixError::DimensionMismatch(format!(
                "Input vector length {} != BSE dimension {}",
                x.len(),
                self.dimension()
            )));
        }

        // Step 1: Diagonal contribution (QP energy differences)
        // y = Delta_qp * x (element-wise multiplication)
        let mut y = &self.delta_qp * x;

        // Step 2: Exchange kernel contribution
        // y += alpha_x * K^x(x)
        let alpha_x = self.config.spin_type.exchange_prefactor();
        if alpha_x.abs() > f64::EPSILON {
            let y_kx = self.apply_exchange_kernel(x)?;
            y = y + alpha_x * &y_kx;
        }

        // Step 3: Direct kernel contribution
        // y += K^d(x) -- note: K^d already returns -(ij|W|ab)*x
        let y_kd = self.apply_direct_kernel(x)?;
        y += &y_kd;

        // Validate output
        if y.iter().any(|&v| !v.is_finite()) {
            return Err(QuasixError::NumericalError(
                "BSE Hamiltonian application produced non-finite values".to_string(),
            ));
        }

        Ok(y)
    }

    /// Apply BSE-TDA Hamiltonian to multiple trial vectors
    ///
    /// More efficient than applying to each vector separately due to
    /// potential BLAS-3 batching in the exchange kernel.
    ///
    /// # Arguments
    ///
    /// * `x_batch` - Trial vectors, shape [nocc * nvirt, n_vectors]
    ///
    /// # Returns
    ///
    /// Result vectors Y = H * X, shape [nocc * nvirt, n_vectors]
    pub fn apply_tda_hamiltonian_batch(&self, x_batch: &Array2<f64>) -> Result<Array2<f64>> {
        let (n_trans, n_vectors) = x_batch.dim();

        if n_trans != self.dimension() {
            return Err(QuasixError::DimensionMismatch(format!(
                "Input matrix rows {} != BSE dimension {}",
                n_trans,
                self.dimension()
            )));
        }

        // Step 1: Diagonal contribution (broadcasted)
        // Y = Delta_qp[:, None] * X
        let delta_qp_col = self.delta_qp.clone().insert_axis(Axis(1));
        let mut y = &delta_qp_col * x_batch;

        // Step 2: Exchange kernel (batched BLAS-3)
        let alpha_x = self.config.spin_type.exchange_prefactor();
        if alpha_x.abs() > f64::EPSILON {
            // Z = B.T @ X  =>  [naux, n_vectors]
            let z = self.df_ia.t().dot(x_batch);
            // Y_kx = B @ Z  =>  [n_trans, n_vectors]
            let y_kx = self.df_ia.dot(&z);
            y = y + alpha_x * &y_kx;
        }

        // Step 3: Direct kernel (apply to each vector)
        // TODO: Optimize with BLAS-3 batched version
        let y_kd_results: Result<Vec<Array1<f64>>> = (0..n_vectors)
            .into_par_iter()
            .map(|v| {
                let x = x_batch.column(v).to_owned();
                self.apply_direct_kernel(&x)
            })
            .collect();

        let y_kd_vecs = y_kd_results?;

        // Combine direct kernel results
        for (v, y_kd) in y_kd_vecs.into_iter().enumerate() {
            for (i, &val) in y_kd.iter().enumerate() {
                y[[i, v]] += val;
            }
        }

        Ok(y)
    }

    /// Get reference to QP energy differences
    #[must_use]
    pub fn delta_qp(&self) -> &Array1<f64> {
        &self.delta_qp
    }

    /// Get reference to the DF tensor (ia|P)
    #[must_use]
    pub fn df_ia(&self) -> &Array2<f64> {
        &self.df_ia
    }

    /// Get reference to the DF tensor (ij|P)
    #[must_use]
    pub fn df_ij(&self) -> &Array2<f64> {
        &self.df_ij
    }

    /// Get reference to the DF tensor (ab|P)
    #[must_use]
    pub fn df_ab(&self) -> &Array2<f64> {
        &self.df_ab
    }

    /// Get reference to the static screened interaction W(0)
    #[must_use]
    pub fn w0(&self) -> &Array2<f64> {
        &self.w0
    }

    /// Solve BSE-TDA eigenvalue problem using Davidson iteration
    ///
    /// Finds the lowest n_roots excitation energies and corresponding
    /// eigenvectors using the matrix-free Davidson algorithm.
    ///
    /// # Arguments
    ///
    /// * `config` - Davidson solver configuration
    ///
    /// # Returns
    ///
    /// `DavidsonResult` containing eigenvalues, eigenvectors, and convergence info
    ///
    /// # Example
    ///
    /// ```ignore
    /// use quasix_core::bse::{BSETDAKernel, DavidsonConfig};
    ///
    /// let kernel = // ... create kernel ...
    /// let config = DavidsonConfig::with_n_roots(5);
    /// let result = kernel.solve_davidson(&config)?;
    ///
    /// println!("Lowest excitation energy: {} Ha", result.eigenvalues[0]);
    /// ```
    pub fn solve_davidson(&self, config: &super::DavidsonConfig) -> Result<super::DavidsonResult> {
        // Create closure that applies H to a batch of vectors
        let apply_h =
            |v: &Array2<f64>| -> Result<Array2<f64>> { self.apply_tda_hamiltonian_batch(v) };

        super::davidson_bse(apply_h, self.delta_qp.view(), config)
    }

    /// Solve BSE-TDA with default configuration
    ///
    /// Convenience method using default DavidsonConfig.
    ///
    /// # Arguments
    ///
    /// * `n_roots` - Number of lowest eigenvalues to compute
    ///
    /// # Returns
    ///
    /// `DavidsonResult` containing eigenvalues, eigenvectors, and convergence info
    pub fn solve(&self, n_roots: usize) -> Result<super::DavidsonResult> {
        let config = super::DavidsonConfig::with_n_roots(n_roots);
        self.solve_davidson(&config)
    }
}

/// Validate all BSE input tensors
///
/// Checks:
/// - Non-zero dimensions
/// - Consistent shapes
/// - Finite values
/// - W0 symmetry
fn validate_bse_inputs(
    nocc: usize,
    nvirt: usize,
    naux: usize,
    delta_qp: &Array1<f64>,
    df_ia: &Array2<f64>,
    df_ij: &Array2<f64>,
    df_ab: &Array2<f64>,
    w0: &Array2<f64>,
) -> Result<()> {
    let n_trans = nocc * nvirt;

    // Check non-zero dimensions
    if nocc == 0 {
        return Err(QuasixError::InvalidInput(
            "Number of occupied orbitals must be > 0".to_string(),
        ));
    }
    if nvirt == 0 {
        return Err(QuasixError::InvalidInput(
            "Number of virtual orbitals must be > 0".to_string(),
        ));
    }
    if naux == 0 {
        return Err(QuasixError::InvalidInput(
            "Number of auxiliary functions must be > 0".to_string(),
        ));
    }

    // Check delta_qp
    if delta_qp.len() != n_trans {
        return Err(QuasixError::DimensionMismatch(format!(
            "delta_qp length {} != nocc*nvirt = {}",
            delta_qp.len(),
            n_trans
        )));
    }

    // Check for finite values in delta_qp
    if delta_qp.iter().any(|&v| !v.is_finite()) {
        return Err(QuasixError::NumericalError(
            "delta_qp contains non-finite values".to_string(),
        ));
    }

    // Check df_ia dimensions
    let (df_ia_rows, df_ia_cols) = df_ia.dim();
    if df_ia_rows != n_trans || df_ia_cols != naux {
        return Err(QuasixError::DimensionMismatch(format!(
            "df_ia shape ({}, {}) != expected ({}, {})",
            df_ia_rows, df_ia_cols, n_trans, naux
        )));
    }

    // Check df_ij dimensions
    let (df_ij_rows, df_ij_cols) = df_ij.dim();
    if df_ij_rows != nocc * nocc || df_ij_cols != naux {
        return Err(QuasixError::DimensionMismatch(format!(
            "df_ij shape ({}, {}) != expected ({}, {})",
            df_ij_rows,
            df_ij_cols,
            nocc * nocc,
            naux
        )));
    }

    // Check df_ab dimensions
    let (df_ab_rows, df_ab_cols) = df_ab.dim();
    if df_ab_rows != nvirt * nvirt || df_ab_cols != naux {
        return Err(QuasixError::DimensionMismatch(format!(
            "df_ab shape ({}, {}) != expected ({}, {})",
            df_ab_rows,
            df_ab_cols,
            nvirt * nvirt,
            naux
        )));
    }

    // Check W0 dimensions
    let (w_rows, w_cols) = w0.dim();
    if w_rows != naux || w_cols != naux {
        return Err(QuasixError::DimensionMismatch(format!(
            "W0 shape ({}, {}) != expected ({}, {})",
            w_rows, w_cols, naux, naux
        )));
    }

    // Check W0 symmetry
    for i in 0..naux {
        for j in (i + 1)..naux {
            let diff = (w0[[i, j]] - w0[[j, i]]).abs();
            if diff > SYMMETRY_THRESHOLD {
                return Err(QuasixError::NumericalError(format!(
                    "W0 not symmetric: W[{},{}]={:.6e} != W[{},{}]={:.6e}, diff={:.6e}",
                    i,
                    j,
                    w0[[i, j]],
                    j,
                    i,
                    w0[[j, i]],
                    diff
                )));
            }
        }
    }

    // Check for finite values in DF tensors
    if df_ia.iter().any(|&v| !v.is_finite()) {
        return Err(QuasixError::NumericalError(
            "df_ia contains non-finite values".to_string(),
        ));
    }
    if df_ij.iter().any(|&v| !v.is_finite()) {
        return Err(QuasixError::NumericalError(
            "df_ij contains non-finite values".to_string(),
        ));
    }
    if df_ab.iter().any(|&v| !v.is_finite()) {
        return Err(QuasixError::NumericalError(
            "df_ab contains non-finite values".to_string(),
        ));
    }
    if w0.iter().any(|&v| !v.is_finite()) {
        return Err(QuasixError::NumericalError(
            "W0 contains non-finite values".to_string(),
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::{arr1, arr2, Array1, Array2};

    /// Create a simple test kernel for minimal system
    fn create_test_kernel() -> BSETDAKernel {
        let nocc = 1;
        let nvirt = 1;
        let naux = 2;

        // Simple QP gaps
        let delta_qp = arr1(&[0.5]); // 0.5 Ha gap

        // Simple DF tensors
        let df_ia = arr2(&[[0.5, 0.3]]);
        let df_ij = arr2(&[[0.4, 0.2]]);
        let df_ab = arr2(&[[0.6, 0.1]]);

        // Simple W0 (symmetric, identity-like)
        let w0 = arr2(&[[1.0, 0.1], [0.1, 1.0]]);

        let config = BSEKernelConfig::default();

        BSETDAKernel::new(nocc, nvirt, naux, delta_qp, df_ia, df_ij, df_ab, w0, config)
            .expect("Failed to create test kernel")
    }

    /// Create a larger test kernel for more comprehensive tests
    fn create_larger_test_kernel() -> BSETDAKernel {
        let nocc = 2;
        let nvirt = 3;
        let naux = 4;
        let n_trans = nocc * nvirt;

        // QP energy differences
        let delta_qp = Array1::from_vec(vec![0.3, 0.4, 0.5, 0.35, 0.45, 0.55]);

        // DF tensor for occ-virt: [6, 4]
        let mut df_ia = Array2::zeros((n_trans, naux));
        for ia in 0..n_trans {
            for p in 0..naux {
                df_ia[[ia, p]] = 0.1 * ((ia + p + 1) as f64).sqrt();
            }
        }

        // DF tensor for occ-occ: [4, 4]
        let mut df_ij = Array2::zeros((nocc * nocc, naux));
        for ij in 0..(nocc * nocc) {
            for p in 0..naux {
                df_ij[[ij, p]] = 0.15 * ((ij + p + 2) as f64).sqrt();
            }
        }

        // DF tensor for virt-virt: [9, 4]
        let mut df_ab = Array2::zeros((nvirt * nvirt, naux));
        for ab in 0..(nvirt * nvirt) {
            for p in 0..naux {
                df_ab[[ab, p]] = 0.12 * ((ab + p + 1) as f64).sqrt();
            }
        }

        // Symmetric W0
        let mut w0 = Array2::eye(naux);
        for i in 0..naux {
            for j in 0..naux {
                if i != j {
                    w0[[i, j]] = 0.05 / ((i as f64 - j as f64).abs() + 1.0);
                    w0[[j, i]] = w0[[i, j]];
                }
            }
        }

        let config = BSEKernelConfig::default();

        BSETDAKernel::new(nocc, nvirt, naux, delta_qp, df_ia, df_ij, df_ab, w0, config)
            .expect("Failed to create larger test kernel")
    }

    #[test]
    fn test_spin_type_prefactor() {
        assert_relative_eq!(SpinType::Singlet.exchange_prefactor(), 2.0);
        assert_relative_eq!(SpinType::Triplet.exchange_prefactor(), 0.0);
    }

    #[test]
    fn test_spin_type_default() {
        assert_eq!(SpinType::default(), SpinType::Singlet);
    }

    #[test]
    fn test_config_default() {
        let config = BSEKernelConfig::default();
        assert_eq!(config.spin_type, SpinType::Singlet);
        assert!(config.use_tda);
        assert!(config.cache_intermediates);
        assert_eq!(config.block_size, 256);
    }

    #[test]
    fn test_config_small_system() {
        let config = BSEKernelConfig::small_system();
        assert_eq!(config.block_size, 64);
        assert!(!config.cache_intermediates);
    }

    #[test]
    fn test_config_large_system() {
        let config = BSEKernelConfig::large_system();
        assert_eq!(config.block_size, 512);
        assert!(config.cache_intermediates);
    }

    #[test]
    fn test_config_triplet() {
        let config = BSEKernelConfig::triplet();
        assert_eq!(config.spin_type, SpinType::Triplet);
    }

    #[test]
    fn test_kernel_creation() {
        let kernel = create_test_kernel();
        assert_eq!(kernel.nocc, 1);
        assert_eq!(kernel.nvirt, 1);
        assert_eq!(kernel.naux, 2);
        assert_eq!(kernel.dimension(), 1);
    }

    #[test]
    fn test_larger_kernel_creation() {
        let kernel = create_larger_test_kernel();
        assert_eq!(kernel.nocc, 2);
        assert_eq!(kernel.nvirt, 3);
        assert_eq!(kernel.dimension(), 6);
    }

    #[test]
    fn test_exchange_kernel_simple() {
        let kernel = create_test_kernel();
        let x = arr1(&[1.0]);

        let y_kx = kernel.apply_exchange_kernel(&x).unwrap();

        // K^x = B_ia @ B_ia.T @ x
        // B_ia = [[0.5, 0.3]], B_ia @ B_ia.T = [[0.5*0.5 + 0.3*0.3]] = [[0.34]]
        // y = 0.34 * 1.0 = 0.34
        assert_relative_eq!(y_kx[0], 0.34, epsilon = 1e-10);
    }

    #[test]
    fn test_direct_kernel_simple() {
        let kernel = create_test_kernel();
        let x = arr1(&[1.0]);

        let y_kd = kernel.apply_direct_kernel(&x).unwrap();

        // For 1x1 system with nocc=1, nvirt=1:
        // K^d_{0,0} = -(00|W|00) = -sum_{PQ} B_{00,P} W_{PQ} B_{00,Q}
        // df_ij = [[0.4, 0.2]], df_ab = [[0.6, 0.1]]
        // W0 = [[1.0, 0.1], [0.1, 1.0]]
        //
        // B_ij @ W = [0.4, 0.2] @ [[1.0, 0.1], [0.1, 1.0]]
        //          = [0.4*1.0 + 0.2*0.1, 0.4*0.1 + 0.2*1.0]
        //          = [0.42, 0.24]
        //
        // (00|W|00) = [0.42, 0.24] @ [0.6, 0.1]
        //           = 0.42*0.6 + 0.24*0.1
        //           = 0.252 + 0.024 = 0.276
        //
        // K^d = -0.276 * 1.0 = -0.276
        assert_relative_eq!(y_kd[0], -0.276, epsilon = 1e-10);
    }

    #[test]
    fn test_hamiltonian_application() {
        let kernel = create_test_kernel();
        let x = arr1(&[1.0]);

        let y = kernel.apply_tda_hamiltonian(&x).unwrap();

        // For singlet:
        // y = delta_qp * x + 2 * K^x(x) + K^d(x)
        // y = 0.5 * 1.0 + 2 * 0.34 + (-0.276)
        // y = 0.5 + 0.68 - 0.276 = 0.904
        assert_relative_eq!(y[0], 0.904, epsilon = 1e-10);
        assert!(y[0].is_finite());
    }

    #[test]
    fn test_triplet_no_exchange() {
        let nocc = 1;
        let nvirt = 1;
        let naux = 2;

        let delta_qp = arr1(&[0.5]);
        let df_ia = arr2(&[[0.5, 0.3]]);
        let df_ij = arr2(&[[0.4, 0.2]]);
        let df_ab = arr2(&[[0.6, 0.1]]);
        let w0 = arr2(&[[1.0, 0.1], [0.1, 1.0]]);

        let config = BSEKernelConfig::triplet();
        let kernel =
            BSETDAKernel::new(nocc, nvirt, naux, delta_qp, df_ia, df_ij, df_ab, w0, config)
                .unwrap();

        let x = arr1(&[1.0]);
        let y = kernel.apply_tda_hamiltonian(&x).unwrap();

        // For triplet: y = delta_qp * x + 0 * K^x(x) + K^d(x)
        // y = 0.5 + 0 + (-0.276) = 0.224
        assert_relative_eq!(y[0], 0.224, epsilon = 1e-10);
    }

    #[test]
    fn test_exchange_kernel_symmetry() {
        let kernel = create_larger_test_kernel();
        let n_trans = kernel.dimension();

        // Apply to unit vectors and check symmetry
        // K^x_{ia,jb} should equal K^x_{jb,ia}
        for ia in 0..n_trans {
            for jb in ia..n_trans {
                let mut x_ia = Array1::zeros(n_trans);
                let mut x_jb = Array1::zeros(n_trans);
                x_ia[ia] = 1.0;
                x_jb[jb] = 1.0;

                let y_ia = kernel.apply_exchange_kernel(&x_ia).unwrap();
                let y_jb = kernel.apply_exchange_kernel(&x_jb).unwrap();

                // K^x_{ia,jb} = y_jb[ia] should equal K^x_{jb,ia} = y_ia[jb]
                assert_relative_eq!(y_jb[ia], y_ia[jb], epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_dimension_validation_delta_qp() {
        let nocc = 2;
        let nvirt = 3;
        let naux = 4;

        // Wrong delta_qp length
        let bad_delta_qp = Array1::zeros(5); // Should be 6
        let df_ia = Array2::zeros((6, 4));
        let df_ij = Array2::zeros((4, 4));
        let df_ab = Array2::zeros((9, 4));
        let w0 = Array2::eye(4);

        let result = BSETDAKernel::new(
            nocc,
            nvirt,
            naux,
            bad_delta_qp,
            df_ia,
            df_ij,
            df_ab,
            w0,
            BSEKernelConfig::default(),
        );

        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, QuasixError::DimensionMismatch(_)));
    }

    #[test]
    fn test_dimension_validation_df_ia() {
        let nocc = 2;
        let nvirt = 3;
        let naux = 4;

        let delta_qp = Array1::zeros(6);
        let bad_df_ia = Array2::zeros((5, 4)); // Should be (6, 4)
        let df_ij = Array2::zeros((4, 4));
        let df_ab = Array2::zeros((9, 4));
        let w0 = Array2::eye(4);

        let result = BSETDAKernel::new(
            nocc,
            nvirt,
            naux,
            delta_qp,
            bad_df_ia,
            df_ij,
            df_ab,
            w0,
            BSEKernelConfig::default(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_dimension_validation_w0() {
        let nocc = 2;
        let nvirt = 3;
        let naux = 4;

        let delta_qp = Array1::zeros(6);
        let df_ia = Array2::zeros((6, 4));
        let df_ij = Array2::zeros((4, 4));
        let df_ab = Array2::zeros((9, 4));
        let bad_w0 = Array2::eye(3); // Should be (4, 4)

        let result = BSETDAKernel::new(
            nocc,
            nvirt,
            naux,
            delta_qp,
            df_ia,
            df_ij,
            df_ab,
            bad_w0,
            BSEKernelConfig::default(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_w0_symmetry_validation() {
        let nocc = 1;
        let nvirt = 1;
        let naux = 2;

        let delta_qp = arr1(&[0.5]);
        let df_ia = arr2(&[[0.5, 0.3]]);
        let df_ij = arr2(&[[0.4, 0.2]]);
        let df_ab = arr2(&[[0.6, 0.1]]);

        // Non-symmetric W0
        let bad_w0 = arr2(&[[1.0, 0.5], [0.1, 1.0]]); // 0.5 != 0.1

        let result = BSETDAKernel::new(
            nocc,
            nvirt,
            naux,
            delta_qp,
            df_ia,
            df_ij,
            df_ab,
            bad_w0,
            BSEKernelConfig::default(),
        );

        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, QuasixError::NumericalError(_)));
    }

    #[test]
    fn test_zero_dimension_validation() {
        let delta_qp = Array1::zeros(0);
        let df_ia = Array2::zeros((0, 2));
        let df_ij = Array2::zeros((0, 2));
        let df_ab = Array2::zeros((0, 2));
        let w0 = Array2::eye(2);

        let result = BSETDAKernel::new(
            0,
            1,
            2,
            delta_qp,
            df_ia,
            df_ij,
            df_ab,
            w0,
            BSEKernelConfig::default(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_hamiltonian_batch_application() {
        let kernel = create_larger_test_kernel();
        let n_trans = kernel.dimension();

        // Create batch of trial vectors
        let mut x_batch = Array2::zeros((n_trans, 3));
        for v in 0..3 {
            for i in 0..n_trans {
                x_batch[[i, v]] = ((i + v + 1) as f64) * 0.1;
            }
        }

        let y_batch = kernel.apply_tda_hamiltonian_batch(&x_batch).unwrap();

        // Verify against individual applications
        for v in 0..3 {
            let x_single = x_batch.column(v).to_owned();
            let y_single = kernel.apply_tda_hamiltonian(&x_single).unwrap();

            for i in 0..n_trans {
                assert_relative_eq!(y_batch[[i, v]], y_single[i], epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_input_vector_dimension_mismatch() {
        let kernel = create_larger_test_kernel();
        let wrong_size_x = Array1::zeros(10); // Should be 6

        assert!(kernel.apply_exchange_kernel(&wrong_size_x).is_err());
        assert!(kernel.apply_direct_kernel(&wrong_size_x).is_err());
        assert!(kernel.apply_tda_hamiltonian(&wrong_size_x).is_err());
    }

    #[test]
    fn test_non_finite_values_rejected() {
        let nocc = 1;
        let nvirt = 1;
        let naux = 2;

        let bad_delta_qp = arr1(&[f64::NAN]);
        let df_ia = arr2(&[[0.5, 0.3]]);
        let df_ij = arr2(&[[0.4, 0.2]]);
        let df_ab = arr2(&[[0.6, 0.1]]);
        let w0 = arr2(&[[1.0, 0.1], [0.1, 1.0]]);

        let result = BSETDAKernel::new(
            nocc,
            nvirt,
            naux,
            bad_delta_qp,
            df_ia,
            df_ij,
            df_ab,
            w0,
            BSEKernelConfig::default(),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_accessor_methods() {
        let kernel = create_test_kernel();

        assert_eq!(kernel.delta_qp().len(), 1);
        assert_eq!(kernel.df_ia().dim(), (1, 2));
        assert_eq!(kernel.df_ij().dim(), (1, 2));
        assert_eq!(kernel.df_ab().dim(), (1, 2));
        assert_eq!(kernel.w0().dim(), (2, 2));
    }

    #[test]
    fn test_positive_excitation_energies() {
        // For a physical system, the BSE eigenvalues should be positive
        // This test verifies the Hamiltonian produces reasonable values
        let kernel = create_larger_test_kernel();
        let n_trans = kernel.dimension();

        // Random normalized vector
        let mut x: Array1<f64> = (0..n_trans).map(|i| (i as f64 + 1.0).sqrt()).collect();
        let norm: f64 = x.iter().map(|&v| v * v).sum::<f64>().sqrt();
        x.mapv_inplace(|v| v / norm);

        let y = kernel.apply_tda_hamiltonian(&x).unwrap();

        // Rayleigh quotient: <x|H|x> should give an estimate of eigenvalue
        let rayleigh: f64 = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum();

        // Should be positive for physical excitations
        assert!(rayleigh > 0.0, "Rayleigh quotient should be positive");
    }
}
