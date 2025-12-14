//! Davidson eigenvalue solver for BSE-TDA
//!
//! This module implements the Davidson iterative eigenvalue algorithm specialized
//! for the BSE-TDA Hamiltonian. The solver finds the lowest n_roots eigenvalues
//! using a matrix-free approach compatible with the S6-1 BSE-TDA kernel.
//!
//! # Algorithm
//!
//! The Davidson algorithm is an iterative subspace method:
//! 1. Initialize subspace with unit vectors at lowest QP gaps
//! 2. Apply Hamiltonian: AV = H @ V (matrix-free via S6-1)
//! 3. Form projected Hamiltonian: H_eff = V^T @ AV
//! 4. Solve small eigenvalue problem for Ritz values θ and vectors y
//! 5. Compute residuals: r = H x - θ x
//! 6. If converged: return eigenvalues and eigenvectors
//! 7. Expand subspace with preconditioned residuals
//! 8. Restart if subspace exceeds max_space
//!
//! # References
//!
//! - Davidson, E.R. (1975). "The iterative calculation of a few of the lowest
//!   eigenvalues and corresponding eigenvectors of large real-symmetric matrices"
//! - PySCF: pyscf/lib/linalg_helper.py (davidson1 function)
//! - Theory: docs/derivations/s6-2/theory.md

#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::Result;
use ndarray::{s, Array1, Array2, ArrayView1};
use ndarray_linalg::{Eigh, Norm, QR, UPLO};

/// Davidson solver configuration
///
/// # Example
///
/// ```ignore
/// use quasix_core::bse::DavidsonConfig;
///
/// let config = DavidsonConfig {
///     n_roots: 5,
///     max_iter: 50,
///     tol_residual: 1e-6,
///     ..Default::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct DavidsonConfig {
    /// Number of eigenvalues to find
    pub n_roots: usize,

    /// Maximum Davidson iterations
    pub max_iter: usize,

    /// Maximum subspace size before restart
    pub max_space: usize,

    /// Convergence tolerance for residual norm
    pub tol_residual: f64,

    /// Linear dependence threshold for orthogonalization
    pub lindep_threshold: f64,

    /// Level shift for preconditioner (prevents division by zero)
    pub level_shift: f64,
}

impl Default for DavidsonConfig {
    fn default() -> Self {
        Self {
            n_roots: 3,
            max_iter: 100,
            max_space: 60, // 20 * n_roots
            tol_residual: 1e-5,
            lindep_threshold: 1e-12,
            level_shift: 0.001, // ~0.027 eV
        }
    }
}

impl DavidsonConfig {
    /// Create configuration for a specific number of roots
    #[must_use]
    pub fn with_n_roots(n_roots: usize) -> Self {
        Self {
            n_roots,
            max_space: 20 * n_roots,
            ..Self::default()
        }
    }

    /// Create configuration for tight convergence
    #[must_use]
    pub fn tight() -> Self {
        Self {
            tol_residual: 1e-8,
            lindep_threshold: 1e-14,
            ..Self::default()
        }
    }
}

/// Davidson solver result
///
/// Contains the converged eigenvalues, eigenvectors, and convergence information.
#[derive(Debug, Clone)]
pub struct DavidsonResult {
    /// Converged eigenvalues (excitation energies)
    pub eigenvalues: Array1<f64>,

    /// Converged eigenvectors (n_trans × n_roots)
    pub eigenvectors: Array2<f64>,

    /// Convergence status for each root
    pub converged: Vec<bool>,

    /// Number of iterations performed
    pub iterations: usize,

    /// Final residual norms for each root
    pub residual_norms: Vec<f64>,
}

impl DavidsonResult {
    /// Check if all requested roots converged
    #[must_use]
    pub fn all_converged(&self) -> bool {
        self.converged.iter().all(|&c| c)
    }

    /// Count number of converged roots
    #[must_use]
    pub fn n_converged(&self) -> usize {
        self.converged.iter().filter(|&&c| c).count()
    }

    /// Get the maximum residual norm
    #[must_use]
    pub fn max_residual(&self) -> f64 {
        self.residual_norms.iter().fold(0.0_f64, |a, &b| a.max(b))
    }
}

/// Initialize subspace with unit vectors at lowest QP gaps
///
/// Creates unit vectors at positions corresponding to the n_roots lowest
/// QP energy differences. This provides a good starting point for the
/// Davidson iteration.
///
/// # Arguments
///
/// * `qp_gaps` - QP energy differences (diagonal of BSE Hamiltonian)
/// * `n_roots` - Number of starting vectors
///
/// # Returns
///
/// Initial orthonormal subspace V, shape [n_trans, n_roots]
fn initialize_subspace(qp_gaps: ArrayView1<f64>, n_roots: usize) -> Result<Array2<f64>> {
    let n_trans = qp_gaps.len();
    let n_init = n_roots.min(n_trans);

    // Sort gaps and get indices of lowest
    let mut indices: Vec<usize> = (0..n_trans).collect();
    indices.sort_by(|&a, &b| {
        qp_gaps[a]
            .partial_cmp(&qp_gaps[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Create unit vectors at lowest gap positions
    let mut v = Array2::zeros((n_trans, n_init));
    for (i, &idx) in indices.iter().take(n_init).enumerate() {
        v[[idx, i]] = 1.0;
    }

    Ok(v)
}

/// Apply diagonal preconditioner to residual
///
/// Computes: t_i = r_i / (D_i - θ + level_shift)
///
/// where D is the diagonal of the BSE Hamiltonian (QP gaps).
///
/// # Arguments
///
/// * `residual` - Residual vector r = Hx - θx
/// * `qp_gaps` - QP energy differences (diagonal)
/// * `ritz_value` - Current Ritz value θ
/// * `level_shift` - Shift to prevent division by zero
///
/// # Returns
///
/// Preconditioned and normalized residual
fn apply_preconditioner(
    residual: ArrayView1<f64>,
    qp_gaps: ArrayView1<f64>,
    ritz_value: f64,
    level_shift: f64,
) -> Array1<f64> {
    let n = residual.len();
    let mut precond = Array1::zeros(n);

    for i in 0..n {
        let denom = qp_gaps[i] - ritz_value + level_shift;
        // Clamp to avoid extreme values
        let denom_safe = if denom.abs() < 1e-8 {
            if denom >= 0.0 {
                1e-8
            } else {
                -1e-8
            }
        } else {
            denom
        };
        precond[i] = residual[i] / denom_safe;
    }

    // Normalize
    let norm = precond.norm_l2();
    if norm > 1e-14 {
        precond /= norm;
    }

    precond
}

/// Orthogonalize new vectors against existing subspace
///
/// Uses modified Gram-Schmidt with re-orthogonalization for numerical stability.
///
/// # Arguments
///
/// * `v` - Current subspace V
/// * `new_vec` - New vector to orthogonalize
/// * `lindep_threshold` - Threshold for linear dependence
///
/// # Returns
///
/// Orthonormalized vector, or None if linearly dependent
fn orthogonalize_against_subspace(
    v: &Array2<f64>,
    new_vec: Array1<f64>,
    lindep_threshold: f64,
) -> Option<Array1<f64>> {
    let mut t = new_vec;

    // First orthogonalization pass
    for j in 0..v.ncols() {
        let v_j = v.column(j);
        let proj: f64 = t.iter().zip(v_j.iter()).map(|(&ti, &vj)| ti * vj).sum();
        t = &t - proj * &v_j;
    }

    // Second orthogonalization pass (for numerical stability)
    for j in 0..v.ncols() {
        let v_j = v.column(j);
        let proj: f64 = t.iter().zip(v_j.iter()).map(|(&ti, &vj)| ti * vj).sum();
        t = &t - proj * &v_j;
    }

    // Check norm (linear dependence test)
    let norm = t.norm_l2();
    if norm < lindep_threshold {
        return None;
    }

    // Normalize
    Some(&t / norm)
}

/// Expand subspace with new trial vectors
///
/// Orthogonalizes new vectors against the current subspace and appends them.
/// If subspace exceeds max_space, performs restart with Ritz vectors.
///
/// # Arguments
///
/// * `v` - Current subspace
/// * `new_vectors` - New trial vectors to add
/// * `ritz_vectors` - Current best approximations (for restart)
/// * `max_space` - Maximum subspace size
/// * `lindep_threshold` - Linear dependence threshold
///
/// # Returns
///
/// Expanded (or restarted) subspace
fn expand_subspace(
    v: Array2<f64>,
    new_vectors: Vec<Array1<f64>>,
    ritz_vectors: &Array2<f64>,
    n_roots: usize,
    max_space: usize,
    lindep_threshold: f64,
) -> Result<Array2<f64>> {
    let current_space = v.ncols();
    let n_new = new_vectors.len();

    // Check if restart needed
    if current_space + n_new > max_space {
        // Restart: keep only Ritz vectors corresponding to requested roots
        let n_keep = n_roots.min(ritz_vectors.ncols());
        let restart_vectors = ritz_vectors.slice(s![.., ..n_keep]).to_owned();

        // QR orthonormalize the restart vectors
        let (q, _) = restart_vectors.qr()?;
        return Ok(q);
    }

    // Orthogonalize and add new vectors
    let mut v_expanded = v;
    for new_vec in new_vectors {
        if let Some(orthonorm) =
            orthogonalize_against_subspace(&v_expanded, new_vec, lindep_threshold)
        {
            // Append new column
            let n_rows = v_expanded.nrows();
            let n_cols = v_expanded.ncols();
            let mut v_new = Array2::zeros((n_rows, n_cols + 1));
            v_new.slice_mut(s![.., ..n_cols]).assign(&v_expanded);
            v_new.column_mut(n_cols).assign(&orthonorm);
            v_expanded = v_new;
        }
    }

    Ok(v_expanded)
}

/// Solve small symmetric eigenvalue problem
///
/// Uses LAPACK DSYEVR for the projected Hamiltonian.
fn solve_small_eigenproblem(h_eff: &Array2<f64>) -> Result<(Array1<f64>, Array2<f64>)> {
    // Symmetrize to ensure numerical symmetry
    let n = h_eff.nrows();
    let mut h_sym = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            h_sym[[i, j]] = 0.5 * (h_eff[[i, j]] + h_eff[[j, i]]);
        }
    }

    // Solve eigenvalue problem
    let (eigenvalues, eigenvectors) = h_sym.eigh(UPLO::Lower)?;

    Ok((eigenvalues, eigenvectors))
}

/// Main Davidson solver for BSE-TDA eigenvalue problem
///
/// Finds the lowest n_roots eigenvalues and eigenvectors of the BSE-TDA
/// Hamiltonian using the Davidson iterative algorithm with a matrix-free
/// Hamiltonian application.
///
/// # Arguments
///
/// * `apply_h` - Function that applies H to a batch of vectors: Y = H @ X
/// * `qp_gaps` - QP energy differences (diagonal of H)
/// * `config` - Davidson configuration
///
/// # Returns
///
/// DavidsonResult containing eigenvalues, eigenvectors, and convergence info
///
/// # Algorithm
///
/// 1. Initialize subspace with unit vectors at lowest QP gaps
/// 2. Apply Hamiltonian to subspace
/// 3. Form and solve projected eigenvalue problem (Rayleigh-Ritz)
/// 4. Compute residuals for each root
/// 5. Expand subspace with preconditioned residuals
/// 6. Restart if subspace exceeds max_space
/// 7. Repeat until convergence or max_iter
pub fn davidson_bse<F>(
    apply_h: F,
    qp_gaps: ArrayView1<f64>,
    config: &DavidsonConfig,
) -> Result<DavidsonResult>
where
    F: Fn(&Array2<f64>) -> Result<Array2<f64>>,
{
    let n_trans = qp_gaps.len();
    let n_roots = config.n_roots.min(n_trans);
    let max_space = config.max_space.min(n_trans);

    // Handle trivial case: n_trans <= n_roots
    if n_trans <= n_roots {
        return solve_trivial_case(&apply_h, qp_gaps, n_trans);
    }

    // Initialize subspace with unit vectors at lowest gaps
    let mut v = initialize_subspace(qp_gaps, n_roots)?;

    // Storage for results
    let mut eigenvalues = Array1::zeros(n_roots);
    let mut eigenvectors = Array2::zeros((n_trans, n_roots));
    let mut converged = vec![false; n_roots];
    let mut residual_norms = vec![f64::MAX; n_roots];

    for iteration in 0..config.max_iter {
        // 1. Apply Hamiltonian to subspace: AV = H @ V
        let av = apply_h(&v)?;

        // 2. Form projected Hamiltonian: H_eff = V^T @ AV
        let h_eff = v.t().dot(&av);

        // 3. Solve small eigenvalue problem
        let (theta, y) = solve_small_eigenproblem(&h_eff)?;

        // 4. Compute Ritz vectors: X = V @ Y
        let n_ritz = n_roots.min(y.ncols());
        let y_slice = y.slice(s![.., ..n_ritz]).to_owned();
        let x = v.dot(&y_slice);

        // Also compute AX = AV @ Y for residual
        let ax = av.dot(&y_slice);

        // 5. Check convergence and compute preconditioned residuals
        let mut new_vectors = Vec::new();

        for j in 0..n_roots {
            if j >= n_ritz {
                // Not enough Ritz vectors
                continue;
            }

            if converged[j] {
                continue;
            }

            // Residual: r_j = Ax_j - θ_j * x_j
            let x_j = x.column(j);
            let ax_j = ax.column(j);
            let r_j: Array1<f64> = ax_j.to_owned() - theta[j] * &x_j;

            let res_norm = r_j.norm_l2();
            residual_norms[j] = res_norm;

            if res_norm < config.tol_residual {
                converged[j] = true;
                eigenvalues[j] = theta[j];
                eigenvectors.column_mut(j).assign(&x_j);
                continue;
            }

            // Preconditioned residual: t = (D - θ)^{-1} r
            let t_j = apply_preconditioner(r_j.view(), qp_gaps, theta[j], config.level_shift);

            if t_j.norm_l2() > config.lindep_threshold {
                new_vectors.push(t_j);
            }
        }

        // 6. Check if all converged
        if converged.iter().all(|&c| c) {
            return Ok(DavidsonResult {
                eigenvalues,
                eigenvectors,
                converged,
                iterations: iteration + 1,
                residual_norms,
            });
        }

        // 7. Expand or restart subspace
        if new_vectors.is_empty() {
            // No new vectors to add - might be stuck
            // Try continuing with current subspace
            continue;
        }

        v = expand_subspace(
            v,
            new_vectors,
            &x,
            n_roots,
            max_space,
            config.lindep_threshold,
        )?;
    }

    // Did not fully converge - return best results
    // Make sure to store final eigenvalues even if not converged
    for j in 0..n_roots {
        if !converged[j] && eigenvectors.column(j).norm_l2() < 0.1 {
            // Copy latest Ritz approximation if not stored yet
            let av = apply_h(&v)?;
            let h_eff = v.t().dot(&av);
            let (theta, y) = solve_small_eigenproblem(&h_eff)?;
            let n_ritz = n_roots.min(y.ncols());
            if j < n_ritz {
                let y_slice = y.slice(s![.., ..n_ritz]).to_owned();
                let x = v.dot(&y_slice);
                eigenvalues[j] = theta[j];
                eigenvectors.column_mut(j).assign(&x.column(j));
            }
        }
    }

    Ok(DavidsonResult {
        eigenvalues,
        eigenvectors,
        converged,
        iterations: config.max_iter,
        residual_norms,
    })
}

/// Solve trivial case when n_trans <= n_roots
///
/// For very small systems, directly build and diagonalize the explicit matrix.
fn solve_trivial_case<F>(
    apply_h: &F,
    _qp_gaps: ArrayView1<f64>,
    n_trans: usize,
) -> Result<DavidsonResult>
where
    F: Fn(&Array2<f64>) -> Result<Array2<f64>>,
{
    // Build full matrix by applying H to identity
    let identity = Array2::eye(n_trans);
    let h_full = apply_h(&identity)?;

    // Solve eigenvalue problem
    let (eigenvalues, eigenvectors) = solve_small_eigenproblem(&h_full)?;

    // Compute residual norms (should be ~0 for exact solution)
    let n_roots = n_trans;
    let mut residual_norms = Vec::with_capacity(n_roots);
    let hx = apply_h(&eigenvectors)?;
    for j in 0..n_roots {
        let r = &hx.column(j) - eigenvalues[j] * &eigenvectors.column(j);
        residual_norms.push(r.norm_l2());
    }

    Ok(DavidsonResult {
        eigenvalues,
        eigenvectors,
        converged: vec![true; n_roots],
        iterations: 1,
        residual_norms,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::arr1;

    /// Create a simple diagonal matrix test case
    fn create_diagonal_apply_h(
        diagonal: Array1<f64>,
    ) -> impl Fn(&Array2<f64>) -> Result<Array2<f64>> {
        move |v: &Array2<f64>| {
            let mut result = v.clone();
            for i in 0..result.nrows() {
                for j in 0..result.ncols() {
                    result[[i, j]] *= diagonal[i];
                }
            }
            Ok(result)
        }
    }

    #[test]
    fn test_davidson_config_default() {
        let config = DavidsonConfig::default();
        assert_eq!(config.n_roots, 3);
        assert_eq!(config.max_iter, 100);
        assert_eq!(config.max_space, 60);
        assert_relative_eq!(config.tol_residual, 1e-5);
    }

    #[test]
    fn test_davidson_config_with_n_roots() {
        let config = DavidsonConfig::with_n_roots(5);
        assert_eq!(config.n_roots, 5);
        assert_eq!(config.max_space, 100); // 20 * 5
    }

    #[test]
    fn test_davidson_config_tight() {
        let config = DavidsonConfig::tight();
        assert_relative_eq!(config.tol_residual, 1e-8);
    }

    #[test]
    fn test_davidson_result_all_converged() {
        let result = DavidsonResult {
            eigenvalues: Array1::zeros(3),
            eigenvectors: Array2::zeros((5, 3)),
            converged: vec![true, true, true],
            iterations: 10,
            residual_norms: vec![1e-8, 1e-8, 1e-8],
        };
        assert!(result.all_converged());
        assert_eq!(result.n_converged(), 3);
    }

    #[test]
    fn test_davidson_result_partial_converged() {
        let result = DavidsonResult {
            eigenvalues: Array1::zeros(3),
            eigenvectors: Array2::zeros((5, 3)),
            converged: vec![true, false, true],
            iterations: 100,
            residual_norms: vec![1e-8, 1e-3, 1e-8],
        };
        assert!(!result.all_converged());
        assert_eq!(result.n_converged(), 2);
        assert_relative_eq!(result.max_residual(), 1e-3);
    }

    #[test]
    fn test_initialize_subspace() {
        let qp_gaps = arr1(&[0.5, 0.2, 0.8, 0.1, 0.6]); // Index 3 has lowest
        let v = initialize_subspace(qp_gaps.view(), 3).unwrap();

        assert_eq!(v.dim(), (5, 3));

        // First vector should be unit vector at position 3 (lowest gap 0.1)
        assert_relative_eq!(v[[3, 0]], 1.0);
        // Second at position 1 (gap 0.2)
        assert_relative_eq!(v[[1, 1]], 1.0);
        // Third at position 0 (gap 0.5)
        assert_relative_eq!(v[[0, 2]], 1.0);
    }

    #[test]
    fn test_apply_preconditioner() {
        let residual = arr1(&[1.0, 2.0, 3.0]);
        let qp_gaps = arr1(&[0.5, 0.6, 0.7]);
        let ritz_value = 0.4;
        let level_shift = 0.001;

        let precond =
            apply_preconditioner(residual.view(), qp_gaps.view(), ritz_value, level_shift);

        // Should be normalized
        assert_relative_eq!(precond.norm_l2(), 1.0, epsilon = 1e-10);

        // Preconditioned values: r[i] / (D[i] - theta + shift)
        // D - theta + shift = [0.101, 0.201, 0.301]
        // raw = [1/0.101, 2/0.201, 3/0.301] = [9.90, 9.95, 9.97] approx
        // After normalization, all should be similar magnitude
    }

    #[test]
    fn test_orthogonalize_against_subspace() {
        // Create subspace with first two unit vectors
        let mut v = Array2::zeros((4, 2));
        v[[0, 0]] = 1.0;
        v[[1, 1]] = 1.0;

        // New vector with components in all directions
        let new_vec = arr1(&[0.5, 0.5, 0.7, 0.2]);

        let orthogonalized = orthogonalize_against_subspace(&v, new_vec, 1e-12)
            .expect("Should not be linearly dependent");

        // Should be normalized
        assert_relative_eq!(orthogonalized.norm_l2(), 1.0, epsilon = 1e-10);

        // Should be orthogonal to both existing vectors
        let dot0: f64 = orthogonalized
            .iter()
            .zip(v.column(0).iter())
            .map(|(a, b)| a * b)
            .sum();
        let dot1: f64 = orthogonalized
            .iter()
            .zip(v.column(1).iter())
            .map(|(a, b)| a * b)
            .sum();
        assert_relative_eq!(dot0, 0.0, epsilon = 1e-10);
        assert_relative_eq!(dot1, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_orthogonalize_linearly_dependent() {
        // Create subspace spanning e1, e2
        let mut v = Array2::zeros((4, 2));
        v[[0, 0]] = 1.0;
        v[[1, 1]] = 1.0;

        // New vector in span of existing subspace
        let new_vec = arr1(&[0.6, 0.8, 0.0, 0.0]);

        let result = orthogonalize_against_subspace(&v, new_vec, 1e-12);
        assert!(result.is_none());
    }

    #[test]
    fn test_davidson_diagonal_matrix() {
        // Diagonal matrix: eigenvalues = diagonal elements
        let diagonal = arr1(&[0.3, 0.1, 0.5, 0.2, 0.4]); // Sorted: 0.1, 0.2, 0.3, 0.4, 0.5
        let apply_h = create_diagonal_apply_h(diagonal.clone());

        let config = DavidsonConfig {
            n_roots: 3,
            max_iter: 50,
            tol_residual: 1e-8,
            ..Default::default()
        };

        let result = davidson_bse(apply_h, diagonal.view(), &config).unwrap();

        // Should converge in very few iterations
        assert!(result.iterations <= 5, "Diagonal should converge quickly");
        assert!(result.all_converged());

        // Eigenvalues should be the 3 smallest diagonal elements
        assert_relative_eq!(result.eigenvalues[0], 0.1, epsilon = 1e-8);
        assert_relative_eq!(result.eigenvalues[1], 0.2, epsilon = 1e-8);
        assert_relative_eq!(result.eigenvalues[2], 0.3, epsilon = 1e-8);
    }

    #[test]
    fn test_davidson_trivial_case() {
        // n_trans <= n_roots should trigger trivial case
        let diagonal = arr1(&[0.3, 0.1]);
        let apply_h = create_diagonal_apply_h(diagonal.clone());

        let config = DavidsonConfig {
            n_roots: 5, // More roots than transitions
            tol_residual: 1e-8,
            ..Default::default()
        };

        let result = davidson_bse(apply_h, diagonal.view(), &config).unwrap();

        assert_eq!(result.iterations, 1);
        assert!(result.all_converged());
        assert_eq!(result.eigenvalues.len(), 2);
        assert_relative_eq!(result.eigenvalues[0], 0.1, epsilon = 1e-10);
        assert_relative_eq!(result.eigenvalues[1], 0.3, epsilon = 1e-10);
    }

    #[test]
    fn test_davidson_symmetric_matrix() {
        // Create a small symmetric matrix and verify Davidson matches eigh
        let n = 6;

        // Build a symmetric positive definite matrix
        let mut h = Array2::zeros((n, n));
        for i in 0..n {
            h[[i, i]] = (i + 1) as f64 * 0.2;
            for j in (i + 1)..n {
                let val = 0.02 * ((i + j) as f64);
                h[[i, j]] = val;
                h[[j, i]] = val;
            }
        }

        let h_for_closure = h.clone();
        let apply_h = move |v: &Array2<f64>| -> Result<Array2<f64>> { Ok(h_for_closure.dot(v)) };

        // Extract diagonal for QP gaps
        let qp_gaps: Array1<f64> = (0..n).map(|i| h[[i, i]]).collect();

        let config = DavidsonConfig {
            n_roots: 3,
            max_iter: 50,
            tol_residual: 1e-8,
            ..Default::default()
        };

        let result = davidson_bse(apply_h, qp_gaps.view(), &config).unwrap();

        // Compare with direct diagonalization
        let (ref_eigs, _) = h.eigh(UPLO::Lower).unwrap();

        assert!(result.all_converged());
        for i in 0..3 {
            assert_relative_eq!(result.eigenvalues[i], ref_eigs[i], epsilon = 1e-8);
        }

        // Check eigenvector orthonormality
        let xtx = result.eigenvectors.t().dot(&result.eigenvectors);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_relative_eq!(xtx[[i, j]], 1.0, epsilon = 1e-10);
                } else {
                    assert_relative_eq!(xtx[[i, j]], 0.0, epsilon = 1e-10);
                }
            }
        }
    }

    #[test]
    fn test_davidson_convergence_iterations() {
        // Test that solver converges within max_iter
        let n = 20;

        // Build a symmetric matrix with some coupling
        let mut h = Array2::zeros((n, n));
        for i in 0..n {
            h[[i, i]] = 0.1 + 0.05 * (i as f64);
            for j in (i + 1)..n {
                let val = 0.01 / ((i as f64 - j as f64).abs() + 1.0);
                h[[i, j]] = val;
                h[[j, i]] = val;
            }
        }

        let h_for_closure = h.clone();
        let apply_h = move |v: &Array2<f64>| -> Result<Array2<f64>> { Ok(h_for_closure.dot(v)) };

        let qp_gaps: Array1<f64> = (0..n).map(|i| h[[i, i]]).collect();

        let config = DavidsonConfig {
            n_roots: 5,
            max_iter: 100,
            tol_residual: 1e-6,
            ..Default::default()
        };

        let result = davidson_bse(apply_h, qp_gaps.view(), &config).unwrap();

        assert!(result.all_converged());
        assert!(
            result.iterations < 100,
            "Should converge in < 100 iterations, got {}",
            result.iterations
        );

        // Verify eigenvalues are correct
        let (ref_eigs, _) = h.eigh(UPLO::Lower).unwrap();
        for i in 0..5 {
            assert_relative_eq!(result.eigenvalues[i], ref_eigs[i], epsilon = 1e-6);
        }
    }
}
