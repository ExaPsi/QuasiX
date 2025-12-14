//! Linear dependency handling for density fitting metric matrices
//!
//! This module provides robust algorithms for detecting and removing linear
//! dependencies in auxiliary basis sets, which is crucial for numerical
//! stability in GW/BSE calculations.

use crate::common::{QuasixError, Result};
use ndarray::{s, Array1, Array2, Array3};
use std::fmt::Write;
use tracing::{debug, info, warn};

use ndarray_linalg::{Eigh, UPLO};

/// Configuration for linear dependency handling
#[derive(Debug, Clone)]
pub struct LinearDependencyConfig {
    /// Primary threshold for eigenvalue cutoff (default: 1e-10)
    pub eigenvalue_threshold: f64,

    /// Relative threshold based on maximum eigenvalue (default: 1e-12)
    pub relative_threshold: f64,

    /// Cholesky diagonal threshold for pivoting (default: 1e-8)
    pub cholesky_threshold: f64,

    /// Maximum condition number before regularization (default: 1e12)
    pub max_condition_number: f64,

    /// Method for handling dependencies
    pub method: DependencyMethod,

    /// Whether to print detailed analysis
    pub verbose: bool,
}

impl Default for LinearDependencyConfig {
    fn default() -> Self {
        Self {
            eigenvalue_threshold: 1e-10,
            relative_threshold: 1e-12,
            cholesky_threshold: 1e-8,
            max_condition_number: 1e12,
            method: DependencyMethod::Automatic,
            verbose: false,
        }
    }
}

impl LinearDependencyConfig {
    /// Conservative configuration for high-accuracy calculations
    pub fn conservative() -> Self {
        Self {
            eigenvalue_threshold: 1e-12,
            relative_threshold: 1e-14,
            cholesky_threshold: 1e-10,
            max_condition_number: 1e14,
            ..Default::default()
        }
    }

    /// Aggressive configuration for large systems with expected dependencies
    pub fn aggressive() -> Self {
        Self {
            eigenvalue_threshold: 1e-8,
            relative_threshold: 1e-10,
            cholesky_threshold: 1e-6,
            max_condition_number: 1e10,
            ..Default::default()
        }
    }

    /// Configuration for periodic systems where dependencies are common
    pub fn periodic_systems() -> Self {
        Self {
            eigenvalue_threshold: 1e-9,
            relative_threshold: 1e-11,
            method: DependencyMethod::EigenvalueDecomposition,
            ..Default::default()
        }
    }
}

/// Method for handling linear dependencies
#[derive(Debug, Clone, PartialEq)]
pub enum DependencyMethod {
    /// Eigenvalue decomposition with truncation
    EigenvalueDecomposition,

    /// Pivoted Cholesky with rank revealing
    PivotedCholesky,

    /// Canonical orthogonalization (X = S^{-1/2})
    CanonicalOrthogonalization,

    /// Symmetric orthogonalization (X = S^{-1/2} U)
    SymmetricOrthogonalization,

    /// Automatic selection based on condition number
    Automatic,
}

/// Analysis results from linear dependency detection
#[derive(Debug, Clone)]
pub struct DependencyAnalysis {
    /// All eigenvalues
    pub eigenvals: Array1<f64>,

    /// Eigenvectors (if computed)
    pub eigenvecs: Option<Array2<f64>>,

    /// Indices of dependent basis functions
    pub dependent_indices: Vec<usize>,

    /// Number of negative eigenvalues
    pub negative_count: usize,

    /// Number of near-zero eigenvalues
    pub near_zero_count: usize,

    /// Condition number of the matrix
    pub condition_number: f64,

    /// Effective threshold used
    pub effective_threshold: f64,

    /// Effective rank after removing dependencies
    pub rank: usize,

    /// Original dimension
    pub n_total: usize,
}

impl DependencyAnalysis {
    /// Generate a human-readable report of the analysis
    pub fn generate_report(&self) -> String {
        let mut report = String::new();

        report.push_str("Linear Dependency Analysis:\n");
        writeln!(&mut report, "  Original dimension: {}", self.n_total).unwrap();

        if self.negative_count > 0 {
            let min_eigenval = self.eigenvals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
            writeln!(
                &mut report,
                "  WARNING: Found {} negative eigenvalues (min: {:.2e})",
                self.negative_count, min_eigenval
            )
            .unwrap();
        }

        if self.near_zero_count > 0 {
            writeln!(
                &mut report,
                "  INFO: Found {} near-zero eigenvalues (|λ| < {:.2e})",
                self.near_zero_count, self.effective_threshold
            )
            .unwrap();
        }

        if self.condition_number > 1e10 {
            writeln!(
                &mut report,
                "  WARNING: Matrix is ill-conditioned (κ = {:.2e})",
                self.condition_number
            )
            .unwrap();
        } else {
            writeln!(
                &mut report,
                "  Condition number: {:.2e}",
                self.condition_number
            )
            .unwrap();
        }

        let retention_rate = 100.0 * self.rank as f64 / self.n_total as f64;
        writeln!(
            &mut report,
            "  Effective rank: {}/{} ({:.1}% retained)",
            self.rank, self.n_total, retention_rate
        )
        .unwrap();

        if retention_rate < 50.0 {
            report.push_str("  WARNING: More than half of auxiliary functions removed\n");
        }

        report
    }

    /// Check if the matrix has significant linear dependencies
    pub fn has_dependencies(&self) -> bool {
        !self.dependent_indices.is_empty()
    }

    /// Check if the matrix should be regularized
    pub fn needs_regularization(&self) -> bool {
        self.condition_number > 1e10 || self.negative_count > 0
    }
}

/// Transformation matrix for projecting to/from reduced space
#[derive(Debug, Clone)]
pub struct TransformationMatrix {
    /// Forward transformation: full -> reduced (n_full × n_reduced)
    pub forward: Array2<f64>,

    /// Backward transformation: reduced -> full (n_reduced × n_full)
    pub backward: Array2<f64>,

    /// Rank of the reduced space
    pub rank: usize,
}

/// Detect linear dependencies in a metric matrix
pub fn detect_linear_dependencies(
    metric: &Array2<f64>,
    config: &LinearDependencyConfig,
) -> Result<DependencyAnalysis> {
    let n = metric.nrows();

    if n != metric.ncols() {
        return Err(QuasixError::InvalidInput(format!(
            "Metric matrix must be square, got {}×{}",
            n,
            metric.ncols()
        )));
    }

    info!("Detecting linear dependencies in {}×{} metric matrix", n, n);

    // Compute eigenvalue decomposition
    let (eigenvals, eigenvecs) = {
        let metric_owned = metric.to_owned();
        match metric_owned.eigh(UPLO::Lower) {
            Ok((vals, vecs)) => (vals, Some(vecs)),
            Err(e) => {
                return Err(QuasixError::NumericalError(format!(
                    "Eigenvalue decomposition failed: {}",
                    e
                )));
            }
        }
    };

    // Analyze spectrum
    let lambda_max = eigenvals.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let lambda_min_positive = eigenvals
        .iter()
        .filter(|&&x| x > 0.0)
        .fold(f64::INFINITY, |a, &b| a.min(b));

    // Determine effective threshold
    let machine_eps_threshold = f64::EPSILON * lambda_max * (n as f64).sqrt();
    let effective_threshold = config
        .eigenvalue_threshold
        .max(config.relative_threshold * lambda_max)
        .max(machine_eps_threshold);

    if config.verbose {
        debug!(
            "Eigenvalue range: [{:.2e}, {:.2e}]",
            eigenvals.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
            lambda_max
        );
        debug!("Effective threshold: {:.2e}", effective_threshold);
    }

    // Identify dependent functions
    let mut dependent_indices = Vec::new();
    let mut near_zero_count = 0;
    let mut negative_count = 0;

    for (i, &lambda) in eigenvals.iter().enumerate() {
        if lambda < -effective_threshold {
            negative_count += 1;
            dependent_indices.push(i);
        } else if lambda.abs() < effective_threshold {
            near_zero_count += 1;
            dependent_indices.push(i);
        }
    }

    // Compute condition number
    let condition_number = if lambda_min_positive > 0.0 {
        lambda_max / lambda_min_positive
    } else {
        f64::INFINITY
    };

    let rank = n - dependent_indices.len();

    let analysis = DependencyAnalysis {
        eigenvals,
        eigenvecs,
        dependent_indices,
        negative_count,
        near_zero_count,
        condition_number,
        effective_threshold,
        rank,
        n_total: n,
    };

    if config.verbose {
        info!("{}", analysis.generate_report());
    }

    Ok(analysis)
}

/// Remove linear dependencies using eigenvalue truncation
pub fn remove_dependencies_eigenvalue(
    metric: &Array2<f64>,
    analysis: &DependencyAnalysis,
) -> Result<(Array2<f64>, TransformationMatrix)> {
    let n = metric.nrows();
    let rank = analysis.rank;

    if analysis.eigenvecs.is_none() {
        return Err(QuasixError::InvalidInput(
            "Eigenvectors required for eigenvalue truncation".to_string(),
        ));
    }

    let eigenvecs = analysis.eigenvecs.as_ref().unwrap();

    info!(
        "Removing {} linear dependencies using eigenvalue truncation",
        n - rank
    );

    // Keep only eigenvectors with eigenvalues above threshold
    let mut kept_indices = Vec::new();
    let mut kept_eigenvals = Vec::new();

    for (i, &lambda) in analysis.eigenvals.iter().enumerate() {
        if lambda > analysis.effective_threshold {
            kept_indices.push(i);
            kept_eigenvals.push(lambda);
        }
    }

    // Build transformation matrix X = U_kept * Λ_kept^{-1/2}
    let mut x_matrix = Array2::zeros((n, rank));
    for (j, &i) in kept_indices.iter().enumerate() {
        let scale = 1.0 / kept_eigenvals[j].sqrt();
        x_matrix
            .column_mut(j)
            .assign(&(&eigenvecs.column(i) * scale));
    }

    // Compute regularized metric in reduced space (should be identity)
    let metric_reg = x_matrix.t().dot(metric).dot(&x_matrix);

    // Check that regularized metric is well-conditioned
    let diag_min = (0..rank)
        .map(|i| metric_reg[[i, i]])
        .fold(f64::INFINITY, f64::min);
    let diag_max = (0..rank)
        .map(|i| metric_reg[[i, i]])
        .fold(f64::NEG_INFINITY, f64::max);

    let reg_condition = diag_max / diag_min;
    if reg_condition > 1e8 {
        warn!(
            "Regularized metric still has high condition number: {:.2e}",
            reg_condition
        );
    }

    Ok((
        metric_reg,
        TransformationMatrix {
            forward: x_matrix.clone(),
            backward: x_matrix.t().to_owned(),
            rank,
        },
    ))
}

/// Remove linear dependencies using pivoted Cholesky decomposition
pub fn remove_dependencies_cholesky(
    metric: &Array2<f64>,
    config: &LinearDependencyConfig,
) -> Result<(Vec<usize>, usize)> {
    let n = metric.nrows();
    let mut kept_functions = Vec::new();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut diag = metric.diag().to_owned();
    let mut working_matrix = metric.to_owned();

    info!("Removing linear dependencies using pivoted Cholesky");

    for k in 0..n {
        // Find pivot (maximum diagonal)
        let (pivot_idx, &pivot_val) = diag
            .slice(s![k..])
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap();

        let pivot_idx = pivot_idx + k;

        // Check for linear dependency
        if pivot_val < config.cholesky_threshold {
            info!(
                "Truncating at rank {} (pivot value {:.2e} < threshold {:.2e})",
                k, pivot_val, config.cholesky_threshold
            );
            break;
        }

        // Swap if needed
        if pivot_idx != k {
            perm.swap(k, pivot_idx);

            // Swap in working matrix
            for i in 0..n {
                working_matrix.swap([i, k], [i, pivot_idx]);
            }
            for i in 0..n {
                working_matrix.swap([k, i], [pivot_idx, i]);
            }

            // Swap in diagonal
            diag.swap(k, pivot_idx);
        }

        kept_functions.push(perm[k]);

        // Update remaining diagonal elements for next iteration
        let pivot_sqrt = pivot_val.sqrt();
        for i in k + 1..n {
            let l_ik = working_matrix[[i, k]] / pivot_sqrt;
            diag[i] -= l_ik * l_ik;
        }
    }

    let rank = kept_functions.len();
    info!("Pivoted Cholesky: kept {}/{} auxiliary functions", rank, n);

    Ok((kept_functions, rank))
}

/// Transform DF tensors to reduced auxiliary basis
pub fn transform_df_tensors(
    j3c_ao: &Array3<f64>,
    transform: &TransformationMatrix,
) -> Result<Array3<f64>> {
    let (n_ao1, n_ao2, n_aux) = j3c_ao.dim();
    let rank = transform.rank;

    if n_aux != transform.forward.nrows() {
        return Err(QuasixError::InvalidInput(format!(
            "DF tensor auxiliary dimension {} doesn't match transformation {}",
            n_aux,
            transform.forward.nrows()
        )));
    }

    info!(
        "Transforming DF tensors from {} to {} auxiliary functions",
        n_aux, rank
    );

    // Transform to reduced auxiliary basis: (μν|P') = (μν|P) X_{PP'}
    let mut j3c_reduced = Array3::zeros((n_ao1, n_ao2, rank));

    // Process in chunks for better memory access
    for i in 0..n_ao1 {
        for j in 0..n_ao2 {
            let j3c_slice = j3c_ao.slice(s![i, j, ..]);
            let j3c_transformed = j3c_slice.dot(&transform.forward);
            j3c_reduced.slice_mut(s![i, j, ..]).assign(&j3c_transformed);
        }
    }

    Ok(j3c_reduced)
}

/// Transform MO DF tensors to reduced auxiliary basis
pub fn transform_mo_df_tensors(
    j3c_ia: &Array2<f64>,
    transform: &TransformationMatrix,
) -> Result<Array2<f64>> {
    let (_n_trans, n_aux) = j3c_ia.dim();

    if n_aux != transform.forward.nrows() {
        return Err(QuasixError::InvalidInput(format!(
            "MO DF tensor auxiliary dimension {} doesn't match transformation {}",
            n_aux,
            transform.forward.nrows()
        )));
    }

    info!(
        "Transforming MO DF tensors from {} to {} auxiliary functions",
        n_aux, transform.rank
    );

    // (ia|P') = (ia|P) X_{PP'}
    Ok(j3c_ia.dot(&transform.forward))
}

/// Apply regularization to a metric matrix
pub fn regularize_metric(
    metric: &Array2<f64>,
    config: &LinearDependencyConfig,
) -> Result<Array2<f64>> {
    let analysis = detect_linear_dependencies(metric, config)?;

    if !analysis.has_dependencies() && analysis.condition_number < config.max_condition_number {
        info!("No regularization needed");
        return Ok(metric.clone());
    }

    match config.method {
        DependencyMethod::EigenvalueDecomposition | DependencyMethod::Automatic => {
            let (metric_reg, _transform) = remove_dependencies_eigenvalue(metric, &analysis)?;
            Ok(metric_reg)
        }
        DependencyMethod::PivotedCholesky => {
            // For pivoted Cholesky, we return a reduced metric
            let (kept_indices, rank) = remove_dependencies_cholesky(metric, config)?;

            let mut metric_reduced = Array2::zeros((rank, rank));
            for (i, &idx_i) in kept_indices.iter().enumerate() {
                for (j, &idx_j) in kept_indices.iter().enumerate() {
                    metric_reduced[[i, j]] = metric[[idx_i, idx_j]];
                }
            }

            Ok(metric_reduced)
        }
        _ => Err(QuasixError::NotImplemented(format!(
            "Regularization method {:?} not yet implemented",
            config.method
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    fn create_rank_deficient_matrix(n: usize, rank: usize) -> Array2<f64> {
        assert!(rank <= n);

        // Create a matrix with specific rank
        let mut matrix = Array2::zeros((n, n));

        // Set first 'rank' diagonal elements to non-zero
        for i in 0..rank {
            matrix[[i, i]] = 1.0 + 0.1 * i as f64;

            // Add some off-diagonal elements
            for j in i + 1..n {
                let val = 0.01 * ((i + j) as f64).sin();
                matrix[[i, j]] = val;
                matrix[[j, i]] = val;
            }
        }

        // Last n-rank rows/columns are linear combinations
        for i in rank..n {
            for j in 0..n {
                let val = if j < rank {
                    matrix[[0, j]] * 0.5 + matrix[[1.min(rank - 1), j]] * 0.5
                } else {
                    1e-15 // Near-zero
                };
                matrix[[i, j]] = val;
                matrix[[j, i]] = val;
            }
        }

        matrix
    }

    #[test]
    fn test_detect_linear_dependencies() {
        let matrix = create_rank_deficient_matrix(10, 7);
        let config = LinearDependencyConfig::default();

        let analysis = detect_linear_dependencies(&matrix, &config).unwrap();

        assert!(analysis.has_dependencies());
        assert_eq!(analysis.rank, 7);
        assert_eq!(analysis.dependent_indices.len(), 3);
    }

    #[test]
    fn test_remove_dependencies_eigenvalue() {
        let matrix = create_rank_deficient_matrix(8, 5);
        let config = LinearDependencyConfig {
            method: DependencyMethod::EigenvalueDecomposition,
            ..Default::default()
        };

        let _analysis = detect_linear_dependencies(&matrix, &config).unwrap();

        // Note: This test requires eigenvalue decomposition functionality
        // which is available through ndarray-linalg
        let (metric_reg, transform) = remove_dependencies_eigenvalue(&matrix, &_analysis).unwrap();

        assert_eq!(metric_reg.nrows(), 5);
        assert_eq!(metric_reg.ncols(), 5);
        assert_eq!(transform.rank, 5);
    }

    #[test]
    fn test_diffuse_basis_dependency() {
        // Create a metric with actual linear dependency
        let n = 20;
        let mut metric = Array2::zeros((n, n));

        // Build a rank-deficient metric
        // First create a matrix with rank n-2
        for i in 0..n - 2 {
            for j in 0..n - 2 {
                let val = if i == j {
                    1.0 + 0.1 * i as f64
                } else {
                    0.01 * ((i + j) as f64).sin()
                };
                metric[[i, j]] = val;
            }
        }

        // Make last two rows/columns linear combinations
        for i in 0..n {
            metric[[n - 2, i]] = 0.5 * metric[[0, i]] + 0.3 * metric[[1, i]];
            metric[[i, n - 2]] = metric[[n - 2, i]];
            metric[[n - 1, i]] = 0.7 * metric[[0, i]] - 0.2 * metric[[1, i]];
            metric[[i, n - 1]] = metric[[n - 1, i]];
        }

        // Add small perturbation to make it numerically rank deficient
        metric[[n - 2, n - 2]] += 1e-11;
        metric[[n - 1, n - 1]] += 1e-11;

        let config = LinearDependencyConfig::aggressive();
        let analysis = detect_linear_dependencies(&metric, &config).unwrap();

        // Should detect dependencies or high condition number
        assert!(
            analysis.has_dependencies() || analysis.condition_number > 1e8,
            "Neither dependencies nor ill-conditioning detected: κ = {:.2e}, rank = {}/{}",
            analysis.condition_number,
            analysis.rank,
            n
        );
    }

    #[test]
    fn test_negative_eigenvalues() {
        let mut matrix = Array2::eye(5);
        matrix[[2, 2]] = -0.1; // Introduce negative eigenvalue

        let config = LinearDependencyConfig::default();
        let analysis = detect_linear_dependencies(&matrix, &config).unwrap();

        assert_eq!(analysis.negative_count, 1);
        assert!(analysis.needs_regularization());
    }

    #[test]
    fn test_condition_number_detection() {
        let mut matrix = Array2::eye(5);
        matrix[[0, 0]] = 1e12;
        matrix[[4, 4]] = 1e-12;

        let config = LinearDependencyConfig::default();
        let analysis = detect_linear_dependencies(&matrix, &config).unwrap();

        assert!(analysis.condition_number > 1e20);
        assert!(analysis.needs_regularization());
    }

    #[test]
    fn test_transform_preservation() {
        // Test that transformation preserves important properties
        let n = 6;
        let matrix = {
            let mut m = Array2::eye(n);
            // Add structure while maintaining positive definiteness
            for i in 0..n {
                for j in i + 1..n {
                    let val = 0.1 * f64::from((i as i32 - j as i32).abs()).recip();
                    m[[i, j]] = val;
                    m[[j, i]] = val;
                }
            }
            m
        };

        let config = LinearDependencyConfig::conservative();
        let analysis = detect_linear_dependencies(&matrix, &config).unwrap();

        // Verify the metric is well-conditioned
        assert!(!analysis.has_dependencies());
        assert!(analysis.condition_number < 100.0);
    }
}
