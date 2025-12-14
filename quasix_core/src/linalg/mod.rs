/// Linear algebra utilities for QuasiX
///
/// This module provides numerical linear algebra operations
/// with a focus on stability and efficiency for quantum chemistry.
pub mod conditioning;
pub mod metric_sqrt;

pub use conditioning::{
    analyze_condition, apply_adaptive_regularization, apply_tikhonov_regularization,
    estimate_condition_power_iteration, is_well_conditioned, ConditionAnalysis, ConditionConfig,
};
pub use metric_sqrt::{LinalgError, MetricSqrt};

use crate::common::Result;
use ndarray::{Array1, Array2};

/// Solve linear system Ax = b using Gaussian elimination
///
/// This is a simple implementation for small systems.
/// For production, would use LAPACK.
pub fn solve_linear_system(a: &Array2<f64>, b: &Array1<f64>) -> Result<Array1<f64>> {
    let n = a.nrows();
    if a.ncols() != n || b.len() != n {
        return Err(crate::common::QuasixError::InvalidInput(
            "Matrix dimensions incompatible".to_string(),
        ));
    }

    // Create augmented matrix
    let mut aug = a.clone();
    let mut rhs = b.clone();

    // Forward elimination
    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        for k in (i + 1)..n {
            if aug[[k, i]].abs() > aug[[max_row, i]].abs() {
                max_row = k;
            }
        }

        // Swap rows
        for j in 0..n {
            let temp = aug[[i, j]];
            aug[[i, j]] = aug[[max_row, j]];
            aug[[max_row, j]] = temp;
        }
        let temp = rhs[i];
        rhs[i] = rhs[max_row];
        rhs[max_row] = temp;

        // Make all rows below this one 0 in current column
        for k in (i + 1)..n {
            let factor = aug[[k, i]] / aug[[i, i]];
            for j in i..n {
                aug[[k, j]] -= factor * aug[[i, j]];
            }
            rhs[k] -= factor * rhs[i];
        }
    }

    // Back substitution
    let mut x = Array1::<f64>::zeros(n);
    for i in (0..n).rev() {
        x[i] = rhs[i];
        for j in (i + 1)..n {
            x[i] -= aug[[i, j]] * x[j];
        }
        x[i] /= aug[[i, i]];
    }

    Ok(x)
}

// Re-export moved to avoid conflicts
