//! Integration tests for Cholesky factorization of auxiliary metric v_{PQ}

use ndarray::{arr1, arr2, Array2};
use quasix_core::df::{compute_cholesky_v, incremental_cholesky_update, CholeskyMethod};

/// Helper function to assert two arrays are approximately equal
fn assert_array_eq(a: &Array2<f64>, b: &Array2<f64>, tol: f64) {
    assert_eq!(a.dim(), b.dim(), "Arrays have different dimensions");
    let max_diff = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_diff < tol,
        "Arrays differ by {:.2e} (tolerance: {:.2e})",
        max_diff,
        tol
    );
}

/// Create a symmetric positive definite matrix for testing
fn create_spd_matrix(n: usize, condition: f64) -> Array2<f64> {
    // Create a matrix with known eigenvalues
    // Since Q is identity, the result is just the diagonal matrix D
    // This makes the condition number exactly as specified
    let mut matrix = Array2::zeros((n, n));
    for i in 0..n {
        // Eigenvalues ranging from 1.0 to condition
        matrix[[i, i]] = 1.0 + (condition - 1.0) * (i as f64) / ((n - 1) as f64);
    }

    matrix
}

#[test]
fn test_cholesky_simple_3x3() {
    // Simple 3x3 positive definite matrix
    let metric = arr2(&[[4.0, 2.0, 1.0], [2.0, 5.0, 2.0], [1.0, 2.0, 6.0]]);

    let result = compute_cholesky_v(&metric, None).unwrap();

    // Check basic properties
    assert_eq!(result.naux, 3);
    assert_eq!(result.method, CholeskyMethod::Standard);
    assert!(result.rank.is_none());
    assert!(result.permutation.is_none());

    // Validate structure
    result.validate().unwrap();

    // Check reconstruction
    let reconstructed = result.reconstruct_metric();
    assert_array_eq(&reconstructed, &metric, 1e-10);

    // Check condition number is reasonable
    assert!(result.condition_number < 10.0);
}

#[test]
fn test_cholesky_well_conditioned() {
    // Well-conditioned 10x10 matrix
    let metric = create_spd_matrix(10, 10.0);

    let result = compute_cholesky_v(&metric, None).unwrap();

    // Should use standard Cholesky
    assert_eq!(result.method, CholeskyMethod::Standard);
    assert_eq!(result.naux, 10);

    // Check reconstruction error
    let recon_error = result.reconstruction_error(&metric);
    assert!(
        recon_error < 1e-10,
        "Reconstruction error: {:.2e}",
        recon_error
    );

    // Validate structure
    result.validate().unwrap();

    // Check that L is lower triangular
    let l = &result.l_matrix;
    for i in 0..10 {
        for j in i + 1..10 {
            assert!(l[[i, j]].abs() < 1e-14);
        }
    }
}

#[test]
fn test_cholesky_ill_conditioned() {
    // Ill-conditioned 8x8 matrix
    let metric = create_spd_matrix(8, 1e12);

    let result = compute_cholesky_v(&metric, None).unwrap();

    // Should use pivoted Cholesky for extremely ill-conditioned (but not rank-deficient) matrices
    // Since condition number is 1e12 > 1e10 and min eigenvalue is 1.0 (not near zero)
    assert!(
        matches!(result.method, CholeskyMethod::Pivoted { .. }),
        "Expected Pivoted Cholesky for ill-conditioned matrix, got {:?}",
        result.method
    );
    assert_eq!(result.naux, 8);
    assert!(result.rank.is_some());

    // Check reconstruction error is acceptable for ill-conditioned matrix
    // With eigenvalue decomposition, we should get reasonable reconstruction
    let recon_error = result.reconstruction_error(&metric);
    assert!(
        recon_error < 1e-3,
        "Reconstruction error: {:.2e}",
        recon_error
    );

    // Validate structure
    result.validate().unwrap();
}

#[test]
fn test_cholesky_rank_deficient() {
    // Create a truly rank-deficient matrix
    let mut metric = Array2::eye(5);
    metric[[4, 4]] = 1e-15; // Nearly singular (below machine epsilon threshold)

    let result = compute_cholesky_v(&metric, Some(1e-10)).unwrap();

    // Should use eigenvalue decomposition for rank-deficient matrices
    // since min eigenvalue is 1e-15 < 1e-14
    assert!(
        matches!(
            result.method,
            CholeskyMethod::EigenvalueDecomposition { .. }
        ),
        "Expected EigenvalueDecomposition for rank-deficient matrix, got {:?}",
        result.method
    );

    if let Some(rank) = result.rank {
        assert!(rank <= 4, "Expected rank <= 4, got {}", rank);
    }
}

#[test]
fn test_cholesky_not_symmetric() {
    let metric = arr2(&[
        [1.0, 2.0],
        [2.1, 4.0], // Not symmetric
    ]);

    let result = compute_cholesky_v(&metric, None);
    assert!(result.is_err());

    if let Err(e) = result {
        let msg = e.to_string();
        assert!(msg.contains("symmetric"));
    }
}

#[test]
fn test_cholesky_not_positive_definite() {
    let metric = arr2(&[
        [1.0, 2.0],
        [2.0, 1.0], // Eigenvalues: 3, -1 (not positive definite)
    ]);

    let result = compute_cholesky_v(&metric, None);
    assert!(result.is_err());

    if let Err(e) = result {
        let msg = e.to_string();
        assert!(
            msg.contains("positive definite") || msg.contains("eigenvalue"),
            "Expected error about positive definiteness, got: {}",
            msg
        );
    }
}

#[test]
fn test_incremental_update() {
    // Start with 3x3 matrix
    let metric = arr2(&[[4.0, 2.0, 1.0], [2.0, 5.0, 2.0], [1.0, 2.0, 6.0]]);

    let cholesky = compute_cholesky_v(&metric, None).unwrap();

    // Add new row/column
    let v_new = arr1(&[0.5, 1.5, 2.5]);
    let v_diag = 8.0;

    let updated = incremental_cholesky_update(&cholesky, &v_new, v_diag).unwrap();

    assert_eq!(updated.naux, 4);
    assert_eq!(updated.method, CholeskyMethod::Standard);

    // Construct the full 4x4 matrix for comparison
    let mut full_metric = Array2::zeros((4, 4));
    full_metric.slice_mut(ndarray::s![..3, ..3]).assign(&metric);
    for i in 0..3 {
        full_metric[[3, i]] = v_new[i];
        full_metric[[i, 3]] = v_new[i];
    }
    full_metric[[3, 3]] = v_diag;

    // Compare with direct computation
    let direct = compute_cholesky_v(&full_metric, None).unwrap();

    // Check that the results match
    assert_array_eq(&updated.l_matrix, &direct.l_matrix, 1e-10);
}

#[test]
fn test_incremental_update_fails_not_pd() {
    // Start with 2x2 matrix
    let metric = arr2(&[[1.0, 0.5], [0.5, 1.0]]);

    let cholesky = compute_cholesky_v(&metric, None).unwrap();

    // Try to add a row/column that makes it not positive definite
    let v_new = arr1(&[10.0, 10.0]); // Too large correlation
    let v_diag = 1.0; // Too small diagonal

    let result = incremental_cholesky_update(&cholesky, &v_new, v_diag);
    assert!(result.is_err());
}

#[test]
fn test_large_matrix_performance() {
    // Test with a reasonably large matrix
    let n = 200;
    let metric = create_spd_matrix(n, 100.0);

    let start = std::time::Instant::now();
    let result = compute_cholesky_v(&metric, None).unwrap();
    let duration = start.elapsed();

    println!(
        "Cholesky factorization of {}x{} matrix took {:?}",
        n, n, duration
    );

    assert_eq!(result.naux, n);
    result.validate().unwrap();

    // Check reconstruction error
    let recon_error = result.reconstruction_error(&metric);
    assert!(
        recon_error < 1e-8,
        "Reconstruction error: {:.2e}",
        recon_error
    );

    // Performance expectation: should complete in reasonable time
    assert!(duration.as_secs() < 5, "Factorization took too long");
}

#[test]
fn test_diagonal_dominance() {
    // Create a diagonally dominant matrix (always positive definite)
    let n = 10;
    let mut metric = Array2::zeros((n, n));

    for i in 0..n {
        for j in 0..n {
            if i == j {
                metric[[i, j]] = 10.0;
            } else {
                metric[[i, j]] = 0.1 * ((i + j) as f64).sin();
            }
        }
    }

    let result = compute_cholesky_v(&metric, None).unwrap();

    assert_eq!(result.method, CholeskyMethod::Standard);
    assert!(result.condition_number < 100.0);

    // Check reconstruction
    let recon_error = result.reconstruction_error(&metric);
    assert!(recon_error < 1e-10);
}

#[test]
fn test_tolerance_parameter() {
    // Test different tolerance values with a moderately ill-conditioned matrix
    // Use condition number 1e8 which is ill-conditioned but not extreme
    let metric = create_spd_matrix(5, 1e8);

    // Default tolerance
    let result1 = compute_cholesky_v(&metric, None).unwrap();

    // Strict tolerance
    let result2 = compute_cholesky_v(&metric, Some(1e-12)).unwrap();

    // Loose tolerance
    let result3 = compute_cholesky_v(&metric, Some(1e-6)).unwrap();

    // All should succeed but may use different methods
    result1.validate().unwrap();
    result2.validate().unwrap();
    result3.validate().unwrap();

    // Check reconstruction errors
    // Standard Cholesky should work well with diagonal matrices
    // even when ill-conditioned, since there's no fill-in
    assert!(
        result1.reconstruction_error(&metric) < 1e-6,
        "Result1 reconstruction error: {:.2e} (method: {:?})",
        result1.reconstruction_error(&metric),
        result1.method
    );
    assert!(
        result2.reconstruction_error(&metric) < 1e-6,
        "Result2 reconstruction error: {:.2e} (method: {:?})",
        result2.reconstruction_error(&metric),
        result2.method
    );
    assert!(
        result3.reconstruction_error(&metric) < 1e-6,
        "Result3 reconstruction error: {:.2e} (method: {:?})",
        result3.reconstruction_error(&metric),
        result3.method
    );
}

#[test]
fn test_validation_checks() {
    // Create a valid Cholesky result
    let metric = arr2(&[[4.0, 2.0], [2.0, 5.0]]);

    let result = compute_cholesky_v(&metric, None).unwrap();

    // Validation should pass
    result.validate().unwrap();

    // Check reconstruction error method
    let error = result.reconstruction_error(&metric);
    assert!(error < 1e-10);
}
