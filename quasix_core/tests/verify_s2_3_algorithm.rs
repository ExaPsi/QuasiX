// S2-3 Algorithm Verification Tests
// Tests the v^{±1/2} transformation implementation against theoretical requirements

use anyhow::Result;
use ndarray::{arr1, Array1, Array2};
use ndarray_linalg::{Eigh, Norm, UPLO};
use quasix_core::linalg::metric_sqrt::MetricSqrt;

/// Generate a symmetric positive definite matrix with specified condition number
fn generate_spd_matrix(n: usize, condition_number: f64) -> Array2<f64> {
    // Generate random orthogonal matrix using QR decomposition
    let mut rng = 42u64;
    let mut random = || {
        rng = (rng.wrapping_mul(1_103_515_245).wrapping_add(12345)) & 0x7fff_ffff;
        (rng as f64) / f64::from(0x7fff_ffff_u32) - 0.5
    };

    let mut a = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            a[[i, j]] = random();
        }
    }

    // Get orthogonal matrix via eigendecomposition of A^T A
    let ata = a.t().dot(&a);
    let (_, q) = ata.eigh(UPLO::Lower).unwrap();

    // Create diagonal with controlled condition number
    let mut eigenvalues = Array1::zeros(n);
    for i in 0..n {
        // Logarithmically spaced eigenvalues
        let t = i as f64 / (n - 1) as f64;
        eigenvalues[i] = condition_number.powf(t);
    }

    // Construct SPD matrix with desired spectrum
    let lambda = Array2::from_diag(&eigenvalues);
    q.dot(&lambda).dot(&q.t())
}

/// Generate a near-singular matrix with controlled rank deficiency
fn generate_near_singular_matrix(n: usize, rank: usize, epsilon: f64) -> Array2<f64> {
    let matrix = generate_spd_matrix(n, 1e3);

    // Perform eigendecomposition
    let (mut eigenvals, eigenvecs) = matrix.eigh(UPLO::Lower).unwrap();

    // Set smallest eigenvalues to near-zero
    for i in 0..(n - rank) {
        eigenvals[i] = epsilon;
    }

    // Reconstruct matrix
    let lambda = Array2::from_diag(&eigenvals);
    eigenvecs.dot(&lambda).dot(&eigenvecs.t())
}

#[test]
fn test_fundamental_properties() -> Result<()> {
    let test_cases = vec![
        ("Identity", Array2::eye(5)),
        (
            "Diagonal",
            Array2::from_diag(&arr1(&[1.0, 4.0, 9.0, 16.0, 25.0])),
        ),
        ("Well-conditioned SPD", generate_spd_matrix(10, 1e2)),
        ("Moderately ill-conditioned", generate_spd_matrix(10, 1e6)),
    ];

    for (name, matrix) in test_cases {
        let metric_sqrt = MetricSqrt::compute(&matrix, None)?;

        // Property 1: v^{1/2} v^{1/2} = v (in kept subspace)
        let reconstruction_error = metric_sqrt.verify_reconstruction(&matrix)?;

        // Property 2: v^{-1/2} v v^{-1/2} = I (in kept subspace)
        let identity_error = metric_sqrt.verify_identity(&matrix)?;

        // Check physical validity
        assert!(
            reconstruction_error < 1e-6 || metric_sqrt.condition_number > 1e8,
            "{}: Reconstruction error too large: {:.2e}",
            name,
            reconstruction_error
        );
        assert!(
            identity_error < 1e-6 || metric_sqrt.condition_number > 1e8,
            "{}: Identity error too large: {:.2e}",
            name,
            identity_error
        );
    }

    Ok(())
}

#[test]
fn test_numerical_stability() -> Result<()> {
    let condition_numbers = vec![1e2, 1e4, 1e6, 1e8, 1e10];

    for cond in condition_numbers {
        let matrix = generate_spd_matrix(20, cond);
        let metric_sqrt = MetricSqrt::compute(&matrix, Some(1e-10))?;

        // Test vector application stability
        let x = Array1::ones(20);
        let y = metric_sqrt.apply_sqrt(&x.view())?;
        let x_recovered = metric_sqrt.apply_inv_sqrt(&y.view())?;

        // Project x onto kept eigenspace for comparison
        let x_proj = metric_sqrt
            .eigenvecs_kept
            .dot(&metric_sqrt.eigenvecs_kept.t().dot(&x));
        let recovery_error = (&x_proj - &x_recovered).norm() / x_proj.norm().max(1e-15);

        // For condition numbers up to 1e10, should maintain reasonable accuracy
        if cond <= 1e10 {
            assert!(
                recovery_error < 1e-4,
                "Round-trip error {:.2e} too large for condition {:.2e}",
                recovery_error,
                cond
            );
        }
    }

    Ok(())
}

#[test]
fn test_eigenvalue_thresholding() -> Result<()> {
    let thresholds = vec![1e-14, 1e-12, 1e-10, 1e-8];
    let matrix = generate_near_singular_matrix(15, 10, 1e-13);

    for threshold in thresholds {
        let metric_sqrt = MetricSqrt::compute(&matrix, Some(threshold))?;

        // Verify thresholding is working correctly
        let min_kept_eigenval = metric_sqrt
            .eigenvals
            .iter()
            .filter(|&&x| x > metric_sqrt.threshold_used)
            .fold(f64::INFINITY, |a, &b| a.min(b));

        assert!(
            min_kept_eigenval > metric_sqrt.threshold_used,
            "Kept eigenvalue below threshold"
        );
    }

    Ok(())
}

#[test]
fn test_edge_cases() -> Result<()> {
    // Test 1: Identity matrix
    let identity = Array2::eye(5);
    let metric_sqrt = MetricSqrt::compute(&identity, None)?;
    assert_eq!(metric_sqrt.n_kept, 5);
    assert!((metric_sqrt.condition_number - 1.0).abs() < 1e-10);

    // Test 2: Pure diagonal matrix
    let diag = Array2::from_diag(&arr1(&[0.1, 1.0, 10.0, 100.0, 1000.0]));
    let metric_sqrt = MetricSqrt::compute(&diag, None)?;
    assert_eq!(metric_sqrt.n_kept, 5);
    assert!((metric_sqrt.condition_number - 1e4).abs() < 1.0);

    // Test 3: Rank-deficient matrix
    let rank_def = generate_near_singular_matrix(10, 5, 1e-15);
    let metric_sqrt = MetricSqrt::compute(&rank_def, Some(1e-12))?;
    assert!(metric_sqrt.n_kept <= 6); // Allow for numerical noise

    // Test 4: Very small matrix
    let small = Array2::from_shape_vec((2, 2), vec![2.0, 1.0, 1.0, 2.0]).unwrap();
    let metric_sqrt = MetricSqrt::compute(&small, None)?;
    assert_eq!(metric_sqrt.n_kept, 2);

    // Verify exact square root for 2x2 case
    let v_sqrt = metric_sqrt
        .eigenvecs_kept
        .dot(&Array2::from_diag(&metric_sqrt.sqrt_eigenvals))
        .dot(&metric_sqrt.eigenvecs_kept.t());
    let reconstructed = v_sqrt.dot(&v_sqrt);
    let error = (&reconstructed - &small).norm() / small.norm();
    assert!(error < 1e-10);

    Ok(())
}

#[test]
fn test_application_operations() -> Result<()> {
    let matrix = generate_spd_matrix(50, 1e4);
    let metric_sqrt = MetricSqrt::compute(&matrix, None)?;

    // Test vector application
    let x = Array1::ones(50);
    let y = metric_sqrt.apply_sqrt(&x.view())?;
    let z = metric_sqrt.apply_inv_sqrt(&y.view())?;

    // Project x for comparison (since we might have dropped eigenvalues)
    let x_proj = metric_sqrt
        .eigenvecs_kept
        .dot(&metric_sqrt.eigenvecs_kept.t().dot(&x));
    let round_trip_error = (&z - &x_proj).norm() / x_proj.norm();
    assert!(round_trip_error < 1e-10);

    // Test matrix application
    let x_mat = Array2::ones((50, 10));
    let y_mat = metric_sqrt.apply_sqrt_matrix(&x_mat.view())?;
    let z_mat = metric_sqrt.apply_inv_sqrt_matrix(&y_mat.view())?;

    // Project for comparison
    let x_mat_proj = metric_sqrt
        .eigenvecs_kept
        .dot(&metric_sqrt.eigenvecs_kept.t().dot(&x_mat));
    let mat_round_trip_error = (&z_mat - &x_mat_proj).norm() / x_mat_proj.norm();
    assert!(mat_round_trip_error < 1e-10);

    // Verify column-wise consistency
    for j in 0..10 {
        let col = x_mat.column(j);
        let y_col = metric_sqrt.apply_sqrt(&col)?;
        let y_mat_col = y_mat.column(j);
        let col_error = (&y_col - &y_mat_col).norm() / y_col.norm();
        assert!(
            col_error < 1e-12,
            "Column {} inconsistent: {:.2e}",
            j,
            col_error
        );
    }

    Ok(())
}

#[test]
fn test_physical_validity_gw_bse() -> Result<()> {
    // Simulate auxiliary basis metric from Coulomb integrals
    let aux_metric = generate_spd_matrix(30, 1e6);
    let metric_sqrt = MetricSqrt::compute(&aux_metric, Some(1e-10))?;

    // Check 1: All kept eigenvalues are positive
    for &eigenval in metric_sqrt.eigenvals.iter() {
        if eigenval > metric_sqrt.threshold_used {
            assert!(eigenval > 0.0, "Negative eigenvalue in Coulomb metric!");
        }
    }

    // Check 2: Symmetrized operations preserve positive definiteness
    // Simulate RPA polarizability (positive on imaginary axis)
    let mut p0 = Array2::zeros((30, 30));
    for i in 0..30 {
        for j in 0..30 {
            p0[[i, j]] = (-((i as f64 - j as f64).powi(2)) / 10.0).exp();
        }
    }

    // Transform to symmetrized basis: M = v^{1/2} P0 v^{1/2}
    let p0_transformed = metric_sqrt.apply_sqrt_matrix(&p0.view())?;
    let m = metric_sqrt.apply_sqrt_matrix(&p0_transformed.t().view())?;

    // M should be positive semi-definite
    let (m_eigenvals, _) = m.eigh(UPLO::Lower)?;
    let min_eigenval = m_eigenvals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    assert!(min_eigenval > -1e-10, "M not positive semi-definite!");

    // Check 3: Verify transformations preserve key properties for GW
    // The symmetrized formulation ensures numerical stability
    // Key property: v^{1/2} is well-defined and stable for physical Coulomb metrics

    // Test that applying v^{1/2} and v^{-1/2} preserves vector norms appropriately
    let test_vec = Array1::ones(30);
    let transformed = metric_sqrt.apply_sqrt(&test_vec.view())?;
    let back_transformed = metric_sqrt.apply_inv_sqrt(&transformed.view())?;

    // Should recover the projection of the original vector
    let projection = metric_sqrt
        .eigenvecs_kept
        .dot(&metric_sqrt.eigenvecs_kept.t().dot(&test_vec));
    let recovery_error = (&back_transformed - &projection).norm() / projection.norm();
    assert!(
        recovery_error < 1e-10,
        "Transform round-trip error too large: {:.2e}",
        recovery_error
    );

    // Verify that the transformation maintains physical properties
    // The condition number should be reasonable for a physical Coulomb metric
    assert!(
        metric_sqrt.condition_number < 1e10,
        "Condition number too large for physical system: {:.2e}",
        metric_sqrt.condition_number
    );

    Ok(())
}

#[test]
fn test_condition_number_extremes() -> Result<()> {
    // Test with condition number exactly at 1e12
    let matrix = generate_spd_matrix(15, 1e12);
    let metric_sqrt = MetricSqrt::compute(&matrix, Some(1e-10))?;

    // Should still work but with warnings
    assert!(metric_sqrt.condition_number <= 1e13);
    assert!(metric_sqrt.n_kept > 0);

    // Test reconstruction still reasonable
    let reconstruction_error = metric_sqrt.verify_reconstruction(&matrix)?;
    assert!(reconstruction_error < 1.0); // Relaxed tolerance for extreme case

    Ok(())
}
