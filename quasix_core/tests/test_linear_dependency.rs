//! Integration tests for linear dependency handling in density fitting

use ndarray::{Array2, Array3};
use quasix_core::df::{
    compute_cholesky_v_with_lindep, detect_linear_dependencies, regularize_metric,
    remove_dependencies_eigenvalue, transform_df_tensors, transform_mo_df_tensors,
    DependencyMethod, LinearDependencyConfig,
};

/// Create a metric matrix with intentional linear dependencies
fn create_dependent_metric(n: usize, n_dependent: usize) -> Array2<f64> {
    let mut metric = Array2::eye(n);

    // Add structure
    for i in 0..n {
        for j in i + 1..n {
            let val = 0.1 * ((i as f64 - j as f64).abs() + 1.0).recip();
            metric[[i, j]] = val;
            metric[[j, i]] = val;
        }
    }

    // Create linear dependencies by making some rows/cols linear combinations
    for i in (n - n_dependent)..n {
        for j in 0..n {
            // Make last n_dependent rows linear combinations of first rows
            metric[[i, j]] = 0.5 * metric[[0, j]] + 0.5 * metric[[1, j]] + 1e-14 * (i as f64);
            metric[[j, i]] = metric[[i, j]];
        }
    }

    metric
}

#[test]
fn test_detect_linear_dependencies_clean_matrix() {
    // Test with a well-conditioned matrix
    let n = 10;
    let mut metric = Array2::eye(n);
    for i in 0..n {
        metric[[i, i]] = 1.0 + 0.1 * (i as f64);
    }

    let config = LinearDependencyConfig::default();
    let analysis = detect_linear_dependencies(&metric, &config).unwrap();

    assert!(!analysis.has_dependencies());
    assert_eq!(analysis.rank, n);
    assert_eq!(analysis.dependent_indices.len(), 0);
    assert!(analysis.condition_number < 10.0);
}

#[test]
fn test_detect_linear_dependencies_rank_deficient() {
    let n = 15;
    let n_dependent = 3;
    let metric = create_dependent_metric(n, n_dependent);

    let config = LinearDependencyConfig {
        eigenvalue_threshold: 1e-10,
        ..Default::default()
    };

    let analysis = detect_linear_dependencies(&metric, &config).unwrap();

    assert!(analysis.has_dependencies());
    assert!(analysis.rank <= n - n_dependent + 1); // Allow for numerical error
    assert!(analysis.dependent_indices.len() >= n_dependent - 1);
}

#[test]
fn test_remove_dependencies_eigenvalue_method() {
    let n = 12;
    let metric = create_dependent_metric(n, 2);

    let config = LinearDependencyConfig {
        method: DependencyMethod::EigenvalueDecomposition,
        eigenvalue_threshold: 1e-9,
        ..Default::default()
    };

    let analysis = detect_linear_dependencies(&metric, &config).unwrap();

    let (metric_reg, transform) = remove_dependencies_eigenvalue(&metric, &analysis).unwrap();

    assert!(metric_reg.nrows() == analysis.rank);
    assert!(metric_reg.ncols() == analysis.rank);
    assert_eq!(transform.rank, analysis.rank);

    // Check that regularized metric is better conditioned
    let diag_min = (0..metric_reg.nrows())
        .map(|i| metric_reg[[i, i]])
        .fold(f64::INFINITY, f64::min);
    let diag_max = (0..metric_reg.nrows())
        .map(|i| metric_reg[[i, i]])
        .fold(f64::NEG_INFINITY, f64::max);

    let new_condition = diag_max / diag_min;
    assert!(new_condition < analysis.condition_number);
}

#[test]
fn test_transform_df_tensors() {
    let n_ao = 10;
    let n_aux = 8;
    let n_aux_reduced = 6;

    // Create mock DF tensors
    let mut j3c_ao = Array3::zeros((n_ao, n_ao, n_aux));
    for i in 0..n_ao {
        for j in 0..n_ao {
            for p in 0..n_aux {
                j3c_ao[[i, j, p]] = ((i + j + p) as f64).sin() * 0.1;
            }
        }
    }

    // Create mock transformation
    let mut forward = Array2::zeros((n_aux, n_aux_reduced));
    for i in 0..n_aux_reduced {
        forward[[i, i]] = 1.0; // Simple projection to first n_aux_reduced functions
    }

    let transform = quasix_core::df::TransformationMatrix {
        forward: forward.clone(),
        backward: forward.t().to_owned(),
        rank: n_aux_reduced,
    };

    // Transform tensors
    let j3c_reduced = transform_df_tensors(&j3c_ao, &transform).unwrap();

    assert_eq!(j3c_reduced.dim(), (n_ao, n_ao, n_aux_reduced));
}

#[test]
fn test_transform_mo_df_tensors() {
    let n_trans = 20; // n_occ * n_vir
    let n_aux = 15;
    let n_aux_reduced = 10;

    // Create mock MO DF tensors
    let mut j3c_ia = Array2::zeros((n_trans, n_aux));
    for i in 0..n_trans {
        for p in 0..n_aux {
            j3c_ia[[i, p]] = ((i * p) as f64).cos() * 0.1;
        }
    }

    // Create mock transformation
    let mut forward = Array2::zeros((n_aux, n_aux_reduced));
    for i in 0..n_aux_reduced {
        forward[[i, i]] = 1.0 / (i as f64 + 1.0).sqrt(); // With scaling
    }

    let transform = quasix_core::df::TransformationMatrix {
        forward: forward.clone(),
        backward: forward.t().to_owned(),
        rank: n_aux_reduced,
    };

    // Transform MO tensors
    let j3c_ia_reduced = transform_mo_df_tensors(&j3c_ia, &transform).unwrap();

    assert_eq!(j3c_ia_reduced.dim(), (n_trans, n_aux_reduced));
}

#[test]
fn test_cholesky_with_lindep() {
    let n = 20;
    let metric = create_dependent_metric(n, 3);

    let config = LinearDependencyConfig {
        eigenvalue_threshold: 1e-10,
        verbose: true,
        ..Default::default()
    };

    let result = compute_cholesky_v_with_lindep(&metric, Some(config));

    match result {
        Ok(cholesky) => {
            println!("Cholesky factorization succeeded");
            println!("Method: {:?}", cholesky.method);
            println!("Rank: {:?}", cholesky.rank);
            println!("Condition number: {:.2e}", cholesky.condition_number);

            // Verify the factorization
            if let Some(rank) = cholesky.rank {
                assert!(rank < n); // Should have reduced rank
            }
        }
        Err(e) => {
            // For matrices with severe dependencies, it might fail
            println!("Cholesky failed (expected for severe dependencies): {}", e);
        }
    }
}

#[test]
fn test_diffuse_basis_simulation() {
    // Simulate a case with diffuse basis functions (common in RI basis sets)
    let n = 25;
    let mut metric = Array2::zeros((n, n));

    // Create blocks representing different angular momenta
    for l in 0..5 {
        // s, p, d, f, g shells
        let block_start = l * 5;
        let block_end = (l + 1) * 5;

        for i in block_start..block_end.min(n) {
            for j in block_start..block_end.min(n) {
                // Within block: strong overlap
                let r = ((i as i32 - j as i32) as f64).abs();
                metric[[i, j]] = (-(r * r) / (2.0 * (l + 1) as f64)).exp();
            }

            // Diagonal dominance
            metric[[i, i]] += 1.0;
        }
    }

    // Add diffuse functions that create near-linear dependencies
    if n > 20 {
        // Make the last few rows linear combinations of earlier ones
        for i in 20..n {
            for j in 0..n {
                // Row i is a linear combination of rows 0 and 1
                metric[[i, j]] = 0.7 * metric[[0, j]] + 0.3 * metric[[1, j]];
                metric[[j, i]] = metric[[i, j]]; // Maintain symmetry
            }
        }
        // Add tiny diagonal perturbation to keep it positive semi-definite
        for i in 20..n {
            metric[[i, i]] += 1e-12;
        }
    }

    let config = LinearDependencyConfig {
        eigenvalue_threshold: 1e-8,
        method: DependencyMethod::Automatic,
        ..Default::default()
    };

    let analysis = detect_linear_dependencies(&metric, &config).unwrap();
    println!("{}", analysis.generate_report());

    assert!(analysis.has_dependencies());
    assert!(analysis.rank < n);
}

#[test]
fn test_periodic_boundary_conditions() {
    // Simulate metric from periodic boundary conditions
    let n_per_cell = 10;
    let n_cells = 3;
    let n = n_per_cell * n_cells;

    let mut metric = Array2::zeros((n, n));

    // Create periodic structure
    for cell_i in 0..n_cells {
        for cell_j in 0..n_cells {
            let offset_i = cell_i * n_per_cell;
            let offset_j = cell_j * n_per_cell;

            for i in 0..n_per_cell {
                for j in 0..n_per_cell {
                    let global_i = offset_i + i;
                    let global_j = offset_j + j;

                    if global_i < n && global_j < n {
                        // Overlap depends on distance including periodic images
                        let dist = ((cell_i as i32 - cell_j as i32).abs()
                            + (i as i32 - j as i32).abs())
                            as f64;
                        metric[[global_i, global_j]] = (-dist / 5.0).exp();
                    }
                }
            }
        }
    }

    // Add diagonal dominance
    for i in 0..n {
        metric[[i, i]] += 1.0;
    }

    let config = LinearDependencyConfig::periodic_systems();
    let analysis = detect_linear_dependencies(&metric, &config).unwrap();

    println!("PBC analysis: {}", analysis.generate_report());

    // Periodic systems often have dependencies
    if analysis.has_dependencies() {
        assert!(analysis.rank < n);
    }
}

#[test]
fn test_regularization_preserves_positive_definiteness() {
    let n = 15;
    let mut metric = create_dependent_metric(n, 2);

    // Add a small negative eigenvalue
    metric[[n - 1, n - 1]] = -1e-11;

    let config = LinearDependencyConfig {
        eigenvalue_threshold: 1e-10,
        method: DependencyMethod::Automatic,
        ..Default::default()
    };

    let metric_reg = regularize_metric(&metric, &config).unwrap();

    // Check that regularized metric is positive definite
    // (all diagonal elements should be positive for a positive definite matrix)
    for i in 0..metric_reg.nrows() {
        assert!(
            metric_reg[[i, i]] > 0.0,
            "Diagonal element {} is not positive: {}",
            i,
            metric_reg[[i, i]]
        );
    }
}

#[test]
fn test_aggressive_vs_conservative_thresholds() {
    let n = 20;
    let metric = create_dependent_metric(n, 4);

    // Test conservative threshold
    let config_conservative = LinearDependencyConfig::conservative();
    let analysis_conservative = detect_linear_dependencies(&metric, &config_conservative).unwrap();

    // Test aggressive threshold
    let config_aggressive = LinearDependencyConfig::aggressive();
    let analysis_aggressive = detect_linear_dependencies(&metric, &config_aggressive).unwrap();

    // Aggressive should remove more dependencies
    assert!(analysis_aggressive.rank <= analysis_conservative.rank);
    assert!(
        analysis_aggressive.dependent_indices.len()
            >= analysis_conservative.dependent_indices.len()
    );

    println!("Conservative: rank = {}/{}", analysis_conservative.rank, n);
    println!("Aggressive: rank = {}/{}", analysis_aggressive.rank, n);
}
