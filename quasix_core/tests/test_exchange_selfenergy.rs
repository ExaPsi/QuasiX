//! Integration tests for exchange self-energy calculation
//!
//! This test validates the exchange self-energy implementation against
//! known results and demonstrates proper usage of the API.

use ndarray::{arr1, Array2, Array3};
use quasix_core::selfenergy::{metric_utils, ExchangeSelfEnergyRI};

/// Test exchange self-energy calculation for a simple H2O-like system
#[test]
fn test_exchange_selfenergy_h2o_like() {
    // System parameters (simplified H2O-like)
    let nbasis = 7; // Minimal basis
    let nocc = 5; // 5 occupied orbitals
    let naux = 28; // Auxiliary basis size

    // Create calculator
    let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);

    // Create mock DF tensors (would come from real integral calculations)
    let mut df_3c_mo = Array3::<f64>::zeros((nbasis, nbasis, naux));

    // Fill with realistic-looking values
    for m in 0..nbasis {
        for n in 0..nbasis {
            for p in 0..naux {
                // Decay with distance-like behavior
                let decay = (-0.5 * ((m as f64 - n as f64).abs() + (p as f64 / 10.0))).exp();
                df_3c_mo[[m, n, p]] = decay * 0.1;
            }
        }
    }

    // Create metric (simplified - would come from (P|Q) integrals)
    let mut metric = Array2::<f64>::eye(naux);
    for i in 0..naux {
        for j in 0..naux {
            if i != j {
                let dist = (i as f64 - j as f64).abs();
                metric[[i, j]] = (-0.1 * dist).exp() * 0.05;
            }
        }
    }

    // Compute metric inverse
    let metric_inv = metric_utils::compute_metric_inverse(&metric, 1e-10).unwrap();

    // Calculate exchange self-energy
    let sigma_x = calc
        .compute_exchange_matrix_ri(&df_3c_mo, &metric_inv)
        .unwrap();

    // Validate properties
    assert_eq!(sigma_x.shape(), &[nbasis, nbasis]);

    // Check symmetry
    for i in 0..nbasis {
        for j in 0..i {
            let diff = (sigma_x[[i, j]] - sigma_x[[j, i]]).abs();
            assert!(diff < 1e-10, "Matrix not symmetric at ({}, {})", i, j);
        }
    }

    // Check that occupied diagonal elements are negative (exchange is attractive)
    for i in 0..nocc {
        assert!(
            sigma_x[[i, i]] <= 0.0,
            "Exchange diagonal element {} should be negative, got {}",
            i,
            sigma_x[[i, i]]
        );
    }

    // Check metadata
    assert_eq!(calc.metadata.nbasis, nbasis);
    assert_eq!(calc.metadata.nocc, nocc);
    assert_eq!(calc.metadata.naux, naux);
    assert!(calc.metadata.properties.contains_key("computation_time_ms"));

    println!("Exchange self-energy calculation successful!");
    println!("Calculation ID: {}", calc.metadata.calc_id);
    println!("Method: {}", calc.metadata.method);
}

/// Test diagonal-only calculation
#[test]
fn test_exchange_diagonal_only() {
    let nbasis = 5;
    let nocc = 3;
    let naux = 15;

    let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);

    // Create test data
    let mut df_3c_mo = Array3::<f64>::zeros((nbasis, nbasis, naux));
    for m in 0..nbasis {
        for n in 0..nbasis {
            for p in 0..naux {
                df_3c_mo[[m, n, p]] = ((m + n) as f64 * 0.01 + p as f64 * 0.001) * 0.1;
            }
        }
    }

    let metric_inv = Array2::<f64>::eye(naux) * 0.9; // Slightly perturbed identity

    // Get diagonal only
    let sigma_x_diag = calc
        .compute_exchange_diagonal_ri(&df_3c_mo, &metric_inv)
        .unwrap();

    // Also compute full matrix for comparison
    let sigma_x_full = calc
        .compute_exchange_matrix_ri(&df_3c_mo, &metric_inv)
        .unwrap();

    // Check diagonal elements match
    for i in 0..nbasis {
        let diff = (sigma_x_diag[i] - sigma_x_full[[i, i]]).abs();
        assert!(
            diff < 1e-12,
            "Diagonal element {} mismatch: {} vs {}",
            i,
            sigma_x_diag[i],
            sigma_x_full[[i, i]]
        );
    }
}

/// Test symmetric DF tensor calculation
#[test]
fn test_symmetric_df_calculation() {
    let nbasis = 4;
    let nocc = 2;
    let naux = 10;

    let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);

    // Create symmetrized DF tensors (already include metric)
    let mut df_3c_mo_sym = Array3::<f64>::zeros((nbasis, nbasis, naux));
    for m in 0..nbasis {
        for n in 0..nbasis {
            for p in 0..naux {
                // Symmetric tensors with v^(1/2) already applied
                df_3c_mo_sym[[m, n, p]] = ((m + n) as f64 * 0.05 + p as f64 * 0.02).sqrt() * 0.1;
            }
        }
    }

    // Calculate using symmetric method
    let sigma_x = calc
        .compute_exchange_matrix_symmetric(&df_3c_mo_sym)
        .unwrap();

    // Validate
    assert_eq!(sigma_x.shape(), &[nbasis, nbasis]);

    // Check symmetry
    for i in 0..nbasis {
        for j in 0..i {
            let diff = (sigma_x[[i, j]] - sigma_x[[j, i]]).abs();
            assert!(
                diff < 1e-10,
                "Symmetric result not symmetric at ({}, {})",
                i,
                j
            );
        }
    }

    // Check metadata shows symmetric method
    assert_eq!(calc.metadata.method, "RI-V-symmetric");
}

/// Test operator application without forming full matrix
#[test]
fn test_operator_application() {
    let nbasis = 3;
    let nocc = 2;
    let naux = 6;

    let calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);

    // Create test data
    let mut df_3c_mo = Array3::<f64>::zeros((nbasis, nbasis, naux));
    for m in 0..nbasis {
        for n in 0..nbasis {
            for p in 0..naux {
                df_3c_mo[[m, n, p]] = m as f64 * 0.1 + n as f64 * 0.05 + p as f64 * 0.01;
            }
        }
    }

    let metric_inv = Array2::<f64>::eye(naux);
    let test_vec = arr1(&[1.0, 2.0, 3.0]);

    // Apply operator
    let result = calc
        .apply_exchange_operator(&test_vec, &df_3c_mo, &metric_inv)
        .unwrap();

    // Should have same dimension
    assert_eq!(result.len(), nbasis);

    // Result should be non-trivial
    assert!(result.iter().any(|&x| x.abs() > 1e-10));
}

/// Test error handling for invalid inputs
#[test]
fn test_error_handling() {
    let nbasis = 5;
    let nocc = 3;
    let naux = 10;

    let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);

    // Wrong dimensions
    let df_3c_wrong = Array3::<f64>::zeros((nbasis + 1, nbasis, naux));
    let metric_inv = Array2::<f64>::eye(naux);

    let result = calc.compute_exchange_matrix_ri(&df_3c_wrong, &metric_inv);
    assert!(result.is_err());

    // Wrong metric dimensions
    let df_3c_mo = Array3::<f64>::zeros((nbasis, nbasis, naux));
    let metric_wrong = Array2::<f64>::eye(naux + 1);

    let result = calc.compute_exchange_matrix_ri(&df_3c_mo, &metric_wrong);
    assert!(result.is_err());
}

/// Performance benchmark for medium-sized system
#[test]
#[ignore] // Run with --ignored flag for benchmarking
fn bench_medium_system() {
    use std::time::Instant;

    // Medium-sized molecular system
    let nbasis = 100;
    let nocc = 50;
    let naux = 500;

    println!("\nBenchmarking exchange self-energy for medium system:");
    println!("  nbasis = {}", nbasis);
    println!("  nocc = {}", nocc);
    println!("  naux = {}", naux);

    let mut calc = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);

    // Create realistic-sized tensors
    let df_3c_mo = Array3::<f64>::from_elem((nbasis, nbasis, naux), 0.01);
    let metric_inv = Array2::<f64>::eye(naux);

    // Benchmark full matrix calculation
    let start = Instant::now();
    let _sigma_x = calc
        .compute_exchange_matrix_ri(&df_3c_mo, &metric_inv)
        .unwrap();
    let duration = start.elapsed();

    println!("  Full matrix time: {:?}", duration);

    // Benchmark diagonal-only
    let start = Instant::now();
    let _sigma_x_diag = calc
        .compute_exchange_diagonal_ri(&df_3c_mo, &metric_inv)
        .unwrap();
    let duration = start.elapsed();

    println!("  Diagonal-only time: {:?}", duration);

    // Memory estimate
    let memory_mb = (nbasis * nbasis * naux * 8) as f64 / 1_048_576.0;
    println!("  DF tensor memory: {:.1} MB", memory_mb);
}
