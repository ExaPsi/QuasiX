//! Benchmark for integral screening in exchange self-energy calculations
//!
//! This example demonstrates the performance improvements achieved through
//! Schwarz-based integral screening for sparse systems.

use ndarray::{Array2, Array3};
use quasix_core::selfenergy::{screening::ScreeningConfig, ExchangeSelfEnergyRI};
use std::time::Instant;

fn main() {
    println!("\n=== Exchange Self-Energy Screening Benchmark ===\n");

    // Test different system sizes
    let test_cases = [
        (20, 5, 40, "Small molecule (H2O-like)"),
        (50, 10, 100, "Medium molecule (benzene-like)"),
        (100, 20, 200, "Large molecule (naphthalene-like)"),
        (200, 40, 400, "Very large molecule"),
    ];

    for (nbasis, nocc, naux, description) in test_cases {
        println!("System: {}", description);
        println!("  nbasis={}, nocc={}, naux={}", nbasis, nocc, naux);

        // Create sparse DF tensor (only ~20% significant elements)
        let df_3c_mo = create_sparse_df_tensor(nbasis, nocc, naux);
        let df_metric_inv = Array2::<f64>::eye(naux);

        // Benchmark without screening
        let mut calc_no_screen = ExchangeSelfEnergyRI::new(nbasis, nocc, naux);
        calc_no_screen.set_screening(ScreeningConfig::disabled());

        let start = Instant::now();
        let _sigma_x_no_screen = calc_no_screen
            .compute_exchange_matrix_screened(&df_3c_mo, &df_metric_inv)
            .expect("Failed to compute exchange matrix without screening");
        let time_no_screen = start.elapsed();

        // Benchmark with screening (threshold 1e-12)
        let mut calc_with_screen = ExchangeSelfEnergyRI::with_screening(
            nbasis,
            nocc,
            naux,
            ScreeningConfig::with_threshold(1e-12),
        );

        let start = Instant::now();
        let _sigma_x_with_screen = calc_with_screen
            .compute_exchange_matrix_screened(&df_3c_mo, &df_metric_inv)
            .expect("Failed to compute exchange matrix with screening");
        let time_with_screen = start.elapsed();

        // Get screening statistics
        let n_computed = calc_with_screen
            .metadata
            .properties
            .get("pairs_computed")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let n_skipped = calc_with_screen
            .metadata
            .properties
            .get("pairs_skipped")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let total_pairs = n_computed + n_skipped;
        let screening_ratio = if total_pairs > 0 {
            (n_skipped as f64 / total_pairs as f64) * 100.0
        } else {
            0.0
        };

        // Calculate speedup
        let speedup = time_no_screen.as_secs_f64() / time_with_screen.as_secs_f64();

        println!("  Without screening: {:?}", time_no_screen);
        println!("  With screening:    {:?}", time_with_screen);
        println!("  Speedup:           {:.2}x", speedup);
        println!(
            "  Pairs computed:    {}/{} ({:.1}% screened)",
            n_computed, total_pairs, screening_ratio
        );
        println!();

        // Additional benchmark with different thresholds
        if nbasis <= 50 {
            // Only for smaller systems to save time
            println!("  Threshold sensitivity:");
            for threshold in &[1e-6, 1e-8, 1e-10, 1e-12, 1e-14] {
                let mut calc = ExchangeSelfEnergyRI::with_screening(
                    nbasis,
                    nocc,
                    naux,
                    ScreeningConfig::with_threshold(*threshold),
                );

                let start = Instant::now();
                let sigma_x = calc
                    .compute_exchange_matrix_screened(&df_3c_mo, &df_metric_inv)
                    .expect("Failed to compute exchange matrix");
                let time = start.elapsed();

                // Check accuracy (Frobenius norm of matrix)
                let norm: f64 = sigma_x.iter().map(|x| x * x).sum::<f64>().sqrt();

                let n_skipped = calc
                    .metadata
                    .properties
                    .get("pairs_skipped")
                    .and_then(|s| s.parse::<usize>().ok())
                    .unwrap_or(0);

                println!(
                    "    Threshold {:.0e}: time={:?}, norm={:.6}, skipped={}",
                    threshold, time, norm, n_skipped
                );
            }
            println!();
        }
    }

    println!("\n=== Screening Effectiveness for Different Sparsity Patterns ===\n");

    // Test different sparsity patterns
    let nbasis = 50;
    let nocc = 10;
    let naux = 100;

    for (sparsity, pattern_name) in &[
        (0.05, "Very sparse (5% significant)"),
        (0.20, "Sparse (20% significant)"),
        (0.50, "Moderate (50% significant)"),
        (0.80, "Dense (80% significant)"),
        (1.00, "Full (100% significant)"),
    ] {
        println!("Pattern: {}", pattern_name);

        let df_3c_mo = create_df_tensor_with_sparsity(nbasis, nocc, naux, *sparsity);
        let df_metric_inv = Array2::<f64>::eye(naux);

        let mut calc = ExchangeSelfEnergyRI::with_screening(
            nbasis,
            nocc,
            naux,
            ScreeningConfig::with_threshold(1e-10),
        );

        let start = Instant::now();
        let _sigma_x = calc
            .compute_exchange_matrix_screened(&df_3c_mo, &df_metric_inv)
            .expect("Failed to compute exchange matrix");
        let time = start.elapsed();

        let screening_ratio = calc
            .metadata
            .properties
            .get("screening_ratio")
            .unwrap_or(&"0%".to_string())
            .clone();

        println!("  Time: {:?}, Screening: {}", time, screening_ratio);
    }

    println!("\n=== Benchmark Complete ===\n");
}

/// Create a sparse DF tensor with only a fraction of significant elements
fn create_sparse_df_tensor(nbasis: usize, nocc: usize, naux: usize) -> Array3<f64> {
    let mut df_3c_mo = Array3::<f64>::zeros((nbasis, nbasis, naux));

    // Create block structure: significant elements only in lower indices
    let significant_block = nbasis / 5; // 20% of orbitals have significant integrals

    for m in 0..nbasis {
        for n in 0..nbasis {
            for p in 0..naux {
                if m < significant_block && n < significant_block {
                    // Significant elements
                    df_3c_mo[[m, n, p]] = 1.0 + ((m + n + p) as f64 * 0.1);
                } else if m < nocc || n < nocc {
                    // Small but non-zero elements in occupied block
                    df_3c_mo[[m, n, p]] = 1e-8;
                } else {
                    // Very small elements elsewhere
                    df_3c_mo[[m, n, p]] = 1e-15;
                }
            }
        }
    }

    df_3c_mo
}

/// Create DF tensor with specified sparsity level
fn create_df_tensor_with_sparsity(
    nbasis: usize,
    _nocc: usize,
    naux: usize,
    sparsity: f64,
) -> Array3<f64> {
    let mut df_3c_mo = Array3::<f64>::zeros((nbasis, nbasis, naux));
    let threshold_idx = ((nbasis as f64) * sparsity) as usize;

    for m in 0..nbasis {
        for n in 0..nbasis {
            for p in 0..naux {
                if m < threshold_idx && n < threshold_idx {
                    // Significant elements
                    df_3c_mo[[m, n, p]] = 1.0 + ((m + n + p) as f64 * 0.01);
                } else {
                    // Small elements
                    df_3c_mo[[m, n, p]] = 1e-14;
                }
            }
        }
    }

    df_3c_mo
}
