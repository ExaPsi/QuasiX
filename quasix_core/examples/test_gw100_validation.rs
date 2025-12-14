//! Example demonstrating GW100 validation functionality

use quasix_core::validation::gw100::{compute_validation_stats, BenchmarkResult, SimdStatistics};

fn main() {
    println!("=== GW100 Validation Module Test ===\n");

    // Test 1: SIMD MAD computation
    println!("Test 1: SIMD-accelerated MAD computation");
    let calculated = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let reference = vec![1.1, 1.9, 3.2, 3.8, 5.1, 5.9, 7.2, 7.8];

    match SimdStatistics::compute_mad_simd(&calculated, &reference) {
        Ok(mad) => println!("  MAD = {:.6} eV", mad),
        Err(e) => println!("  Error: {}", e),
    }

    // Test 2: Correlation computation
    println!("\nTest 2: Correlation coefficient");
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![1.05, 2.05, 3.05, 4.05, 5.05];

    match SimdStatistics::compute_correlation(&x, &y) {
        Ok(r2) => println!("  R² = {:.6}", r2),
        Err(e) => println!("  Error: {}", e),
    }

    // Test 3: RMSE computation
    println!("\nTest 3: Root Mean Square Error");
    match SimdStatistics::compute_rmse(&calculated, &reference) {
        Ok(rmse) => println!("  RMSE = {:.6} eV", rmse),
        Err(e) => println!("  Error: {}", e),
    }

    // Test 4: Outlier detection
    println!("\nTest 4: Outlier detection");
    let values = vec![0.1, 0.12, 0.09, 0.11, 0.10, 1.5, 0.08, 0.13, 0.095];
    let outliers = SimdStatistics::detect_outliers(&values, 3.5);
    println!(
        "  Found {} outlier(s) at indices: {:?}",
        outliers.len(),
        outliers
    );

    // Test 5: Benchmark result
    println!("\nTest 5: Benchmark result creation");
    let result = BenchmarkResult::new(
        "H2O".to_string(),
        "def2-tzvp".to_string(),
        "def2-tzvp-jkfit".to_string(),
        vec![-15.0, -10.0, -5.0, 1.0, 2.0],
        vec![0.8, 0.85, 0.9, 0.92, 0.95],
        vec![-15.1, -10.05, -5.02, 1.05, 2.1],
        2, // HOMO
        3, // LUMO
        10.5,
    );

    println!("  Molecule: {}", result.molecule);
    println!("  IP error: {:.3} eV", result.ip_error);
    println!("  EA error: {:.3} eV", result.ea_error);
    println!("  MAD: {:.3} eV", result.mad());
    println!(
        "  Passes 0.1 eV threshold: {}",
        result.passes_threshold(0.1)
    );

    // Test 6: Validation statistics
    println!("\nTest 6: Computing validation statistics");
    let results = vec![
        BenchmarkResult::new(
            "H2O".to_string(),
            "def2-tzvp".to_string(),
            "def2-tzvp-jkfit".to_string(),
            vec![-15.0, -10.0, -5.0],
            vec![0.8, 0.85, 0.9],
            vec![-15.1, -10.1, -5.05],
            1,
            2,
            10.0,
        ),
        BenchmarkResult::new(
            "NH3".to_string(),
            "def2-tzvp".to_string(),
            "def2-tzvp-jkfit".to_string(),
            vec![-12.0, -8.0, -3.0],
            vec![0.82, 0.88, 0.91],
            vec![-12.15, -8.1, -3.08],
            1,
            2,
            12.0,
        ),
    ];

    match compute_validation_stats(&results) {
        Ok(stats) => {
            println!("  {}", stats.summary());
            println!(
                "  Passes validation (MAD < 0.2, R² > 0.9): {}",
                stats.passes_validation(0.2, 0.9)
            );
        }
        Err(e) => println!("  Error computing stats: {}", e),
    }

    println!("\n=== All tests completed ===");
}
