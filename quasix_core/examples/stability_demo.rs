//! Demonstration of numerical stability enhancements in frequency grids
//!
//! Run with: cargo run --example stability_demo --release

use quasix_core::freq::{
    stability_validation::{benchmark_stability_improvements, StabilityValidator},
    FrequencyGrid, GridType, KahanSum,
};

fn main() {
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║   QuasiX S3-1: Numerical Stability Enhancements Demo      ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");

    // 1. Demonstrate catastrophic cancellation prevention
    demonstrate_cancellation_prevention();

    // 2. Show condition number monitoring
    demonstrate_condition_monitoring();

    // 3. Compare standard vs stabilized integration
    compare_integration_methods();

    // 4. Run comprehensive validation
    run_validation_suite();

    // 5. Run performance benchmarks
    benchmark_stability_improvements();
}

fn demonstrate_cancellation_prevention() {
    println!("\n┌─────────────────────────────────────────────────────────────┐");
    println!("│ 1. CATASTROPHIC CANCELLATION PREVENTION                    │");
    println!("└─────────────────────────────────────────────────────────────┘\n");

    // Example: Large values that cancel
    let values = vec![1e20, 1.23456789, -1e20, 2.34567890, 1e-10];

    println!("Test values: [1e20, 1.23456789, -1e20, 2.34567890, 1e-10]");
    println!("Expected sum: 3.58024679 + 1e-10\n");

    // Standard summation
    let standard_sum: f64 = values.iter().sum();
    println!("Standard summation:     {:.15}", standard_sum);

    // Kahan summation
    let mut kahan = KahanSum::new();
    for &v in &values {
        kahan.add(v);
    }
    let kahan_sum = kahan.sum();
    println!("Kahan summation:        {:.15}", kahan_sum);

    // Error comparison
    let expected = 3.58024679 + 1e-10;
    println!("\nError comparison:");
    println!("  Standard error: {:.2e}", (standard_sum - expected).abs());
    println!("  Kahan error:    {:.2e}", (kahan_sum - expected).abs());

    let improvement = (standard_sum - expected).abs() / ((kahan_sum - expected).abs() + 1e-20);
    println!("  Improvement:    {:.0}x better accuracy", improvement);
}

fn demonstrate_condition_monitoring() {
    println!("\n┌─────────────────────────────────────────────────────────────┐");
    println!("│ 2. CONDITION NUMBER MONITORING                             │");
    println!("└─────────────────────────────────────────────────────────────┘\n");

    let omega_max_values = vec![10.0, 100.0, 1000.0, 10000.0];

    println!("Grid condition numbers for different omega_max:\n");
    println!("  omega_max  | n_points | condition number | status");
    println!("  -----------|----------|------------------|--------");

    for omega_max in omega_max_values {
        let grid = FrequencyGrid::new(50, GridType::ModifiedGaussLegendre { omega_max }).unwrap();

        let status = if let Some(cond) = grid.condition_number {
            if cond < 1e6 {
                "Good    "
            } else if cond < 1e8 {
                "Warning "
            } else {
                "Critical"
            }
        } else {
            "Unknown "
        };

        println!(
            "  {:9.1} |    50    | {:15.2e} | {}",
            omega_max,
            grid.condition_number.unwrap_or(0.0),
            status
        );
    }

    println!("\nNote: Condition numbers > 1e8 indicate potential accuracy loss");
}

fn compare_integration_methods() {
    println!("\n┌─────────────────────────────────────────────────────────────┐");
    println!("│ 3. INTEGRATION ACCURACY COMPARISON                         │");
    println!("└─────────────────────────────────────────────────────────────┘\n");

    // Create a challenging integrand with cancellation
    let grid = FrequencyGrid::new(64, GridType::GaussLegendre).unwrap();

    println!("Test: Integrate (1e8 * sin(x) + x^2) over [-1, 1]");
    println!("      Large oscillatory term should cancel, leaving 2/3\n");

    // Standard integration (without Kahan)
    let standard_result: f64 = grid
        .points
        .iter()
        .zip(grid.weights.iter())
        .map(|(&x, &w)| w * (1e8 * x.sin() + x * x))
        .sum();

    // Stabilized integration (with Kahan)
    let stabilized_result = grid.integrate(|x| 1e8 * x.sin() + x * x);

    let exact_value = 2.0 / 3.0; // ∫_{-1}^{1} x^2 dx = 2/3

    println!("Results:");
    println!("  Exact value:              {:.15}", exact_value);
    println!("  Standard integration:     {:.15}", standard_result);
    println!("  Stabilized integration:   {:.15}", stabilized_result);
    println!("\nErrors:");
    println!(
        "  Standard error:   {:.2e}",
        (standard_result - exact_value).abs()
    );
    println!(
        "  Stabilized error: {:.2e}",
        (stabilized_result - exact_value).abs()
    );

    let improvement =
        (standard_result - exact_value).abs() / ((stabilized_result - exact_value).abs() + 1e-20);
    println!("  Improvement:      {:.0}x better accuracy", improvement);
}

fn run_validation_suite() {
    println!("\n┌─────────────────────────────────────────────────────────────┐");
    println!("│ 4. COMPREHENSIVE STABILITY VALIDATION                      │");
    println!("└─────────────────────────────────────────────────────────────┘\n");

    let validator = StabilityValidator::new().with_verbose(false);

    match validator.run_full_validation() {
        Ok(summary) => {
            println!("Validation Results:");
            println!("  Total tests:  {}", summary.total_tests());
            println!("  Passed:       {}", summary.passed_tests());
            println!(
                "  Failed:       {}",
                summary.total_tests() - summary.passed_tests()
            );

            if summary.all_passed() {
                println!("\n✓ All stability tests PASSED");
            } else {
                println!("\n✗ Some tests failed - review needed");
            }
        }
        Err(e) => {
            println!("Error running validation: {}", e);
        }
    }
}

/// Additional test: Demonstrate stability at extreme frequencies
#[allow(dead_code)]
fn test_extreme_frequencies() {
    println!("\n┌─────────────────────────────────────────────────────────────┐");
    println!("│ EXTREME FREQUENCY STABILITY TEST                           │");
    println!("└─────────────────────────────────────────────────────────────┘\n");

    let omega_max = 1e6; // Very large frequency
    let n_points = 100;

    println!("Creating grid with omega_max = {:.2e}...", omega_max);

    match FrequencyGrid::new(n_points, GridType::ModifiedGaussLegendre { omega_max }) {
        Ok(grid) => {
            // Check that all values are finite
            let all_finite = grid.points.iter().all(|&x| x.is_finite())
                && grid.weights.iter().all(|&w| w.is_finite());

            if all_finite {
                println!("✓ Successfully created stable grid");

                if let Some(cutoff) = grid.omega_stable {
                    println!("  Maximum stable frequency: {:.2e}", cutoff);
                }

                if let Some(cond) = grid.condition_number {
                    println!("  Condition number: {:.2e}", cond);
                }

                // Test integration stability
                let result = grid.integrate(|omega| omega.exp() * (-omega));
                println!("  Test integral converged: {}", result.is_finite());
            } else {
                println!("✗ Grid contains non-finite values");
            }
        }
        Err(e) => {
            println!("✗ Failed to create grid: {}", e);
        }
    }
}
