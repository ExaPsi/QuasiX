//! Performance verification tests for condition monitoring
//!
//! This test ensures that the optimized condition monitoring
//! meets the < 5% overhead requirement for production use.

use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::linalg::conditioning::{
    estimate_condition_power_iteration, apply_adaptive_regularization,
};
use quasix_core::linalg::conditioning_optimized::{
    OptimizedConditionMonitor, BatchConditionPool, MonitorConfig,
};
use rand::prelude::*;
use std::time::{Duration, Instant};

/// Generate test matrix with specified properties
fn generate_test_matrix(n: usize, condition_number: f64, hermitian: bool) -> Array2<Complex64> {
    let mut rng = StdRng::seed_from_u64(42);

    if hermitian {
        // Generate Hermitian matrix with controlled condition number
        let mut matrix = Array2::zeros((n, n));

        // Set eigenvalues to control condition number
        for i in 0..n {
            let eigenvalue = if i == 0 {
                condition_number
            } else if i == n - 1 {
                1.0
            } else {
                1.0 + rng.gen::<f64>() * (condition_number - 1.0).log10()
            };
            matrix[[i, i]] = Complex64::new(eigenvalue, 0.0);
        }

        // Add small off-diagonal elements
        for i in 0..n {
            for j in i + 1..n {
                let val = Complex64::new(
                    rng.gen::<f64>() * 0.01,
                    rng.gen::<f64>() * 0.01,
                );
                matrix[[i, j]] = val;
                matrix[[j, i]] = val.conj();
            }
        }

        matrix
    } else {
        // Generate general complex matrix
        let mut matrix = Array2::zeros((n, n));

        for i in 0..n {
            for j in 0..n {
                matrix[[i, j]] = Complex64::new(
                    rng.gen::<f64>() - 0.5,
                    rng.gen::<f64>() - 0.5,
                );
            }
        }

        // Scale to achieve approximate condition number
        let scale = condition_number.sqrt();
        matrix[[0, 0]] *= scale;

        matrix
    }
}

/// Simulate a typical GW computation workload
fn simulate_gw_workload(matrix: &Array2<Complex64>) -> f64 {
    let n = matrix.nrows();
    let mut result = 0.0;

    // Simulate matrix operations typical in GW
    for _ in 0..10 {
        // Matrix-vector multiplication (dominant operation)
        let v = Array1::from_elem(n, Complex64::new(1.0, 0.0));
        let w = matrix.dot(&v);
        result += w[0].re;
    }

    result
}

#[test]
fn test_overhead_well_conditioned() {
    println!("\n=== Testing Overhead for Well-Conditioned Matrices ===");

    for &n in &[50, 100, 200, 500, 1000] {
        let matrix = generate_test_matrix(n, 10.0, true);
        let monitor = OptimizedConditionMonitor::new(n);

        // Measure baseline workload time
        let start = Instant::now();
        let baseline_result = simulate_gw_workload(&matrix);
        let baseline_time = start.elapsed();

        // Measure workload with condition monitoring
        let mut matrix_copy = matrix.clone();
        let start = Instant::now();
        let condition = monitor.estimate_condition(&matrix_copy);
        if condition > 1e12 {
            monitor.regularize_inplace(&mut matrix_copy, condition);
        }
        let monitored_result = simulate_gw_workload(&matrix_copy);
        let monitored_time = start.elapsed();

        // Calculate overhead
        let overhead_percent = 100.0 * (monitored_time.as_secs_f64() - baseline_time.as_secs_f64())
            / baseline_time.as_secs_f64();

        println!(
            "N={:4}: baseline={:6.3}ms, monitored={:6.3}ms, overhead={:5.2}%, condition={:.2e}",
            n,
            baseline_time.as_secs_f64() * 1000.0,
            monitored_time.as_secs_f64() * 1000.0,
            overhead_percent,
            condition
        );

        // Verify overhead is < 5% for well-conditioned matrices
        assert!(
            overhead_percent < 5.0,
            "Overhead {}% exceeds 5% limit for N={}",
            overhead_percent,
            n
        );

        // Results should be similar
        assert!((baseline_result - monitored_result).abs() < 1e-10);
    }
}

#[test]
fn test_overhead_ill_conditioned() {
    println!("\n=== Testing Overhead for Ill-Conditioned Matrices ===");

    for &n in &[50, 100, 200, 500] {
        let matrix = generate_test_matrix(n, 1e10, true);
        let monitor = OptimizedConditionMonitor::new(n);

        // Measure baseline workload time
        let start = Instant::now();
        let baseline_result = simulate_gw_workload(&matrix);
        let baseline_time = start.elapsed();

        // Measure workload with condition monitoring
        let mut matrix_copy = matrix.clone();
        let start = Instant::now();
        let condition = monitor.estimate_condition(&matrix_copy);
        if condition > 1e12 {
            monitor.regularize_inplace(&mut matrix_copy, condition);
        }
        let monitored_result = simulate_gw_workload(&matrix_copy);
        let monitored_time = start.elapsed();

        // Calculate overhead
        let overhead_percent = 100.0 * (monitored_time.as_secs_f64() - baseline_time.as_secs_f64())
            / baseline_time.as_secs_f64();

        println!(
            "N={:4}: baseline={:6.3}ms, monitored={:6.3}ms, overhead={:5.2}%, condition={:.2e}",
            n,
            baseline_time.as_secs_f64() * 1000.0,
            monitored_time.as_secs_f64() * 1000.0,
            overhead_percent,
            condition
        );

        // For ill-conditioned matrices, overhead can be higher but should still be reasonable
        assert!(
            overhead_percent < 10.0,
            "Overhead {}% exceeds 10% limit for ill-conditioned N={}",
            overhead_percent,
            n
        );
    }
}

#[test]
fn test_batch_performance() {
    println!("\n=== Testing Batch Processing Performance ===");

    let n = 200;
    let n_freq = 50;

    // Generate batch of matrices with varying condition numbers
    let mut matrices: Vec<_> = (0..n_freq)
        .map(|i| {
            let condition = 10.0_f64.powf(2.0 + (i as f64) * 0.1);
            generate_test_matrix(n, condition, true)
        })
        .collect();

    // Test with different thread counts
    for &n_threads in &[1, 2, 4, 8] {
        let pool = BatchConditionPool::new(n, n_threads);

        let start = Instant::now();
        let conditions = pool.process_batch(&mut matrices);
        let elapsed = start.elapsed();

        let throughput = n_freq as f64 / elapsed.as_secs_f64();

        println!(
            "Threads={}: time={:6.3}ms, throughput={:5.1} matrices/sec",
            n_threads,
            elapsed.as_secs_f64() * 1000.0,
            throughput
        );

        // Verify all matrices were processed
        assert_eq!(conditions.len(), n_freq);

        // Verify reasonable condition estimates
        for (i, &cond) in conditions.iter().enumerate() {
            let expected = 10.0_f64.powf(2.0 + (i as f64) * 0.1);
            assert!(
                cond > expected * 0.1 && cond < expected * 10.0,
                "Condition estimate {} far from expected {} at index {}",
                cond, expected, i
            );
        }
    }
}

#[test]
fn test_fast_path_effectiveness() {
    println!("\n=== Testing Fast Path Effectiveness ===");

    let n = 500;
    let monitor = OptimizedConditionMonitor::new(n);

    // Test 1: Strongly diagonal dominant (should hit fast path)
    let mut matrix = Array2::zeros((n, n));
    for i in 0..n {
        matrix[[i, i]] = Complex64::new(10.0, 0.0);
        for j in 0..n {
            if i != j {
                matrix[[i, j]] = Complex64::new(0.001, 0.001);
            }
        }
    }

    let start = Instant::now();
    let cond1 = monitor.estimate_condition(&matrix);
    let time1 = start.elapsed();

    println!(
        "Diagonal dominant: condition={:.2e}, time={:.3}μs (fast path)",
        cond1,
        time1.as_nanos() as f64 / 1000.0
    );

    // Test 2: Dense matrix (should not hit fast path)
    for i in 0..n {
        for j in 0..n {
            if i != j {
                matrix[[i, j]] = Complex64::new(0.5, 0.5);
            }
        }
    }

    let start = Instant::now();
    let cond2 = monitor.estimate_condition(&matrix);
    let time2 = start.elapsed();

    println!(
        "Dense matrix: condition={:.2e}, time={:.3}μs (full estimate)",
        cond2,
        time2.as_nanos() as f64 / 1000.0
    );

    // Fast path should be at least 10x faster
    assert!(
        time1 < time2 / 10,
        "Fast path not sufficiently faster: {:?} vs {:?}",
        time1, time2
    );

    // Check statistics
    let stats = monitor.get_stats();
    assert_eq!(stats.fast_path_hits, 1);
    assert_eq!(stats.full_estimates, 1);
}

#[test]
fn test_accuracy_vs_speed_tradeoff() {
    println!("\n=== Testing Accuracy vs Speed Tradeoff ===");

    let n = 100;

    // Generate test matrix with known condition number
    let true_condition = 1e6;
    let matrix = generate_test_matrix(n, true_condition, true);

    // Test original implementation
    let start = Instant::now();
    let original_estimate = estimate_condition_power_iteration(&matrix, 30, 1e-6);
    let original_time = start.elapsed();

    // Test optimized implementation
    let monitor = OptimizedConditionMonitor::new(n);
    let start = Instant::now();
    let optimized_estimate = monitor.estimate_condition(&matrix);
    let optimized_time = start.elapsed();

    // Calculate relative error
    let relative_error = (optimized_estimate - original_estimate).abs() / original_estimate;
    let speedup = original_time.as_secs_f64() / optimized_time.as_secs_f64();

    println!(
        "Original: estimate={:.2e}, time={:.3}ms",
        original_estimate,
        original_time.as_secs_f64() * 1000.0
    );
    println!(
        "Optimized: estimate={:.2e}, time={:.3}ms",
        optimized_estimate,
        optimized_time.as_secs_f64() * 1000.0
    );
    println!("Relative error: {:.2}%", relative_error * 100.0);
    println!("Speedup: {:.1}x", speedup);

    // Optimized should be within 50% accuracy
    assert!(
        relative_error < 0.5,
        "Accuracy loss too high: {:.2}%",
        relative_error * 100.0
    );

    // Optimized should be at least 2x faster
    assert!(
        speedup > 2.0,
        "Insufficient speedup: {:.1}x",
        speedup
    );
}

#[test]
fn test_large_matrix_performance() {
    println!("\n=== Testing Large Matrix Performance (N=5000) ===");

    let n = 5000;

    // Use custom config for large matrices
    let mut config = MonitorConfig::default();
    config.max_power_iterations = 10; // Fewer iterations for speed
    config.power_tolerance = 1e-2; // Relaxed tolerance
    config.block_size = 128; // Larger blocks for big matrices

    let monitor = OptimizedConditionMonitor::with_config(n, config);

    // Generate sparse-ish matrix (typical for QuasiX)
    let mut matrix = Array2::zeros((n, n));
    let mut rng = StdRng::seed_from_u64(42);

    // Diagonal elements
    for i in 0..n {
        matrix[[i, i]] = Complex64::new(1.0 + rng.gen::<f64>(), 0.0);
    }

    // Add some off-diagonal elements (sparse pattern)
    for _ in 0..n * 10 {
        let i = rng.gen_range(0..n);
        let j = rng.gen_range(0..n);
        if i != j {
            let val = Complex64::new(
                rng.gen::<f64>() * 0.01,
                rng.gen::<f64>() * 0.01,
            );
            matrix[[i, j]] = val;
            matrix[[j, i]] = val.conj();
        }
    }

    let start = Instant::now();
    let condition = monitor.estimate_condition(&matrix);
    let elapsed = start.elapsed();

    println!(
        "N=5000: condition={:.2e}, time={:.3}ms",
        condition,
        elapsed.as_secs_f64() * 1000.0
    );

    // Should complete in < 1 second
    assert!(
        elapsed < Duration::from_secs(1),
        "Large matrix estimation too slow: {:?}",
        elapsed
    );

    // Check that result is reasonable
    assert!(condition > 1.0 && condition < 1e10);
}

fn main() {
    // Run all tests with output
    test_overhead_well_conditioned();
    test_overhead_ill_conditioned();
    test_batch_performance();
    test_fast_path_effectiveness();
    test_accuracy_vs_speed_tradeoff();
    test_large_matrix_performance();

    println!("\n✅ All performance tests passed!");
}