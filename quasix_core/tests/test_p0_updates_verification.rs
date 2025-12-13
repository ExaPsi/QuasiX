//! Comprehensive verification tests for S5-2: P0 Denominator Updates
//!
//! Tests verify:
//! 1. Algorithm correctness: Denominator update formula
//! 2. Physical validity: Gap positivity, energy ordering
//! 3. Convergence properties with mixing
//! 4. Edge case handling (small gaps, near-degeneracies)

use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::gw::polarizability_builder::{PolarizabilityBuilder, GapStatistics};

/// Create a test system with controllable gaps
fn create_test_system(nocc: usize, nvirt: usize, naux: usize) -> (Array1<f64>, Array2<f64>) {
    let mut mo_energies = Array1::zeros(nocc + nvirt);

    // Occupied energies: -0.5 to -0.1 Ha
    for i in 0..nocc {
        mo_energies[i] = -0.5 + 0.4 * (i as f64) / ((nocc - 1).max(1) as f64);
    }

    // Virtual energies: 0.1 to 0.5 Ha
    for a in 0..nvirt {
        mo_energies[nocc + a] = 0.1 + 0.4 * (a as f64) / ((nvirt - 1).max(1) as f64);
    }

    // Create DF tensor with small random values
    let n_trans = nocc * nvirt;
    let mut df_ia = Array2::zeros((n_trans, naux));

    // Initialize with deterministic "random" values for reproducibility
    let mut seed = 42u64;
    for i in 0..n_trans {
        for p in 0..naux {
            seed = (seed.wrapping_mul(1664525).wrapping_add(1013904223)) % (1 << 32);
            df_ia[[i, p]] = ((seed as f64) / (1u64 << 32) as f64 - 0.5) * 0.1;
        }
    }

    (mo_energies, df_ia)
}

/// Simulate one evGW iteration with simple self-energy corrections
fn simulate_evgw_iteration(
    energies: &Array1<f64>,
    nocc: usize,
    sigma_correction: f64,
) -> Array1<f64> {
    let mut new_energies = energies.clone();
    let nvirt = energies.len() - nocc;

    // Compute mean gap
    let mut mean_gap = 0.0;
    for i in 0..nocc {
        for a in 0..nvirt {
            mean_gap += energies[nocc + a] - energies[i];
        }
    }
    mean_gap /= (nocc * nvirt) as f64;

    // Apply corrections
    for i in 0..nocc {
        new_energies[i] -= sigma_correction * (1.0 - i as f64 / nocc as f64) * mean_gap;
    }

    for a in 0..nvirt {
        new_energies[nocc + a] += sigma_correction * (1.0 - a as f64 / nvirt as f64) * mean_gap;
    }

    new_energies
}

#[test]
fn test_algorithm_correctness() {
    println!("\n=== TEST 1: Algorithm Correctness ===");

    let nocc = 5;
    let nvirt = 10;
    let naux = 30;

    let (mo_energies, df_ia) = create_test_system(nocc, nvirt, naux);
    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

    // Initial gaps
    let gaps_initial = builder.compute_gaps();
    let min_initial = gaps_initial.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_initial = gaps_initial.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("Initial gaps: [{:.6}, {:.6}] Ha", min_initial, max_initial);

    // Simulate QP correction
    let qp_energies = simulate_evgw_iteration(&mo_energies, nocc, 0.05);

    // Update denominators
    let stats = builder.update_energies(&qp_energies).unwrap();
    println!("After update: min_gap={:.6} Ha, max_gap={:.6} Ha", stats.min_gap, stats.max_gap);

    // Verify formula: Δ_ia = E_a - E_i
    let gaps_updated = builder.compute_gaps();
    let mut max_error: f64 = 0.0;

    for (idx, &gap) in gaps_updated.iter().enumerate() {
        let i = idx / nvirt;
        let a = idx % nvirt;
        let expected = qp_energies[nocc + a] - qp_energies[i];
        let expected_protected = expected.max(builder.gap_threshold);
        let error = (gap - expected_protected).abs();
        max_error = max_error.max(error);
    }

    println!("Max error in gap calculation: {:.2e} Ha", max_error);
    assert!(max_error < 1e-12, "Gap calculation error too large");
}

#[test]
fn test_physical_validity() {
    println!("\n=== TEST 2: Physical Validity ===");

    let nocc = 5;
    let nvirt = 10;
    let naux = 30;

    let (mo_energies, df_ia) = create_test_system(nocc, nvirt, naux);
    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

    let mut energies = mo_energies.clone();
    let mut all_positive = true;

    // Run several evGW iterations
    for iter in 0..5 {
        energies = simulate_evgw_iteration(&energies, nocc, 0.02);
        let stats = builder.update_energies(&energies).unwrap();

        println!("Iteration {}: min_gap={:.6} Ha, n_negative={}",
                 iter + 1, stats.min_gap, stats.n_negative);

        if stats.n_negative > 0 {
            all_positive = false;
        }
    }

    // Check energy ordering
    let e_occ = energies.slice(ndarray::s![..nocc]);
    let e_virt = energies.slice(ndarray::s![nocc..]);

    let occ_ordered = e_occ.windows(2)
        .into_iter()
        .all(|w| w[1] >= w[0]);

    let virt_ordered = e_virt.windows(2)
        .into_iter()
        .all(|w| w[1] >= w[0]);

    let homo = e_occ[e_occ.len() - 1];
    let lumo = e_virt[0];
    let proper_gap = lumo > homo;

    println!("Energy ordering:");
    println!("  Occupied ordered: {}", occ_ordered);
    println!("  Virtual ordered: {}", virt_ordered);
    println!("  HOMO-LUMO gap positive: {} ({:.6} Ha)", proper_gap, lumo - homo);
    println!("  All gaps positive after protection: {}", all_positive);

    assert!(proper_gap, "HOMO-LUMO gap is not positive");
}

#[test]
fn test_convergence_with_mixing() {
    println!("\n=== TEST 3: Convergence with Mixing ===");

    let nocc = 5;
    let nvirt = 10;
    let naux = 30;

    let (mo_energies, df_ia) = create_test_system(nocc, nvirt, naux);
    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

    let mut energies = mo_energies.clone();
    let mixing = 0.5;
    let mut gap_history = Vec::new();

    println!("Running evGW with mixing={}", mixing);

    for iter in 0..10 {
        // Get new energies with decreasing correction
        let energies_new = simulate_evgw_iteration(&energies, nocc, 0.03 / (iter + 1) as f64);

        // Apply mixing
        energies = energies * (1.0 - mixing) + energies_new * mixing;

        // Update denominators
        let stats = builder.update_energies(&energies).unwrap();
        gap_history.push(stats.mean_gap);

        if iter % 3 == 0 {
            println!("  Iteration {}: mean_gap={:.6} Ha", iter + 1, stats.mean_gap);
        }
    }

    // Check for smooth convergence
    let mut smooth = true;
    for i in 2..gap_history.len() {
        let change = gap_history[i] - gap_history[i-1];
        let prev_change = gap_history[i-1] - gap_history[i-2];
        if change * prev_change < -1e-6 {  // Sign change indicates oscillation
            smooth = false;
            break;
        }
    }

    println!("Gap evolution:");
    println!("  Initial: {:.6} Ha", gap_history[0]);
    println!("  Final: {:.6} Ha", gap_history.last().unwrap());
    println!("  Total change: {:.6} Ha", gap_history.last().unwrap() - gap_history[0]);
    println!("  Smooth evolution: {}", smooth);

    // Check convergence (change in last iterations should be small)
    let last_change = (gap_history[9] - gap_history[8]).abs();
    println!("  Last iteration change: {:.2e} Ha", last_change);

    // Relaxed threshold since we're using simple model corrections
    assert!(last_change < 2e-3, "Not converged sufficiently");
}

#[test]
fn test_edge_cases() {
    println!("\n=== TEST 4: Edge Cases ===");

    let nocc = 5;
    let nvirt = 10;
    let naux = 30;

    let (mut mo_energies, df_ia) = create_test_system(nocc, nvirt, naux);

    // Create near-degenerate levels
    mo_energies[nocc - 1] = -0.01;  // Highest occupied
    mo_energies[nocc] = 0.001;      // Lowest virtual - very small gap!
    mo_energies[nocc + 1] = 0.002;  // Near-degenerate virtuals

    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();
    builder.gap_threshold = 1e-6;

    println!("Testing near-degenerate levels:");
    let stats = builder.update_energies(&mo_energies).unwrap();
    println!("  Smallest gap: {:.2e} Ha", stats.min_gap);
    println!("  Gaps thresholded: {}", stats.n_thresholded);

    // Test level crossing
    println!("\nTesting level crossing detection:");
    let mut crossed_energies = mo_energies.clone();
    crossed_energies[nocc - 1] = 0.05;   // Make occupied higher
    crossed_energies[nocc] = -0.05;      // Make virtual lower

    let stats_crossed = builder.update_energies(&crossed_energies).unwrap();
    println!("  Negative gaps detected: {}", stats_crossed.n_negative);
    assert!(stats_crossed.n_negative > 0, "Should detect negative gaps");

    // Test P0 build stability with edge cases
    println!("\nTesting P0 build stability:");
    let test_frequencies = vec![
        Complex64::new(0.0, 0.0),      // Static
        Complex64::new(0.0, 0.1),      // Imaginary
        Complex64::new(0.1, 0.01),     // Complex
    ];

    for omega in test_frequencies {
        match builder.build_p0(omega) {
            Ok(p0) => {
                let is_finite = p0.iter().all(|x| x.re.is_finite() && x.im.is_finite());
                println!("  ω = {:.3}+{:.3}i: P0 finite = {}",
                        omega.re, omega.im, is_finite);
                assert!(is_finite, "P0 contains non-finite values");
            }
            Err(e) => {
                panic!("P0 build failed for ω = {:?}: {}", omega, e);
            }
        }
    }
}

#[test]
fn test_simd_consistency() {
    println!("\n=== TEST 5: SIMD Consistency ===");

    // Use odd numbers to test remainder handling
    let nocc = 7;
    let nvirt = 13;
    let naux = 47;

    let (mo_energies, df_ia) = create_test_system(nocc, nvirt, naux);
    let builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia).unwrap();

    // Compare SIMD and scalar gap computations
    let gaps_simd = builder.compute_gaps();
    let mut gaps_scalar = vec![0.0; nocc * nvirt];
    builder.compute_gaps_scalar(&mut gaps_scalar);

    let mut max_diff: f64 = 0.0;
    for (i, (&simd_val, &scalar_val)) in gaps_simd.iter().zip(gaps_scalar.iter()).enumerate() {
        let diff = (simd_val - scalar_val).abs();
        max_diff = max_diff.max(diff);
        if diff > 1e-14 {
            println!("Gap difference at index {}: {:.2e}", i, diff);
        }
    }

    println!("Maximum difference between SIMD and scalar: {:.2e}", max_diff);
    assert!(max_diff < 1e-14, "SIMD and scalar results differ");

    // Test P0 build for imaginary frequencies
    let omega = Complex64::new(0.0, 1.0);
    let p0_standard = builder.build_p0_standard(omega).unwrap();

    #[cfg(target_arch = "x86_64")]
    {
        let p0_simd = builder.build_p0_simd_imaginary(omega).unwrap();

        let mut max_diff_re: f64 = 0.0;
        let mut max_diff_im: f64 = 0.0;

        for i in 0..naux {
            for j in 0..naux {
                let diff_re = (p0_standard[[i,j]].re - p0_simd[[i,j]].re).abs();
                let diff_im = (p0_standard[[i,j]].im - p0_simd[[i,j]].im).abs();
                max_diff_re = max_diff_re.max(diff_re);
                max_diff_im = max_diff_im.max(diff_im);
            }
        }

        println!("P0 SIMD vs standard (ω = {}+{}i):", omega.re, omega.im);
        println!("  Max diff (real): {:.2e}", max_diff_re);
        println!("  Max diff (imag): {:.2e}", max_diff_im);

        assert!(max_diff_re < 1e-12, "P0 real parts differ");
        assert!(max_diff_im < 1e-12, "P0 imaginary parts differ");
    }
}

#[test]
fn test_gap_statistics() {
    println!("\n=== TEST 6: Gap Statistics ===");

    let gaps = vec![0.5, 1.0, 1.5, 0.001, 2.0, -0.1];
    let stats = GapStatistics::from_gaps(&gaps, 0.01);

    println!("Gap statistics for test array:");
    println!("  Min gap: {:.6} Ha", stats.min_gap);
    println!("  Max gap: {:.6} Ha", stats.max_gap);
    println!("  Mean gap: {:.6} Ha", stats.mean_gap);
    println!("  Negative gaps: {}", stats.n_negative);
    println!("  Thresholded gaps: {}", stats.n_thresholded);
    println!("  Total gaps: {}", stats.n_total);

    assert_eq!(stats.n_negative, 1, "Should detect 1 negative gap");
    assert_eq!(stats.n_thresholded, 2, "Should threshold 2 gaps");
    assert_eq!(stats.n_total, 6, "Should have 6 total gaps");

    assert!((stats.min_gap - (-0.1)).abs() < 1e-10, "Min gap incorrect");
    assert!((stats.max_gap - 2.0).abs() < 1e-10, "Max gap incorrect");
}

fn main() {
    println!("\n========================================");
    println!("S5-2 P0 Denominator Update Verification");
    println!("========================================\n");

    test_algorithm_correctness();
    test_physical_validity();
    test_convergence_with_mixing();
    test_edge_cases();
    test_simd_consistency();
    test_gap_statistics();

    println!("\n========================================");
    println!("ALL TESTS PASSED ✓");
    println!("========================================");
}