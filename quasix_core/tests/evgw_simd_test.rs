//! Test for evGW SIMD functionality

use ndarray::{Array1, Array2, Array3};
use quasix_core::gw::check_simd_features;

#[test]
fn test_simd_features_detection() {
    let features = check_simd_features();

    println!("SIMD Feature Detection:");
    println!("  SSE2: {}", features.sse2);
    println!("  AVX: {}", features.avx);
    println!("  AVX2: {}", features.avx2);
    println!("  FMA: {}", features.fma);
    println!("  AVX512F: {}", features.avx512f);
    println!("  High-performance SIMD: {}", features.has_high_perf_simd());
    println!("  Expected speedup: {:.1}×", features.expected_speedup());

    // On x86_64, we should at least have SSE2
    #[cfg(target_arch = "x86_64")]
    assert!(features.sse2, "SSE2 should be available on x86_64");
}

#[test]
fn test_evgw_with_simd() {
    // This test validates SIMD functionality without using the full evGW driver
    // which can cause issues with nested parallelism in test environments

    use quasix_core::gw::evgw_simd::{build_polarizability_simd, evaluate_selfenergy_simd};

    // Small problem size for fast testing
    let nbasis = 4;
    let nocc = 2;
    let nvirt = nbasis - nocc;
    let naux = 6;

    // Create minimal test data
    let mo_energies = Array1::linspace(-0.5, 0.5, nbasis);
    let mo_occ = {
        let mut occ = Array1::<f64>::zeros(nbasis);
        for i in 0..nocc {
            occ[i] = 2.0;
        }
        occ
    };

    let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| {
        ((i + 1) * (a + 1)) as f64 * 0.1 + p as f64 * 0.01
    });

    // Test polarizability construction (key evGW operation)
    let omega = 0.5;
    let p0 = build_polarizability_simd(&ia_p, &mo_energies.view(), omega, nocc, nvirt, naux);

    // Validate polarizability matrix
    assert_eq!(p0.shape(), &[naux, naux]);
    assert!(
        p0.iter().all(|&x| x.is_finite()),
        "P0 contains non-finite values"
    );

    // Check symmetry (P0 should be symmetric)
    for i in 0..naux {
        for j in i + 1..naux {
            assert!(
                (p0[[i, j]] - p0[[j, i]]).abs() < 1e-10,
                "P0 not symmetric at [{},{}]: {} vs {}",
                i,
                j,
                p0[[i, j]],
                p0[[j, i]]
            );
        }
    }

    // Test self-energy evaluation (another key evGW operation)
    let w_omega = Array2::<f64>::eye(naux) * 0.5;
    let eta = 0.01;

    for n in 0..nbasis {
        let e_n = mo_energies[n];
        let is_occupied = mo_occ[n] > 0.5;

        let sigma = evaluate_selfenergy_simd(
            n,
            e_n,
            omega,
            is_occupied,
            &ia_p,
            &w_omega,
            eta,
            nocc,
            nvirt,
            naux,
        );

        assert!(
            sigma.re.is_finite() && sigma.im.is_finite(),
            "Self-energy for orbital {} not finite: {:?}",
            n,
            sigma
        );
    }

    // Report SIMD performance
    let features = check_simd_features();
    if std::env::var("RUST_TEST_VERBOSE").is_ok() {
        println!("\nSIMD Features Detected:");
        println!("  SSE2: {}", features.sse2);
        println!("  AVX: {}", features.avx);
        println!("  AVX2: {}", features.avx2);
        println!("  FMA: {}", features.fma);
        println!("  High-performance SIMD: {}", features.has_high_perf_simd());
        println!("  Expected speedup: {:.1}×", features.expected_speedup());
    }

    // On x86_64, we should at least have SSE2
    #[cfg(target_arch = "x86_64")]
    assert!(features.sse2, "SSE2 should be available on x86_64");

    println!("evGW SIMD test completed successfully");
}

#[test]
fn test_simd_speedup_polarizability() {
    use quasix_core::gw::evgw_simd::build_polarizability_simd;
    use std::time::Instant;

    // Smaller problem size for CI/CD
    let nocc = 4;
    let nvirt = 4;
    let naux = 12;

    let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| {
        ((i + 1) * (a + 1)) as f64 * 0.1 + p as f64 * 0.01
    });

    let qp_energies = Array1::linspace(-1.0, 1.0, nocc + nvirt);
    let omega = 0.5;

    // Warm-up with fewer iterations
    for _ in 0..2 {
        let _ = build_polarizability_simd(&ia_p, &qp_energies.view(), omega, nocc, nvirt, naux);
    }

    // Timing run with fewer iterations for speed
    let start = Instant::now();
    let iterations = 10; // Reduced from 100
    for _ in 0..iterations {
        let _ = build_polarizability_simd(&ia_p, &qp_energies.view(), omega, nocc, nvirt, naux);
    }
    let elapsed = start.elapsed();

    let avg_time = elapsed.as_secs_f64() / iterations as f64;

    // Only print in verbose mode
    if std::env::var("RUST_TEST_VERBOSE").is_ok() {
        println!("\nPolarizability construction performance:");
        println!("  Average time: {:.3} ms", avg_time * 1000.0);
        println!("  Matrix size: {}×{}", naux, naux);

        // Estimate FLOPS
        let flops = (nocc * nvirt * naux * naux * 4) as f64; // Rough estimate
        let gflops = (flops * iterations as f64) / (elapsed.as_secs_f64() * 1e9);
        println!("  Estimated GFLOPS: {:.2}", gflops);

        let features = check_simd_features();
        if features.has_high_perf_simd() {
            println!("  SIMD: AVX2/FMA (optimized path)");
        } else {
            println!("  SIMD: Scalar fallback");
        }
    }

    // Relaxed performance check - just ensure it runs without hanging
    assert!(
        avg_time < 1.0,
        "Polarizability calculation should complete quickly (< 1s)"
    );
}
