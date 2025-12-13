//! Simple test for evGW SIMD functionality without full driver

use ndarray::{Array1, Array2, Array3};
use quasix_core::gw::{
    check_simd_features,
    evgw_simd::{build_polarizability_simd, evaluate_selfenergy_simd},
};

#[test]
fn test_simd_features() {
    let features = check_simd_features();

    println!("SIMD Features:");
    println!("  SSE2: {}", features.sse2);
    println!("  AVX: {}", features.avx);
    println!("  AVX2: {}", features.avx2);
    println!("  FMA: {}", features.fma);

    // On x86_64, SSE2 should always be available
    #[cfg(target_arch = "x86_64")]
    assert!(features.sse2, "SSE2 should be available on x86_64");
}

#[test]
fn test_polarizability_simd() {
    // Small test problem
    let nocc = 2;
    let nvirt = 2;
    let naux = 4;

    let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| (i + a + p) as f64 * 0.1);

    let qp_energies = Array1::from(vec![-0.5, -0.3, 0.3, 0.5]);
    let omega = 0.1;

    // Test polarizability construction
    let p0 = build_polarizability_simd(&ia_p, &qp_energies.view(), omega, nocc, nvirt, naux);

    // Check dimensions
    assert_eq!(p0.shape(), &[naux, naux]);

    // Check values are finite
    assert!(p0.iter().all(|&x| x.is_finite()));

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
}

#[test]
fn test_selfenergy_simd() {
    // Small test problem
    let nocc = 2;
    let nvirt = 2;
    let naux = 4;

    let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| (i + a + p) as f64 * 0.05);

    let w_omega = Array2::eye(naux) * 0.5;

    // Test for occupied orbital
    let n = 0;
    let e_n = -0.5;
    let omega = 0.1;
    let is_occupied = true;
    let eta = 0.01;

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

    // Check result is finite
    assert!(
        sigma.re.is_finite() && sigma.im.is_finite(),
        "Self-energy not finite: {:?}",
        sigma
    );

    // For small eta, imaginary part should be small
    assert!(
        sigma.im.abs() < 1.0,
        "Imaginary part too large: {}",
        sigma.im
    );
}

#[test]
fn test_simd_performance_baseline() {
    use std::time::Instant;

    // Moderate size for baseline performance
    let nocc = 5;
    let nvirt = 5;
    let naux = 10;

    let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| {
        ((i + 1) * (a + 1)) as f64 * 0.1 + p as f64 * 0.01
    });

    let qp_energies = Array1::linspace(-1.0, 1.0, nocc + nvirt);
    let omega = 0.5;

    // Time polarizability construction
    let start = Instant::now();
    for _ in 0..5 {
        let _ = build_polarizability_simd(&ia_p, &qp_energies.view(), omega, nocc, nvirt, naux);
    }
    let elapsed = start.elapsed();

    // Should complete quickly
    assert!(
        elapsed.as_secs() < 2,
        "Polarizability took too long: {:?}",
        elapsed
    );

    if std::env::var("RUST_TEST_VERBOSE").is_ok() {
        println!(
            "Polarizability ({}x{}) completed in {:?}",
            naux, naux, elapsed
        );

        let features = check_simd_features();
        if features.has_high_perf_simd() {
            println!("Using AVX2/FMA optimization");
        } else {
            println!("Using scalar fallback");
        }
    }
}
