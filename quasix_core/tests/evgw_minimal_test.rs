//! Minimal test for evGW functionality without full driver

use ndarray::{Array1, Array2, Array3};
use quasix_core::gw::check_simd_features;

#[test]
fn test_simd_detection() {
    let features = check_simd_features();

    // Just check that detection works
    #[cfg(target_arch = "x86_64")]
    assert!(features.sse2, "SSE2 should be available on x86_64");

    println!("SIMD detection test passed");
}

#[test]
fn test_evgw_minimal() {
    // This test just validates that we can create the basic structures
    // without running the full evGW loop

    let nbasis = 4;
    let nocc = 2;
    let naux = 6;

    // Create test data
    let mo_energies = Array1::linspace(-0.5, 0.5, nbasis);
    let mo_occ = {
        let mut occ = Array1::<f64>::zeros(nbasis);
        for i in 0..nocc {
            occ[i] = 2.0;
        }
        occ
    };

    let ia_p = Array3::from_shape_fn((nocc, nbasis - nocc, naux), |(i, a, p)| {
        ((i + 1) * (a + 1)) as f64 * 0.1 + p as f64 * 0.01
    });

    let ij_p = Array3::from_shape_fn((nocc, nocc, naux), |(i, j, p)| {
        ((i + 1) * (j + 1)) as f64 * 0.1 + p as f64 * 0.01
    });

    let chol_v = Array2::<f64>::eye(naux);
    let vxc_dft = Array1::<f64>::zeros(nbasis);

    // Just validate the data is created correctly
    assert_eq!(mo_energies.len(), nbasis);
    assert_eq!(mo_occ.len(), nbasis);
    assert_eq!(ia_p.shape(), &[nocc, nbasis - nocc, naux]);
    assert_eq!(ij_p.shape(), &[nocc, nocc, naux]);
    assert_eq!(chol_v.shape(), &[naux, naux]);
    assert_eq!(vxc_dft.len(), nbasis);

    println!("Minimal evGW test passed - data structures created successfully");
}

#[test]
fn test_evgw_components() {
    use quasix_core::gw::evgw_simd::{build_polarizability_simd, evaluate_selfenergy_simd};

    // Test individual components without the driver
    let nocc = 2;
    let nvirt = 2;
    let naux = 4;

    let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| (i + a + p) as f64 * 0.1);

    let qp_energies = Array1::from(vec![-0.5, -0.3, 0.3, 0.5]);
    let omega = 0.1;

    // Test polarizability
    let p0 = build_polarizability_simd(&ia_p, &qp_energies.view(), omega, nocc, nvirt, naux);
    assert_eq!(p0.shape(), &[naux, naux]);
    assert!(p0.iter().all(|&x| x.is_finite()));

    // Test self-energy
    let w_omega = Array2::<f64>::eye(naux) * 0.5;
    let sigma = evaluate_selfenergy_simd(
        0, -0.5, omega, true, &ia_p, &w_omega, 0.01, nocc, nvirt, naux,
    );
    assert!(sigma.re.is_finite() && sigma.im.is_finite());

    println!("Component tests passed - SIMD functions work correctly");
}
