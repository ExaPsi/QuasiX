//! Integration test for polarizability module S3-2

use approx::assert_relative_eq;
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::dielectric::{
    DielectricFunction, PolarizabilityConfig, PolarizabilityRI, ScreenedCoulomb,
};
use quasix_core::freq::{FrequencyGrid, GridType};

/// Test full workflow: DF tensor -> P0(ω) -> ε -> W
#[test]
fn test_polarizability_workflow() {
    // System parameters (H2O-like)
    let nocc = 5;
    let nvirt = 10;
    let naux = 30;
    let n_trans = nocc * nvirt;

    // 1. Generate mock DF tensor (ia|P)
    let mut df_ia = Array2::<f64>::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            // Create some reasonable values with structure
            df_ia[[i, p]] = (-(i as f64 - p as f64).powi(2) / 10.0).exp();
        }
    }

    // 2. Generate mock orbital energies
    let e_occ = Array1::linspace(-1.0, -0.5, nocc);
    let e_virt = Array1::linspace(0.1, 1.0, nvirt);

    // 3. Create polarizability calculator
    let config = PolarizabilityConfig {
        eta: 1e-4,
        threshold: 1e-12,
        enforce_hermiticity: true,
        batch_size: 4,
        max_threads: Some(2),
    };
    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux).with_config(config);

    // 4. Generate frequency grid
    let grid = FrequencyGrid::new(10, GridType::ModifiedGaussLegendre { omega_max: 10.0 }).unwrap();

    // 5. Compute P0 at multiple frequencies
    let p0_batch = pol_calc
        .compute_p0_with_grid(&grid, &df_ia, &e_occ, &e_virt)
        .unwrap();

    assert_eq!(p0_batch.len(), 10);

    // Check properties of P0 at each frequency
    for p0 in &p0_batch {
        assert_eq!(p0.dim(), (naux, naux));

        // Check Hermiticity
        for i in 0..naux {
            for j in 0..naux {
                assert_relative_eq!(p0[[i, j]].re, p0[[j, i]].re, epsilon = 1e-12);
                assert_relative_eq!(p0[[i, j]].im, -p0[[j, i]].im, epsilon = 1e-12);
            }
        }
    }

    // 6. Test symmetrization with mock V^(1/2)
    let vsqrt = Array2::<f64>::eye(naux) * 0.5;
    let omega = Complex64::new(0.0, 1.0);
    let p0_single = pol_calc.compute_p0(omega, &df_ia, &e_occ, &e_virt).unwrap();
    let m_matrix = pol_calc.symmetrize_p0(&p0_single, &vsqrt).unwrap();

    // 7. Compute dielectric matrix
    let diel_func = DielectricFunction::new(naux);
    let epsilon = diel_func.compute_epsilon(&m_matrix).unwrap();

    // Check that ε = I - M
    for i in 0..naux {
        for j in 0..naux {
            let expected = if i == j {
                Complex64::new(1.0, 0.0) - m_matrix[[i, j]]
            } else {
                -m_matrix[[i, j]]
            };
            assert_relative_eq!(epsilon[[i, j]].re, expected.re, epsilon = 1e-12);
            assert_relative_eq!(epsilon[[i, j]].im, expected.im, epsilon = 1e-12);
        }
    }

    // 8. Compute screened interaction W
    let epsilon_inv = diel_func.compute_epsilon_inv(&epsilon).unwrap();
    let screened = ScreenedCoulomb::new(naux);
    let w = screened.compute_w(&epsilon_inv, &vsqrt).unwrap();

    // W should be Hermitian
    for i in 0..naux {
        for j in 0..naux {
            assert_relative_eq!(
                w[[i, j]].re,
                w[[j, i]].re,
                epsilon = 1e-10 // Slightly relaxed due to numerical errors
            );
            assert_relative_eq!(w[[i, j]].im, -w[[j, i]].im, epsilon = 1e-10);
        }
    }

    // 9. Compute W - V for correlation
    let v = Array2::<f64>::eye(naux); // Mock bare Coulomb
    let w_minus_v = screened.compute_w_minus_v(&w, &v).unwrap();

    // Check dimensions
    assert_eq!(w_minus_v.dim(), (naux, naux));
}

/// Test static polarizability calculation
#[test]
fn test_static_polarizability() {
    use quasix_core::dielectric::StaticPolarizability;

    let nocc = 5;
    let nvirt = 10;
    let naux = 20;
    let n_trans = nocc * nvirt;

    // Generate mock data
    let mut df_ia = Array2::<f64>::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            df_ia[[i, p]] = (-(i as f64 - p as f64).powi(2) / 20.0).exp();
        }
    }

    let e_occ = Array1::linspace(-1.0, -0.5, nocc);
    let e_virt = Array1::linspace(0.1, 1.0, nvirt);

    // Compute static polarizability
    let static_calc = StaticPolarizability::new(nocc, nvirt, naux);
    let p0_static = static_calc.compute(&df_ia, &e_occ, &e_virt).unwrap();

    // Static polarizability should be real and symmetric
    assert_eq!(p0_static.dim(), (naux, naux));

    for i in 0..naux {
        for j in 0..naux {
            assert_relative_eq!(p0_static[[i, j]], p0_static[[j, i]], epsilon = 1e-12);
            // Check values are reasonable (not infinite or NaN)
            assert!(p0_static[[i, j]].is_finite());
        }
    }

    // Static polarizability should be positive semi-definite
    // (all eigenvalues should be non-negative for physical systems)
    // This is a simplified check - just ensure diagonal elements are positive
    for i in 0..naux {
        assert!(
            p0_static[[i, i]] >= 0.0,
            "Diagonal element {} is negative: {}",
            i,
            p0_static[[i, i]]
        );
    }
}

/// Test parallel batch processing
#[test]
fn test_batch_processing() {
    let nocc = 4;
    let nvirt = 8;
    let naux = 15;
    let n_trans = nocc * nvirt;

    // Setup test data
    let df_ia = Array2::<f64>::ones((n_trans, naux)) * 0.1;
    let e_occ = Array1::linspace(-1.0, -0.6, nocc);
    let e_virt = Array1::linspace(0.2, 0.8, nvirt);

    // Create calculator with specific batch size
    let config = PolarizabilityConfig {
        batch_size: 3,
        max_threads: Some(2),
        ..Default::default()
    };

    let pol_calc = PolarizabilityRI::new(nocc, nvirt, naux).with_config(config);

    // Test with multiple frequencies
    let frequencies = vec![
        Complex64::new(0.0, 0.5),
        Complex64::new(0.0, 1.0),
        Complex64::new(0.0, 1.5),
        Complex64::new(0.0, 2.0),
        Complex64::new(0.0, 2.5),
        Complex64::new(0.0, 3.0),
    ];

    let p0_batch = pol_calc
        .compute_p0_batch(&frequencies, &df_ia, &e_occ, &e_virt)
        .unwrap();

    assert_eq!(p0_batch.len(), 6);

    // Each P0 should be different due to different frequencies
    for i in 0..5 {
        let diff_norm: f64 = (&p0_batch[i] - &p0_batch[i + 1])
            .iter()
            .map(|x| x.norm_sqr())
            .sum::<f64>()
            .sqrt();
        assert!(
            diff_norm > 1e-6,
            "P0 matrices at different frequencies should differ"
        );
    }
}

/// Test error handling for invalid inputs
#[test]
fn test_invalid_inputs() {
    let pol_calc = PolarizabilityRI::new(5, 10, 30);

    // Wrong DF tensor dimensions
    let df_wrong = Array2::<f64>::zeros((40, 30)); // Should be 50x30
    let e_occ = Array1::zeros(5);
    let e_virt = Array1::zeros(10);
    let omega = Complex64::new(0.0, 1.0);

    let result = pol_calc.compute_p0(omega, &df_wrong, &e_occ, &e_virt);
    assert!(result.is_err());

    // Wrong orbital energy dimensions
    let df_correct = Array2::<f64>::zeros((50, 30));
    let e_occ_wrong = Array1::zeros(4); // Should be 5
    let result = pol_calc.compute_p0(omega, &df_correct, &e_occ_wrong, &e_virt);
    assert!(result.is_err());

    let e_virt_wrong = Array1::zeros(9); // Should be 10
    let result = pol_calc.compute_p0(omega, &df_correct, &e_occ, &e_virt_wrong);
    assert!(result.is_err());
}
