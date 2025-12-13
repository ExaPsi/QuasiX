//! Verification tests for SIMD polarizability implementation
//!
//! This test suite ensures that the SIMD-optimized implementation
//! produces results identical to the standard implementation within
//! numerical tolerance.

use approx::assert_relative_eq;
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::dielectric::{PolarizabilityRI, PolarizabilitySIMDCore};

/// Generate realistic test system
fn generate_test_system(
    nocc: usize,
    nvirt: usize,
    naux: usize,
) -> (Array2<f64>, Array1<f64>, Array1<f64>) {
    let n_trans = nocc * nvirt;

    // Generate realistic DF tensor with Gaussian-like structure
    let mut df_ia = Array2::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            let r_sq = ((i as f64 - n_trans as f64 / 2.0).powi(2)
                + (p as f64 - naux as f64 / 2.0).powi(2))
                / (naux as f64).powi(2);
            df_ia[[i, p]] = (-r_sq).exp();
        }
    }

    // Realistic orbital energies (Hartree units)
    let e_occ = Array1::linspace(-1.5, -0.3, nocc);
    let e_virt = Array1::linspace(0.05, 1.5, nvirt);

    (df_ia, e_occ, e_virt)
}

#[test]
fn test_simd_vs_standard_accuracy() {
    // Test different system sizes
    let test_cases = [
        ("small", 5, 10, 30),
        ("medium", 10, 50, 100),
        ("large", 15, 75, 150),
    ];

    for (name, nocc, nvirt, naux) in &test_cases {
        println!(
            "Testing {} system: nocc={}, nvirt={}, naux={}",
            name, nocc, nvirt, naux
        );

        let (df_ia, e_occ, e_virt) = generate_test_system(*nocc, *nvirt, *naux);

        // Create calculators
        let standard_calc = PolarizabilityRI::new(*nocc, *nvirt, *naux);
        let simd_calc = PolarizabilitySIMDCore::new(*nocc, *nvirt, *naux);

        // Test at different frequency types
        let test_frequencies = [
            ("static", Complex64::new(0.0, 1e-4)),
            ("imaginary", Complex64::new(0.0, 1.0)),
            ("real", Complex64::new(0.5, 1e-4)),
            ("complex", Complex64::new(0.5, 0.5)),
        ];

        for (freq_type, omega) in &test_frequencies {
            println!(
                "  Testing {} frequency: ω = {:.3}+{:.3}i",
                freq_type, omega.re, omega.im
            );

            // Compute with both implementations
            let p0_standard = standard_calc
                .compute_p0(*omega, &df_ia, &e_occ, &e_virt)
                .unwrap();
            let p0_simd = simd_calc
                .compute_p0_simd(*omega, &df_ia, &e_occ, &e_virt)
                .unwrap();

            // Check dimensions
            assert_eq!(p0_standard.dim(), p0_simd.dim());

            // Check accuracy - should be identical within machine precision
            let mut max_error: f64 = 0.0;
            for i in 0..*naux {
                for j in 0..*naux {
                    let diff = (p0_standard[[i, j]] - p0_simd[[i, j]]).norm();
                    max_error = max_error.max(diff);

                    // Verify within tight tolerance
                    assert_relative_eq!(
                        p0_standard[[i, j]].re,
                        p0_simd[[i, j]].re,
                        epsilon = 1e-10,
                        max_relative = 1e-10
                    );
                    assert_relative_eq!(
                        p0_standard[[i, j]].im,
                        p0_simd[[i, j]].im,
                        epsilon = 1e-10,
                        max_relative = 1e-10
                    );
                }
            }

            println!("    Max error: {:.2e}", max_error);
            assert!(
                max_error < 1e-10,
                "SIMD implementation deviates from standard"
            );

            // Check Hermiticity
            for i in 0..*naux {
                // Diagonal must be real
                assert!(p0_simd[[i, i]].im.abs() < 1e-12);

                for j in i + 1..*naux {
                    // Off-diagonal Hermiticity
                    assert_relative_eq!(p0_simd[[i, j]].re, p0_simd[[j, i]].re, epsilon = 1e-12);
                    assert_relative_eq!(p0_simd[[i, j]].im, -p0_simd[[j, i]].im, epsilon = 1e-12);
                }
            }
        }
    }
}

#[test]
fn test_simd_batch_accuracy() {
    let (nocc, nvirt, naux) = (10, 50, 80);
    let (df_ia, e_occ, e_virt) = generate_test_system(nocc, nvirt, naux);

    // Create frequency grid
    let frequencies: Vec<Complex64> = (0..8)
        .map(|i| Complex64::new(0.0, 0.1 * (i + 1) as f64))
        .collect();

    // Compute with both implementations
    let standard_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let simd_calc = PolarizabilitySIMDCore::new(nocc, nvirt, naux);

    let p0_standard_batch = standard_calc
        .compute_p0_batch(&frequencies, &df_ia, &e_occ, &e_virt)
        .unwrap();
    let p0_simd_batch = simd_calc
        .compute_p0_batch_simd(&frequencies, &df_ia, &e_occ, &e_virt)
        .unwrap();

    // Check each frequency point
    assert_eq!(p0_standard_batch.len(), p0_simd_batch.len());

    for (idx, (p0_std, p0_simd)) in p0_standard_batch
        .iter()
        .zip(p0_simd_batch.iter())
        .enumerate()
    {
        println!("Checking frequency point {}", idx);

        for i in 0..naux {
            for j in 0..naux {
                assert_relative_eq!(
                    p0_std[[i, j]].re,
                    p0_simd[[i, j]].re,
                    epsilon = 1e-10,
                    max_relative = 1e-10
                );
                assert_relative_eq!(
                    p0_std[[i, j]].im,
                    p0_simd[[i, j]].im,
                    epsilon = 1e-10,
                    max_relative = 1e-10
                );
            }
        }
    }
}

#[test]
fn test_simd_cpu_features() {
    let calc = PolarizabilitySIMDCore::new(5, 10, 30);
    let features = calc.cpu_features();

    println!("SIMD CPU features detected: {}", features);

    // Should detect at least portable SIMD
    assert!(!features.is_empty());
    assert!(
        features.contains("Portable SIMD")
            || features.contains("AVX2")
            || features.contains("AVX-512")
    );
}

#[test]
fn test_simd_hermiticity_enforcement() {
    let (nocc, nvirt, naux) = (8, 16, 40);
    let (df_ia, e_occ, e_virt) = generate_test_system(nocc, nvirt, naux);

    let calc = PolarizabilitySIMDCore::new(nocc, nvirt, naux);

    // Test at complex frequency where Hermiticity is non-trivial
    let omega = Complex64::new(0.3, 0.2);
    let p0 = calc
        .compute_p0_simd(omega, &df_ia, &e_occ, &e_virt)
        .unwrap();

    // Verify Hermiticity
    let mut max_hermiticity_error: f64 = 0.0;

    for i in 0..naux {
        // Diagonal elements must be real
        let diag_im_error = p0[[i, i]].im.abs();
        max_hermiticity_error = max_hermiticity_error.max(diag_im_error);
        assert!(
            diag_im_error < 1e-12,
            "Diagonal element [{},{}] not real: {:.2e}",
            i,
            i,
            diag_im_error
        );

        for j in i + 1..naux {
            // Off-diagonal: P[i,j] = conj(P[j,i])
            let re_error = (p0[[i, j]].re - p0[[j, i]].re).abs();
            let im_error = (p0[[i, j]].im + p0[[j, i]].im).abs();

            max_hermiticity_error = max_hermiticity_error.max(re_error);
            max_hermiticity_error = max_hermiticity_error.max(im_error);

            assert!(
                re_error < 1e-12,
                "Hermiticity violated for real part [{},{}]: {:.2e}",
                i,
                j,
                re_error
            );
            assert!(
                im_error < 1e-12,
                "Hermiticity violated for imag part [{},{}]: {:.2e}",
                i,
                j,
                im_error
            );
        }
    }

    println!("Maximum Hermiticity error: {:.2e}", max_hermiticity_error);
}

#[test]
fn test_simd_edge_cases() {
    // Test with very small system
    let (df_ia, e_occ, e_virt) = generate_test_system(1, 1, 2);
    let calc = PolarizabilitySIMDCore::new(1, 1, 2);
    let omega = Complex64::new(0.0, 1.0);
    let p0 = calc.compute_p0_simd(omega, &df_ia, &e_occ, &e_virt);
    assert!(p0.is_ok());

    // Test with non-SIMD-aligned dimensions
    let (df_ia, e_occ, e_virt) = generate_test_system(3, 7, 11);
    let calc = PolarizabilitySIMDCore::new(3, 7, 11);
    let p0 = calc.compute_p0_simd(omega, &df_ia, &e_occ, &e_virt);
    assert!(p0.is_ok());

    // Test with prime number dimensions
    let (df_ia, e_occ, e_virt) = generate_test_system(7, 13, 17);
    let calc = PolarizabilitySIMDCore::new(7, 13, 17);
    let p0 = calc.compute_p0_simd(omega, &df_ia, &e_occ, &e_virt);
    assert!(p0.is_ok());
}
