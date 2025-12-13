//! Integration tests for P0 denominator updates in evGW loop
//!
//! This test validates that the polarizability denominators are correctly
//! updated during evGW iterations using current quasiparticle energies.

use ndarray::{Array1, Array2, Array3};
use quasix_core::gw::{PolarizabilityBuilder, GapStatistics};

/// Create a small test system with well-defined energy structure
fn create_test_system() -> (Array1<f64>, Array3<f64>) {
    let nocc = 3;
    let nvirt = 4;
    let naux = 10;

    // Create initial MO energies with clear HOMO-LUMO gap
    let mut mo_energies = Array1::zeros(nocc + nvirt);

    // Occupied energies: -0.8, -0.6, -0.4 Ha
    mo_energies[0] = -0.8;
    mo_energies[1] = -0.6;
    mo_energies[2] = -0.4;

    // Virtual energies: 0.2, 0.4, 0.6, 0.8 Ha
    mo_energies[3] = 0.2;
    mo_energies[4] = 0.4;
    mo_energies[5] = 0.6;
    mo_energies[6] = 0.8;

    // Create DF tensor (ia|P) with physically reasonable values
    let mut df_ia = Array3::zeros((nocc, nvirt, naux));
    for i in 0..nocc {
        for a in 0..nvirt {
            for p in 0..naux {
                // Simple Gaussian-like overlap
                df_ia[[i, a, p]] = (-0.1 * ((i + a + p) as f64)).exp();
            }
        }
    }

    (mo_energies, df_ia)
}

#[test]
fn test_gap_statistics_computation() {
    let gaps = vec![0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    let stats = GapStatistics::from_gaps(&gaps, 0.01);

    assert_eq!(stats.n_total, 6);
    assert_eq!(stats.n_negative, 0);
    assert_eq!(stats.n_thresholded, 0);
    assert!((stats.min_gap - 0.5).abs() < 1e-10);
    assert!((stats.max_gap - 1.0).abs() < 1e-10);
    assert!((stats.mean_gap - 0.75).abs() < 1e-10);
}

#[test]
fn test_gap_thresholding() {
    // Include some very small and negative gaps
    let gaps = vec![-0.1, 0.0001, 0.005, 0.5, 1.0];
    let threshold = 0.01;
    let stats = GapStatistics::from_gaps(&gaps, threshold);

    assert_eq!(stats.n_total, 5);
    assert_eq!(stats.n_negative, 1); // The -0.1 gap
    assert_eq!(stats.n_thresholded, 3); // -0.1, 0.0001, and 0.005
}

#[test]
fn test_polarizability_builder_creation() {
    let (mo_energies, df_ia) = create_test_system();

    // Convert 3D to 2D for builder
    let nocc = 3;
    let nvirt = 4;
    let naux = 10;
    let n_trans = nocc * nvirt;

    let mut df_ia_2d = Array2::zeros((n_trans, naux));
    let mut idx = 0;
    for i in 0..nocc {
        for a in 0..nvirt {
            for p in 0..naux {
                df_ia_2d[[idx, p]] = df_ia[[i, a, p]];
            }
            idx += 1;
        }
    }

    let builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia_2d);
    assert!(builder.is_ok());

    let builder = builder.unwrap();
    assert_eq!(builder.nocc, 3);
    assert_eq!(builder.nvirt, 4);
    assert_eq!(builder.naux, 10);
}

#[test]
fn test_energy_update_increases_gaps() {
    let (mo_energies, df_ia) = create_test_system();

    let nocc = 3;
    let nvirt = 4;
    let naux = 10;
    let n_trans = nocc * nvirt;

    let mut df_ia_2d = Array2::zeros((n_trans, naux));
    let mut idx = 0;
    for i in 0..nocc {
        for a in 0..nvirt {
            for p in 0..naux {
                df_ia_2d[[idx, p]] = df_ia[[i, a, p]];
            }
            idx += 1;
        }
    }

    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia_2d).unwrap();

    // Simulate evGW update: lower occupied, raise virtual energies
    let mut qp_energies = mo_energies.clone();
    for i in 0..nocc {
        qp_energies[i] -= 0.05; // Lower occupied energies
    }
    for i in nocc..(nocc + nvirt) {
        qp_energies[i] += 0.05; // Raise virtual energies
    }

    let stats = builder.update_energies(&qp_energies).unwrap();

    // All gaps should have increased
    assert!(stats.min_gap > 0.6); // Original min was ~0.6 (0.2 - (-0.4))
    assert_eq!(stats.n_negative, 0);
    assert_eq!(stats.n_total, 12); // 3 occ * 4 virt
}

#[test]
fn test_gap_evolution_during_convergence() {
    let (mo_energies, df_ia) = create_test_system();

    let nocc = 3;
    let nvirt = 4;
    let naux = 10;
    let n_trans = nocc * nvirt;

    let mut df_ia_2d = Array2::zeros((n_trans, naux));
    let mut idx = 0;
    for i in 0..nocc {
        for a in 0..nvirt {
            for p in 0..naux {
                df_ia_2d[[idx, p]] = df_ia[[i, a, p]];
            }
            idx += 1;
        }
    }

    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia_2d.clone()).unwrap();

    // Simulate multiple evGW iterations with damping
    let mut qp_energies = mo_energies.clone();
    let mut gap_history = Vec::new();

    for iter in 0..5 {
        // Apply decreasing corrections (simulating convergence)
        let damping = 0.5_f64.powi(iter);
        for i in 0..nocc {
            qp_energies[i] -= 0.02 * damping;
        }
        for i in nocc..(nocc + nvirt) {
            qp_energies[i] += 0.02 * damping;
        }

        let stats = builder.update_energies(&qp_energies).unwrap();
        gap_history.push(stats.mean_gap);
    }

    // Check that gaps increase monotonically (or stay constant)
    for i in 1..gap_history.len() {
        assert!(gap_history[i] >= gap_history[i-1] - 1e-10);
    }
}

#[test]
fn test_negative_gap_detection() {
    let (mo_energies, _) = create_test_system();

    let nocc = 3;
    let nvirt = 4;
    let naux = 10;

    let df_ia_2d = Array2::zeros((nocc * nvirt, naux));
    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia_2d).unwrap();

    // Create pathological QP energies with negative gap
    let mut qp_energies = mo_energies.clone();
    // Make HOMO higher than LUMO (unphysical but tests detection)
    qp_energies[2] = 0.3;  // HOMO now above original LUMO
    qp_energies[3] = -0.1; // LUMO now below original HOMO

    let stats = builder.update_energies(&qp_energies).unwrap();

    // Should have detected negative gaps
    // The gap between HOMO (index 2) and LUMO (index 3) is -0.1 - 0.3 = -0.4
    assert!(stats.n_negative > 0);
    assert!(stats.n_thresholded > 0); // Negative gaps also get thresholded
}

#[test]
fn test_p0_build_with_updated_denominators() {
    let (mo_energies, df_ia) = create_test_system();

    let nocc = 3;
    let nvirt = 4;
    let naux = 10;
    let n_trans = nocc * nvirt;

    // Create non-zero DF tensor for meaningful test
    let mut df_ia_2d = Array2::zeros((n_trans, naux));
    let mut idx = 0;
    for i in 0..nocc {
        for a in 0..nvirt {
            for p in 0..naux {
                // Non-zero values to see actual changes
                df_ia_2d[[idx, p]] = 0.1 * (1.0 + idx as f64) * (1.0 + p as f64).sqrt();
            }
            idx += 1;
        }
    }

    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia_2d).unwrap();

    // Build P0 at imaginary frequency with initial energies
    let omega = num_complex::Complex64::new(0.0, 1.0);
    let p0_initial = builder.build_p0(omega).unwrap();

    // Update energies (increase gaps)
    let mut qp_energies = mo_energies.clone();
    for i in 0..nocc {
        qp_energies[i] -= 0.1;
    }
    for i in nocc..(nocc + nvirt) {
        qp_energies[i] += 0.1;
    }

    builder.update_energies(&qp_energies).unwrap();

    // Build P0 with updated energies
    let p0_updated = builder.build_p0(omega).unwrap();

    // With larger gaps, P0 should generally decrease in magnitude
    let norm_initial: f64 = p0_initial.iter().map(|x| x.norm_sqr()).sum();
    let norm_updated: f64 = p0_updated.iter().map(|x| x.norm_sqr()).sum();

    // This assertion depends on the specific frequency and gap changes
    // For imaginary frequency and increased gaps, P0 typically decreases
    println!("P0 norm initial: {}, updated: {}", norm_initial, norm_updated);

    // At minimum, they should be different
    assert!((norm_initial - norm_updated).abs() > 1e-10);
}

/// Test that the gap statistics properly track the minimum gap
#[test]
fn test_minimum_gap_tracking() {
    let (mo_energies, df_ia) = create_test_system();

    let nocc = 3;
    let nvirt = 4;
    let naux = 10;
    let n_trans = nocc * nvirt;

    let mut df_ia_2d = Array2::zeros((n_trans, naux));
    let mut idx = 0;
    for i in 0..nocc {
        for a in 0..nvirt {
            for p in 0..naux {
                df_ia_2d[[idx, p]] = df_ia[[i, a, p]];
            }
            idx += 1;
        }
    }

    let mut builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia_2d).unwrap();

    let stats = builder.update_energies(&mo_energies).unwrap();

    // The minimum gap should be LUMO - HOMO = 0.2 - (-0.4) = 0.6
    assert!((stats.min_gap - 0.6).abs() < 1e-10);

    // The maximum gap should be highest virtual - lowest occupied = 0.8 - (-0.8) = 1.6
    assert!((stats.max_gap - 1.6).abs() < 1e-10);
}