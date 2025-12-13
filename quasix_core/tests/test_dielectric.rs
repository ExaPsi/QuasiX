//! Integration tests for dielectric function and polarizability computation

use approx::assert_relative_eq;
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::dielectric::{DielectricMatrix, PolarizabilityRI};

/// Setup test system (H2O-like)
fn setup_test_system() -> (Array2<f64>, Array1<f64>, Array2<f64>, usize, usize) {
    let nocc = 5;
    let nvirt = 3;
    let naux = 18;

    // Generate mock DF tensor (ia|P)
    let df_tensor = Array2::from_shape_fn((nocc * nvirt, naux), |(ia, p)| {
        ((ia + 1) as f64 * (p + 1) as f64).sin() * 0.1
    });

    // Generate orbital energies
    let mut orbital_energies = Array1::zeros(nocc + nvirt);
    for i in 0..nocc {
        orbital_energies[i] = -0.5 - 0.1 * i as f64; // Occupied
    }
    for a in 0..nvirt {
        orbital_energies[nocc + a] = 0.1 + 0.2 * a as f64; // Virtual
    }

    // Generate positive definite metric matrix
    let mut metric = Array2::eye(naux);
    for i in 0..naux {
        for j in 0..naux {
            if i != j {
                metric[[i, j]] = 0.05 * ((i + j) as f64).sin();
                metric[[j, i]] = metric[[i, j]];
            }
        }
        metric[[i, i]] = 2.0 + 0.1 * i as f64;
    }

    (df_tensor, orbital_energies, metric, nocc, nvirt)
}

#[test]
fn test_p0_computation() {
    let (df_tensor, orbital_energies, _, nocc, nvirt) = setup_test_system();

    // Use PolarizabilityRI which is the actual type
    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);

    // Test at imaginary frequency
    let omega = Complex64::new(0.0, 1.0);

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    let p0 = p0_calc
        .compute_p0(omega, &df_tensor, &e_occ, &e_virt)
        .unwrap();

    assert_eq!(p0.dim(), (18, 18));

    // P⁰ should be Hermitian
    for i in 0..18 {
        // Diagonal should be real
        assert!(p0[[i, i]].im.abs() < 1e-10);

        // Off-diagonal Hermiticity
        for j in i + 1..18 {
            let diff = (p0[[i, j]] - p0[[j, i]].conj()).norm();
            assert!(
                diff < 1e-10,
                "P⁰ not Hermitian at [{},{}]: diff = {}",
                i,
                j,
                diff
            );
        }
    }
}

#[test]
fn test_dielectric_computation() {
    let (df_tensor, orbital_energies, metric, nocc, nvirt) = setup_test_system();

    // Compute P⁰
    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);
    let omega = Complex64::new(0.0, 1.0);

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    let p0 = p0_calc
        .compute_p0(omega, &df_tensor, &e_occ, &e_virt)
        .unwrap();

    // Compute dielectric
    let dielectric = DielectricMatrix::new(&metric).unwrap();
    let m = dielectric.compute_m_matrix(&p0).unwrap();
    let epsilon = dielectric.compute_epsilon(&m).unwrap();

    assert_eq!(epsilon.dim(), (18, 18));

    // Diagonal should be close to 1 for small P⁰
    for i in 0..18 {
        assert!(
            epsilon[[i, i]].re > 0.5,
            "Diagonal element {} = {}",
            i,
            epsilon[[i, i]].re
        );
        assert!(
            epsilon[[i, i]].re < 2.0,
            "Diagonal element {} = {}",
            i,
            epsilon[[i, i]].re
        );
    }
}

#[test]
fn test_screening_computation() {
    let (df_tensor, orbital_energies, metric, nocc, nvirt) = setup_test_system();

    // Compute P⁰
    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);
    let omega = Complex64::new(0.0, 1.0);

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    let p0 = p0_calc
        .compute_p0(omega, &df_tensor, &e_occ, &e_virt)
        .unwrap();

    // Compute screening
    let dielectric = DielectricMatrix::new(&metric).unwrap();
    let m = dielectric.compute_m_matrix(&p0).unwrap();
    let epsilon = dielectric.compute_epsilon(&m).unwrap();
    let epsilon_inv = dielectric.compute_epsilon_inverse(&epsilon).unwrap();
    let w = dielectric.compute_screened_coulomb(&epsilon_inv).unwrap();

    assert_eq!(w.dim(), (18, 18));

    // W should be Hermitian
    for i in 0..18 {
        for j in i..18 {
            let diff = (w[[i, j]] - w[[j, i]].conj()).norm();
            assert!(
                diff < 1e-8,
                "W not Hermitian at [{},{}]: diff = {}",
                i,
                j,
                diff
            );
        }
    }
}

#[test]
fn test_frequency_grid() {
    let (df_tensor, orbital_energies, _, nocc, nvirt) = setup_test_system();

    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);

    // Test multiple frequencies
    let frequencies = [
        Complex64::new(0.0, 1.0),
        Complex64::new(0.0, 2.0),
        Complex64::new(0.0, 3.0),
    ];

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    let mut p0_grid = Vec::new();
    for omega in frequencies.iter() {
        let p0 = p0_calc
            .compute_p0(*omega, &df_tensor, &e_occ, &e_virt)
            .unwrap();
        p0_grid.push(p0);
    }

    assert_eq!(p0_grid.len(), 3);
    for p0 in p0_grid.iter() {
        assert_eq!(p0.dim(), (18, 18));
    }
}

#[test]
fn test_p0_reproducibility() {
    let (df_tensor, orbital_energies, _, nocc, nvirt) = setup_test_system();

    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);

    let omega = Complex64::new(0.0, 1.0);

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    // First computation
    let p0_1 = p0_calc
        .compute_p0(omega, &df_tensor, &e_occ, &e_virt)
        .unwrap();

    // Second computation
    let p0_2 = p0_calc
        .compute_p0(omega, &df_tensor, &e_occ, &e_virt)
        .unwrap();

    // Results should be identical
    for i in 0..18 {
        for j in 0..18 {
            assert_eq!(p0_1[[i, j]], p0_2[[i, j]]);
        }
    }
}

#[test]
fn test_transition_energies() {
    let (_, orbital_energies, _, nocc, nvirt) = setup_test_system();

    // Calculate transition energies manually
    let mut trans_energies = Vec::new();
    for i in 0..nocc {
        for a in 0..nvirt {
            let e_occ = orbital_energies[i];
            let e_virt = orbital_energies[nocc + a];
            trans_energies.push(e_virt - e_occ);
        }
    }

    assert_eq!(trans_energies.len(), nocc * nvirt);

    // All transition energies should be positive (HOMO-LUMO gap)
    for energy in trans_energies.iter() {
        assert!(*energy > 0.0, "Negative transition energy: {}", energy);
    }

    // Check specific values
    let e_homo = orbital_energies[nocc - 1];
    let e_lumo = orbital_energies[nocc];
    let homo_lumo_gap = e_lumo - e_homo;

    // Last transition should be HOMO->LUMO
    assert_relative_eq!(
        trans_energies[(nocc - 1) * nvirt],
        homo_lumo_gap,
        epsilon = 1e-10
    );
}

#[test]
fn test_numerical_stability() {
    let (df_tensor, orbital_energies, _, nocc, nvirt) = setup_test_system();

    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    // Compute at transition energy (small denominator)
    let trans_energy = e_virt[0] - e_occ[nocc - 1];
    let omega = Complex64::new(trans_energy, 1e-10); // Small imaginary part for stability

    let p0 = p0_calc
        .compute_p0(omega, &df_tensor, &e_occ, &e_virt)
        .unwrap();

    // Should still get finite results
    for val in p0.iter() {
        assert!(val.is_finite(), "Non-finite value in P⁰");
    }
}

#[test]
fn test_parallel_computation() {
    let (df_tensor, orbital_energies, _, nocc, nvirt) = setup_test_system();

    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);

    let frequencies = vec![
        Complex64::new(0.0, 1.0),
        Complex64::new(0.0, 2.0),
        Complex64::new(0.0, 3.0),
        Complex64::new(0.0, 4.0),
    ];

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    // Compute P⁰ for each frequency using rayon for parallelism
    use rayon::prelude::*;
    let p0_batch: Vec<_> = frequencies
        .par_iter()
        .map(|&omega| {
            p0_calc
                .compute_p0(omega, &df_tensor, &e_occ, &e_virt)
                .unwrap()
        })
        .collect();

    assert_eq!(p0_batch.len(), 4);
    for p0 in p0_batch.iter() {
        assert_eq!(p0.dim(), (18, 18));
    }
}

#[test]
fn test_dielectric_batch_processing() {
    let (df_tensor, orbital_energies, metric, nocc, nvirt) = setup_test_system();

    let p0_calc = PolarizabilityRI::new(nocc, nvirt, 18);
    let dielectric = DielectricMatrix::new(&metric).unwrap();

    let frequencies = [Complex64::new(0.0, 1.0), Complex64::new(0.0, 2.0)];

    // Extract occupied and virtual energies
    let e_occ = orbital_energies.slice(ndarray::s![..nocc]).to_owned();
    let e_virt = orbital_energies.slice(ndarray::s![nocc..]).to_owned();

    for omega in frequencies.iter() {
        let p0 = p0_calc
            .compute_p0(*omega, &df_tensor, &e_occ, &e_virt)
            .unwrap();
        let m = dielectric.compute_m_matrix(&p0).unwrap();
        let epsilon = dielectric.compute_epsilon(&m).unwrap();

        // Check that epsilon is valid
        assert_eq!(epsilon.dim(), (18, 18));
        for i in 0..18 {
            assert!(epsilon[[i, i]].re > 0.0, "Non-positive diagonal element");
        }
    }
}
