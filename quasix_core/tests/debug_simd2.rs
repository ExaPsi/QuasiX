//! Debug test focusing on matrix structure

use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::dielectric::{PolarizabilityRI, PolarizabilitySIMDCore};

#[test]
fn debug_matrix_structure() {
    // Small test where we can inspect every element
    let nocc = 2;
    let nvirt = 3;
    let naux = 5;

    let n_trans = nocc * nvirt;

    // Create distinct DF values to trace them
    let mut df_ia = Array2::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            // Make each element unique
            df_ia[[i, p]] = (i * naux + p + 1) as f64;
        }
    }

    // Simple energies
    let e_occ = Array1::from_vec(vec![-1.0, -0.5]);
    let e_virt = Array1::from_vec(vec![0.5, 1.0, 1.5]);

    // Static frequency for simplicity
    let omega = Complex64::new(0.0, 1e-4);

    // Compute with both
    let standard_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let simd_calc = PolarizabilitySIMDCore::new(nocc, nvirt, naux);

    let p0_standard = standard_calc
        .compute_p0(omega, &df_ia, &e_occ, &e_virt)
        .unwrap();
    let p0_simd = simd_calc
        .compute_p0_simd(omega, &df_ia, &e_occ, &e_virt)
        .unwrap();

    println!("DF tensor shape: {} x {}", n_trans, naux);

    // Print a few key elements to see the pattern
    println!("\nDiagonal elements:");
    for i in 0..naux.min(5) {
        println!("  Standard[{},{}] = {:.4}", i, i, p0_standard[[i, i]].re);
        println!("  SIMD[{},{}] = {:.4}", i, i, p0_simd[[i, i]].re);
        println!(
            "  Ratio = {:.4}",
            p0_simd[[i, i]].re / p0_standard[[i, i]].re
        );
        println!();
    }

    // Check symmetry
    println!("Off-diagonal check (should be Hermitian):");
    for i in 0..2 {
        for j in i + 1..3 {
            let std_ij = p0_standard[[i, j]];
            let std_ji = p0_standard[[j, i]];
            let simd_ij = p0_simd[[i, j]];
            let simd_ji = p0_simd[[j, i]];

            println!(
                "  Standard[{},{}] = {:.4}, [{},{}] = {:.4}",
                i, j, std_ij.re, j, i, std_ji.re
            );
            println!(
                "  SIMD[{},{}] = {:.4}, [{},{}] = {:.4}",
                i, j, simd_ij.re, j, i, simd_ji.re
            );

            // Check Hermiticity
            assert!(
                (std_ij - std_ji.conj()).norm() < 1e-10,
                "Standard not Hermitian"
            );
            assert!(
                (simd_ij - simd_ji.conj()).norm() < 1e-10,
                "SIMD not Hermitian"
            );
        }
    }

    // Compute norms
    let std_norm: f64 = p0_standard.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt();
    let simd_norm: f64 = p0_simd.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt();

    println!("\nFrobenius norms:");
    println!("  Standard: {:.4}", std_norm);
    println!("  SIMD: {:.4}", simd_norm);
    println!("  Ratio: {:.4}", simd_norm / std_norm);
}
