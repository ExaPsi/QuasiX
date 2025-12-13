//! Debug test to find SIMD discrepancy

use ndarray::{Array1, Array2};
use num_complex::Complex64;
use quasix_core::dielectric::{PolarizabilityRI, PolarizabilitySIMDCore};

#[test]
fn debug_simd_issue() {
    // Test case matching the failing test
    let nocc = 5;
    let nvirt = 10;
    let naux = 30;

    // Generate DF tensor matching the test
    let n_trans = nocc * nvirt;
    let mut df_ia = Array2::zeros((n_trans, naux));
    for i in 0..n_trans {
        for p in 0..naux {
            let r_sq = ((i as f64 - n_trans as f64 / 2.0).powi(2)
                + (p as f64 - naux as f64 / 2.0).powi(2))
                / (naux as f64).powi(2);
            df_ia[[i, p]] = (-r_sq).exp();
        }
    }

    // Realistic energies
    let e_occ = Array1::linspace(-1.5, -0.3, nocc);
    let e_virt = Array1::linspace(0.05, 1.5, nvirt);

    // Static frequency
    let omega = Complex64::new(0.0, 1e-4);

    // Compute with both
    let standard_calc = PolarizabilityRI::new(nocc, nvirt, naux);
    let simd_calc = PolarizabilitySIMDCore::new(nocc, nvirt, naux);
    println!("Using SIMD: {}", simd_calc.cpu_features());

    let p0_standard = standard_calc
        .compute_p0(omega, &df_ia, &e_occ, &e_virt)
        .unwrap();
    let p0_simd = simd_calc
        .compute_p0_simd(omega, &df_ia, &e_occ, &e_virt)
        .unwrap();

    println!("First few elements:");
    println!(
        "Standard P0[0,0] = {:.6}+{:.6}i",
        p0_standard[[0, 0]].re,
        p0_standard[[0, 0]].im
    );
    println!(
        "SIMD P0[0,0] = {:.6}+{:.6}i",
        p0_simd[[0, 0]].re,
        p0_simd[[0, 0]].im
    );
    println!("Ratio = {:.6}", p0_simd[[0, 0]].re / p0_standard[[0, 0]].re);

    // Check for systematic difference
    let mut sum_standard = 0.0;
    let mut sum_simd = 0.0;
    for i in 0..naux {
        for j in 0..naux {
            sum_standard += p0_standard[[i, j]].re;
            sum_simd += p0_simd[[i, j]].re;
        }
    }
    println!("\nTotal sum:");
    println!("Standard: {:.6}", sum_standard);
    println!("SIMD: {:.6}", sum_simd);
    println!("Ratio: {:.6}", sum_simd / sum_standard);
}
