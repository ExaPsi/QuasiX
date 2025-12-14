// ! PySCF Contour Deformation (CD) implementation for correlation self-energy
//!
//! This module implements the exact PySCF CD formula from gw_cd.py
//! lines 169-178, including target-dependent Green's function.

use crate::common::{QuasixError, Result};
use ndarray::{s, Array1, Array2, Array3};
use ndarray_linalg::Inverse;
use num_complex::Complex;
use std::f64::consts::PI;

/// Compute correlation self-energy using PySCF CD formula
///
/// This implements the contour deformation Green's function:
///   emo = omega - 1j*eta*sign(epsilon_n - epsilon_m) - epsilon_m
///   G_m^(n)(iω) = weight * emo / (emo² + ω²)
///
/// Key difference from AC: Green's function depends on TARGET orbital!
///
/// # Formula
///
/// For each target orbital n:
///   Σᶜ_n(ω=0) = -1/π Σ_w Σ_m (n|W(iω_w)|m) * G_m^(n)(iω_w)
///
/// where:
/// - W(iω) = Σ_PQ (P|nm) [ε^{-1}_PQ - δ_PQ] (Q|mn)
/// - sign = sign(ε_n - ε_m)  (CRITICAL: depends on target!)
/// - emo = -1j*η*sign - ε_m  (for ω=0 static self-energy)
/// - G_m^(n)(iω) = weight * emo / (emo² + ω²)
///
/// # Arguments
///
/// * `mo_energy` - Orbital energies [nmo] (Hartree)
/// * `eta` - Broadening parameter (Ha, typically 0.01)
/// * `p0_iw` - Polarizability P₀(iω) [nfreq, naux, naux]
/// * `lpq` - DF tensor (P|pq) [naux, nmo, nmo]
/// * `freqs` - Imaginary frequency grid [nfreq] (Ha)
/// * `weights` - Quadrature weights [nfreq]
///
/// # Returns
///
/// Complex Σᶜ diagonal [nmo] (should be nearly real)
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_cd.py::get_sigmaI_diag()` lines 169-178
///
/// # Examples
///
/// ```ignore
/// let sigma_c = compute_sigma_c_pyscf_cd(
///     &mo_energy,
///     0.01,  // eta
///     &p0_iw,
///     &lpq,
///     &freqs,
///     &weights,
/// )?;
/// ```
pub fn compute_sigma_c_pyscf_cd(
    mo_energy: &Array1<f64>,
    eta: f64,
    p0_iw: &Array3<f64>,
    lpq: &Array3<f64>,
    freqs: &Array1<f64>,
    weights: &Array1<f64>,
) -> Result<Array1<Complex<f64>>> {
    let nmo = mo_energy.len();
    let (nfreq, naux_p1, naux_p2) = p0_iw.dim();
    let (naux_l, nmo_l1, nmo_l2) = lpq.dim();

    // Validate dimensions
    if naux_p1 != naux_p2 {
        return Err(QuasixError::DimensionMismatch(format!(
            "P₀ not square: {} x {}",
            naux_p1, naux_p2
        )));
    }

    if naux_l != naux_p1 {
        return Err(QuasixError::DimensionMismatch(format!(
            "Lpq naux {} != P₀ naux {}",
            naux_l, naux_p1
        )));
    }

    if nmo_l1 != nmo || nmo_l2 != nmo {
        return Err(QuasixError::DimensionMismatch(format!(
            "Lpq MO dimensions ({}, {}) != nmo {}",
            nmo_l1, nmo_l2, nmo
        )));
    }

    if weights.len() != nfreq {
        return Err(QuasixError::DimensionMismatch(format!(
            "Weights length {} != nfreq {}",
            weights.len(),
            nfreq
        )));
    }

    if freqs.len() != nfreq {
        return Err(QuasixError::DimensionMismatch(format!(
            "Frequencies length {} != nfreq {}",
            freqs.len(),
            nfreq
        )));
    }

    let naux = naux_p1;

    eprintln!("\n=== PySCF CD Σᶜ CALCULATION ===");
    eprintln!("Parameters:");
    eprintln!("  nmo = {}, naux = {}, nfreq = {}", nmo, naux, nfreq);
    eprintln!("  eta = {:.4e} Ha", eta);
    eprintln!("  Orbital energies:");
    for n in 0..(nmo.min(5)) {
        eprintln!("    ε[{}] = {:.6e} Ha", n, mo_energy[n]);
    }

    // Allocate result
    let mut sigma_c = Array1::<Complex<f64>>::zeros(nmo);

    // Pre-compute W(iω) for all frequencies (same for all targets)
    // This is expensive but done once
    let mut w_iw = Array3::<f64>::zeros((nfreq, naux, naux));

    eprintln!("\nPre-computing W(iω) for all frequencies...");
    for ifreq in 0..nfreq {
        let p0_freq = p0_iw.slice(s![ifreq, .., ..]);

        // Compute epsilon^{-1} - I = (I - P₀)^{-1} - I
        let mut epsilon = Array2::<f64>::eye(naux);
        for i in 0..naux {
            for j in 0..naux {
                epsilon[[i, j]] -= p0_freq[[i, j]];
            }
        }

        let epsilon_inv = match epsilon.inv() {
            Ok(inv) => inv,
            Err(e) => {
                eprintln!("ERROR: Epsilon singular at freq {}: {}", ifreq, e);
                return Err(QuasixError::NumericalError(format!(
                    "Epsilon singular at frequency {}",
                    ifreq
                )));
            }
        };

        // Pi_inv = epsilon_inv - I
        let mut pi_inv = epsilon_inv.clone();
        for i in 0..naux {
            pi_inv[[i, i]] -= 1.0;
        }

        // Store Pi_inv for this frequency
        w_iw.slice_mut(s![ifreq, .., ..]).assign(&pi_inv);
    }

    // Now compute Σᶜ_n for each target orbital
    eprintln!("\nComputing Σᶜ for each target orbital...");

    for n_target in 0..nmo {
        let epsilon_n = mo_energy[n_target];

        if n_target < 3 {
            eprintln!(
                "\n=== TARGET ORBITAL n={} (ε_n={:.6e} Ha) ===",
                n_target, epsilon_n
            );
        }

        let mut sigma_n = Complex::new(0.0, 0.0);

        // Loop over frequencies
        for ifreq in 0..nfreq {
            let omega_freq = freqs[ifreq];
            let weight = weights[ifreq];

            // Get Pi_inv for this frequency
            let pi_inv = w_iw.slice(s![ifreq, .., ..]);

            // Compute W_mn = Σ_PQ (n|P) Pi_inv[P,Q] (Q|m)
            let mut w_mn = Array2::<f64>::zeros((nmo, nmo));

            for m in 0..nmo {
                for n in 0..nmo {
                    let mut w_val = 0.0;
                    for p in 0..naux {
                        for q in 0..naux {
                            w_val += lpq[[p, n, m]] * pi_inv[[p, q]] * lpq[[q, m, n]];
                        }
                    }
                    w_mn[[m, n]] = w_val; // W[m,n] (PySCF convention)
                }
            }

            // Compute Green's function FOR THIS TARGET at this frequency
            // PySCF CD formula (lines 174-175):
            //   emo = omega - 1j*eta*sign - mo_energy[m]
            //   g0 = wts[None,:]*emo[:,None] / ((emo**2)[:,None]+(freqs**2)[None,:])
            //
            // For static self-energy (omega=0):
            //   emo = -1j*eta*sign(epsilon_n - epsilon_m) - epsilon_m

            let omega = 0.0; // Static self-energy

            for m in 0..nmo {
                let epsilon_m = mo_energy[m];

                // CRITICAL: sign factor depends on target n
                let sign = if epsilon_n > epsilon_m { 1.0 } else { -1.0 };

                // emo = omega - 1j*eta*sign - epsilon_m
                let emo = Complex::new(
                    omega - epsilon_m, // Real part
                    -eta * sign,       // Imaginary part (contour shift)
                );

                // G_m^(n)(iω) = weight * emo / (emo² + ω²)
                let emo_squared = emo * emo;
                let freq_squared = Complex::new(omega_freq * omega_freq, 0.0);
                let denominator = emo_squared + freq_squared;

                let g_m = Complex::new(weight, 0.0) * emo / denominator;

                // Accumulate: Σᶜ_n += W[m,n] * G_m^(n)
                sigma_n += w_mn[[m, n_target]] * g_m;

                // Debug first target, first frequency, first two m
                if n_target == 0 && ifreq == 0 && m < 2 {
                    eprintln!("    m={}: ε_m={:.6e}, sign={}", m, epsilon_m, sign);
                    eprintln!("      emo = {:.6e} + {:.6e}i", emo.re, emo.im);
                    eprintln!("      G_m = {:.6e} + {:.6e}i", g_m.re, g_m.im);
                    eprintln!("      W[m,n] = {:.6e}", w_mn[[m, n_target]]);
                    eprintln!(
                        "      term = {:.6e} + {:.6e}i",
                        (w_mn[[m, n_target]] * g_m).re,
                        (w_mn[[m, n_target]] * g_m).im
                    );
                }
            }
        }

        // Apply -1/π normalization
        sigma_c[n_target] = sigma_n * (-1.0 / PI);

        if n_target < 3 {
            eprintln!(
                "  Σᶜ[{}] = {:.8e} + {:.8e}i Ha",
                n_target, sigma_c[n_target].re, sigma_c[n_target].im
            );
            eprintln!(
                "         = {:.6e} eV (real part)",
                sigma_c[n_target].re * 27.211386
            );
        }
    }

    eprintln!("\n=== FINAL RESULTS ===");
    for n in 0..(nmo.min(5)) {
        eprintln!(
            "  Σᶜ[{}] = {:.8e} Ha = {:.6e} eV (Im: {:.3e})",
            n,
            sigma_c[n].re,
            sigma_c[n].re * 27.211386,
            sigma_c[n].im
        );
    }

    // Sign check
    eprintln!("\n=== SIGN CHECK ===");
    if nmo >= 2 {
        let homo_sign = if sigma_c[0].re < 0.0 {
            "NEGATIVE ✓"
        } else {
            "POSITIVE ✗"
        };
        let lumo_sign = if sigma_c[1].re > 0.0 {
            "POSITIVE ✓"
        } else {
            "NEGATIVE ✗"
        };
        eprintln!("  HOMO (n=0): Σᶜ={:.6e} Ha → {}", sigma_c[0].re, homo_sign);
        eprintln!("  LUMO (n=1): Σᶜ={:.6e} Ha → {}", sigma_c[1].re, lumo_sign);
    }
    eprintln!("=== PySCF CD Σᶜ COMPLETE ===\n");

    // Validate result
    if sigma_c.iter().any(|z| !z.is_finite()) {
        return Err(QuasixError::NumericalError(
            "Σᶜ contains NaN or Inf".to_string(),
        ));
    }

    Ok(sigma_c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::array;

    #[test]
    fn test_green_function_sign_factor() {
        // Test that sign factor is computed correctly
        let mo_energy = array![-0.578554, 0.671143]; // H2/STO-3G
        let eta = 0.01;

        // For HOMO (n=0, ε_n = -0.578):
        //   m=0 (ε_m = -0.578): sign = 0 (ε_n = ε_m, use +1 convention)
        //   m=1 (ε_m = +0.671): sign = -1 (ε_n < ε_m)

        let epsilon_n = mo_energy[0];
        let epsilon_m = mo_energy[1];
        let sign = if epsilon_n > epsilon_m { 1.0 } else { -1.0 };

        assert_eq!(sign, -1.0, "HOMO→LUMO should have sign=-1");

        // emo = -ε_m - 1j*η*sign (for ω=0)
        let emo = Complex::new(-epsilon_m, -eta * sign);

        // Expected: emo = -0.671 + 0.01i (positive imaginary part!)
        assert_abs_diff_eq!(emo.re, -0.671143, epsilon = 1e-6);
        assert_abs_diff_eq!(emo.im, 0.01, epsilon = 1e-6);
    }

    #[test]
    fn test_sigma_c_cd_dimensions() {
        // Test dimension validation
        let nmo = 2;
        let naux = 3;
        let nfreq = 4;

        let mo_energy = Array1::<f64>::zeros(nmo);
        let p0 = Array3::<f64>::zeros((nfreq, naux, naux));
        let lpq = Array3::<f64>::zeros((naux, nmo, nmo));
        let freqs = Array1::<f64>::ones(nfreq);
        let weights = Array1::<f64>::ones(nfreq);

        let result = compute_sigma_c_pyscf_cd(&mo_energy, 0.01, &p0, &lpq, &freqs, &weights);
        assert!(result.is_ok(), "Valid dimensions should work");

        // Invalid: non-square P0
        let p0_bad = Array3::<f64>::zeros((nfreq, naux, naux + 1));
        let result = compute_sigma_c_pyscf_cd(&mo_energy, 0.01, &p0_bad, &lpq, &freqs, &weights);
        assert!(result.is_err(), "Should reject non-square P0");
    }
}
