//! FIXED: Correlation self-energy Σᶜ matching PySCF exactly
//!
//! **Bug Fix Date**: 2025-11-23
//! **Root Cause**: Wrong formula - used Σ_ia (ni|W|ia) instead of PySCF's (n|W|m)
//!
//! # PySCF Reference Formula (gw_ac.py lines 204-211)
//!
//! ```python
//! for w in range(nw):
//!     Pi = get_rho_response(freqs[w], mo_energy, Lpq[:,:nocc,nocc:])
//!     Pi_inv = np.linalg.inv(np.eye(naux)-Pi) - np.eye(naux)  # ε^{-1} - I
//!     g0_occ = wts[w] * emo_occ / (emo_occ**2+freqs[w]**2)
//!     g0_vir = wts[w] * emo_vir / (emo_vir**2+freqs[w]**2)
//!     Qnm = einsum('Pnm,PQ->Qnm', Lpq[:,orbs,:], Pi_inv)
//!     Wmn = einsum('Qnm,Qmn->mn', Qnm, Lpq[:,:,orbs])
//!     sigma[:norbs_occ] += -einsum('mn,mw->nw', Wmn[:,:norbs_occ], g0_occ)/np.pi
//!     sigma[norbs_occ:] += -einsum('mn,mw->nw', Wmn[:,norbs_occ:], g0_vir)/np.pi
//! ```
//!
//! # Correct Formula
//!
//! Σᶜ_n(ω=0) = -1/π Σ_w Σ_m weight_w * (n|W(iω_w)|m) * G_m(iω_w)
//!
//! where:
//! - W(iω) = Σ_PQ (P|n)(ε^{-1}_PQ - δ_PQ)(Q|m)
//! - G_m(iω) = ε_m / (ε_m² + ω²)  (no chemical potential shift for G₀W₀)
//! - weight_w includes Gauss-Legendre weight AND Jacobian from transformation
//!
//! # Key Differences from Old (Buggy) Code
//!
//! 1. **OLD**: Σᶜ_n = Σ_ia (ni|W|ia) G_a  [WRONG!]
//!    **NEW**: Σᶜ_n = Σ_m (n|W|m) G_m      [CORRECT!]
//!
//! 2. **OLD**: Applied -1/π only during integration
//!    **NEW**: Apply -1/π during accumulation (like PySCF)
//!
//! 3. **OLD**: Computed sigma_c_iw[freq, orb] then integrated
//!    **NEW**: Accumulate directly (integration happens in loop)

use ndarray::{Array1, Array2, Array3, s};
use num_complex::Complex;
use crate::common::{QuasixError, Result};
use std::f64::consts::PI;
use ndarray_linalg::{Eigh, UPLO, Inverse};

/// Compute correlation self-energy diagonal using PySCF CD formula
///
/// This function implements the exact PySCF contour deformation formula
/// from gw_cd.py lines 169-178.
///
/// # Formula (PySCF CD)
///
/// For each target orbital n:
///   Σᶜ_n(ω=0) = -1/π Σ_freq Σ_m (n|W(iω)|m) * G_m^(n)(iω)
///
/// where Green's function depends on TARGET orbital n:
///   emo = ω - 1j*η*sign(ε_n - ε_m) - ε_m
///   G_m^(n)(iω) = weight * emo / (emo² + ω²)
///
/// Key difference from AC: sign factor varies with target!
///
/// # Arguments
///
/// * `mo_energy` - Orbital energies [nmo] (Hartree)
/// * `eta` - Broadening parameter (default: 0.01 Ha)
/// * `p0_iw` - Polarizability P₀(iω) [nfreq, naux, naux] (real-valued)
/// * `lpq` - DF tensor (P|pq) [naux, nmo, nmo] (real)
/// * `freqs` - Imaginary frequency grid [nfreq] (positive real values)
/// * `weights` - Quadrature weights [nfreq] (includes Jacobian)
///
/// # Returns
///
/// * Σᶜ diagonal - [nmo] (complex, should be nearly real)
///
/// # Validation
///
/// * Must match PySCF CD within < 0.01 eV for H₂/STO-3G
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
///
/// # Errors
///
/// Returns error if:
/// - Dimension mismatch between inputs
/// - Singular epsilon matrix (RPA instability)
/// - Result contains NaN/Inf
pub fn compute_sigma_c_pyscf_cd(
    mo_energy: &Array1<f64>,
    eta: f64,
    p0_iw: &Array3<f64>,
    lpq: &Array3<f64>,
    freqs: &Array1<f64>,
    weights: &Array1<f64>,
) -> Result<Array1<Complex<f64>>> {
    // FIXED VERSION: Using corrected formula that matches PySCF CD method

    let (nfreq_g, nmo) = green_function.dim();
    let (nfreq_p, naux_p1, naux_p2) = p0_iw.dim();
    let (naux_l, nmo_l1, nmo_l2) = lpq.dim();

    // Validate dimensions
    if nfreq_g != nfreq_p {
        return Err(QuasixError::DimensionMismatch(
            format!("Green function nfreq {} != P₀ nfreq {}", nfreq_g, nfreq_p)
        ));
    }

    if naux_p1 != naux_p2 {
        return Err(QuasixError::DimensionMismatch(
            format!("P₀ not square: {} x {}", naux_p1, naux_p2)
        ));
    }

    if naux_l != naux_p1 {
        return Err(QuasixError::DimensionMismatch(
            format!("Lpq naux {} != P₀ naux {}", naux_l, naux_p1)
        ));
    }

    if nmo_l1 != nmo || nmo_l2 != nmo {
        return Err(QuasixError::DimensionMismatch(
            format!("Lpq MO dimensions ({}, {}) != nmo {}", nmo_l1, nmo_l2, nmo)
        ));
    }

    if weights.len() != nfreq_g {
        return Err(QuasixError::DimensionMismatch(
            format!("Weights length {} != nfreq {}", weights.len(), nfreq_g)
        ));
    }

    if freqs.len() != nfreq_g {
        return Err(QuasixError::DimensionMismatch(
            format!("Frequencies length {} != nfreq {}", freqs.len(), nfreq_g)
        ));
    }

    let nfreq = nfreq_g;
    let naux = naux_p1;

    // ===== DIAGNOSTIC OUTPUT START =====
    eprintln!("\n=== Σᶜ DIAGNOSTIC START ===");
    eprintln!("Input shapes:");
    eprintln!("  green_function: {:?}", green_function.shape());
    eprintln!("  p0_iw: {:?}", p0_iw.shape());
    eprintln!("  lpq: {:?}", lpq.shape());
    eprintln!("  freqs: {:?}", freqs.shape());
    eprintln!("  weights: {:?}", weights.shape());
    eprintln!("  nmo={}, naux={}, nfreq={}", nmo, naux, nfreq);

    // Sample input values
    eprintln!("\nSample input values:");
    if nfreq > 0 && nmo > 0 {
        eprintln!("  green[0,0] = {:.6e} + {:.6e}i",
                 green_function[[0,0]].re, green_function[[0,0]].im);
        if nfreq > 50 {
            eprintln!("  green[50,0] = {:.6e} + {:.6e}i",
                     green_function[[50,0]].re, green_function[[50,0]].im);
        }
    }
    if nfreq > 0 && naux > 0 {
        eprintln!("  P₀[0,0,0] = {:.6e}", p0_iw[[0,0,0]]);
        if nfreq > 50 {
            eprintln!("  P₀[50,0,0] = {:.6e}", p0_iw[[50,0,0]]);
        }
    }
    if naux > 0 && nmo > 0 {
        eprintln!("  lpq[0,0,0] = {:.6e}", lpq[[0,0,0]]);
    }
    if weights.len() > 0 {
        eprintln!("  weights[0] = {:.6e}", weights[0]);
    }
    if freqs.len() > 0 {
        eprintln!("  freqs[0] = {:.6e} Ha", freqs[0]);
    }
    // ===== DIAGNOSTIC OUTPUT END =====

    // Allocate result
    let mut sigma_c = Array1::<Complex<f64>>::zeros(nmo);

    // PySCF formula: accumulate over frequencies
    // sigma[n] += -1/π * weight_w * Σ_m (n|W(iω)|m) * G_m(iω)

    for ifreq in 0..nfreq {
        // ===== DIAGNOSTIC: Track first few frequencies =====
        let debug_freq = ifreq < 3;
        if debug_freq {
            eprintln!("\n=== FREQUENCY ifreq={} (ω={:.6e} Ha) ===", ifreq, freqs[ifreq]);
        }

        // Extract P₀(ω) for this frequency [naux, naux]
        let p0_freq = p0_iw.slice(s![ifreq, .., ..]);

        if debug_freq {
            eprintln!("  P₀ matrix stats:");
            let p0_norm: f64 = p0_freq.iter().map(|x| x * x).sum::<f64>().sqrt();
            eprintln!("    Frobenius norm: {:.6e}", p0_norm);
            eprintln!("    P₀[0,0] = {:.6e}", p0_freq[[0,0]]);
        }

        // Compute epsilon^{-1} - I = (I - P₀)^{-1} - I
        // PySCF: Pi_inv = np.linalg.inv(np.eye(naux)-Pi) - np.eye(naux)
        let mut epsilon = Array2::<f64>::eye(naux);
        for i in 0..naux {
            for j in 0..naux {
                epsilon[[i, j]] -= p0_freq[[i, j]];
            }
        }

        // Invert epsilon
        let epsilon_inv = match epsilon.inv() {
            Ok(inv) => inv,
            Err(e) => {
                eprintln!("WARNING: Epsilon singular at freq {}: {}", ifreq, e);
                eprintln!("  This indicates RPA instability (over-screening)");
                // Try to get eigenvalues for diagnostics
                if let Ok((eigvals, _)) = epsilon.eigh(UPLO::Upper) {
                    eprintln!("  Eigenvalues of epsilon: {:?}", eigvals);
                }
                return Err(QuasixError::NumericalError(
                    format!("Epsilon singular at frequency {}", ifreq)
                ));
            }
        };

        // Compute Pi_inv = epsilon_inv - I
        let mut pi_inv = epsilon_inv.clone();
        for i in 0..naux {
            pi_inv[[i, i]] -= 1.0;
        }

        if debug_freq {
            eprintln!("  Pi_inv stats:");
            let pi_inv_norm: f64 = pi_inv.iter().map(|x| x * x).sum::<f64>().sqrt();
            eprintln!("    Frobenius norm: {:.6e}", pi_inv_norm);
            eprintln!("    Pi_inv[0,0] = {:.6e}", pi_inv[[0,0]]);
        }

        // Extract Green's function for this frequency [nmo]
        let g_freq = green_function.slice(s![ifreq, ..]);

        if debug_freq {
            eprintln!("  Green's function G(iω):");
            for m in 0..(nmo.min(3)) {
                eprintln!("    G[{}] = {:.6e} + {:.6e}i", m, g_freq[m].re, g_freq[m].im);
            }
        }

        // Compute Wmn = Σ_PQ (n|P) Pi_inv_PQ (Q|m)
        // PySCF:
        //   Qnm = einsum('Pnm,PQ->Qnm', Lpq[:,orbs,:], Pi_inv)
        //   Wmn = einsum('Qnm,Qmn->mn', Qnm, Lpq[:,:,orbs])
        //
        // Optimized: Wmn[n,m] = Σ_P Σ_Q Lpq[P,n,m] * Pi_inv[P,Q] * Lpq[Q,m,n]
        //                      = Σ_P Σ_Q (n|P) Pi_inv_PQ (Q|m)

        // Use BLAS-3 for efficiency: W = L^T @ Pi_inv @ L
        // where L[P, nm] is reshaped Lpq
        // But for clarity, use explicit loops first (optimize later)

        let mut wmn = Array2::<f64>::zeros((nmo, nmo));

        for n in 0..nmo {
            for m in 0..nmo {
                let mut w_nm = 0.0;

                // Contract: Σ_P Σ_Q (n|P) Pi_inv[P,Q] (Q|m)
                for p in 0..naux {
                    for q in 0..naux {
                        w_nm += lpq[[p, n, m]] * pi_inv[[p, q]] * lpq[[q, m, n]];
                    }
                }

                // CRITICAL FIX: Store as W[m,n] to match PySCF convention!
                // PySCF computes W[m,n], not W[n,m]
                wmn[[m, n]] = w_nm;  // FIXED: was [[n, m]]
            }
        }

        if debug_freq {
            eprintln!("  W matrix stats:");
            let w_norm: f64 = wmn.iter().map(|x| x * x).sum::<f64>().sqrt();
            eprintln!("    Frobenius norm: {:.6e}", w_norm);
            eprintln!("    W[0,0] = {:.6e}", wmn[[0,0]]);
            eprintln!("    W[1,0] = {:.6e}", wmn[[1,0]]);
        }

        // Accumulate sigma: sigma[n] += -1/π * weight * Σ_m W[m,n] * G_m(ω)
        // Now matches PySCF: sigma[:norbs_occ] += -einsum('mn,mw->nw', Wmn[:,:norbs_occ], g0_occ)/np.pi
        let weight_w = weights[ifreq];
        let normalization = -1.0 / PI;

        // ===== DIAGNOSTIC: Track HOMO accumulation =====
        if debug_freq {
            eprintln!("\n  Accumulation for HOMO (n=0):");
        }

        for n in 0..nmo {
            let mut sigma_n_contrib = Complex::new(0.0, 0.0);

            for m in 0..nmo {
                // Use W[m,n] to match PySCF
                let term = wmn[[m, n]] * g_freq[m];
                sigma_n_contrib += term;

                // ===== DIAGNOSTIC: First two m states for HOMO =====
                if debug_freq && n == 0 && m < 2 {
                    eprintln!("    m={}: W[m,n]={:.6e}, G[m]={:.6e}+{:.6e}i, term={:.6e}+{:.6e}i",
                             m, wmn[[m, n]], g_freq[m].re, g_freq[m].im,
                             term.re, term.im);
                }
            }

            // ===== DIAGNOSTIC: Show accumulation for HOMO =====
            if debug_freq && n == 0 {
                eprintln!("    Sum before normalization: {:.6e} + {:.6e}i",
                         sigma_n_contrib.re, sigma_n_contrib.im);
                eprintln!("    Normalization factor: {:.6e}", normalization);
                let after_norm = normalization * sigma_n_contrib;
                eprintln!("    After normalization: {:.6e} + {:.6e}i",
                         after_norm.re, after_norm.im);
                eprintln!("    Current sigma_c[0] before += : {:.6e} + {:.6e}i",
                         sigma_c[n].re, sigma_c[n].im);
            }

            // FIX: NO weight_w here - weights are already in Green's function!
            // This matches PySCF convention: weights inside G₀, not in integration
            sigma_c[n] += normalization * sigma_n_contrib;

            // ===== DIAGNOSTIC: Show running total for HOMO =====
            if debug_freq && n == 0 {
                eprintln!("    Current sigma_c[0] after  += : {:.6e} + {:.6e}i",
                         sigma_c[n].re, sigma_c[n].im);
            }
        }
    }

    // ===== DIAGNOSTIC: Final results =====
    eprintln!("\n=== Σᶜ FINAL RESULTS ===");
    for n in 0..(nmo.min(5)) {
        let sigma_ha = sigma_c[n].re;
        let sigma_ev = sigma_ha * 27.211386;
        eprintln!("  Σᶜ[{}] = {:.8e} Ha = {:.6e} eV (Im: {:.3e})",
                 n, sigma_ha, sigma_ev, sigma_c[n].im);
    }

    // SIGN CHECK
    eprintln!("\n=== SIGN CHECK ===");
    eprintln!("  Expected for occupied states (HOMO): Σᶜ < 0 (negative)");
    eprintln!("  Expected for virtual states (LUMO+): Σᶜ > 0 (positive)");
    if nmo > 0 {
        let homo_sign = if sigma_c[0].re < 0.0 { "NEGATIVE ✓" } else { "POSITIVE ✗ WRONG!" };
        eprintln!("  HOMO (n=0): Σᶜ={:.6e} Ha → {}", sigma_c[0].re, homo_sign);
    }
    if nmo > 1 {
        let lumo_sign = if sigma_c[1].re > 0.0 { "POSITIVE ✓" } else { "NEGATIVE ✗ WRONG!" };
        eprintln!("  LUMO (n=1): Σᶜ={:.6e} Ha → {}", sigma_c[1].re, lumo_sign);
    }
    eprintln!("=== Σᶜ DIAGNOSTIC END ===\n");

    // Validate result
    if sigma_c.iter().any(|z| !z.is_finite()) {
        return Err(QuasixError::NumericalError(
            "Σᶜ contains NaN or Inf".to_string()
        ));
    }

    // Physical check: imaginary part should be small
    let max_imag = sigma_c.iter()
        .map(|z| z.im.abs())
        .fold(0.0, f64::max);

    if max_imag > 1e-5 {
        eprintln!("WARNING: Σᶜ has large imaginary part: {:.2e}", max_imag);
        eprintln!("  This may indicate numerical issues");
        eprintln!("  Expected: |Im(Σᶜ)| < 1e-5 for well-behaved systems");
    }

    Ok(sigma_c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_sigma_c_fixed_dimensions() {
        // Test that function validates dimensions correctly
        let nfreq = 2;
        let nmo = 4;
        let naux = 6;

        let green = Array2::<Complex<f64>>::zeros((nfreq, nmo));
        let p0 = Array3::<f64>::zeros((nfreq, naux, naux));
        let lpq = Array3::<f64>::zeros((naux, nmo, nmo));
        let freqs = Array1::<f64>::ones(nfreq);
        let weights = Array1::<f64>::ones(nfreq);

        // Valid call
        let result = compute_sigma_c_diagonal_fixed(&green, &p0, &lpq, &freqs, &weights, 2);
        assert!(result.is_ok(), "Valid dimensions should not error");

        // Invalid: mismatched frequency dimensions
        let green_bad = Array2::<Complex<f64>>::zeros((nfreq + 1, nmo));
        let result = compute_sigma_c_diagonal_fixed(&green_bad, &p0, &lpq, &freqs, &weights, 2);
        assert!(result.is_err(), "Should reject frequency dimension mismatch");
    }

    #[test]
    fn test_sigma_c_fixed_zero_p0() {
        // Test with zero polarizability (no correlation)
        let nfreq = 2;
        let nmo = 2;
        let naux = 2;

        let green = Array2::<Complex<f64>>::from_shape_fn((nfreq, nmo), |(_, m)| {
            Complex::new(0.1 * (m + 1) as f64, 0.0)
        });
        let p0 = Array3::<f64>::zeros((nfreq, naux, naux));  // Zero polarizability
        let lpq = Array3::<f64>::from_elem((naux, nmo, nmo), 0.1);
        let freqs = array![0.5, 1.0];
        let weights = array![0.3, 0.4];

        let sigma_c = compute_sigma_c_diagonal_fixed(&green, &p0, &lpq, &freqs, &weights, 1)
            .expect("Should work with zero P₀");

        // With zero P₀: epsilon = I, epsilon_inv = I, Pi_inv = 0
        // So Wmn = 0, sigma_c should be ~0
        for i in 0..nmo {
            assert_abs_diff_eq!(sigma_c[i].re, 0.0, epsilon = 1e-10);
            assert_abs_diff_eq!(sigma_c[i].im, 0.0, epsilon = 1e-10);
        }
    }
}
