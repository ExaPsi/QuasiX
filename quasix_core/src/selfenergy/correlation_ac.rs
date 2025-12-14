//! G0W0 Correlation Self-Energy via Analytic Continuation (AC) Method
//!
//! This module implements the correlation self-energy calculation matching
//! PySCF's gw_ac.py exactly for validation < 1e-8 tolerance.
//!
//! # Physical Background
//!
//! The correlation self-energy is computed on the imaginary frequency axis
//! and analytically continued to real frequencies using Pade approximants.
//!
//! # PySCF Reference
//!
//! - File: `pyscf/gw/gw_ac.py`
//! - Functions: `get_sigma_diag()`, `thiele()`, `pade_thiele()`
//! - Reference: T. Zhu and G.K.-L. Chan, arXiv:2007.03148 (2020)
//!
//! # Key Formulas
//!
//! ## Polarizability (Eq. W2 from derivation)
//!
//! P0_PQ(iw) = 4 * sum_ia L_Pia * (e_ia / (w^2 + e_ia^2)) * L_Qia
//!
//! where factor 4 = 2 (spin) * 2 (particle-hole symmetry)
//!
//! ## Screened Interaction
//!
//! W_c = (I - P0)^{-1} - I  (this is epsilon^{-1} - I)
//!
//! ## Self-Energy on Imaginary Axis (Eq. W4a/W4b)
//!
//! For occupied orbital n:
//!   Sigma_c(iw_k) = -1/pi * sum_w wts_w * sum_m W_mn(iw_w) * G_m(-iw_k, iw_w)
//!
//! For virtual orbital n:
//!   Sigma_c(iw_k) = -1/pi * sum_w wts_w * sum_m W_mn(iw_w) * G_m(+iw_k, iw_w)
//!
//! where G_m(iw', iw) = (iw' + E_F - e_m) / ((iw' + E_F - e_m)^2 + w^2)
//!
//! # Validation Target
//!
//! H2/STO-3G reference:
//! - Sigma_c HOMO (real) = -0.040900545031 Ha
//! - Sigma_c LUMO (real) = +0.040900545031 Ha

use crate::common::{QuasixError, Result};
use ndarray::{s, Array1, Array2, Array3, Axis};
use ndarray_linalg::{Inverse, InverseC};
use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

/// Configuration for AC self-energy calculation
#[derive(Debug, Clone)]
pub struct ACConfig {
    /// Number of frequency points for W integration
    pub nw: usize,
    /// Cutoff for imaginary frequencies in self-energy grid (Hartree)
    pub iw_cutoff: f64,
    /// Frequency scaling parameter x0 (PySCF default: 0.5)
    pub x0: f64,
    /// Broadening parameter eta (unused on imaginary axis, default 0)
    pub eta: f64,
}

impl Default for ACConfig {
    fn default() -> Self {
        Self {
            nw: 100,
            iw_cutoff: 5.0,
            x0: 0.5,
            eta: 0.0,
        }
    }
}

/// Generate scaled Legendre-Gauss quadrature for [0, infinity)
///
/// Maps GL nodes from [-1, 1] to [0, infinity) using the transformation:
///   omega = x0 * (1 + t) / (1 - t)
///   weights = w * 2 * x0 / (1 - t)^2
///
/// # Arguments
///
/// * `nw` - Number of frequency points
/// * `x0` - Scaling parameter (default 0.5)
///
/// # Returns
///
/// (frequencies, weights) both of length nw
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_ac.py::_get_scaled_legendre_roots()`
pub fn get_scaled_legendre_roots(nw: usize, x0: f64) -> (Array1<f64>, Array1<f64>) {
    // Get standard Gauss-Legendre nodes and weights on [-1, 1]
    let (gl_nodes, gl_weights) = gauss_legendre_nodes_weights(nw);

    // Transform to [0, infinity)
    let mut freqs = Array1::<f64>::zeros(nw);
    let mut wts = Array1::<f64>::zeros(nw);

    for i in 0..nw {
        let t = gl_nodes[i];
        let one_minus_t = 1.0 - t;

        // omega = x0 * (1 + t) / (1 - t)
        freqs[i] = x0 * (1.0 + t) / one_minus_t;

        // weight = w * 2 * x0 / (1 - t)^2
        wts[i] = gl_weights[i] * 2.0 * x0 / (one_minus_t * one_minus_t);
    }

    (freqs, wts)
}

/// Compute Gauss-Legendre nodes and weights on [-1, 1]
///
/// Uses the Golub-Welsch algorithm via eigenvalue decomposition
/// of the symmetric tridiagonal Jacobi matrix.
///
/// Now uses LAPACK DSYEV for reliable eigenvalue computation.
fn gauss_legendre_nodes_weights(n: usize) -> (Array1<f64>, Array1<f64>) {
    use ndarray_linalg::Eigh;

    if n == 0 {
        return (Array1::zeros(0), Array1::zeros(0));
    }

    if n == 1 {
        return (Array1::from_vec(vec![0.0]), Array1::from_vec(vec![2.0]));
    }

    // Build tridiagonal Jacobi matrix for Legendre polynomials
    // J is symmetric tridiagonal with zeros on diagonal
    // and beta_i = i / sqrt((2i-1)(2i+1)) on off-diagonals
    let mut jacobi = Array2::<f64>::zeros((n, n));

    for i in 1..n {
        let i_f64 = i as f64;
        let beta = i_f64 / ((2.0 * i_f64 - 1.0) * (2.0 * i_f64 + 1.0)).sqrt();
        jacobi[[i - 1, i]] = beta;
        jacobi[[i, i - 1]] = beta;
    }

    // Use LAPACK DSYEV for symmetric eigenvalue problem
    let (eigenvalues, eigenvectors) = jacobi
        .eigh(ndarray_linalg::UPLO::Upper)
        .expect("LAPACK DSYEV failed for Jacobi matrix");

    // Nodes are eigenvalues (already sorted by LAPACK)
    let nodes = eigenvalues;

    // Weights are 2 * (first component of eigenvector)^2
    let mut weights = Array1::<f64>::zeros(n);
    for i in 0..n {
        weights[i] = 2.0 * eigenvectors[[0, i]].powi(2);
    }

    (nodes, weights)
}

/// Compute polarizability P0(iw) in auxiliary basis
///
/// P0_PQ(iw) = 4 * sum_ia L_Pia * (e_ia / (w^2 + e_ia^2)) * L_Qia
///
/// # Arguments
///
/// * `omega` - Imaginary frequency (positive real value)
/// * `mo_energy` - MO energies [nmo]
/// * `lpq_ov` - DF tensor occupied-virtual block [naux, nocc, nvir]
/// * `nocc` - Number of occupied orbitals
///
/// # Returns
///
/// P0(iw) [naux, naux] (real symmetric matrix)
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_ac.py::get_rho_response()`
///
/// # Performance
///
/// Optimized using broadcast operations and BLAS-3 matrix multiply.
pub fn get_rho_response(omega: f64, mo_energy: &Array1<f64>, lpq_ov: &Array3<f64>) -> Array2<f64> {
    let (naux, nocc, nvir) = lpq_ov.dim();
    let omega_sq = omega * omega;

    // Compute chi_ia = e_ia / (omega^2 + e_ia^2) as 2D array [nocc, nvir]
    let chi = Array2::from_shape_fn((nocc, nvir), |(i, a)| {
        let e_ia = mo_energy[i] - mo_energy[nocc + a];
        e_ia / (omega_sq + e_ia * e_ia)
    });

    // Reshape Lpq_ov to [naux, nocc*nvir] for efficient operations
    let lpq_flat = lpq_ov
        .view()
        .into_shape_with_order((naux, nocc * nvir))
        .unwrap();

    // Flatten chi to [nocc*nvir]
    let chi_flat = chi.into_shape_with_order(nocc * nvir).unwrap();

    // Compute Pia_flat = Lpq_flat * chi_flat (broadcast multiply along axis 1)
    // Using ndarray broadcast multiplication
    let pia_flat = &lpq_flat * &chi_flat;

    // Pi = 4 * Pia @ Lpq^T using BLAS dgemm via ndarray
    let pi = 4.0 * pia_flat.dot(&lpq_flat.t());

    pi
}

/// Compute correlation self-energy on imaginary axis using AC method
///
/// This is the main entry point for computing Sigma_c(iw) matching PySCF exactly.
///
/// # Arguments
///
/// * `mo_energy` - MO energies [nmo]
/// * `lpq` - Full DF tensor [naux, nmo, nmo]
/// * `orbs` - Orbital indices to compute self-energy for
/// * `nocc` - Number of occupied orbitals (CRITICAL: must be passed from Python,
///            NOT inferred from eigenvalue signs which fails for DFT with negative LUMO)
/// * `config` - AC configuration
///
/// # Returns
///
/// * `sigma` - Sigma_c(iw) [norbs, nw_sigma] (complex)
/// * `omega` - Frequency points [norbs, nw_sigma] (complex)
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_ac.py::get_sigma_diag()`
///
/// # Note on nocc parameter
///
/// Previously this function inferred nocc by counting negative eigenvalues.
/// This is WRONG for DFT calculations where LUMO can have negative energy
/// (e.g., H2O/PBE has LUMO = -0.0009 Ha, CO/PBE has LUMO = -0.07 Ha).
/// The fix is to pass nocc explicitly from Python based on electron count.
pub fn get_sigma_diag(
    mo_energy: &Array1<f64>,
    lpq: &Array3<f64>,
    orbs: &[usize],
    nocc: usize,
    config: &ACConfig,
) -> Result<(Array2<Complex64>, Array2<Complex64>)> {
    let (naux, nmo, nmo2) = lpq.dim();
    if nmo != nmo2 {
        return Err(QuasixError::DimensionMismatch(format!(
            "Lpq not square in MO dimensions: {} x {}",
            nmo, nmo2
        )));
    }

    // Validate nocc parameter
    if nocc == 0 || nocc >= nmo {
        return Err(QuasixError::InvalidInput(format!(
            "Invalid nocc={}: must be in range [1, {})",
            nocc, nmo
        )));
    }

    let norbs = orbs.len();
    let _norbs_occ = orbs.iter().filter(|&&o| o < nocc).count();

    // Get frequency grid
    let (freqs, wts) = get_scaled_legendre_roots(config.nw, config.x0);

    // Determine frequency indices below cutoff (bounds-safe approach)
    // This fixes segfault for molecules like CO where many frequencies are below cutoff
    let freq_indices: Vec<usize> = freqs
        .iter()
        .enumerate()
        .filter(|(_, &f)| f < config.iw_cutoff)
        .map(|(idx, _)| idx)
        .collect();

    // nw_sigma = number of frequencies below cutoff + 1 (for omega=0)
    // Cap at config.nw to ensure we never exceed array bounds
    let nw_sigma = (freq_indices.len() + 1).min(config.nw + 1);

    // Fermi level at midpoint of HOMO-LUMO gap
    let ef = (mo_energy[nocc - 1] + mo_energy[nocc]) / 2.0;

    // Build frequency arrays for occupied and virtual orbitals
    // Occupied: omega = -i*freqs (negative imaginary)
    // Virtual:  omega = +i*freqs (positive imaginary)
    let mut omega_occ = Array1::<Complex64>::zeros(nw_sigma);
    let mut omega_vir = Array1::<Complex64>::zeros(nw_sigma);

    omega_occ[0] = Complex64::new(0.0, 0.0); // omega=0
    omega_vir[0] = Complex64::new(0.0, 0.0); // omega=0

    // Use the actual frequency indices to avoid out-of-bounds access
    for (k, &freq_idx) in freq_indices.iter().enumerate() {
        let k_sigma = k + 1; // +1 because k=0 is omega=0
        if k_sigma < nw_sigma {
            omega_occ[k_sigma] = Complex64::new(0.0, -freqs[freq_idx]); // -i*freq
            omega_vir[k_sigma] = Complex64::new(0.0, freqs[freq_idx]); // +i*freq
        }
    }

    // Compute emo = omega + ef - mo_energy for Green's function
    // Shape: [nmo, nw_sigma]
    let mut emo_occ = Array2::<Complex64>::zeros((nmo, nw_sigma));
    let mut emo_vir = Array2::<Complex64>::zeros((nmo, nw_sigma));

    for m in 0..nmo {
        for k in 0..nw_sigma {
            emo_occ[[m, k]] = omega_occ[k] + ef - mo_energy[m];
            emo_vir[[m, k]] = omega_vir[k] + ef - mo_energy[m];
        }
    }

    // Extract occupied-virtual block of Lpq
    let lpq_ov = lpq.slice(s![.., ..nocc, nocc..]).to_owned();

    // Allocate output
    let mut sigma = Array2::<Complex64>::zeros((norbs, nw_sigma));
    let mut omega_out = Array2::<Complex64>::zeros((norbs, nw_sigma));

    // Store omega values for each orbital
    for (idx, &orb) in orbs.iter().enumerate() {
        if orb < nocc {
            omega_out.row_mut(idx).assign(&omega_occ);
        } else {
            omega_out.row_mut(idx).assign(&omega_vir);
        }
    }

    // Pre-extract orbital slices for efficient access (done once, outside freq loop)
    // lpq_orbs[n_idx] = Lpq[:, orbs[n_idx], :] with shape [naux, nmo]
    // lpq_orbs_t[n_idx] = Lpq[:, :, orbs[n_idx]] with shape [naux, nmo]
    let lpq_orbs: Vec<Array2<f64>> = orbs
        .iter()
        .map(|&o| lpq.slice(s![.., o, ..]).to_owned())
        .collect();
    let lpq_orbs_t: Vec<Array2<f64>> = orbs
        .iter()
        .map(|&o| lpq.slice(s![.., .., o]).to_owned())
        .collect();

    // Parallelize over frequency points (100 iterations)
    // Each frequency's computation is independent
    let norm_re = -1.0 / PI;

    // PARALLELISM STRATEGY (2025-12-13): Hybrid Rayon + BLAS threading
    //
    // CRITICAL FIX: Previous approach used BlasThreadGuard::new() which set
    // OPENBLAS_NUM_THREADS=1, causing ALL BLAS operations to run single-threaded.
    // This resulted in ~1x speedup at 32 threads (3.1% efficiency)!
    //
    // New approach: Let each Rayon task use BLAS with its default threading.
    // OpenBLAS is thread-safe for independent matrix operations (different data).
    // Each frequency computation is completely independent - no shared state.
    //
    // For optimal performance on multi-socket systems:
    // - Set RAYON_NUM_THREADS = num_sockets (e.g., 2 for dual-socket)
    // - Let BLAS use cores_per_socket threads per task
    // - This minimizes NUMA cross-socket memory traffic
    //
    // Alternative: Use RAYON_NUM_THREADS = nw (frequency parallelism)
    // with OPENBLAS_NUM_THREADS = 1 for many small independent tasks.
    // Best when nw >> num_cores.

    let all_freq_contribs: Vec<Array2<Complex64>> = (0..config.nw)
        .into_par_iter()
        .map(|w| {
            let freq_w = freqs[w];
            let wt_w = wts[w];

            // Compute P0(iw)
            let pi = get_rho_response(freq_w, mo_energy, &lpq_ov);

            // Compute Pi_inv = (I - Pi)^{-1} - I
            let mut epsilon = Array2::<f64>::eye(naux);
            epsilon = &epsilon - &pi;

            // OPTIMIZATION (2025-12-03): Use Cholesky inversion for SPD matrices
            // On imaginary axis, epsilon = I - Pi is symmetric positive definite.
            // Cholesky inversion (LAPACK DPOTRI) is ~1.5-2x faster than LU (DGETRI).
            let epsilon_inv = match epsilon.invc() {
                Ok(inv) => inv,
                Err(chol_err) => {
                    // Cholesky failed - matrix may not be strictly SPD
                    // Fall back to general LU decomposition
                    match epsilon.inv() {
                        Ok(inv) => inv,
                        Err(lu_err) => {
                            eprintln!(
                                "WARNING: Matrix inversion failed at freq {}: Chol={:?}, LU={:?}",
                                w, chol_err, lu_err
                            );
                            // Return identity as fallback (no screening correction)
                            Array2::<f64>::eye(naux)
                        }
                    }
                }
            };
            let pi_inv = &epsilon_inv - &Array2::<f64>::eye(naux);
            let pi_inv_t = pi_inv.t();

            // Compute Green's function factors with weights
            let freq_w_sq = freq_w * freq_w;
            let freq_w_sq_c = Complex64::new(freq_w_sq, 0.0);
            let wt_w_c = Complex64::new(wt_w, 0.0);

            let g0_occ = emo_occ.mapv(|emo| wt_w_c * emo / (emo * emo + freq_w_sq_c));
            let g0_vir = emo_vir.mapv(|emo| wt_w_c * emo / (emo * emo + freq_w_sq_c));

            // Compute sigma contribution for all orbitals at this frequency
            let mut sigma_w = Array2::<Complex64>::zeros((norbs, nw_sigma));

            for (n_idx, &orb_n) in orbs.iter().enumerate() {
                // Contract: Q_nm = Pi_inv.T @ Lpq[:, orb_n, :]
                // PySCF: Qnm[Q,n,m] = sum_P Lpq[P,n,m] * Pi_inv[P,Q]
                let q_nm = pi_inv_t.dot(&lpq_orbs[n_idx]);

                // Contract: W_m = sum_Q Q_nm[Q,m] * Lpq[Q,m,orb_n]
                let wmn_col = (&q_nm * &lpq_orbs_t[n_idx]).sum_axis(Axis(0));

                // Select appropriate Green's function based on orbital occupancy
                let g0 = if orb_n < nocc { &g0_occ } else { &g0_vir };

                // Compute sigma contribution
                for k in 0..nw_sigma {
                    let mut sum_re = 0.0;
                    let mut sum_im = 0.0;
                    for m in 0..nmo {
                        let wmn_m = wmn_col[m];
                        let g0_mk = g0[[m, k]];
                        sum_re += wmn_m * g0_mk.re;
                        sum_im += wmn_m * g0_mk.im;
                    }
                    sigma_w[[n_idx, k]] = Complex64::new(norm_re * sum_re, norm_re * sum_im);
                }
            }
            sigma_w
        })
        .collect();

    // Accumulate contributions from all frequencies
    for sigma_w in all_freq_contribs {
        sigma = &sigma + &sigma_w;
    }

    Ok((sigma, omega_out))
}

/// Thiele's continued fraction coefficients via reciprocal differences
///
/// Given data points (z_n, f_n), compute coefficients a_n for the
/// continued fraction representation:
///
/// f(z) = a_0 + (z - z_0) / (a_1 + (z - z_1) / (a_2 + ...))
///
/// # Arguments
///
/// * `fn_vals` - Function values f(z_n) [nfit]
/// * `zn` - Sample points z_n [nfit]
///
/// # Returns
///
/// Continued fraction coefficients [nfit]
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_ac.py::thiele()`
pub fn thiele(fn_vals: &Array1<Complex64>, zn: &Array1<Complex64>) -> Array1<Complex64> {
    let nfit = fn_vals.len();
    if nfit == 0 {
        return Array1::zeros(0);
    }

    // Build reciprocal difference table
    // g[i,0] = f_i
    // g[i,j] = (z_i - z_{j-1}) / (g[i,j-1] - g[j-1,j-1]) + g[j-2,j-1]
    let mut g = Array2::<Complex64>::zeros((nfit, nfit));

    // Initialize first column with function values
    for i in 0..nfit {
        g[[i, 0]] = fn_vals[i];
    }

    // Compute reciprocal differences
    // PySCF formula: g[i:,i] = (g[i-1,i-1] - g[i:,i-1]) / ((zn[i:] - zn[i-1]) * g[i:,i-1])
    // Element-wise: g[i,j] = (g[j-1,j-1] - g[i,j-1]) / ((zn[i] - zn[j-1]) * g[i,j-1])
    for j in 1..nfit {
        for i in j..nfit {
            let g_ij_minus_1 = g[[i, j - 1]];
            let denom = (zn[i] - zn[j - 1]) * g_ij_minus_1;
            if denom.norm() < 1e-14 {
                // Avoid division by zero
                g[[i, j]] = Complex64::new(1e10, 0.0);
            } else {
                g[[i, j]] = (g[[j - 1, j - 1]] - g_ij_minus_1) / denom;
            }
        }
    }

    // Extract diagonal as continued fraction coefficients
    let mut a = Array1::<Complex64>::zeros(nfit);
    for i in 0..nfit {
        a[i] = g[[i, i]];
    }

    a
}

/// Evaluate Pade approximant using Thiele's continued fraction
///
/// Given continued fraction coefficients a_n and sample points z_n,
/// evaluate at frequency omega using backward recursion.
///
/// # Arguments
///
/// * `omega` - Frequency to evaluate at (real, relative to Fermi level)
/// * `zn` - Sample points used in fitting [nfit]
/// * `coeff` - Thiele coefficients [nfit]
///
/// # Returns
///
/// Sigma_c(omega) as complex value
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_ac.py::pade_thiele()`
#[inline]
pub fn pade_thiele(omega: f64, zn: &Array1<Complex64>, coeff: &Array1<Complex64>) -> Complex64 {
    let nfit = coeff.len();
    if nfit == 0 {
        return Complex64::new(0.0, 0.0);
    }

    if nfit == 1 {
        return coeff[0];
    }

    // Convert omega to complex
    let z = Complex64::new(omega, 0.0);

    // Backward recursion for continued fraction evaluation
    // X = a_{n-1} * (z - z_{n-2})
    // for i = n-2, ..., 0:
    //   X = a_i * (z - z_{i-1}) / (1 + X)
    // result = a_0 / (1 + X)

    let mut x = coeff[nfit - 1] * (z - zn[nfit - 2]);

    for i in (1..nfit - 1).rev() {
        let denom = Complex64::new(1.0, 0.0) + x;
        if denom.norm() < 1e-14 {
            x = Complex64::new(1e10, 0.0);
        } else {
            x = coeff[i] * (z - zn[i - 1]) / denom;
        }
    }

    let denom = Complex64::new(1.0, 0.0) + x;
    if denom.norm() < 1e-14 {
        coeff[0]
    } else {
        coeff[0] / denom
    }
}

/// Perform analytic continuation using Pade approximant
///
/// Select subset of frequency points and compute Thiele coefficients
/// for analytic continuation to real axis.
///
/// # Arguments
///
/// * `sigma` - Self-energy on imaginary axis [norbs, nw_sigma]
/// * `omega` - Frequency points [norbs, nw_sigma]
///
/// # Returns
///
/// (coefficients, omega_fit) for each orbital
/// - coefficients: [npade*2, norbs]
/// - omega_fit: [norbs, npade*2]
///
/// # PySCF Reference
///
/// `pyscf/gw/gw_ac.py::AC_pade_thiele_diag()`
pub fn ac_pade_thiele_diag(
    sigma: &Array2<Complex64>,
    omega: &Array2<Complex64>,
) -> (Array2<Complex64>, Array2<Complex64>) {
    let (norbs, nw) = sigma.dim();

    // PySCF point selection:
    // idx = range(1, 40, 6)  -> [1, 7, 13, 19, 25, 31, 37]
    // sigma1 = sigma[:, idx]
    // sigma2 = sigma[:, (idx[-1]+4)::4]  -> [41, 45, 49, ...]
    // sigma = hstack((sigma1, sigma2))

    let mut idx1: Vec<usize> = Vec::new();
    let mut i = 1;
    while i < 40.min(nw) {
        idx1.push(i);
        i += 6;
    }

    let mut idx2: Vec<usize> = Vec::new();
    if !idx1.is_empty() {
        let start = idx1[idx1.len() - 1] + 4;
        let mut j = start;
        while j < nw {
            idx2.push(j);
            j += 4;
        }
    }

    // Combine indices
    let mut all_idx: Vec<usize> = idx1;
    all_idx.extend(idx2);

    // Limit to even number for Pade
    let npade = all_idx.len() / 2;
    let nfit = npade * 2;

    if nfit == 0 {
        // Not enough points for Pade
        return (Array2::zeros((0, norbs)), Array2::zeros((norbs, 0)));
    }

    // Extract selected points
    let mut sigma_fit = Array2::<Complex64>::zeros((norbs, nfit));
    let mut omega_fit = Array2::<Complex64>::zeros((norbs, nfit));

    for (k, &idx) in all_idx.iter().take(nfit).enumerate() {
        for p in 0..norbs {
            sigma_fit[[p, k]] = sigma[[p, idx]];
            omega_fit[[p, k]] = omega[[p, idx]];
        }
    }

    // Compute Thiele coefficients for each orbital
    let mut coeff = Array2::<Complex64>::zeros((nfit, norbs));

    for p in 0..norbs {
        let sigma_p = sigma_fit.row(p).to_owned();
        let omega_p = omega_fit.row(p).to_owned();
        let coeff_p = thiele(&sigma_p, &omega_p);

        for k in 0..nfit {
            coeff[[k, p]] = coeff_p[k];
        }
    }

    (coeff, omega_fit)
}

/// Compute linearized QP energy correction
///
/// E_n^QP = e_n + Z_n * (Sigma_x_nn + Re[Sigma_c(e_n)] - V_xc_nn)
///
/// where Z_n = 1 / (1 - dSigma/dw |_{w=e_n})
///
/// # Arguments
///
/// * `mo_energy_hf` - HF orbital energies [nmo]
/// * `sigma_x_diag` - Exchange self-energy diagonal [nmo]
/// * `v_xc_diag` - DFT XC potential diagonal [nmo]
/// * `sigma_coeff` - Pade coefficients [nfit, norbs]
/// * `omega_fit` - Fitted frequency points [norbs, nfit]
/// * `orbs` - Orbital indices
/// * `nocc` - Number of occupied orbitals
///
/// # Returns
///
/// QP energies [norbs]
pub fn compute_qp_energies_linearized(
    mo_energy_hf: &Array1<f64>,
    sigma_x_diag: &Array1<f64>,
    v_xc_diag: &Array1<f64>,
    sigma_coeff: &Array2<Complex64>,
    omega_fit: &Array2<Complex64>,
    orbs: &[usize],
    nocc: usize,
) -> Array1<f64> {
    let norbs = orbs.len();
    let ef = (mo_energy_hf[nocc - 1] + mo_energy_hf[nocc]) / 2.0;

    let mut qp_energy = Array1::<f64>::zeros(norbs);

    let de = 1e-6; // Numerical derivative step

    for (idx, &p) in orbs.iter().enumerate() {
        let ep = mo_energy_hf[p];

        // Evaluate Sigma_c at ep (relative to Fermi level)
        let omega_p = omega_fit.row(idx).to_owned();
        let coeff_p = sigma_coeff.column(idx).to_owned();

        let sigma_r = pade_thiele(ep - ef, &omega_p, &coeff_p).re;

        // Numerical derivative for Z-factor
        let sigma_r_plus = pade_thiele(ep - ef + de, &omega_p, &coeff_p).re;
        let dsigma = (sigma_r_plus - sigma_r) / de;

        // Z-factor
        let zn = 1.0 / (1.0 - dsigma);

        // QP correction
        let correction = sigma_x_diag[p] + sigma_r - v_xc_diag[p];

        qp_energy[idx] = ep + zn * correction;
    }

    qp_energy
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_scaled_legendre_roots() {
        // Test that quadrature integrates polynomial correctly
        let (freqs, wts) = get_scaled_legendre_roots(20, 0.5);

        assert_eq!(freqs.len(), 20);
        assert_eq!(wts.len(), 20);

        // All frequencies should be positive
        assert!(freqs.iter().all(|&f| f >= 0.0));

        // All weights should be positive
        assert!(wts.iter().all(|&w| w > 0.0));

        // Test that frequencies are increasing
        for i in 1..freqs.len() {
            assert!(freqs[i] > freqs[i - 1]);
        }
    }

    #[test]
    fn test_gauss_legendre_nodes() {
        // Test 2-point Gauss-Legendre
        let (nodes, weights) = gauss_legendre_nodes_weights(2);

        // Known values: nodes = +/- 1/sqrt(3), weights = 1
        let expected_node = 1.0 / 3.0_f64.sqrt();
        assert_abs_diff_eq!(nodes[0], -expected_node, epsilon = 1e-10);
        assert_abs_diff_eq!(nodes[1], expected_node, epsilon = 1e-10);
        assert_abs_diff_eq!(weights[0], 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(weights[1], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_thiele_simple() {
        // Test Thiele on simple polynomial f(z) = 1
        let zn = Array1::from_vec(vec![
            Complex64::new(0.0, 0.1),
            Complex64::new(0.0, 0.2),
            Complex64::new(0.0, 0.3),
        ]);
        let fn_vals = Array1::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
        ]);

        let coeff = thiele(&fn_vals, &zn);

        // For constant function, a_0 = 1, rest should be ~0 or very large
        assert_abs_diff_eq!(coeff[0].re, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_pade_thiele_evaluation() {
        // Test Pade approximant for analytic continuation
        // pade_thiele evaluates at REAL frequencies (omega) using fit on imaginary axis
        let zn = Array1::from_vec(vec![
            Complex64::new(0.0, 0.5),
            Complex64::new(0.0, 1.0),
            Complex64::new(0.0, 2.0),
            Complex64::new(0.0, 4.0),
        ]);

        // Simple function: f(z) = 1/(1+z)
        let fn_vals = zn.mapv(|z| Complex64::new(1.0, 0.0) / (Complex64::new(1.0, 0.0) + z));

        let coeff = thiele(&fn_vals, &zn);

        // Test at REAL frequencies - this is what pade_thiele is designed for
        // At omega = 0.1 (real), expected = 1/(1+0.1) = 0.909...
        let omega_test = 0.1;
        let z_real = Complex64::new(omega_test, 0.0);
        let expected = Complex64::new(1.0, 0.0) / (Complex64::new(1.0, 0.0) + z_real);
        let result = pade_thiele(omega_test, &zn, &coeff);

        // Pade should recover the analytic function well for small real omega
        // when fitted on imaginary axis
        assert_abs_diff_eq!(result.re, expected.re, epsilon = 0.01);
        assert!(
            result.im.abs() < 0.1,
            "Imaginary part should be small: {}",
            result.im
        );

        // Also test at omega = 0 where all terms should contribute equally
        let result_zero = pade_thiele(0.0, &zn, &coeff);
        let expected_zero = Complex64::new(1.0, 0.0); // f(0) = 1/(1+0) = 1
        assert_abs_diff_eq!(result_zero.re, expected_zero.re, epsilon = 0.01);
    }

    #[test]
    fn test_rho_response_shape() {
        // Test that polarizability has correct shape
        let naux = 5;
        let nocc = 2;
        let nvir = 3;
        let nmo = nocc + nvir;

        let lpq_ov = Array3::<f64>::from_elem((naux, nocc, nvir), 0.1);
        let mut mo_energy = Array1::<f64>::zeros(nmo);
        mo_energy[0] = -0.5;
        mo_energy[1] = -0.3;
        mo_energy[2] = 0.2;
        mo_energy[3] = 0.4;
        mo_energy[4] = 0.6;

        let pi = get_rho_response(0.5, &mo_energy, &lpq_ov);

        assert_eq!(pi.dim(), (naux, naux));

        // Should be symmetric
        for i in 0..naux {
            for j in 0..naux {
                assert_abs_diff_eq!(pi[[i, j]], pi[[j, i]], epsilon = 1e-12);
            }
        }
    }
}
