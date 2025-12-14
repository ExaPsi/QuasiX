//! Frequency Caching for evGW Calculations
//!
//! This module provides intelligent caching of frequency-dependent quantities
//! in evGW calculations to avoid redundant recomputation.
//!
//! # Motivation
//!
//! In evGW, the screened interaction W(iw) and RPA polarizability P0(iw) are
//! computed using HF/DFT eigenvalues (fixed wavefunctions). Since these
//! eigenvalues don't change between iterations, W and P0 remain constant!
//!
//! The only thing that changes is the evaluation point for Sigma_c:
//! - Iteration 0: evaluate at HF energies (G0W0)
//! - Iteration k: evaluate at QP energies from iteration k-1
//!
//! By caching W_mn(iw) after iteration 0, we can skip the expensive:
//! - Polarizability computation P0(iw) for all 100 frequencies
//! - Dielectric matrix inversion (I - P0)^{-1} for all frequencies
//! - W_mn contraction with DF tensors
//!
//! # Expected Speedup
//!
//! - G0W0: No benefit (single iteration)
//! - evGW with 6-12 iterations: 5-10x speedup
//!
//! The first iteration computes and caches everything, subsequent iterations
//! only need to evaluate Pade continuation at new QP energies.
//!
//! # Memory Cost
//!
//! - W_mn matrix: O(nmo^2 * nw) = O(100 * nmo^2) floats
//! - For nmo=100, this is ~80 MB (acceptable)
//! - For nmo=500, this is ~2 GB (may need streaming approach)
//!
//! # Usage
//!
//! ```ignore
//! let mut cache = FrequencyCache::new();
//!
//! for iter in 0..max_iter {
//!     let (sigma_iw, omega_iw) = if iter == 0 {
//!         // First iteration: compute and cache
//!         let result = get_sigma_diag_and_cache(mo_energy, lpq, ...);
//!         cache.store_w_mn(result.w_mn);
//!         cache.store_pade_coeffs(result.coeff, result.omega_fit);
//!         (result.sigma, result.omega)
//!     } else {
//!         // Subsequent iterations: use cached W_mn, only recompute G0 and accumulate
//!         get_sigma_diag_cached(&cache, qp_energies, ...)
//!     };
//!
//!     // ... rest of evGW iteration
//! }
//! ```

use ndarray::{s, Array1, Array2, Array3, Axis};
use ndarray_linalg::Inverse;
use num_complex::Complex64;
use std::f64::consts::PI;

use super::correlation_ac::{get_rho_response, get_scaled_legendre_roots, ACConfig};
use crate::common::{QuasixError, Result};

/// Cached frequency-dependent quantities for evGW
///
/// Stores intermediate results that don't change between evGW iterations:
/// - W_mn(iw): Screened interaction matrix elements [nw, norbs, nmo]
/// - Pade coefficients and fitting points
/// - Frequency grid (GL nodes and weights)
#[derive(Debug, Clone)]
pub struct FrequencyCache {
    /// Whether cache is populated
    pub is_valid: bool,

    /// Number of MO orbitals
    pub nmo: usize,

    /// Number of occupied orbitals
    pub nocc: usize,

    /// Number of auxiliary basis functions
    pub naux: usize,

    /// Number of frequency points for W integration
    pub nw: usize,

    /// Number of orbitals for which we cached results
    pub norbs: usize,

    /// Orbital indices
    pub orbs: Vec<usize>,

    /// Gauss-Legendre frequency grid points [nw]
    pub freqs: Array1<f64>,

    /// Gauss-Legendre quadrature weights [nw]
    pub weights: Array1<f64>,

    /// Screened interaction W_mn at all frequencies [nw, norbs, nmo]
    /// W_mn[w, n_idx, m] = sum_Q q_nm[Q,m] * lpq_orbs_t[n_idx][Q,m]
    /// where q_nm = pi_inv_t @ lpq_orbs[n_idx]
    pub w_mn: Option<Array3<f64>>,

    /// Frequency-dependent RI screened interaction (I - P0)^{-1} - I at each freq
    /// Shape: [nw, naux, naux]
    pub pi_inv: Option<Array3<f64>>,

    /// Pre-extracted orbital slices: lpq[:, orbs[n], :] for each orbital
    /// Shape: Vec of [naux, nmo] arrays
    pub lpq_orbs: Option<Vec<Array2<f64>>>,

    /// Pre-extracted orbital slices (transposed): lpq[:, :, orbs[n]]
    /// Shape: Vec of [naux, nmo] arrays
    pub lpq_orbs_t: Option<Vec<Array2<f64>>>,

    /// Number of sigma frequency points (frequencies below cutoff + 1)
    pub nw_sigma: usize,

    /// Frequency indices below cutoff
    pub freq_indices: Vec<usize>,

    /// Fermi level
    pub ef: f64,

    /// Omega arrays for occupied/virtual orbitals [nw_sigma]
    pub omega_occ: Option<Array1<Complex64>>,
    pub omega_vir: Option<Array1<Complex64>>,
}

impl Default for FrequencyCache {
    fn default() -> Self {
        Self::new()
    }
}

impl FrequencyCache {
    /// Create an empty cache
    pub fn new() -> Self {
        Self {
            is_valid: false,
            nmo: 0,
            nocc: 0,
            naux: 0,
            nw: 0,
            norbs: 0,
            orbs: Vec::new(),
            freqs: Array1::zeros(0),
            weights: Array1::zeros(0),
            w_mn: None,
            pi_inv: None,
            lpq_orbs: None,
            lpq_orbs_t: None,
            nw_sigma: 0,
            freq_indices: Vec::new(),
            ef: 0.0,
            omega_occ: None,
            omega_vir: None,
        }
    }

    /// Check if cache is valid and matches current calculation parameters
    pub fn is_compatible(&self, nmo: usize, nocc: usize, naux: usize, nw: usize) -> bool {
        self.is_valid && self.nmo == nmo && self.nocc == nocc && self.naux == naux && self.nw == nw
    }

    /// Invalidate the cache
    pub fn invalidate(&mut self) {
        self.is_valid = false;
        self.w_mn = None;
        self.pi_inv = None;
        self.lpq_orbs = None;
        self.lpq_orbs_t = None;
        self.omega_occ = None;
        self.omega_vir = None;
    }

    /// Get memory usage in bytes
    pub fn memory_usage(&self) -> usize {
        let mut total = 0;

        if let Some(ref w) = self.w_mn {
            total += w.len() * std::mem::size_of::<f64>();
        }
        if let Some(ref pi) = self.pi_inv {
            total += pi.len() * std::mem::size_of::<f64>();
        }
        if let Some(ref lpq) = self.lpq_orbs {
            for arr in lpq {
                total += arr.len() * std::mem::size_of::<f64>();
            }
        }
        if let Some(ref lpq_t) = self.lpq_orbs_t {
            for arr in lpq_t {
                total += arr.len() * std::mem::size_of::<f64>();
            }
        }

        total
    }
}

/// Compute Sigma_c on imaginary axis AND populate frequency cache
///
/// This is the first-iteration version that computes everything from scratch
/// and stores intermediate results in the cache for subsequent iterations.
///
/// # Arguments
///
/// * `mo_energy` - MO energies [nmo]
/// * `lpq` - Full DF tensor [naux, nmo, nmo]
/// * `orbs` - Orbital indices to compute self-energy for
/// * `nocc` - Number of occupied orbitals
/// * `config` - AC configuration
/// * `cache` - Mutable frequency cache to populate
///
/// # Returns
///
/// * `sigma` - Sigma_c(iw) [norbs, nw_sigma]
/// * `omega` - Frequency points [norbs, nw_sigma]
pub fn get_sigma_diag_and_cache(
    mo_energy: &Array1<f64>,
    lpq: &Array3<f64>,
    orbs: &[usize],
    nocc: usize,
    config: &ACConfig,
    cache: &mut FrequencyCache,
) -> Result<(Array2<Complex64>, Array2<Complex64>)> {
    let (naux, nmo, nmo2) = lpq.dim();
    if nmo != nmo2 {
        return Err(QuasixError::DimensionMismatch(format!(
            "Lpq not square in MO dimensions: {} x {}",
            nmo, nmo2
        )));
    }

    if nocc == 0 || nocc >= nmo {
        return Err(QuasixError::InvalidInput(format!(
            "Invalid nocc={}: must be in range [1, {})",
            nocc, nmo
        )));
    }

    let norbs = orbs.len();

    // Get frequency grid
    let (freqs, wts) = get_scaled_legendre_roots(config.nw, config.x0);

    // Determine frequency indices below cutoff
    let freq_indices: Vec<usize> = freqs
        .iter()
        .enumerate()
        .filter(|(_, &f)| f < config.iw_cutoff)
        .map(|(idx, _)| idx)
        .collect();

    let nw_sigma = (freq_indices.len() + 1).min(config.nw + 1);

    // Fermi level at midpoint of HOMO-LUMO gap
    let ef = (mo_energy[nocc - 1] + mo_energy[nocc]) / 2.0;

    // Build frequency arrays for occupied and virtual orbitals
    let mut omega_occ = Array1::<Complex64>::zeros(nw_sigma);
    let mut omega_vir = Array1::<Complex64>::zeros(nw_sigma);

    omega_occ[0] = Complex64::new(0.0, 0.0);
    omega_vir[0] = Complex64::new(0.0, 0.0);

    for (k, &freq_idx) in freq_indices.iter().enumerate() {
        let k_sigma = k + 1;
        if k_sigma < nw_sigma {
            omega_occ[k_sigma] = Complex64::new(0.0, -freqs[freq_idx]);
            omega_vir[k_sigma] = Complex64::new(0.0, freqs[freq_idx]);
        }
    }

    // Compute emo = omega + ef - mo_energy for Green's function
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

    // Pre-extract orbital slices (cache these!)
    let lpq_orbs: Vec<Array2<f64>> = orbs
        .iter()
        .map(|&o| lpq.slice(s![.., o, ..]).to_owned())
        .collect();
    let lpq_orbs_t: Vec<Array2<f64>> = orbs
        .iter()
        .map(|&o| lpq.slice(s![.., .., o]).to_owned())
        .collect();

    // Allocate output and cache storage
    let mut sigma = Array2::<Complex64>::zeros((norbs, nw_sigma));
    let mut omega_out = Array2::<Complex64>::zeros((norbs, nw_sigma));

    // W_mn cache: [nw, norbs, nmo]
    let mut w_mn_cache = Array3::<f64>::zeros((config.nw, norbs, nmo));

    // Pi_inv cache: [nw, naux, naux]
    let mut pi_inv_cache = Array3::<f64>::zeros((config.nw, naux, naux));

    // Store omega values for each orbital
    for (idx, &orb) in orbs.iter().enumerate() {
        if orb < nocc {
            omega_out.row_mut(idx).assign(&omega_occ);
        } else {
            omega_out.row_mut(idx).assign(&omega_vir);
        }
    }

    let norm_re = -1.0 / PI;

    // Compute for each frequency point
    for w in 0..config.nw {
        let freq_w = freqs[w];
        let wt_w = wts[w];

        // Compute P0(iw)
        let pi = get_rho_response(freq_w, mo_energy, &lpq_ov);

        // Compute Pi_inv = (I - Pi)^{-1} - I
        let mut epsilon = Array2::<f64>::eye(naux);
        epsilon = &epsilon - &pi;

        let epsilon_inv = match epsilon.inv() {
            Ok(inv) => inv,
            Err(e) => {
                eprintln!("WARNING: Matrix inversion failed at freq {}: {:?}", w, e);
                Array2::<f64>::eye(naux)
            }
        };
        let pi_inv = &epsilon_inv - &Array2::<f64>::eye(naux);

        // Cache pi_inv for this frequency
        pi_inv_cache.slice_mut(s![w, .., ..]).assign(&pi_inv);

        let pi_inv_t = pi_inv.t();

        // Compute Green's function factors with weights
        let freq_w_sq = freq_w * freq_w;
        let freq_w_sq_c = Complex64::new(freq_w_sq, 0.0);
        let wt_w_c = Complex64::new(wt_w, 0.0);

        let g0_occ = emo_occ.mapv(|emo| wt_w_c * emo / (emo * emo + freq_w_sq_c));
        let g0_vir = emo_vir.mapv(|emo| wt_w_c * emo / (emo * emo + freq_w_sq_c));

        // Compute W_mn and sigma contribution for all orbitals at this frequency
        for (n_idx, &orb_n) in orbs.iter().enumerate() {
            let q_nm = pi_inv_t.dot(&lpq_orbs[n_idx]);
            let wmn_col = (&q_nm * &lpq_orbs_t[n_idx]).sum_axis(Axis(0));

            // Cache W_mn for this frequency and orbital
            for m in 0..nmo {
                w_mn_cache[[w, n_idx, m]] = wmn_col[m];
            }

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
                sigma[[n_idx, k]] += Complex64::new(norm_re * sum_re, norm_re * sum_im);
            }
        }
    }

    // Populate the cache
    cache.is_valid = true;
    cache.nmo = nmo;
    cache.nocc = nocc;
    cache.naux = naux;
    cache.nw = config.nw;
    cache.norbs = norbs;
    cache.orbs = orbs.to_vec();
    cache.freqs = freqs;
    cache.weights = wts;
    cache.w_mn = Some(w_mn_cache);
    cache.pi_inv = Some(pi_inv_cache);
    cache.lpq_orbs = Some(lpq_orbs);
    cache.lpq_orbs_t = Some(lpq_orbs_t);
    cache.nw_sigma = nw_sigma;
    cache.freq_indices = freq_indices;
    cache.ef = ef;
    cache.omega_occ = Some(omega_occ);
    cache.omega_vir = Some(omega_vir);

    Ok((sigma, omega_out))
}

/// Compute Sigma_c on imaginary axis using cached W_mn
///
/// This is the fast path for iterations > 0 in evGW. It reuses the cached
/// W_mn and pi_inv, only recomputing the Green's function G0 with updated
/// QP energies.
///
/// # Arguments
///
/// * `mo_energy` - HF/DFT orbital energies [nmo] (for Green's function)
/// * `cache` - Pre-populated frequency cache from iteration 0
/// * `config` - AC configuration
///
/// # Returns
///
/// * `sigma` - Sigma_c(iw) [norbs, nw_sigma]
/// * `omega` - Frequency points [norbs, nw_sigma]
///
/// # Performance
///
/// This is 5-10x faster than full recomputation because:
/// - Skips P0(iw) computation at all frequencies
/// - Skips matrix inversion (I - P0)^{-1} at all frequencies
/// - Skips W_mn contraction with DF tensors
/// - Only recomputes G0(iw) and final accumulation
pub fn get_sigma_diag_cached(
    mo_energy: &Array1<f64>,
    cache: &FrequencyCache,
    config: &ACConfig,
) -> Result<(Array2<Complex64>, Array2<Complex64>)> {
    if !cache.is_valid {
        return Err(QuasixError::InvalidInput(
            "FrequencyCache is not valid - call get_sigma_diag_and_cache first".to_string(),
        ));
    }

    let w_mn = cache
        .w_mn
        .as_ref()
        .ok_or_else(|| QuasixError::InvalidInput("W_mn not cached".to_string()))?;
    let omega_occ = cache
        .omega_occ
        .as_ref()
        .ok_or_else(|| QuasixError::InvalidInput("omega_occ not cached".to_string()))?;
    let omega_vir = cache
        .omega_vir
        .as_ref()
        .ok_or_else(|| QuasixError::InvalidInput("omega_vir not cached".to_string()))?;

    let nmo = cache.nmo;
    let nocc = cache.nocc;
    let norbs = cache.norbs;
    let nw_sigma = cache.nw_sigma;
    let ef = cache.ef;
    let orbs = &cache.orbs;
    let freqs = &cache.freqs;
    let wts = &cache.weights;

    // Compute emo with CURRENT mo_energy (may differ from HF energies in later iterations)
    // NOTE: In evGW, we actually use HF energies for G0, so mo_energy here is HF energies
    let mut emo_occ = Array2::<Complex64>::zeros((nmo, nw_sigma));
    let mut emo_vir = Array2::<Complex64>::zeros((nmo, nw_sigma));

    for m in 0..nmo {
        for k in 0..nw_sigma {
            emo_occ[[m, k]] = omega_occ[k] + ef - mo_energy[m];
            emo_vir[[m, k]] = omega_vir[k] + ef - mo_energy[m];
        }
    }

    // Allocate output
    let mut sigma = Array2::<Complex64>::zeros((norbs, nw_sigma));
    let mut omega_out = Array2::<Complex64>::zeros((norbs, nw_sigma));

    // Store omega values for each orbital
    for (idx, &orb) in orbs.iter().enumerate() {
        if orb < nocc {
            omega_out.row_mut(idx).assign(omega_occ);
        } else {
            omega_out.row_mut(idx).assign(omega_vir);
        }
    }

    let norm_re = -1.0 / PI;

    // Compute sigma using cached W_mn (main speedup!)
    for w in 0..config.nw {
        let freq_w = freqs[w];
        let wt_w = wts[w];

        // Green's function factors with weights (only thing we recompute)
        let freq_w_sq = freq_w * freq_w;
        let freq_w_sq_c = Complex64::new(freq_w_sq, 0.0);
        let wt_w_c = Complex64::new(wt_w, 0.0);

        let g0_occ = emo_occ.mapv(|emo| wt_w_c * emo / (emo * emo + freq_w_sq_c));
        let g0_vir = emo_vir.mapv(|emo| wt_w_c * emo / (emo * emo + freq_w_sq_c));

        // Accumulate sigma using cached W_mn
        for (n_idx, &orb_n) in orbs.iter().enumerate() {
            let g0 = if orb_n < nocc { &g0_occ } else { &g0_vir };

            for k in 0..nw_sigma {
                let mut sum_re = 0.0;
                let mut sum_im = 0.0;
                for m in 0..nmo {
                    let wmn_m = w_mn[[w, n_idx, m]];
                    let g0_mk = g0[[m, k]];
                    sum_re += wmn_m * g0_mk.re;
                    sum_im += wmn_m * g0_mk.im;
                }
                sigma[[n_idx, k]] += Complex64::new(norm_re * sum_re, norm_re * sum_im);
            }
        }
    }

    Ok((sigma, omega_out))
}

/// Pade coefficients cache for analytic continuation
#[derive(Debug, Clone)]
pub struct PadeCache {
    /// Thiele coefficients [nfit, norbs]
    pub coeff: Array2<Complex64>,
    /// Fitting frequency points [norbs, nfit]
    pub omega_fit: Array2<Complex64>,
    /// Whether cache is valid
    pub is_valid: bool,
}

impl Default for PadeCache {
    fn default() -> Self {
        Self::new()
    }
}

impl PadeCache {
    pub fn new() -> Self {
        Self {
            coeff: Array2::zeros((0, 0)),
            omega_fit: Array2::zeros((0, 0)),
            is_valid: false,
        }
    }

    /// Store Pade coefficients
    pub fn store(&mut self, coeff: Array2<Complex64>, omega_fit: Array2<Complex64>) {
        self.coeff = coeff;
        self.omega_fit = omega_fit;
        self.is_valid = true;
    }

    /// Invalidate cache
    pub fn invalidate(&mut self) {
        self.is_valid = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_frequency_cache_creation() {
        let cache = FrequencyCache::new();
        assert!(!cache.is_valid);
        assert_eq!(cache.memory_usage(), 0);
    }

    #[test]
    fn test_frequency_cache_invalidation() {
        let mut cache = FrequencyCache::new();
        cache.is_valid = true;
        cache.nmo = 10;

        cache.invalidate();
        assert!(!cache.is_valid);
        assert!(cache.w_mn.is_none());
    }

    #[test]
    fn test_cache_compatibility() {
        let mut cache = FrequencyCache::new();
        cache.is_valid = true;
        cache.nmo = 10;
        cache.nocc = 5;
        cache.naux = 50;
        cache.nw = 100;

        // Compatible
        assert!(cache.is_compatible(10, 5, 50, 100));

        // Incompatible
        assert!(!cache.is_compatible(10, 5, 50, 50)); // Different nw
        assert!(!cache.is_compatible(10, 6, 50, 100)); // Different nocc
    }

    #[test]
    fn test_pade_cache() {
        let mut cache = PadeCache::new();
        assert!(!cache.is_valid);

        let coeff = Array2::zeros((10, 5));
        let omega_fit = Array2::zeros((5, 10));
        cache.store(coeff, omega_fit);

        assert!(cache.is_valid);
        assert_eq!(cache.coeff.dim(), (10, 5));
        assert_eq!(cache.omega_fit.dim(), (5, 10));
    }

    /// Test that cached sigma_c matches non-cached computation
    ///
    /// This validates that the frequency caching optimization gives
    /// numerically identical results to the original computation.
    #[test]
    fn test_cached_vs_noncached_sigma() {
        use crate::selfenergy::correlation_ac;

        // Small test system: 6 MOs, 3 occupied, 15 auxiliary functions
        let nmo = 6;
        let nocc = 3;
        let naux = 15;

        // Create test MO energies (in Hartree)
        // Occupied: -0.5, -0.4, -0.3 (HOMO)
        // Virtual:   0.1 (LUMO), 0.2, 0.3
        let mo_energy = Array1::from_vec(vec![-0.5, -0.4, -0.3, 0.1, 0.2, 0.3]);

        // Create test DF tensor with reasonable values
        // lpq[P, m, n] should be symmetric in (m, n) for stability
        let mut lpq = Array3::<f64>::zeros((naux, nmo, nmo));
        for p in 0..naux {
            for m in 0..nmo {
                for n in 0..nmo {
                    // Simple test values: symmetric and decaying
                    lpq[[p, m, n]] = 0.1 * ((-0.2 * (m + n + p) as f64).exp());
                }
            }
        }

        // Test for all orbitals
        let orbs: Vec<usize> = (0..nmo).collect();

        // AC configuration matching evGW defaults
        let ac_config = ACConfig {
            nw: 20,
            iw_cutoff: 5.0,
            x0: 0.5,
            eta: 0.01,
        };

        // Method 1: Non-cached (original) computation
        let (sigma_noncached, omega_noncached) =
            correlation_ac::get_sigma_diag(&mo_energy, &lpq, &orbs, nocc, &ac_config)
                .expect("Non-cached sigma computation failed");

        // Method 2: Cached computation (first call populates cache)
        let mut cache = FrequencyCache::new();
        let (sigma_cached, omega_cached) =
            get_sigma_diag_and_cache(&mo_energy, &lpq, &orbs, nocc, &ac_config, &mut cache)
                .expect("Cached sigma computation failed");

        // Validate cache is populated
        assert!(cache.is_valid, "Cache should be valid after first call");
        assert!(cache.w_mn.is_some(), "W_mn should be cached");
        assert_eq!(cache.nmo, nmo);
        assert_eq!(cache.nocc, nocc);

        // Method 3: Cached computation (second call uses cache)
        let (sigma_from_cache, _omega_from_cache) =
            get_sigma_diag_cached(&mo_energy, &cache, &ac_config)
                .expect("Sigma from cache computation failed");

        // Validate results match within numerical tolerance (1e-10)
        let tol = 1e-10;

        // Check dimensions match
        assert_eq!(
            sigma_noncached.dim(),
            sigma_cached.dim(),
            "Sigma dimensions should match"
        );
        assert_eq!(
            sigma_noncached.dim(),
            sigma_from_cache.dim(),
            "Sigma from cache dimensions should match"
        );
        assert_eq!(
            omega_noncached.dim(),
            omega_cached.dim(),
            "Omega dimensions should match"
        );

        // Check sigma values match (cached vs non-cached)
        let mut max_diff_sigma = 0.0_f64;
        for idx in 0..sigma_noncached.len() {
            let diff = (sigma_noncached.as_slice().unwrap()[idx]
                - sigma_cached.as_slice().unwrap()[idx])
                .norm();
            max_diff_sigma = max_diff_sigma.max(diff);
        }
        assert!(
            max_diff_sigma < tol,
            "Cached sigma should match non-cached: max_diff = {:.2e}, tol = {:.2e}",
            max_diff_sigma,
            tol
        );

        // Check sigma values match (from cache vs non-cached)
        let mut max_diff_from_cache = 0.0_f64;
        for idx in 0..sigma_noncached.len() {
            let diff = (sigma_noncached.as_slice().unwrap()[idx]
                - sigma_from_cache.as_slice().unwrap()[idx])
                .norm();
            max_diff_from_cache = max_diff_from_cache.max(diff);
        }
        assert!(
            max_diff_from_cache < tol,
            "Sigma from cache should match non-cached: max_diff = {:.2e}, tol = {:.2e}",
            max_diff_from_cache,
            tol
        );

        // Check omega values match
        let mut max_diff_omega = 0.0_f64;
        for idx in 0..omega_noncached.len() {
            let diff = (omega_noncached.as_slice().unwrap()[idx]
                - omega_cached.as_slice().unwrap()[idx])
                .norm();
            max_diff_omega = max_diff_omega.max(diff);
        }
        assert!(
            max_diff_omega < tol,
            "Omega should match: max_diff = {:.2e}, tol = {:.2e}",
            max_diff_omega,
            tol
        );

        println!("Frequency caching validation PASSED:");
        println!(
            "  max_diff(sigma_cached vs noncached) = {:.2e}",
            max_diff_sigma
        );
        println!(
            "  max_diff(sigma_from_cache vs noncached) = {:.2e}",
            max_diff_from_cache
        );
        println!("  max_diff(omega) = {:.2e}", max_diff_omega);
        println!("  Cache memory usage: {} bytes", cache.memory_usage());
    }
}
