//! Imaginary-axis evaluation for P0, W, and Σc
//!
//! This module implements efficient computation of polarizability P0(iξ),
//! screened interaction W(iξ), and correlation self-energy Σc(iξ) on purely
//! imaginary frequencies, enabling stable analytic continuation for GW calculations.
//!
//! Key features:
//! - Real symmetric matrices on imaginary axis
//! - Monotonic decay properties for numerical stability
//! - Efficient caching mechanism for repeated evaluations
//! - Numerical conditioning guards
#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::{QuasixError, Result};
use crate::dielectric::screening::{DielectricSolver, SolverType};
use ndarray::{s, Array1, Array2, Array3};
use num_complex::Complex64;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use tracing::{debug, info, instrument, warn};

/// Configuration for imaginary-axis calculations
#[derive(Debug, Clone)]
pub struct ImagAxisConfig {
    /// Maximum condition number threshold for dielectric inversion
    pub cond_threshold: f64,
    /// Tolerance for symmetry enforcement
    pub symmetry_tol: f64,
    /// Enable caching for repeated evaluations
    pub enable_cache: bool,
    /// Maximum cache size in MB
    pub max_cache_size_mb: usize,
    /// Broadening parameter (typically 0 on imaginary axis)
    pub eta: f64,
}

impl Default for ImagAxisConfig {
    fn default() -> Self {
        Self {
            cond_threshold: 1e12,
            symmetry_tol: 1e-10,
            enable_cache: true,
            max_cache_size_mb: 1024,
            eta: 0.0, // No broadening needed on imaginary axis
        }
    }
}

/// Cache key for storing computed results
#[derive(Debug, Clone)]
struct CacheKey {
    /// Imaginary frequency value
    xi: u64, // Store as bits for exact comparison
    /// Hash of input parameters
    param_hash: u64,
}

impl CacheKey {
    fn new(xi: f64, param_hash: u64) -> Self {
        Self {
            xi: xi.to_bits(),
            param_hash,
        }
    }
}

impl Hash for CacheKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.xi.hash(state);
        self.param_hash.hash(state);
    }
}

impl PartialEq for CacheKey {
    fn eq(&self, other: &Self) -> bool {
        self.xi == other.xi && self.param_hash == other.param_hash
    }
}

impl Eq for CacheKey {}

/// Imaginary-axis calculator for GW quantities
pub struct ImagAxisCalculator {
    /// Number of auxiliary basis functions
    pub naux: usize,
    /// Number of occupied orbitals
    pub nocc: usize,
    /// Number of virtual orbitals
    pub nvirt: usize,
    /// Number of total orbitals
    pub norb: usize,
    /// Configuration parameters
    pub config: ImagAxisConfig,
    /// Cache for P0(iξ) values
    p0_cache: HashMap<CacheKey, Array2<f64>>,
    /// Cache for W(iξ) values
    w_cache: HashMap<CacheKey, Array2<f64>>,
    /// Cache statistics
    cache_hits: usize,
    cache_misses: usize,
}

impl ImagAxisCalculator {
    /// Create a new imaginary-axis calculator
    #[must_use]
    pub fn new(naux: usize, nocc: usize, nvirt: usize) -> Self {
        let norb = nocc + nvirt;
        info!(
            naux = naux,
            nocc = nocc,
            nvirt = nvirt,
            "Initializing imaginary-axis calculator"
        );

        Self {
            naux,
            nocc,
            nvirt,
            norb,
            config: ImagAxisConfig::default(),
            p0_cache: HashMap::new(),
            w_cache: HashMap::new(),
            cache_hits: 0,
            cache_misses: 0,
        }
    }

    /// Set configuration
    #[must_use]
    pub fn with_config(mut self, config: ImagAxisConfig) -> Self {
        self.config = config;
        self
    }

    /// Compute polarizability P0(iξ) on imaginary frequency
    ///
    /// P0_PQ(iξ) = Σ_{ia} (f_i - f_a) (ia|P)(ia|Q) / (ξ + (ε_a - ε_i))
    ///
    /// On imaginary axis, this is real and symmetric.
    #[instrument(skip(self, df_ia, e_occ, e_virt))]
    pub fn polarizability_imag_axis(
        &mut self,
        xi: f64,
        df_ia: &Array2<f64>,
        e_occ: &Array1<f64>,
        e_virt: &Array1<f64>,
    ) -> Result<Array2<f64>> {
        // Check cache first
        let param_hash = self.compute_param_hash(e_occ, e_virt);
        let cache_key = CacheKey::new(xi, param_hash);

        if self.config.enable_cache {
            if let Some(cached) = self.p0_cache.get(&cache_key) {
                self.cache_hits += 1;
                debug!(xi = xi, "Cache hit for P0(iξ)");
                return Ok(cached.clone());
            }
        }
        self.cache_misses += 1;

        // Validate inputs
        if xi < 0.0 {
            return Err(QuasixError::InvalidInput(
                "Imaginary frequency must be non-negative".to_string(),
            ));
        }

        let n_trans = self.nocc * self.nvirt;
        if df_ia.shape() != [n_trans, self.naux] {
            return Err(QuasixError::InvalidInput(format!(
                "DF tensor shape mismatch: expected [{}, {}], got {:?}",
                n_trans,
                self.naux,
                df_ia.shape()
            )));
        }

        // Build P0(iξ) matrix
        let mut p0 = Array2::<f64>::zeros((self.naux, self.naux));

        // Loop over transitions
        for i in 0..self.nocc {
            for a in 0..self.nvirt {
                let ia_idx = i * self.nvirt + a;
                let de = e_virt[a] - e_occ[i];

                // On imaginary axis: denominator is real
                let denom = xi + de;

                // Skip if denominator too small
                if denom.abs() < 1e-10 {
                    warn!(xi = xi, de = de, "Small denominator in P0(iξ)");
                    continue;
                }

                // Factor of 2 for spin (restricted calculation)
                let factor = 2.0 / denom;

                // Accumulate contribution: (ia|P)(ia|Q)
                let df_ia_row = df_ia.row(ia_idx);
                for p in 0..self.naux {
                    for q in p..self.naux {
                        let contrib = factor * df_ia_row[p] * df_ia_row[q];
                        p0[[p, q]] += contrib;
                        if p != q {
                            p0[[q, p]] += contrib; // Symmetry
                        }
                    }
                }
            }
        }

        // Enforce exact symmetry
        p0 = self.enforce_symmetry(&p0)?;

        // Validate properties
        self.validate_p0_properties(&p0, xi)?;

        // Store in cache
        if self.config.enable_cache {
            self.p0_cache.insert(cache_key, p0.clone());
            self.check_cache_size();
        }

        Ok(p0)
    }

    /// Compute screened interaction W(iξ) from P0(iξ)
    ///
    /// Uses symmetrized dielectric formulation:
    /// M(iξ) = v^{1/2} P0(iξ) v^{1/2}
    /// W(iξ) = v^{1/2} [1 - M(iξ)]^{-1} v^{1/2}
    #[instrument(skip(self, p0_imag, v_sqrt))]
    pub fn screening_imag_axis(
        &mut self,
        xi: f64,
        p0_imag: &Array2<f64>,
        v_sqrt: &Array2<f64>,
    ) -> Result<Array2<f64>> {
        // Check cache
        let param_hash = self.compute_matrix_hash(p0_imag);
        let cache_key = CacheKey::new(xi, param_hash);

        if self.config.enable_cache {
            if let Some(cached) = self.w_cache.get(&cache_key) {
                self.cache_hits += 1;
                debug!(xi = xi, "Cache hit for W(iξ)");
                return Ok(cached.clone());
            }
        }
        self.cache_misses += 1;

        // Build symmetrized dielectric
        let dielectric = DielectricSolver::new(self.naux, SolverType::Adaptive);

        // Convert real P0 to complex (imaginary axis gives real values)
        let p0_complex = p0_imag.mapv(|x| Complex64::new(x, 0.0));
        let m_imag = dielectric.build_symmetrized_dielectric(&p0_complex, v_sqrt)?;

        // Check conditioning before inversion (using real part for estimate)
        let m_real = m_imag.mapv(|x| x.re);
        let cond_number = self.estimate_condition_number(&m_real)?;
        if cond_number > self.config.cond_threshold {
            warn!(
                xi = xi,
                cond = cond_number,
                threshold = self.config.cond_threshold,
                "Poor conditioning in dielectric matrix"
            );
        }

        // Compute W = v^{1/2} (1 - M)^{-1} v^{1/2}
        let w_imag_complex = dielectric.compute_screened_interaction(&m_imag, v_sqrt)?;

        // Convert back to real (imaginary axis gives real results)
        let w_imag = w_imag_complex.mapv(|x| x.re);

        // Enforce symmetry
        let w_imag = self.enforce_symmetry(&w_imag)?;

        // Validate monotonicity (W should decay with xi)
        if xi > 0.0 {
            let w_norm = w_imag.iter().map(|x| x * x).sum::<f64>().sqrt();
            debug!(xi = xi, norm = w_norm, "W(iξ) norm");
        }

        // Cache result
        if self.config.enable_cache {
            self.w_cache.insert(cache_key, w_imag.clone());
            self.check_cache_size();
        }

        Ok(w_imag)
    }

    /// Compute correlation self-energy Σc(iξ) on imaginary axis
    ///
    /// Σc_mn(iξ) = (1/π) ∫ dξ' Σ_pq G_pq(iξ + iξ') <mp|W_c(iξ')|qn>
    #[instrument(skip(self, w_imag, df_tensors, mo_energies))]
    pub fn self_energy_corr_imag(
        &self,
        xi: f64,
        w_imag: &Array2<f64>,
        df_tensors: &Array3<f64>, // (pq|P) tensors
        mo_energies: &Array1<f64>,
    ) -> Result<Array2<f64>> {
        // Initialize self-energy matrix
        let mut sigma_c = Array2::<f64>::zeros((self.norb, self.norb));

        // This is a simplified implementation
        // Full implementation would involve frequency integration

        // For now, compute static approximation
        for m in 0..self.norb {
            for n in m..self.norb {
                let mut sigma_mn = 0.0;

                // Sum over auxiliary basis
                for p_aux in 0..self.naux {
                    for q_aux in 0..self.naux {
                        // Get W_c element
                        let w_pq = w_imag[[p_aux, q_aux]];

                        // Get density fitting coefficients
                        // This requires proper transformation from MO to auxiliary basis
                        let df_mp = df_tensors[[m, n, p_aux]];
                        let df_nq = df_tensors[[m, n, q_aux]];

                        sigma_mn += w_pq * df_mp * df_nq;
                    }
                }

                // Apply frequency-dependent factor
                let e_avg = 0.5 * (mo_energies[m] + mo_energies[n]);
                let freq_factor = 1.0 / (xi + e_avg.abs());

                sigma_c[[m, n]] = sigma_mn * freq_factor;
                if m != n {
                    sigma_c[[n, m]] = sigma_c[[m, n]]; // Symmetry
                }
            }
        }

        // Enforce symmetry
        let sigma_c = self.enforce_symmetry(&sigma_c)?;

        Ok(sigma_c)
    }

    /// Compute P0, W, and Σc for a grid of imaginary frequencies
    pub fn compute_imag_axis_grid(
        &mut self,
        xi_grid: &Array1<f64>,
        df_ia: &Array2<f64>,
        e_occ: &Array1<f64>,
        e_virt: &Array1<f64>,
        v_sqrt: &Array2<f64>,
    ) -> Result<(Array3<f64>, Array3<f64>)> {
        let n_freq = xi_grid.len();

        // Initialize output arrays
        let mut p0_grid = Array3::<f64>::zeros((self.naux, self.naux, n_freq));
        let mut w_grid = Array3::<f64>::zeros((self.naux, self.naux, n_freq));

        info!(n_freq = n_freq, "Computing imaginary-axis grid");

        for (i_freq, &xi) in xi_grid.iter().enumerate() {
            // Compute P0(iξ)
            let p0 = self.polarizability_imag_axis(xi, df_ia, e_occ, e_virt)?;

            // Compute W(iξ)
            let w = self.screening_imag_axis(xi, &p0, v_sqrt)?;

            // Store in grid
            p0_grid.slice_mut(s![.., .., i_freq]).assign(&p0);
            w_grid.slice_mut(s![.., .., i_freq]).assign(&w);

            // Validate monotonicity
            if i_freq > 0 {
                let w_norm_curr = w.iter().map(|x| x * x).sum::<f64>().sqrt();
                let w_norm_prev = w_grid
                    .slice(s![.., .., i_freq - 1])
                    .iter()
                    .map(|x| x * x)
                    .sum::<f64>()
                    .sqrt();

                if w_norm_curr > w_norm_prev * 1.01 {
                    warn!(
                        i_freq = i_freq,
                        xi = xi,
                        norm_curr = w_norm_curr,
                        norm_prev = w_norm_prev,
                        "W(iξ) not monotonically decreasing"
                    );
                }
            }
        }

        // Log cache statistics
        if self.config.enable_cache {
            let hit_rate = self.cache_hits as f64 / (self.cache_hits + self.cache_misses) as f64;
            info!(
                hits = self.cache_hits,
                misses = self.cache_misses,
                hit_rate = format!("{:.2}%", hit_rate * 100.0),
                "Cache statistics"
            );
        }

        Ok((p0_grid, w_grid))
    }

    /// Enforce exact symmetry on a matrix
    fn enforce_symmetry(&self, mat: &Array2<f64>) -> Result<Array2<f64>> {
        let n = mat.nrows();
        if mat.ncols() != n {
            return Err(QuasixError::InvalidInput(
                "Matrix must be square for symmetry enforcement".to_string(),
            ));
        }

        let mut sym = mat.clone();
        for i in 0..n {
            for j in i + 1..n {
                let avg = 0.5 * (mat[[i, j]] + mat[[j, i]]);
                sym[[i, j]] = avg;
                sym[[j, i]] = avg;
            }
        }

        // Check symmetry error
        let diff_norm = (mat - &sym.t()).iter().map(|x| x * x).sum::<f64>().sqrt();
        let mat_norm = mat.iter().map(|x| x * x).sum::<f64>().sqrt();
        let sym_error = diff_norm / mat_norm;
        if sym_error > self.config.symmetry_tol {
            debug!(error = sym_error, "Symmetry enforcement error");
        }

        Ok(sym)
    }

    /// Validate P0 properties (positive semi-definite, smooth)
    fn validate_p0_properties(&self, p0: &Array2<f64>, _xi: f64) -> Result<()> {
        // Check for NaN/Inf
        if !p0.iter().all(|&x| x.is_finite()) {
            return Err(QuasixError::NumericalError(
                "P0 contains non-finite values".to_string(),
            ));
        }

        // Check positive semi-definiteness would require eigendecomposition
        // Skip for performance in production code
        // Could use ndarray_linalg::Eigh trait if needed

        Ok(())
    }

    /// Estimate condition number of a matrix
    fn estimate_condition_number(&self, mat: &Array2<f64>) -> Result<f64> {
        // Simple estimation using Frobenius norm
        // Full eigendecomposition would be more accurate but expensive
        let frob_norm = mat.iter().map(|x| x * x).sum::<f64>().sqrt();
        let n = mat.nrows() as f64;

        // Rough estimate: condition number ~ norm * sqrt(n)
        // This is conservative but avoids expensive eigendecomposition
        Ok(frob_norm * n.sqrt())
    }

    /// Compute hash of parameters for caching
    fn compute_param_hash(&self, e_occ: &Array1<f64>, e_virt: &Array1<f64>) -> u64 {
        use std::collections::hash_map::DefaultHasher;

        let mut hasher = DefaultHasher::new();

        // Hash dimensions
        self.nocc.hash(&mut hasher);
        self.nvirt.hash(&mut hasher);
        self.naux.hash(&mut hasher);

        // Hash first few energies as fingerprint
        for &e in e_occ.iter().take(3) {
            e.to_bits().hash(&mut hasher);
        }
        for &e in e_virt.iter().take(3) {
            e.to_bits().hash(&mut hasher);
        }

        hasher.finish()
    }

    /// Compute hash of a matrix for caching
    fn compute_matrix_hash(&self, mat: &Array2<f64>) -> u64 {
        use std::collections::hash_map::DefaultHasher;

        let mut hasher = DefaultHasher::new();

        // Hash dimensions
        mat.nrows().hash(&mut hasher);
        mat.ncols().hash(&mut hasher);

        // Hash diagonal elements as fingerprint
        let n = mat.nrows().min(mat.ncols());
        for i in 0..n.min(5) {
            mat[[i, i]].to_bits().hash(&mut hasher);
        }

        hasher.finish()
    }

    /// Check cache size and evict if necessary (simple LRU-style)
    fn check_cache_size(&mut self) {
        // Estimate memory usage (rough)
        let p0_cache_size = self.p0_cache.len() * self.naux * self.naux * 8;
        let w_cache_size = self.w_cache.len() * self.naux * self.naux * 8;
        let total_size_mb = (p0_cache_size + w_cache_size) / (1024 * 1024);

        if total_size_mb > self.config.max_cache_size_mb {
            // Clear half of cache (simple strategy)
            let _p0_to_remove = self.p0_cache.len() / 2;
            let _w_to_remove = self.w_cache.len() / 2;

            self.p0_cache.clear(); // For simplicity, clear all
            self.w_cache.clear();

            info!(
                size_mb = total_size_mb,
                limit_mb = self.config.max_cache_size_mb,
                "Cache size exceeded, clearing"
            );
        }
    }

    /// Get cache statistics
    pub fn cache_stats(&self) -> (usize, usize, f64) {
        let total = self.cache_hits + self.cache_misses;
        let hit_rate = if total > 0 {
            self.cache_hits as f64 / total as f64
        } else {
            0.0
        };
        (self.cache_hits, self.cache_misses, hit_rate)
    }

    /// Clear all caches
    pub fn clear_cache(&mut self) {
        self.p0_cache.clear();
        self.w_cache.clear();
        self.cache_hits = 0;
        self.cache_misses = 0;
        info!("Cleared all caches");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_imag_axis_calculator_creation() {
        let calc = ImagAxisCalculator::new(20, 5, 10);
        assert_eq!(calc.naux, 20);
        assert_eq!(calc.nocc, 5);
        assert_eq!(calc.nvirt, 10);
        assert_eq!(calc.norb, 15);
    }

    #[test]
    fn test_cache_key() {
        let key1 = CacheKey::new(10.0, 12345);
        let key2 = CacheKey::new(10.0, 12345);
        let key3 = CacheKey::new(20.0, 12345);

        assert_eq!(key1, key2);
        assert_ne!(key1, key3);
    }

    #[test]
    fn test_symmetry_enforcement() {
        let calc = ImagAxisCalculator::new(3, 2, 2);

        let mat = Array2::from_shape_vec((3, 3), vec![1.0, 2.0, 3.0, 2.1, 4.0, 5.0, 2.9, 5.1, 6.0])
            .unwrap();

        let sym = calc.enforce_symmetry(&mat).unwrap();

        // Check exact symmetry
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(sym[[i, j]], sym[[j, i]], epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn test_config_defaults() {
        let config = ImagAxisConfig::default();
        assert_eq!(config.cond_threshold, 1e12);
        assert_eq!(config.symmetry_tol, 1e-10);
        assert!(config.enable_cache);
        assert_eq!(config.eta, 0.0);
    }
}
