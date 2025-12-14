//! Exciton density analysis from BSE-TDA eigenvectors
//!
//! This module provides analysis tools for understanding exciton character
//! from BSE-TDA eigenvectors computed by the Davidson solver (S6-2).
//!
//! # Theory Reference
//!
//! See: `docs/derivations/s6-4/theory.md`
//!
//! The exciton wavefunction is encoded in the BSE-TDA eigenvector:
//!
//! ```text
//! |Ψ^S⟩ = Σ_{ia} X^S_{ia} |i → a⟩
//! ```
//!
//! where i denotes occupied and a denotes virtual orbitals.
//!
//! # Features
//!
//! - **Amplitude heatmap**: 2D visualization of |X^S_{ia}|²
//! - **NTO decomposition**: SVD-based Natural Transition Orbitals
//! - **Participation ratio**: Multi-configurational character measure
//! - **Density matrices**: Hole and electron densities in MO basis
//!
//! # Usage
//!
//! ```ignore
//! use quasix_core::bse::exciton::{analyze_exciton, ExcitonAnalysisResult};
//!
//! // After Davidson solver
//! let analysis = analyze_exciton(
//!     eigenvalues.view(),
//!     eigenvectors.view(),
//!     nocc, nvirt,
//!     0  // First excited state
//! )?;
//!
//! println!("Participation ratio: {:.2}", analysis.participation_ratio);
//! println!("Dominant transition: {} -> {}", analysis.dominant_transition.0, analysis.dominant_transition.1);
//! ```
//!
//! # References
//!
//! 1. Martin, R.L. "Natural Transition Orbitals" J. Chem. Phys. 118, 4775 (2003)
//! 2. PySCF: pyscf/tddft/analyze.py (NTO analysis)

use crate::common::{QuasixError, Result};
use ndarray::{Array1, Array2, ArrayView1, ArrayView2};

/// Module version for compatibility tracking
pub const MODULE_VERSION: &str = "0.6.4";

/// Threshold for numerical symmetry checks
pub const SYMMETRY_THRESHOLD: f64 = 1e-10;

/// Result of exciton analysis for a single excited state
#[derive(Debug, Clone)]
pub struct ExcitonAnalysisResult {
    /// State index (0-based)
    pub state_index: usize,
    /// Excitation energy (Hartree)
    pub excitation_energy: f64,
    /// Amplitude matrix |X^S_{ia}|², shape [nocc, nvirt]
    pub amplitude_squared: Array2<f64>,
    /// NTO occupation numbers (λ_n = σ_n²), sorted descending
    pub nto_occupations: Array1<f64>,
    /// Participation ratio: PR = 1 / Σλ_n²
    pub participation_ratio: f64,
    /// Hole density matrix in MO basis, shape [nocc, nocc]
    pub hole_density: Array2<f64>,
    /// Electron density matrix in MO basis, shape [nvirt, nvirt]
    pub electron_density: Array2<f64>,
    /// Dominant transition (i, a, weight)
    pub dominant_transition: (usize, usize, f64),
    /// SVD U matrix (hole NTO coefficients), shape [nocc, min(nocc, nvirt)]
    pub nto_u: Array2<f64>,
    /// SVD singular values
    pub nto_singular_values: Array1<f64>,
    /// SVD Vt matrix (particle NTO coefficients), shape [min(nocc, nvirt), nvirt]
    pub nto_vt: Array2<f64>,
}

/// Configuration for exciton analysis
#[derive(Debug, Clone)]
pub struct ExcitonConfig {
    /// Number of top NTOs to report (default: 5)
    pub n_top_ntos: usize,
    /// Number of top transitions to report (default: 5)
    pub n_top_transitions: usize,
    /// Threshold for significant NTO occupation (default: 1e-4)
    pub nto_threshold: f64,
}

impl Default for ExcitonConfig {
    fn default() -> Self {
        Self {
            n_top_ntos: 5,
            n_top_transitions: 5,
            nto_threshold: 1e-4,
        }
    }
}

/// Compute amplitude heatmap |X^S_{ia}|² for visualization
///
/// Reshapes the eigenvector from [n_trans] to [nocc, nvirt] and computes
/// the squared magnitude for heatmap display.
///
/// # Arguments
///
/// * `eigenvector` - BSE eigenvector X^S, shape [n_trans] where n_trans = nocc * nvirt
/// * `nocc` - Number of occupied orbitals
/// * `nvirt` - Number of virtual orbitals
///
/// # Returns
///
/// Array2 of shape [nocc, nvirt] containing |X^S_{ia}|²
///
/// # Errors
///
/// Returns error if eigenvector length doesn't match nocc * nvirt
pub fn compute_amplitude_heatmap(
    eigenvector: ArrayView1<f64>,
    nocc: usize,
    nvirt: usize,
) -> Result<Array2<f64>> {
    let n_trans = nocc * nvirt;
    if eigenvector.len() != n_trans {
        return Err(QuasixError::DimensionMismatch(format!(
            "Eigenvector length {} != nocc*nvirt = {}",
            eigenvector.len(),
            n_trans
        )));
    }

    // Reshape to [nocc, nvirt] - row-major order (i varies slowest)
    // X_{ia} where compound index = i * nvirt + a
    let mut amplitude_squared = Array2::zeros((nocc, nvirt));
    for i in 0..nocc {
        for a in 0..nvirt {
            let idx = i * nvirt + a;
            amplitude_squared[[i, a]] = eigenvector[idx] * eigenvector[idx];
        }
    }

    Ok(amplitude_squared)
}

/// Reshape eigenvector to amplitude matrix
///
/// Converts flat eigenvector [n_trans] to matrix [nocc, nvirt]
fn reshape_to_amplitude_matrix(
    eigenvector: ArrayView1<f64>,
    nocc: usize,
    nvirt: usize,
) -> Result<Array2<f64>> {
    let n_trans = nocc * nvirt;
    if eigenvector.len() != n_trans {
        return Err(QuasixError::DimensionMismatch(format!(
            "Eigenvector length {} != nocc*nvirt = {}",
            eigenvector.len(),
            n_trans
        )));
    }

    let mut amplitude = Array2::zeros((nocc, nvirt));
    for i in 0..nocc {
        for a in 0..nvirt {
            let idx = i * nvirt + a;
            amplitude[[i, a]] = eigenvector[idx];
        }
    }

    Ok(amplitude)
}

/// Compute NTO decomposition via SVD
///
/// Decomposes the amplitude matrix X^S = U · Σ · V^T where:
/// - U contains hole NTOs in the occupied space
/// - Σ contains singular values (√λ_n)
/// - V^T contains particle NTOs in the virtual space
///
/// NTO occupation numbers are λ_n = σ_n²
///
/// # Arguments
///
/// * `amplitude_matrix` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Tuple (U, singular_values, Vt) where:
/// - U: [nocc, k] hole NTO coefficients (k = min(nocc, nvirt))
/// - singular_values: [k] singular values
/// - Vt: [k, nvirt] particle NTO coefficients
///
/// # Errors
///
/// Returns error if SVD computation fails
pub fn compute_nto_decomposition(
    amplitude_matrix: ArrayView2<f64>,
) -> Result<(Array2<f64>, Array1<f64>, Array2<f64>)> {
    let (m, n) = amplitude_matrix.dim();
    let k = m.min(n);

    // Simple SVD via eigendecomposition of X^T X and X X^T
    // For small matrices (typical BSE), this is sufficient

    if k == 0 {
        return Ok((
            Array2::zeros((m, 0)),
            Array1::zeros(0),
            Array2::zeros((0, n)),
        ));
    }

    // Compute X^T X for right singular vectors (V)
    let xtx = amplitude_matrix.t().dot(&amplitude_matrix);

    // Compute X X^T for left singular vectors (U)
    let xxt = amplitude_matrix.dot(&amplitude_matrix.t());

    // Get eigenvalues and eigenvectors of X^T X
    // λ_i(X^T X) = σ_i² (squared singular values)
    let (eigenvalues_v, eigenvectors_v) = symmetric_eigendecomposition(&xtx)?;

    // Get eigenvectors of X X^T
    let (eigenvalues_u, eigenvectors_u) = symmetric_eigendecomposition(&xxt)?;

    // Singular values are sqrt of eigenvalues, sorted descending
    let mut sorted_indices: Vec<usize> = (0..eigenvalues_v.len()).collect();
    sorted_indices.sort_by(|&a, &b| {
        eigenvalues_v[b]
            .abs()
            .partial_cmp(&eigenvalues_v[a].abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut singular_values = Array1::zeros(k);
    let mut u = Array2::zeros((m, k));
    let mut vt = Array2::zeros((k, n));

    for (new_idx, &old_idx) in sorted_indices.iter().take(k).enumerate() {
        // Singular value = sqrt(eigenvalue)
        let sigma = eigenvalues_v[old_idx].max(0.0).sqrt();
        singular_values[new_idx] = sigma;

        // V columns (transposed to Vt rows)
        if old_idx < eigenvectors_v.ncols() {
            for j in 0..n {
                vt[[new_idx, j]] = eigenvectors_v[[j, old_idx]];
            }
        }

        // U columns - find matching eigenvector from eigenvalues_u
        // Find the eigenvector with closest eigenvalue
        let target_eval = eigenvalues_v[old_idx];
        let u_idx = (0..eigenvalues_u.len())
            .min_by(|&a, &b| {
                (eigenvalues_u[a] - target_eval)
                    .abs()
                    .partial_cmp(&(eigenvalues_u[b] - target_eval).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(0);

        if u_idx < eigenvectors_u.ncols() {
            for i in 0..m {
                u[[i, new_idx]] = eigenvectors_u[[i, u_idx]];
            }
        }

        // Fix sign ambiguity: ensure U * S * Vt reconstructs X
        if sigma > 1e-14 {
            // Compute what u * sigma * vt_row gives
            let reconstructed: f64 = (0..m)
                .map(|i| u[[i, new_idx]] * sigma * (0..n).map(|j| vt[[new_idx, j]]).sum::<f64>())
                .sum();
            let original: f64 = amplitude_matrix.iter().copied().sum();

            // If signs don't match, flip U
            if reconstructed * original < 0.0 {
                for i in 0..m {
                    u[[i, new_idx]] = -u[[i, new_idx]];
                }
            }
        }
    }

    Ok((u, singular_values, vt))
}

/// Simple symmetric eigendecomposition using Jacobi method
///
/// For small matrices typical in BSE analysis, this is adequate.
fn symmetric_eigendecomposition(matrix: &Array2<f64>) -> Result<(Array1<f64>, Array2<f64>)> {
    let n = matrix.nrows();
    if n != matrix.ncols() {
        return Err(QuasixError::DimensionMismatch(
            "Matrix must be square for eigendecomposition".into(),
        ));
    }

    if n == 0 {
        return Ok((Array1::zeros(0), Array2::zeros((0, 0))));
    }

    // For 1x1 case
    if n == 1 {
        return Ok((
            Array1::from_vec(vec![matrix[[0, 0]]]),
            Array2::from_elem((1, 1), 1.0),
        ));
    }

    // Use Jacobi iteration for small symmetric matrices
    let mut a = matrix.to_owned();
    let mut v = Array2::eye(n);
    let max_iter = 100 * n * n;
    let tolerance = 1e-14;

    for _ in 0..max_iter {
        // Find largest off-diagonal element
        let mut max_val = 0.0;
        let mut p = 0;
        let mut q = 1;

        for i in 0..n {
            for j in (i + 1)..n {
                if a[[i, j]].abs() > max_val {
                    max_val = a[[i, j]].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < tolerance {
            break;
        }

        // Compute rotation angle
        let theta = if (a[[q, q]] - a[[p, p]]).abs() < 1e-14 {
            std::f64::consts::PI / 4.0
        } else {
            0.5 * ((2.0 * a[[p, q]]) / (a[[q, q]] - a[[p, p]])).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Apply Jacobi rotation
        let mut a_new = a.clone();

        // Update diagonal elements
        a_new[[p, p]] = c * c * a[[p, p]] - 2.0 * c * s * a[[p, q]] + s * s * a[[q, q]];
        a_new[[q, q]] = s * s * a[[p, p]] + 2.0 * c * s * a[[p, q]] + c * c * a[[q, q]];
        a_new[[p, q]] = 0.0;
        a_new[[q, p]] = 0.0;

        // Update off-diagonal elements
        for i in 0..n {
            if i != p && i != q {
                let aip = c * a[[i, p]] - s * a[[i, q]];
                let aiq = s * a[[i, p]] + c * a[[i, q]];
                a_new[[i, p]] = aip;
                a_new[[p, i]] = aip;
                a_new[[i, q]] = aiq;
                a_new[[q, i]] = aiq;
            }
        }

        // Update eigenvector matrix
        let mut v_new = v.clone();
        for i in 0..n {
            let vip = c * v[[i, p]] - s * v[[i, q]];
            let viq = s * v[[i, p]] + c * v[[i, q]];
            v_new[[i, p]] = vip;
            v_new[[i, q]] = viq;
        }

        a = a_new;
        v = v_new;
    }

    // Extract eigenvalues from diagonal
    let eigenvalues: Array1<f64> = (0..n).map(|i| a[[i, i]]).collect();

    Ok((eigenvalues, v))
}

/// Compute participation ratio from NTO occupations
///
/// The participation ratio quantifies multi-configurational character:
///
/// PR = 1 / Σ_n λ_n²
///
/// where λ_n = σ_n² are the NTO occupation numbers.
///
/// - PR = 1: Single NTO pair dominates (pure single-particle transition)
/// - PR > 1: Multi-configurational excitation
///
/// # Arguments
///
/// * `nto_occupations` - NTO occupation numbers λ_n (must sum to 1)
///
/// # Returns
///
/// Participation ratio value
pub fn compute_participation_ratio(nto_occupations: ArrayView1<f64>) -> f64 {
    let sum_sq: f64 = nto_occupations.iter().map(|&x| x * x).sum();
    if sum_sq > 1e-14 {
        1.0 / sum_sq
    } else {
        1.0
    }
}

/// Compute hole density matrix in MO basis
///
/// The hole density describes where the electron was removed from:
///
/// ρ^h_{ij} = Σ_a X^S_{ia} (X^S_{ja})*
///
/// For real wavefunctions, this is ρ^h = X @ X^T
///
/// # Arguments
///
/// * `amplitude_matrix` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Hole density matrix, shape [nocc, nocc]
pub fn compute_hole_density(amplitude_matrix: ArrayView2<f64>) -> Array2<f64> {
    // ρ^h = X @ X^T -> [nocc, nocc]
    amplitude_matrix.dot(&amplitude_matrix.t())
}

/// Compute electron density matrix in MO basis
///
/// The electron density describes where the electron was placed:
///
/// ρ^e_{ab} = Σ_i (X^S_{ia})* X^S_{ib}
///
/// For real wavefunctions, this is ρ^e = X^T @ X
///
/// # Arguments
///
/// * `amplitude_matrix` - Amplitude matrix X^S, shape [nocc, nvirt]
///
/// # Returns
///
/// Electron density matrix, shape [nvirt, nvirt]
pub fn compute_electron_density(amplitude_matrix: ArrayView2<f64>) -> Array2<f64> {
    // ρ^e = X^T @ X -> [nvirt, nvirt]
    amplitude_matrix.t().dot(&amplitude_matrix)
}

/// Find dominant i→a transition
///
/// Identifies the transition with largest |X^S_{ia}|² contribution.
///
/// # Arguments
///
/// * `amplitude_squared` - Squared amplitudes |X^S_{ia}|², shape [nocc, nvirt]
///
/// # Returns
///
/// Tuple (i, a, weight) for the dominant transition
pub fn find_dominant_transition(amplitude_squared: ArrayView2<f64>) -> (usize, usize, f64) {
    let (nocc, nvirt) = amplitude_squared.dim();
    let mut max_val = 0.0;
    let mut max_i = 0;
    let mut max_a = 0;

    for i in 0..nocc {
        for a in 0..nvirt {
            if amplitude_squared[[i, a]] > max_val {
                max_val = amplitude_squared[[i, a]];
                max_i = i;
                max_a = a;
            }
        }
    }

    (max_i, max_a, max_val)
}

/// Verify SVD roundtrip reconstruction
///
/// Checks that U @ diag(S) @ Vt reconstructs the original matrix
/// within the specified tolerance.
///
/// # Arguments
///
/// * `u` - Left singular vectors, shape [m, k]
/// * `s` - Singular values, shape [k]
/// * `vt` - Right singular vectors (transposed), shape [k, n]
/// * `original` - Original matrix to compare against, shape [m, n]
/// * `tolerance` - Maximum allowed reconstruction error
///
/// # Returns
///
/// true if reconstruction error is within tolerance
pub fn verify_svd_roundtrip(
    u: ArrayView2<f64>,
    s: ArrayView1<f64>,
    vt: ArrayView2<f64>,
    original: ArrayView2<f64>,
    tolerance: f64,
) -> bool {
    let (m, n) = original.dim();
    let k = s.len();

    if u.dim() != (m, k) || vt.dim() != (k, n) {
        return false;
    }

    // Reconstruct: X_reconstructed = U @ diag(S) @ Vt
    let mut reconstructed = Array2::zeros((m, n));

    for i in 0..m {
        for j in 0..n {
            let mut sum = 0.0;
            for l in 0..k {
                sum += u[[i, l]] * s[l] * vt[[l, j]];
            }
            reconstructed[[i, j]] = sum;
        }
    }

    // Check max absolute difference
    let max_diff = original
        .iter()
        .zip(reconstructed.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, |acc, x| acc.max(x));

    max_diff < tolerance
}

/// Full exciton analysis for a single excited state
///
/// Performs complete analysis including:
/// - Amplitude heatmap
/// - NTO decomposition and occupations
/// - Participation ratio
/// - Hole and electron density matrices
/// - Dominant transition identification
///
/// # Arguments
///
/// * `eigenvalues` - Excitation energies [n_roots] (Hartree)
/// * `eigenvectors` - BSE eigenvectors X_{ia,n}, shape [n_trans, n_roots]
/// * `nocc` - Number of occupied orbitals
/// * `nvirt` - Number of virtual orbitals
/// * `state_index` - Which state to analyze (0-indexed)
///
/// # Returns
///
/// ExcitonAnalysisResult containing all analysis data
///
/// # Errors
///
/// Returns error if dimensions are inconsistent or state_index is out of range
pub fn analyze_exciton(
    eigenvalues: ArrayView1<f64>,
    eigenvectors: ArrayView2<f64>,
    nocc: usize,
    nvirt: usize,
    state_index: usize,
) -> Result<ExcitonAnalysisResult> {
    let n_roots = eigenvalues.len();
    let (n_trans, n_roots_ev) = eigenvectors.dim();

    // Validate dimensions
    if n_roots_ev != n_roots {
        return Err(QuasixError::DimensionMismatch(format!(
            "Eigenvalues length {} != eigenvectors columns {}",
            n_roots, n_roots_ev
        )));
    }

    if n_trans != nocc * nvirt {
        return Err(QuasixError::DimensionMismatch(format!(
            "n_trans {} != nocc*nvirt = {}",
            n_trans,
            nocc * nvirt
        )));
    }

    if state_index >= n_roots {
        return Err(QuasixError::NumericalError(format!(
            "state_index {} >= n_roots {}",
            state_index, n_roots
        )));
    }

    // Extract eigenvector for this state
    let eigenvector = eigenvectors.column(state_index);
    let excitation_energy = eigenvalues[state_index];

    // Reshape to amplitude matrix [nocc, nvirt]
    let amplitude = reshape_to_amplitude_matrix(eigenvector, nocc, nvirt)?;

    // Compute amplitude squared for heatmap
    let amplitude_squared = amplitude.mapv(|x| x * x);

    // NTO decomposition via SVD
    let (nto_u, nto_singular_values, nto_vt) = compute_nto_decomposition(amplitude.view())?;

    // NTO occupations = σ² (normalized to sum = 1)
    let nto_occupations_raw: Array1<f64> = nto_singular_values.mapv(|s| s * s);
    let occ_sum = nto_occupations_raw.sum();
    let nto_occupations = if occ_sum > 1e-14 {
        nto_occupations_raw.mapv(|x| x / occ_sum)
    } else {
        nto_occupations_raw
    };

    // Participation ratio
    let participation_ratio = compute_participation_ratio(nto_occupations.view());

    // Density matrices
    let hole_density = compute_hole_density(amplitude.view());
    let electron_density = compute_electron_density(amplitude.view());

    // Dominant transition
    let dominant_transition = find_dominant_transition(amplitude_squared.view());

    Ok(ExcitonAnalysisResult {
        state_index,
        excitation_energy,
        amplitude_squared,
        nto_occupations,
        participation_ratio,
        hole_density,
        electron_density,
        dominant_transition,
        nto_u,
        nto_singular_values,
        nto_vt,
    })
}

/// Analyze multiple excited states
///
/// Convenience function to analyze all requested states in one call.
///
/// # Arguments
///
/// * `eigenvalues` - Excitation energies [n_roots] (Hartree)
/// * `eigenvectors` - BSE eigenvectors, shape [n_trans, n_roots]
/// * `nocc` - Number of occupied orbitals
/// * `nvirt` - Number of virtual orbitals
/// * `n_states` - Number of states to analyze (starting from 0)
///
/// # Returns
///
/// Vector of ExcitonAnalysisResult for each state
pub fn analyze_multiple_excitons(
    eigenvalues: ArrayView1<f64>,
    eigenvectors: ArrayView2<f64>,
    nocc: usize,
    nvirt: usize,
    n_states: usize,
) -> Result<Vec<ExcitonAnalysisResult>> {
    let actual_states = n_states.min(eigenvalues.len());
    let mut results = Vec::with_capacity(actual_states);

    for state_idx in 0..actual_states {
        results.push(analyze_exciton(
            eigenvalues,
            eigenvectors,
            nocc,
            nvirt,
            state_idx,
        )?);
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::{arr1, arr2, Array1, Array2};

    #[test]
    fn test_amplitude_heatmap_simple() {
        // 2x2 system: nocc=2, nvirt=2
        // Eigenvector in flat form: [X_00, X_01, X_10, X_11]
        let eigenvector = arr1(&[0.6, 0.0, 0.0, 0.8]);
        let nocc = 2;
        let nvirt = 2;

        let heatmap = compute_amplitude_heatmap(eigenvector.view(), nocc, nvirt).unwrap();

        assert_eq!(heatmap.shape(), &[2, 2]);
        assert_relative_eq!(heatmap[[0, 0]], 0.36, epsilon = 1e-10);
        assert_relative_eq!(heatmap[[0, 1]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(heatmap[[1, 0]], 0.0, epsilon = 1e-10);
        assert_relative_eq!(heatmap[[1, 1]], 0.64, epsilon = 1e-10);
    }

    #[test]
    fn test_amplitude_normalization() {
        // Normalized eigenvector: |X|² sums to 1
        let eigenvector = arr1(&[0.6, 0.0, 0.0, 0.8]); // 0.36 + 0.64 = 1.0
        let nocc = 2;
        let nvirt = 2;

        let heatmap = compute_amplitude_heatmap(eigenvector.view(), nocc, nvirt).unwrap();
        let sum: f64 = heatmap.iter().sum();

        assert_relative_eq!(sum, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_nto_occupation_sum() {
        // NTO occupations should sum to 1 for normalized eigenvector
        let amplitude = arr2(&[[0.6, 0.0], [0.0, 0.8]]);

        let (_, s, _) = compute_nto_decomposition(amplitude.view()).unwrap();

        // Occupations = s²
        let occupations: Array1<f64> = s.mapv(|x| x * x);
        let sum = occupations.sum();

        // Sum should be ||X||_F² = 0.36 + 0.64 = 1.0
        assert_relative_eq!(sum, 1.0, epsilon = 1e-8);
    }

    #[test]
    fn test_nto_occupations_ordered() {
        // Singular values (and thus occupations) should be ordered descending
        let amplitude = arr2(&[[0.3, 0.4], [0.5, 0.6]]);

        let (_, s, _) = compute_nto_decomposition(amplitude.view()).unwrap();

        // Check ordering
        for i in 1..s.len() {
            assert!(s[i - 1] >= s[i], "Singular values should be descending");
        }
    }

    #[test]
    fn test_participation_ratio_single_nto() {
        // If one NTO has λ=1, others have λ=0, then PR = 1
        let occupations = arr1(&[1.0, 0.0, 0.0]);
        let pr = compute_participation_ratio(occupations.view());
        assert_relative_eq!(pr, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_participation_ratio_equal_ntos() {
        // If n NTOs have λ=1/n each, PR = n
        let n = 4;
        let occupations = Array1::from_elem(n, 1.0 / n as f64);
        let pr = compute_participation_ratio(occupations.view());
        assert_relative_eq!(pr, n as f64, epsilon = 1e-10);
    }

    #[test]
    fn test_participation_ratio_two_dominant() {
        // Two equal NTOs: λ₁ = λ₂ = 0.5 → PR = 1/(0.25 + 0.25) = 2
        let occupations = arr1(&[0.5, 0.5]);
        let pr = compute_participation_ratio(occupations.view());
        assert_relative_eq!(pr, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_hole_density_symmetric() {
        let amplitude = arr2(&[[0.6, 0.2], [0.1, 0.8]]);
        let hole_dm = compute_hole_density(amplitude.view());

        // Hole density should be symmetric
        for i in 0..hole_dm.nrows() {
            for j in 0..hole_dm.ncols() {
                assert_relative_eq!(
                    hole_dm[[i, j]],
                    hole_dm[[j, i]],
                    epsilon = SYMMETRY_THRESHOLD
                );
            }
        }
    }

    #[test]
    fn test_electron_density_symmetric() {
        let amplitude = arr2(&[[0.6, 0.2], [0.1, 0.8]]);
        let elec_dm = compute_electron_density(amplitude.view());

        // Electron density should be symmetric
        for i in 0..elec_dm.nrows() {
            for j in 0..elec_dm.ncols() {
                assert_relative_eq!(
                    elec_dm[[i, j]],
                    elec_dm[[j, i]],
                    epsilon = SYMMETRY_THRESHOLD
                );
            }
        }
    }

    #[test]
    fn test_hole_density_trace() {
        // Trace of hole density = Σ_{ij} ρ^h_{ii} = Σ_i Σ_a |X_{ia}|² = 1 (normalized)
        #[allow(clippy::approx_constant)] // Intentionally using numeric constant for test data
        let amplitude = arr2(&[[0.6, 0.2], [0.1, 0.7071067811865476]]); // Normalized
                                                                        // Actually need proper normalization: sqrt(0.36 + 0.04 + 0.01 + 0.5) = sqrt(0.91)
        let norm: f64 = amplitude.iter().map(|x| x * x).sum::<f64>().sqrt();
        let normalized = amplitude.mapv(|x| x / norm);

        let hole_dm = compute_hole_density(normalized.view());
        let trace: f64 = (0..hole_dm.nrows()).map(|i| hole_dm[[i, i]]).sum();

        assert_relative_eq!(trace, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_svd_roundtrip() {
        let amplitude = arr2(&[[0.6, 0.2], [0.1, 0.8]]);
        let (u, s, vt) = compute_nto_decomposition(amplitude.view()).unwrap();

        assert!(verify_svd_roundtrip(
            u.view(),
            s.view(),
            vt.view(),
            amplitude.view(),
            1e-8
        ));
    }

    #[test]
    fn test_svd_orthogonality_u() {
        // U columns should be orthonormal
        let amplitude = arr2(&[[0.6, 0.2], [0.1, 0.8]]);
        let (u, _, _) = compute_nto_decomposition(amplitude.view()).unwrap();

        let utu = u.t().dot(&u);
        let identity = Array2::eye(u.ncols());

        for i in 0..utu.nrows() {
            for j in 0..utu.ncols() {
                assert_relative_eq!(utu[[i, j]], identity[[i, j]], epsilon = 1e-8);
            }
        }
    }

    #[test]
    fn test_dominant_transition() {
        let amplitude_sq = arr2(&[[0.1, 0.6], [0.3, 0.0]]);
        let (i, a, weight) = find_dominant_transition(amplitude_sq.view());

        assert_eq!(i, 0);
        assert_eq!(a, 1);
        assert_relative_eq!(weight, 0.6, epsilon = 1e-10);
    }

    #[test]
    fn test_analyze_exciton_full() {
        // Create test data: 2 occupied, 2 virtual, 1 root
        let nocc = 2;
        let nvirt = 2;
        let n_trans = nocc * nvirt;

        let eigenvalues = arr1(&[0.3]); // 0.3 Ha excitation
        let eigenvector = arr1(&[0.6, 0.0, 0.0, 0.8]); // Normalized: 0.36 + 0.64 = 1
        let eigenvectors = eigenvector.into_shape_with_order((n_trans, 1)).unwrap();

        let result = analyze_exciton(eigenvalues.view(), eigenvectors.view(), nocc, nvirt, 0)
            .expect("Analysis should succeed");

        assert_eq!(result.state_index, 0);
        assert_relative_eq!(result.excitation_energy, 0.3, epsilon = 1e-10);

        // Check amplitude squared
        assert_relative_eq!(result.amplitude_squared[[0, 0]], 0.36, epsilon = 1e-10);
        assert_relative_eq!(result.amplitude_squared[[1, 1]], 0.64, epsilon = 1e-10);

        // NTO occupations should sum to 1
        let occ_sum: f64 = result.nto_occupations.sum();
        assert_relative_eq!(occ_sum, 1.0, epsilon = 1e-8);

        // Participation ratio for diagonal case with 0.36 and 0.64
        // PR = 1/(0.36² + 0.64²) = 1/(0.1296 + 0.4096) = 1/0.5392 ≈ 1.85
        assert!(result.participation_ratio > 1.0 && result.participation_ratio < 2.1);

        // Dominant transition should be (1, 1) with weight 0.64
        assert_eq!(result.dominant_transition.0, 1);
        assert_eq!(result.dominant_transition.1, 1);
        assert_relative_eq!(result.dominant_transition.2, 0.64, epsilon = 1e-10);
    }

    #[test]
    fn test_dimension_mismatch() {
        let eigenvalues = arr1(&[0.3]);
        let eigenvectors = arr2(&[[1.0]]);
        let nocc = 2; // Wrong: expects 2*2=4 elements
        let nvirt = 2;

        let result = analyze_exciton(eigenvalues.view(), eigenvectors.view(), nocc, nvirt, 0);

        assert!(result.is_err());
    }

    #[test]
    fn test_state_index_out_of_range() {
        let eigenvalues = arr1(&[0.3]);
        let eigenvectors = arr2(&[[0.6], [0.0], [0.0], [0.8]]);

        let result = analyze_exciton(eigenvalues.view(), eigenvectors.view(), 2, 2, 5);

        assert!(result.is_err());
    }

    #[test]
    fn test_single_transition_system() {
        // H2-like: 1 occupied, 1 virtual
        let eigenvalues = arr1(&[0.5]);
        let eigenvectors = arr2(&[[1.0]]); // Single transition, fully normalized

        let result =
            analyze_exciton(eigenvalues.view(), eigenvectors.view(), 1, 1, 0).expect("Should work");

        assert_eq!(result.participation_ratio, 1.0);
        assert_eq!(result.dominant_transition, (0, 0, 1.0));
    }

    #[test]
    fn test_multiple_excitons() {
        let eigenvalues = arr1(&[0.3, 0.4, 0.5]);
        let eigenvectors = arr2(&[
            [0.8, 0.5, 0.3],
            [0.2, 0.3, 0.4],
            [0.3, 0.5, 0.6],
            [0.4, 0.6, 0.6],
        ]);

        let results =
            analyze_multiple_excitons(eigenvalues.view(), eigenvectors.view(), 2, 2, 3).unwrap();

        assert_eq!(results.len(), 3);
        assert_eq!(results[0].state_index, 0);
        assert_eq!(results[1].state_index, 1);
        assert_eq!(results[2].state_index, 2);
    }
}
