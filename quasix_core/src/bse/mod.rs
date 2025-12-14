//! Bethe-Salpeter Equation (BSE) solver module
//!
//! This module provides BSE functionality for computing optical excitations:
//!
//! ## Implemented Features (Sprint 6)
//! - **BSE-TDA kernel** (`kernel` module): Matrix-free BSE-TDA kernel for singlet/triplet excitations
//!   - Exchange kernel K^x = (ia|jb) via density fitting
//!   - Direct kernel K^d = -(ij|W|ab) with static screening
//!   - Full Hamiltonian application y = H * x
//! - **Davidson solver** (`davidson` module): Iterative eigenvalue solver for BSE-TDA
//!   - Matrix-free approach using kernel's apply_tda_hamiltonian_batch
//!   - Diagonal preconditioner with level shift
//!   - Subspace expansion with restart capability
//! - **Optical properties** (`optical` module): Oscillator strengths and absorption spectra
//!   - Oscillator strength formula: f_n = (2/3) * E_n * |d_n|^2
//!   - Lorentzian and Gaussian broadening for absorption spectra
//!   - Thomas-Reiche-Kuhn sum rule validation
//!
//! - **Exciton analysis** (`exciton` module): NTO decomposition and density analysis
//!   - Amplitude heatmaps |X^S_{ia}|²
//!   - Natural Transition Orbitals via SVD
//!   - Participation ratio for multi-configurational character
//!   - Hole and electron density matrices
//!
//! ## Planned Features
//! - Full BSE (beyond TDA)
//!
//! ## Usage
//!
//! For matrix-free BSE-TDA calculations, use the `BSETDAKernel`:
//!
//! ```ignore
//! use quasix_core::bse::{BSETDAKernel, BSEKernelConfig, SpinType, DavidsonConfig};
//!
//! // Create kernel with input tensors
//! let kernel = BSETDAKernel::new(
//!     nocc, nvirt, naux,
//!     delta_qp,  // QP energy differences
//!     df_ia,     // DF tensor (ia|P)
//!     df_ij,     // DF tensor (ij|P)
//!     df_ab,     // DF tensor (ab|P)
//!     w0,        // Static screened interaction
//!     BSEKernelConfig::default()
//! )?;
//!
//! // Solve for lowest excitation energies using Davidson
//! let config = DavidsonConfig::with_n_roots(5);
//! let result = kernel.solve_davidson(&config)?;
//!
//! println!("Excitation energies: {:?}", result.eigenvalues);
//! println!("Converged in {} iterations", result.iterations);
//! ```
//!
//! For computing optical properties from BSE-TDA results:
//!
//! ```ignore
//! use quasix_core::bse::{compute_oscillator_strengths, compute_absorption_spectrum, OpticalConfig};
//!
//! // Compute oscillator strengths
//! let optical = compute_oscillator_strengths(eigenvalues, eigenvectors, transition_dipoles)?;
//! println!("Oscillator strengths: {:?}", optical.oscillator_strengths);
//! println!("Sum rule: f_sum = {:.2}", optical.f_sum);
//!
//! // Compute absorption spectrum
//! let config = OpticalConfig::default();
//! let spectrum = compute_absorption_spectrum(eigenvalues, optical.oscillator_strengths.view(), &config)?;
//! ```
//!
//! ## Real Implementations
//! For GW calculations (completed in Sprints 3-5), use:
//! - `quasix_core::gw::evgw` for evGW calculations
//! - `quasix_core::selfenergy` for self-energy evaluations
//! - `quasix_core::qp::solver` for quasiparticle equation solving

#![allow(clippy::many_single_char_names)] // Mathematical notation

// Matrix-free BSE-TDA kernel implementation
pub mod kernel;

// Davidson eigenvalue solver
pub mod davidson;

// Optical properties (oscillator strengths and absorption spectra)
pub mod optical;

// Exciton density analysis (NTOs, participation ratio, density matrices)
pub mod exciton;

// Re-export kernel types at module level
pub use kernel::{BSEKernelConfig, BSETDAKernel, SpinType, MODULE_VERSION, SYMMETRY_THRESHOLD};

// Re-export Davidson types at module level
pub use davidson::{davidson_bse, DavidsonConfig, DavidsonResult};

// Re-export optical types at module level
pub use optical::{
    compute_absorption_spectrum, compute_optical_properties, compute_oscillator_strengths,
    BroadeningType, OpticalConfig, OpticalResult, SpectrumResult,
};

// Re-export exciton analysis types at module level
pub use exciton::{
    analyze_exciton, analyze_multiple_excitons, compute_amplitude_heatmap,
    compute_electron_density, compute_hole_density, compute_nto_decomposition,
    compute_participation_ratio, find_dominant_transition, verify_svd_roundtrip,
    ExcitonAnalysisResult, ExcitonConfig, MODULE_VERSION as EXCITON_MODULE_VERSION,
    SYMMETRY_THRESHOLD as EXCITON_SYMMETRY_THRESHOLD,
};

use crate::common::Result;
use ndarray::{Array1, Array2, Array3};

/// Legacy BSE kernel builder (scaffolding)
///
/// **Note**: For matrix-free BSE-TDA operations, use [`BSETDAKernel`] instead.
/// This struct is kept for backward compatibility with existing scaffolding.
#[deprecated(
    since = "0.6.0",
    note = "Use BSETDAKernel for matrix-free BSE-TDA operations"
)]
pub struct BSEKernelLegacy {
    /// Number of occupied orbitals
    pub nocc: usize,
    /// Number of virtual orbitals
    pub nvirt: usize,
    /// Number of auxiliary basis functions
    pub naux: usize,
    /// Use Tamm-Dancoff approximation
    pub use_tda: bool,
}

#[allow(deprecated)]
impl BSEKernelLegacy {
    /// Create a new BSE kernel builder
    #[must_use]
    pub fn new(nocc: usize, nvirt: usize, naux: usize, use_tda: bool) -> Self {
        Self {
            nocc,
            nvirt,
            naux,
            use_tda,
        }
    }

    /// Build exchange kernel Kx
    pub fn build_exchange_kernel(&self, _df_tensor: &Array3<f64>) -> Result<Array2<f64>> {
        todo!("Implement exchange kernel Kx - use BSETDAKernel::apply_exchange_kernel instead")
    }

    /// Build direct kernel Kd (screened Coulomb)
    pub fn build_direct_kernel(&self, _w: &Array2<f64>) -> Result<Array2<f64>> {
        todo!("Implement direct kernel Kd - use BSETDAKernel::apply_direct_kernel instead")
    }

    /// Apply BSE kernel to vector (matrix-free)
    pub fn apply_kernel(&self, _vec: &Array1<f64>) -> Result<Array1<f64>> {
        todo!("Use BSETDAKernel::apply_tda_hamiltonian instead")
    }

    /// Get kernel dimension
    pub fn dimension(&self) -> usize {
        self.nocc * self.nvirt
    }
}

/// Legacy Davidson solver struct (scaffolding)
///
/// **Deprecated**: Use [`DavidsonConfig`] and [`BSETDAKernel::solve_davidson`] instead.
///
/// This struct is kept for backward compatibility. For new code, use:
/// ```ignore
/// use quasix_core::bse::{BSETDAKernel, DavidsonConfig};
///
/// let config = DavidsonConfig::with_n_roots(5);
/// let result = kernel.solve_davidson(&config)?;
/// ```
#[deprecated(
    since = "0.6.2",
    note = "Use DavidsonConfig and BSETDAKernel::solve_davidson instead"
)]
pub struct DavidsonSolver {
    /// Maximum iterations
    pub max_iter: usize,
    /// Convergence threshold
    pub tolerance: f64,
    /// Number of roots to find
    pub nroots: usize,
    /// Maximum subspace size
    pub max_space: usize,
}

#[allow(deprecated)]
impl DavidsonSolver {
    /// Create a new Davidson solver
    #[must_use]
    pub fn new(nroots: usize) -> Self {
        Self {
            max_iter: 100,
            tolerance: 1e-6,
            nroots,
            max_space: nroots * 10,
        }
    }

    /// Solve BSE eigenvalue problem using the matrix-free kernel
    ///
    /// # Arguments
    ///
    /// * `kernel` - BSE-TDA kernel for matrix-vector products
    ///
    /// # Returns
    ///
    /// DavidsonResult containing eigenvalues, eigenvectors, and convergence info
    pub fn solve(&self, kernel: &BSETDAKernel) -> Result<DavidsonResult> {
        let config = DavidsonConfig {
            n_roots: self.nroots,
            max_iter: self.max_iter,
            max_space: self.max_space,
            tol_residual: self.tolerance,
            ..DavidsonConfig::default()
        };
        kernel.solve_davidson(&config)
    }
}

/// Optical properties calculator
///
/// **Note**: This is scaffolding for Sprint 6-3.
pub struct OpticalProperties {
    /// Transition dipole moments
    pub dipoles: Array2<f64>,
    /// BSE eigenvalues (excitation energies)
    pub energies: Array1<f64>,
    /// BSE eigenvectors
    pub eigenvectors: Array2<f64>,
}

impl OpticalProperties {
    /// Create new optical properties
    #[must_use]
    pub fn new(dipoles: Array2<f64>, energies: Array1<f64>, eigenvectors: Array2<f64>) -> Self {
        Self {
            dipoles,
            energies,
            eigenvectors,
        }
    }

    /// Compute oscillator strengths
    pub fn oscillator_strengths(&self) -> Result<Array1<f64>> {
        todo!("Implement oscillator strength calculation (Sprint 6-3)")
    }

    /// Compute absorption spectrum
    pub fn absorption_spectrum(
        &self,
        _omega_range: (f64, f64),
        _broadening: f64,
    ) -> Result<Array1<f64>> {
        todo!("Implement absorption spectrum (Sprint 6-3)")
    }

    /// Compute exciton density
    pub fn exciton_density(&self, _state: usize) -> Result<Array2<f64>> {
        todo!("Implement exciton density calculation (Sprint 6-3)")
    }
}

/// BSE-TDA driver
///
/// High-level driver that combines the BSE-TDA kernel with the Davidson solver.
///
/// This driver provides a convenient interface for running complete BSE-TDA
/// calculations with automatic eigenvalue solving.
pub struct BSETDADriver {
    /// BSE-TDA kernel for matrix-free operations
    pub kernel: BSETDAKernel,
    /// Davidson solver configuration
    pub config: DavidsonConfig,
}

impl BSETDADriver {
    /// Create a new BSE-TDA driver
    ///
    /// # Arguments
    ///
    /// * `kernel` - Pre-configured BSE-TDA kernel
    /// * `nroots` - Number of eigenvalues to compute
    #[must_use]
    pub fn new(kernel: BSETDAKernel, nroots: usize) -> Self {
        Self {
            kernel,
            config: DavidsonConfig::with_n_roots(nroots),
        }
    }

    /// Create driver with custom Davidson configuration
    #[must_use]
    pub fn with_config(kernel: BSETDAKernel, config: DavidsonConfig) -> Self {
        Self { kernel, config }
    }

    /// Run BSE-TDA calculation
    ///
    /// Computes the lowest excitation energies using the Davidson algorithm.
    ///
    /// # Returns
    ///
    /// DavidsonResult containing excitation energies and exciton wavefunctions
    pub fn run(&self) -> Result<DavidsonResult> {
        let result = self.kernel.solve_davidson(&self.config)?;

        // Validate that excitation energies are physical (positive)
        self.validate_energies(&result.eigenvalues)?;

        Ok(result)
    }

    /// Validate excitation energies
    ///
    /// Checks that all excitation energies are positive (physical constraint).
    pub fn validate_energies(&self, energies: &Array1<f64>) -> Result<()> {
        for &e in energies {
            if e < 0.0 {
                return Err(crate::common::QuasixError::NumericalError(format!(
                    "Negative excitation energy: {}",
                    e
                )));
            }
        }
        Ok(())
    }

    /// Get the number of transitions (BSE dimension)
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.kernel.dimension()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{arr1, arr2};

    /// Create a minimal test kernel
    fn create_test_kernel() -> BSETDAKernel {
        let nocc = 1;
        let nvirt = 1;
        let naux = 2;

        let delta_qp = arr1(&[0.5]);
        let df_ia = arr2(&[[0.5, 0.3]]);
        let df_ij = arr2(&[[0.4, 0.2]]);
        let df_ab = arr2(&[[0.6, 0.1]]);
        let w0 = arr2(&[[1.0, 0.1], [0.1, 1.0]]);

        BSETDAKernel::new(
            nocc,
            nvirt,
            naux,
            delta_qp,
            df_ia,
            df_ij,
            df_ab,
            w0,
            BSEKernelConfig::default(),
        )
        .expect("Failed to create test kernel")
    }

    #[test]
    fn test_bse_tda_kernel_creation() {
        let kernel = create_test_kernel();
        assert_eq!(kernel.nocc, 1);
        assert_eq!(kernel.nvirt, 1);
        assert_eq!(kernel.dimension(), 1);
    }

    #[test]
    fn test_davidson_config_creation() {
        let config = DavidsonConfig::with_n_roots(5);
        assert_eq!(config.n_roots, 5);
        assert_eq!(config.max_space, 100); // 20 * n_roots
    }

    #[test]
    fn test_bse_tda_driver_creation() {
        let kernel = create_test_kernel();
        let driver = BSETDADriver::new(kernel, 3);
        assert_eq!(driver.kernel.nocc, 1);
        assert_eq!(driver.config.n_roots, 3);
    }

    #[test]
    fn test_davidson_solve_trivial() {
        // Test Davidson solver on trivial 1x1 case
        let kernel = create_test_kernel();
        let result = kernel.solve(1).expect("Davidson should converge");

        assert!(result.all_converged());
        assert_eq!(result.iterations, 1);
        assert_eq!(result.eigenvalues.len(), 1);

        // Eigenvalue should equal the only diagonal element of H
        // H = delta_qp + 2*Kx + Kd = 0.5 + 2*0.34 + (-0.276) = 0.904
        approx::assert_relative_eq!(result.eigenvalues[0], 0.904, epsilon = 1e-8);
    }

    #[test]
    fn test_validate_energies_positive() {
        let kernel = create_test_kernel();
        let driver = BSETDADriver::new(kernel, 1);

        let good_energies = arr1(&[0.1, 0.2, 0.3]);
        assert!(driver.validate_energies(&good_energies).is_ok());
    }

    #[test]
    fn test_validate_energies_negative() {
        let kernel = create_test_kernel();
        let driver = BSETDADriver::new(kernel, 1);

        let bad_energies = arr1(&[0.1, -0.2, 0.3]);
        assert!(driver.validate_energies(&bad_energies).is_err());
    }

    #[test]
    fn test_spin_type_export() {
        // Verify SpinType is properly exported
        assert_eq!(SpinType::Singlet.exchange_prefactor(), 2.0);
        assert_eq!(SpinType::Triplet.exchange_prefactor(), 0.0);
    }

    #[test]
    fn test_config_export() {
        // Verify BSEKernelConfig is properly exported
        let config = BSEKernelConfig::default();
        assert_eq!(config.spin_type, SpinType::Singlet);
    }
}
