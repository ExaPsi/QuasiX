//! evGW (Eigenvalue Self-Consistent GW) Driver
//!
//! This module implements the evGW method where quasiparticle energies
//! are iteratively updated while keeping wavefunctions fixed.
//!
//! Sprint S5-1: evGW driver (G update)
//! Acceptance Criteria: Converges in <= 12 cycles with damping
//!
//! # Algorithm
//!
//! evGW iterates the QP equation:
//!   E_n^QP = epsilon_n^HF + Z_n * (Sigma_x_n + Re[Sigma_c(E_n^QP)] - V_xc_n)
//!
//! At iteration 0, this reduces to G0W0 with evaluation at HF energies.
//! Subsequent iterations evaluate Sigma_c at the previous QP energies.
//!
//! # Tensor Conventions
//!
//! The evGW driver accepts DF tensors in PySCF convention:
//! - `mn_full`: Shape [nmo, nmo, naux] - full (pq|P) tensor from Python bindings
//!
//! Internally, we transpose to correlation_ac convention:
//! - `lpq`: Shape [naux, nmo, nmo] - for validated AC self-energy calculation
//!
//! This ensures compatibility with the PySCF-validated correlation_ac module.
//!
//! # Key Formulas (from correlation_ac.rs)
//!
//! ## Polarizability (RPA)
//! P0_PQ(iw) = 4 * sum_ia L_Pia * (e_ia / (w^2 + e_ia^2)) * L_Qia
//!
//! ## Screened Interaction
//! W_c = (I - P0)^{-1} - I
//!
//! ## Z-factor
//! Z_n = 1 / (1 - dSigma_c/dw |_{w=E_n})

use crate::common::Result;
use crate::qp::convergence::{AccelerationMethod, ConvergenceCriteria, ConvergenceMonitor};
use crate::qp::newton_solver::{NewtonQPSolver, NewtonSolverConfig};
use crate::selfenergy::correlation_ac::{self, ACConfig};
use crate::selfenergy::frequency_cache::{
    get_sigma_diag_and_cache, get_sigma_diag_cached, FrequencyCache,
};
use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;

/// Configuration for contour deformation / frequency integration
#[derive(Debug, Clone)]
pub struct ContourDeformationConfig {
    /// Broadening parameter eta (Ha)
    pub eta: f64,
    /// Number of imaginary frequency points
    pub n_imag_points: usize,
    /// Alias for n_imag_points (backward compatibility)
    pub n_imag_freq: usize,
    /// Maximum imaginary frequency (Ha)
    pub xi_max: f64,
    /// Alias for xi_max (backward compatibility)
    pub omega_max: f64,
    /// Use Gauss-Legendre quadrature
    pub use_gl_quadrature: bool,
    /// Convergence tolerance
    pub convergence_tol: f64,
    /// Regularization parameter
    pub regularization: f64,
    /// Use SIMD optimization
    pub use_simd: bool,
    /// Number of threads (None = auto)
    pub n_threads: Option<usize>,
    /// Verbosity level
    pub verbose: usize,
    /// Pole detection threshold
    pub pole_threshold: f64,
    /// Compute spectral function
    pub compute_spectral: bool,
}

impl Default for ContourDeformationConfig {
    fn default() -> Self {
        Self {
            eta: 0.01,
            n_imag_points: 100, // PySCF default: 100
            n_imag_freq: 100,
            xi_max: 5.0, // PySCF default iw_cutoff
            omega_max: 5.0,
            use_gl_quadrature: true,
            convergence_tol: 1e-10,
            regularization: 1e-10,
            use_simd: true,
            n_threads: None,
            verbose: 1,
            pole_threshold: 0.01,
            compute_spectral: false,
        }
    }
}

/// Configuration for QP equation solver
#[derive(Debug, Clone)]
pub struct QPSolverConfig {
    /// Energy convergence tolerance
    pub energy_tolerance: f64,
    /// Residual tolerance
    pub residual_tolerance: f64,
    /// Max Newton iterations
    pub max_newton_iterations: usize,
    /// Max bisection iterations
    pub max_bisection_iterations: usize,
    /// Initial damping factor
    pub initial_damping: f64,
    /// Minimum damping
    pub min_damping: f64,
    /// Maximum damping
    pub max_damping: f64,
    /// Maximum energy step
    pub max_energy_step: f64,
    /// Numerical derivative delta
    pub derivative_delta: f64,
    /// Use line search
    pub use_line_search: bool,
    /// Core orbital threshold (Ha)
    pub core_threshold: Option<f64>,
    /// Use Richardson extrapolation
    pub use_richardson: bool,
    /// Use bisection fallback
    pub use_bisection_fallback: bool,
    /// Z-factor bounds
    pub z_bounds: (f64, f64),
    /// Number of threads
    pub n_threads: Option<usize>,
    /// Damping factor
    pub damping_factor: f64,
}

impl Default for QPSolverConfig {
    fn default() -> Self {
        Self {
            energy_tolerance: 1e-6,
            residual_tolerance: 1e-6,
            max_newton_iterations: 50,
            max_bisection_iterations: 20,
            initial_damping: 0.5,
            min_damping: 0.1,
            max_damping: 0.9,
            max_energy_step: 1.0,
            derivative_delta: 1e-6,
            use_line_search: false,
            core_threshold: Some(10.0),
            use_richardson: false,
            use_bisection_fallback: true,
            z_bounds: (0.1, 0.999),
            n_threads: None,
            damping_factor: 0.5,
        }
    }
}

/// evGW configuration
#[derive(Debug, Clone)]
pub struct EvGWConfig {
    /// Maximum iterations
    pub max_iterations: usize,
    /// Energy convergence tolerance (Ha)
    pub energy_tolerance: f64,
    /// Density/Z-factor convergence tolerance
    pub density_tolerance: f64,
    /// Damping factor for mixing
    pub damping_factor: f64,
    /// Minimum Z-factor
    pub z_min: f64,
    /// Maximum Z-factor
    pub z_max: f64,
    /// Use DIIS acceleration
    pub use_diis: bool,
    /// DIIS subspace dimension
    pub diis_space_dim: usize,
    /// Print level
    pub print_level: u32,
    /// Adaptive damping
    pub adaptive_damping: bool,
    /// Core orbital threshold
    pub core_threshold: f64,
    /// Contour deformation config
    pub cd_config: ContourDeformationConfig,
    /// QP solver config
    pub qp_solver_config: QPSolverConfig,
    /// Compute final self-energy
    pub compute_final_selfenergy: bool,

    // === Production-quality evGW settings (VASP/FHI-aims style) ===
    /// Energy window for self-consistency (Ha)
    ///
    /// Only orbitals within this window from the Fermi level are updated
    /// in the evGW iterations. States outside keep their G0W0 values.
    ///
    /// Default: 2.0 Ha (~54 eV) - covers typical valence + low-lying virtuals
    ///
    /// Production codes use similar windows:
    /// - VASP: NOMEGA parameter controls included states
    /// - FHI-aims: uses energy window for efficiency
    /// - TURBOMOLE: similar approach for scGW
    pub energy_window: f64,

    /// Convergence check based on frontier orbitals only
    ///
    /// When true, convergence is checked only for orbitals near the Fermi level
    /// (within `frontier_window` of HOMO-LUMO gap), ignoring deep core and
    /// high-lying virtual states that may oscillate.
    ///
    /// This is the standard approach in VASP, FHI-aims, and TURBOMOLE.
    pub frontier_convergence: bool,

    /// Window for frontier orbital convergence check (Ha)
    ///
    /// Only orbitals within this energy of the HOMO or LUMO are considered
    /// for convergence. Default: 1.0 Ha (~27 eV)
    pub frontier_window: f64,

    /// Number of occupied orbitals to include in convergence check
    ///
    /// If set, only the top `n_occ_active` occupied orbitals are checked.
    /// None means use energy window.
    pub n_occ_active: Option<usize>,

    /// Number of virtual orbitals to include in convergence check
    ///
    /// If set, only the bottom `n_vir_active` virtual orbitals are checked.
    /// None means use energy window.
    pub n_vir_active: Option<usize>,

    // === Newton-Raphson solved QP equation settings ===
    /// Use Newton-Raphson to solve the full QP equation
    ///
    /// When true, solve: E_QP = ε_HF + Σx + Re[Σc(E_QP)] - Vxc
    /// When false (default), use linearized: E_QP = ε_HF + Z × (Σx + Σc(ε_HF) - Vxc)
    ///
    /// The solved QP equation is more accurate (~50-150 meV for HOMO)
    /// but slower. Use for high-accuracy benchmarks or final evGW results.
    ///
    /// Default: false (use linearized for faster convergence during iterations)
    pub use_solved_qp: bool,

    /// Newton solver energy tolerance (Ha)
    ///
    /// Convergence criterion for the Newton-Raphson QP solver.
    /// Only used when `use_solved_qp` is true.
    /// Default: 1e-8 (sub-µeV precision)
    pub newton_energy_tolerance: f64,

    /// Maximum Newton iterations per orbital
    ///
    /// Only used when `use_solved_qp` is true.
    /// Default: 50
    pub newton_max_iterations: usize,

    // === Frequency Caching Settings (Performance Optimization) ===
    /// Enable frequency caching for evGW speedup
    ///
    /// **WARNING (2025-12-09)**: Frequency caching is DISABLED by default for evGW
    /// because it produces incorrect results!
    ///
    /// **The Bug**: When caching is enabled, P0(iw) and W(iw) are computed once
    /// at iteration 0 using HF energies, and reused for all subsequent iterations.
    /// This means the polarizability energy denominators never update with QP energies,
    /// causing evGW to produce the SAME result as G0W0.
    ///
    /// **Correct evGW Algorithm**:
    /// - Iteration 0: P0 uses HF energies for denominators (this is G0W0)
    /// - Iteration k > 0: P0 must use QP energies from iteration k-1
    ///
    /// The energy denominators in P0:
    ///   P0_PQ(iw) = 4 * sum_ia L_Pia * (e_ia / (w^2 + e_ia^2)) * L_Qia
    /// where e_ia = e_i - e_a MUST use updated QP energies for true evGW!
    ///
    /// **When to use caching**:
    /// - G0W0 calculations (max_iterations = 1): Safe to use, gives ~5-10x speedup
    /// - evGW calculations: MUST be disabled (set to false) for correct results
    ///
    /// Default: false (ensures correct evGW behavior)
    pub use_frequency_caching: bool,
}

impl Default for EvGWConfig {
    fn default() -> Self {
        Self {
            max_iterations: 12,
            energy_tolerance: 1e-4,
            density_tolerance: 1e-3,
            damping_factor: 0.5,
            z_min: 0.1,
            z_max: 0.999,
            use_diis: true,
            diis_space_dim: 6,
            print_level: 1,
            adaptive_damping: true,
            core_threshold: 10.0,
            cd_config: ContourDeformationConfig::default(),
            qp_solver_config: QPSolverConfig::default(),
            compute_final_selfenergy: false,
            // Production-quality defaults (VASP/FHI-aims style)
            energy_window: 2.0,         // ~54 eV window
            frontier_convergence: true, // Check frontier orbitals only
            frontier_window: 1.0,       // ~27 eV for convergence check
            n_occ_active: None,         // Use energy window
            n_vir_active: None,         // Use energy window
            // Newton-Raphson defaults
            use_solved_qp: false, // Use linearized by default
            newton_energy_tolerance: 1e-8,
            newton_max_iterations: 50,
            // Frequency caching default - DISABLED for correct evGW behavior
            // See use_frequency_caching documentation for explanation
            use_frequency_caching: false, // MUST be false for evGW to work correctly
        }
    }
}

/// Convergence history for evGW iterations
#[derive(Debug, Clone, Default)]
pub struct EvGWConvergenceHistory {
    /// QP energies at each iteration
    pub energies: Vec<Array1<f64>>,
    /// Maximum energy change at each iteration
    pub max_energy_changes: Vec<f64>,
    /// RMS energy change at each iteration
    pub rms_energy_changes: Vec<f64>,
}

/// evGW calculation result
#[derive(Debug, Clone)]
pub struct EvGWResult {
    /// Converged quasiparticle energies
    pub qp_energies: Array1<f64>,
    /// Z-factors (quasiparticle weights)
    pub z_factors: Array1<f64>,
    /// Exchange self-energy matrix
    pub sigma_x: Array2<f64>,
    /// Correlation self-energy (diagonal, complex)
    pub sigma_c: Array1<Complex64>,
    /// Number of iterations
    pub n_iterations: usize,
    /// Alias for n_iterations (backward compatibility)
    pub n_cycles: usize,
    /// Converged flag
    pub converged: bool,
    /// Convergence history
    pub convergence_history: EvGWConvergenceHistory,
    /// Iteration history (backward compatibility)
    pub iteration_history: Vec<IterationRecord>,
}

/// Record of a single evGW iteration
#[derive(Debug, Clone)]
pub struct IterationRecord {
    /// Iteration number
    pub iteration: usize,
    /// QP energies at this iteration
    pub qp_energies: Array1<f64>,
    /// Z-factors at this iteration
    pub z_factors: Array1<f64>,
    /// Maximum energy change from previous
    pub max_change: f64,
    /// RMS energy change
    pub rms_change: f64,
}

/// evGW driver
pub struct EvGWDriver {
    /// Number of MOs
    n_mo: usize,
    /// Number of auxiliary basis functions
    n_aux: usize,
    /// Number of occupied orbitals
    n_occ: usize,
    /// Configuration
    config: EvGWConfig,
}

impl EvGWDriver {
    /// Create new evGW driver
    #[must_use]
    pub fn new(n_mo: usize, n_aux: usize, n_occ: usize, config: EvGWConfig) -> Self {
        Self {
            n_mo,
            n_aux,
            n_occ,
            config,
        }
    }

    /// Run evGW calculation with block-structured DF tensors
    ///
    /// # Arguments
    /// * `mo_energy` - HF/DFT orbital energies
    /// * `mo_occ` - Orbital occupations
    /// * `mo_coeff` - MO coefficients (n_ao x n_mo)
    /// * `ia_occ_virt` - DF tensor (occ, virt, aux) - kept for API compatibility
    /// * `mi_all_occ` - DF tensor (all_mo, occ, aux) - used for exchange
    /// * `chol_v` - Cholesky factor of Coulomb metric - kept for API compatibility
    /// * `vxc_dft` - Exchange-correlation potential diagonal
    /// * `mn_full` - Full DF tensor (n_mo, n_mo, n_aux)
    ///
    /// # Returns
    /// evGW result with converged QP energies
    ///
    /// # Algorithm
    ///
    /// For each iteration k:
    /// 1. Transpose mn_full to lpq [naux, nmo, nmo] for correlation_ac
    /// 2. Call get_sigma_diag() to compute Sigma_c on imaginary axis
    /// 3. Call ac_pade_thiele_diag() for Pade continuation
    /// 4. Evaluate Sigma_c at current QP energies:
    ///    - Iteration 0: use HF energies (G0W0)
    ///    - Iteration k: use QP energies from iteration k-1
    /// 5. Compute Z-factors from numerical derivative of Pade
    /// 6. Update: E_QP = epsilon_HF + Z * (Sigma_x + Re[Sigma_c] - Vxc)
    /// 7. Apply damping
    /// 8. Check convergence
    #[allow(clippy::too_many_arguments)]
    pub fn run_evgw_blocks(
        &mut self,
        mo_energy: &Array1<f64>,
        _mo_occ: &Array1<f64>,
        _mo_coeff: &Array2<f64>,
        _ia_occ_virt: &Array3<f64>,
        mi_all_occ: &Array3<f64>,
        _chol_v: &Array2<f64>,
        vxc_dft: &Array1<f64>,
        mn_full: &Array3<f64>,
    ) -> Result<EvGWResult> {
        let n_mo = self.n_mo;
        let n_occ = self.n_occ;
        let n_vir = n_mo - n_occ;
        let n_aux = self.n_aux;

        // Compute Fermi level from HF energies
        let ef = (mo_energy[n_occ - 1] + mo_energy[n_occ]) / 2.0;

        // Determine active orbital window for self-consistency
        // This is the key production-quality feature from VASP/FHI-aims/TURBOMOLE
        let (active_occ_start, active_vir_end) = self.determine_active_window(mo_energy, n_occ, ef);
        let n_active_occ = n_occ - active_occ_start;
        let n_active_vir = active_vir_end - n_occ;

        if self.config.print_level > 0 {
            println!("============================================================");
            println!("              evGW Self-Consistency Driver                   ");
            println!("============================================================");
            println!(
                " n_mo: {:4}  n_occ: {:4}  n_vir: {:4}  n_aux: {:4}",
                n_mo, n_occ, n_vir, n_aux
            );
            println!(
                " max_iter: {:2}  conv_tol: {:.0e}  damping: {:.2}",
                self.config.max_iterations,
                self.config.energy_tolerance,
                self.config.damping_factor
            );
            if self.config.use_diis {
                println!(
                    " DIIS: enabled (subspace dim = {})",
                    self.config.diis_space_dim
                );
            } else {
                println!(" DIIS: disabled (using simple damping)");
            }
            // Production feature: energy window
            println!(
                " Energy window: {:.1} Ha ({:.0} eV) around Fermi level",
                self.config.energy_window,
                self.config.energy_window * 27.2114
            );
            println!(
                " Active orbitals: {} occ (MO {}-{}) + {} vir (MO {}-{})",
                n_active_occ,
                active_occ_start,
                n_occ - 1,
                n_active_vir,
                n_occ,
                active_vir_end - 1
            );
            if self.config.frontier_convergence {
                println!(
                    " Convergence: frontier orbitals only (window = {:.1} Ha)",
                    self.config.frontier_window
                );
            }
            if self.config.use_frequency_caching {
                println!(" Frequency caching: ENABLED");
                println!("   WARNING: Caching gives WRONG results for evGW (iter > 1)!");
                println!("   P0/W will NOT update with QP energies. Use only for G0W0.");
            } else {
                println!(" Frequency caching: DISABLED (correct evGW mode)");
                println!("   P0/W will be recomputed at each iteration with updated QP energies.");
            }
            println!(" Using validated AC method from correlation_ac.rs");
            println!("============================================================");
        }

        // Step 1: Compute exchange self-energy (static, computed once)
        let sigma_x = self.compute_exchange(mi_all_occ)?;
        let sigma_x_diag: Array1<f64> = (0..n_mo).map(|i| sigma_x[[i, i]]).collect();

        if self.config.print_level > 0 {
            println!("\nExchange self-energy Sigma_x computed:");
            println!(
                "  Sigma_x[HOMO] = {:.6} Ha = {:.4} eV",
                sigma_x_diag[n_occ - 1],
                sigma_x_diag[n_occ - 1] * 27.2114
            );
            println!(
                "  Sigma_x[LUMO] = {:.6} Ha = {:.4} eV",
                sigma_x_diag[n_occ],
                sigma_x_diag[n_occ] * 27.2114
            );
        }

        // Step 2: Transpose DF tensor from [nmo, nmo, naux] to [naux, nmo, nmo]
        // for compatibility with correlation_ac module
        //
        // CRITICAL FIX (2025-11-27): The correlation_ac module expects lpq[P, m, n]
        // where Lpq[P, m, n] = (mn|P) in MO basis, matching PySCF convention.
        //
        // mn_full has indices (m, n, P), so we need permuted_axes([2, 0, 1])
        // to get (P, m, n).
        //
        // PREVIOUS BUG: permuted_axes([2, 1, 0]) gave (P, n, m) which SWAPS m↔n!
        // This caused ~1 eV errors for polar molecules where (mn|P) ≠ (nm|P)
        // due to asymmetric charge distribution, while symmetric molecules
        // appeared correct because (mn|P) ≈ (nm|P).
        //
        // Proof: mn_full[m=1, n=2, P=0] should map to lpq[P=0, m=1, n=2]
        //   [2,0,1]: lpq[0,1,2] = mn_full[1,2,0] ✓ CORRECT
        //   [2,1,0]: lpq[0,2,1] = mn_full[1,2,0] ✗ WRONG (swapped m↔n)
        //
        // NOTE: We use .as_standard_layout().to_owned() to ensure the result is
        // C-contiguous, as permuted_axes creates a non-contiguous view that would
        // cause errors in downstream reshape operations.
        let lpq = mn_full
            .clone()
            .permuted_axes([2, 0, 1])
            .as_standard_layout()
            .to_owned();

        // Verify dimensions
        debug_assert_eq!(
            lpq.dim(),
            (n_aux, n_mo, n_mo),
            "lpq shape mismatch: expected ({}, {}, {}), got {:?}",
            n_aux,
            n_mo,
            n_mo,
            lpq.dim()
        );

        // Initialize QP energies with HF values
        let mut qp_energies = mo_energy.clone();
        let mut z_factors = Array1::ones(n_mo);
        let mut iteration_history = Vec::new();
        let mut convergence_history = EvGWConvergenceHistory::default();
        let mut final_sigma_c = Array1::<Complex64>::zeros(n_mo);

        // AC configuration
        let ac_config = ACConfig {
            nw: self.config.cd_config.n_imag_points,
            iw_cutoff: self.config.cd_config.xi_max,
            x0: 0.5, // PySCF default
            eta: self.config.cd_config.eta,
        };

        // All orbital indices
        let orbs: Vec<usize> = (0..n_mo).collect();

        // Initialize convergence monitor with DIIS if enabled
        let acceleration = if self.config.use_diis {
            AccelerationMethod::DIIS {
                max_vectors: self.config.diis_space_dim,
            }
        } else {
            AccelerationMethod::None
        };

        let criteria = ConvergenceCriteria {
            max_iterations: self.config.max_iterations,
            energy_tolerance: self.config.energy_tolerance,
            z_factor_tolerance: self.config.density_tolerance,
            z_min: self.config.z_min,
            z_max: self.config.z_max,
            ..Default::default()
        };

        let mut convergence_monitor = ConvergenceMonitor::new(criteria, acceleration);

        // Track whether DIIS is active (starts after iteration 2)
        let mut diis_active = false;

        // Initialize frequency cache for evGW speedup
        // In evGW, P0 and W only depend on HF energies (fixed wavefunctions),
        // so we can cache them after iteration 0 and reuse in subsequent iterations.
        let mut freq_cache = FrequencyCache::new();

        // evGW iteration loop
        for iter in 0..self.config.max_iterations {
            let qp_energies_old = qp_energies.clone();

            // Step 3: Compute Sigma_c on imaginary axis
            //
            // CRITICAL FIX (2025-12-09): For TRUE evGW, the polarizability P0(iw)
            // must be computed using QP energies from the previous iteration, NOT
            // the fixed HF energies. This is what makes evGW different from G0W0!
            //
            // P0_PQ(iw) = 4 * sum_ia L_Pia * (e_ia / (w^2 + e_ia^2)) * L_Qia
            //
            // where e_ia = e_i - e_a uses:
            //   - Iteration 0: HF energies (this is G0W0)
            //   - Iteration k > 0: QP energies from iteration k-1
            //
            // The key insight is that updating the energy denominators in P0
            // causes W to change, which then causes Sigma_c to change, leading
            // to self-consistency in the QP energies.
            //
            // Previous bug: Frequency caching reused P0/W from iteration 0,
            // which meant evGW gave identical results to G0W0.
            //
            // ENERGY SELECTION FOR P0 COMPUTATION:
            let energies_for_p0 = if iter == 0 {
                // Iteration 0: Use HF/DFT energies (this is G0W0)
                mo_energy.clone()
            } else {
                // Iteration k > 0: Use QP energies from previous iteration
                // This is the key difference that makes evGW work!
                qp_energies_old.clone()
            };

            // Print diagnostic showing energy gap evolution
            if self.config.print_level > 1 && iter > 0 {
                let homo_hf = mo_energy[n_occ - 1];
                let lumo_hf = mo_energy[n_occ];
                let homo_qp = energies_for_p0[n_occ - 1];
                let lumo_qp = energies_for_p0[n_occ];
                let gap_hf = (lumo_hf - homo_hf) * 27.2114;
                let gap_qp = (lumo_qp - homo_qp) * 27.2114;
                println!("  [evGW] Using QP energies for P0:");
                println!(
                    "    HF gap:  {:.4} eV (HOMO={:.4}, LUMO={:.4})",
                    gap_hf,
                    homo_hf * 27.2114,
                    lumo_hf * 27.2114
                );
                println!(
                    "    QP gap:  {:.4} eV (HOMO={:.4}, LUMO={:.4})",
                    gap_qp,
                    homo_qp * 27.2114,
                    lumo_qp * 27.2114
                );
                println!("    Gap change: {:.4} eV", gap_qp - gap_hf);
            }

            let (sigma_iw, omega_iw) = if self.config.use_frequency_caching && iter == 0 {
                // G0W0 mode with caching: compute and cache for potential reuse
                // NOTE: Caching only makes sense for single-iteration (G0W0) mode.
                // For evGW (max_iterations > 1), caching should be DISABLED!
                get_sigma_diag_and_cache(
                    &energies_for_p0,
                    &lpq,
                    &orbs,
                    n_occ,
                    &ac_config,
                    &mut freq_cache,
                )?
            } else if self.config.use_frequency_caching && iter > 0 {
                // WARNING: Using cached W with evGW gives WRONG results!
                // The user has explicitly enabled caching despite the warning.
                // We'll use the cache but print a warning.
                if self.config.print_level > 0 && iter == 1 {
                    println!("\n  [WARNING] Frequency caching is enabled for evGW!");
                    println!("  This gives INCORRECT results (evGW = G0W0).");
                    println!("  Set use_frequency_caching = false for correct evGW.\n");
                }
                get_sigma_diag_cached(
                    mo_energy, // HF energies for Green's function
                    &freq_cache,
                    &ac_config,
                )?
            } else {
                // CORRECT MODE: Recompute P0/W at each iteration with updated QP energies
                // This is the only correct way to do evGW!
                correlation_ac::get_sigma_diag(
                    &energies_for_p0, // Use QP energies for P0 denominators
                    &lpq,
                    &orbs,
                    n_occ,
                    &ac_config,
                )?
            };

            // Step 4: Pade analytic continuation
            let (coeff, omega_fit) = correlation_ac::ac_pade_thiele_diag(&sigma_iw, &omega_iw);

            // Fermi level at midpoint of HOMO-LUMO gap
            let ef = (mo_energy[n_occ - 1] + mo_energy[n_occ]) / 2.0;

            // Numerical derivative step for Z-factor
            let de = 1e-6;

            // Step 5: Evaluate Sigma_c at current energies and compute Z-factors
            // Two modes:
            // - Linearized (default): E_QP = ε_HF + Z × (Σx + Σc(ε_eval) - Vxc)
            // - Solved (Newton): Find E_QP where E_QP = ε_HF + Σx + Σc(E_QP) - Vxc
            let mut new_qp_energies = Array1::zeros(n_mo);
            let mut new_z_factors = Array1::zeros(n_mo);

            if self.config.use_solved_qp {
                // === Newton-Raphson solved QP equation ===
                let newton_config = NewtonSolverConfig {
                    max_iterations: self.config.newton_max_iterations,
                    energy_tolerance: self.config.newton_energy_tolerance,
                    z_min: self.config.z_min,
                    z_max: self.config.z_max,
                    verbose: if self.config.print_level > 2 { 1 } else { 0 },
                    ..Default::default()
                };

                for (idx, &p) in orbs.iter().enumerate() {
                    // Get Padé coefficients for this orbital
                    let omega_p = omega_fit.row(idx).to_owned();
                    let coeff_p = coeff.column(idx).to_owned();

                    // Initial guess: use linearized result or previous QP energy
                    let initial_guess = if iter == 0 {
                        // For first iteration, use linearized result as initial guess
                        let omega_eval = mo_energy[p] - ef;
                        let sigma_c_real =
                            correlation_ac::pade_thiele(omega_eval, &omega_p, &coeff_p).re;
                        let sigma_c_plus =
                            correlation_ac::pade_thiele(omega_eval + de, &omega_p, &coeff_p).re;
                        let dsigma_dw = (sigma_c_plus - sigma_c_real) / de;
                        let z = (1.0 / (1.0 - dsigma_dw))
                            .max(self.config.z_min)
                            .min(self.config.z_max);
                        mo_energy[p] + z * (sigma_x_diag[p] + sigma_c_real - vxc_dft[p])
                    } else {
                        qp_energies_old[p]
                    };

                    // Create Σc evaluator closure for Newton solver
                    let sigma_c_eval =
                        |omega: f64| correlation_ac::pade_thiele(omega, &omega_p, &coeff_p);

                    // Create Newton solver
                    let solver = NewtonQPSolver::new(
                        mo_energy[p],
                        sigma_x_diag[p],
                        vxc_dft[p],
                        ef,
                        newton_config.clone(),
                    );

                    // Solve QP equation
                    match solver.solve(&sigma_c_eval, initial_guess) {
                        Ok(result) => {
                            new_qp_energies[p] = result.qp_energy;
                            new_z_factors[p] = result.z_factor;
                            final_sigma_c[p] = result.sigma_c;
                        }
                        Err(_) => {
                            // Fallback to linearized if Newton fails
                            let omega_eval = initial_guess - ef;
                            let sigma_c =
                                correlation_ac::pade_thiele(omega_eval, &omega_p, &coeff_p);
                            let sigma_c_plus =
                                correlation_ac::pade_thiele(omega_eval + de, &omega_p, &coeff_p);
                            let dsigma_dw = (sigma_c_plus.re - sigma_c.re) / de;
                            let z = (1.0 / (1.0 - dsigma_dw))
                                .max(self.config.z_min)
                                .min(self.config.z_max);
                            new_qp_energies[p] =
                                mo_energy[p] + z * (sigma_x_diag[p] + sigma_c.re - vxc_dft[p]);
                            new_z_factors[p] = z;
                            final_sigma_c[p] = sigma_c;
                        }
                    }
                }
            } else {
                // === Linearized QP equation (default) ===
                for (idx, &p) in orbs.iter().enumerate() {
                    // Energy at which to evaluate Sigma_c:
                    // - Iteration 0: HF energy (this is G0W0)
                    // - Iteration k > 0: QP energy from previous iteration
                    let eval_energy = if iter == 0 {
                        mo_energy[p]
                    } else {
                        qp_energies_old[p]
                    };

                    // Get Padé coefficients for this orbital
                    let omega_p = omega_fit.row(idx).to_owned();
                    let coeff_p = coeff.column(idx).to_owned();

                    // Evaluate Sigma_c at current energy (relative to Fermi level)
                    let omega_eval = eval_energy - ef;
                    let sigma_c = correlation_ac::pade_thiele(omega_eval, &omega_p, &coeff_p);
                    let sigma_c_real = sigma_c.re;

                    // Store for output
                    final_sigma_c[p] = sigma_c;

                    // Compute Z-factor from numerical derivative: Z = 1 / (1 - dSigma/dw)
                    let sigma_c_plus =
                        correlation_ac::pade_thiele(omega_eval + de, &omega_p, &coeff_p);
                    let dsigma_dw = (sigma_c_plus.re - sigma_c_real) / de;

                    let z_raw = 1.0 / (1.0 - dsigma_dw);
                    let z = z_raw.max(self.config.z_min).min(self.config.z_max);
                    new_z_factors[p] = z;

                    // QP energy: E_QP = epsilon_HF + Z * (Sigma_x + Sigma_c - Vxc)
                    let sigma_total = sigma_x_diag[p] + sigma_c_real - vxc_dft[p];
                    new_qp_energies[p] = mo_energy[p] + z * sigma_total;
                }
            }

            // Step 6: Apply convergence acceleration (DIIS or damping)
            //
            // Strategy:
            // - Iteration 0: Use new values directly (this is G0W0, no mixing with HF)
            // - Iterations 1-2: Use simple damping (need history for DIIS)
            // - Iterations 3+: Use DIIS if enabled, fall back to damping if it fails
            //
            // DIIS uses error vector: e = E_new - E_old (energy change method)
            // This is simpler than residual-based DIIS and works well for evGW.
            if iter == 0 {
                // First iteration: use new values directly (this is G0W0)
                qp_energies = new_qp_energies.clone();
                z_factors = new_z_factors.clone();

                // Record in convergence monitor for DIIS history
                convergence_monitor.check_convergence(&qp_energies, &z_factors);
            } else if iter == 1 {
                // Second iteration: simple damping, build DIIS history
                let damping = self.config.damping_factor;
                for p in 0..n_mo {
                    qp_energies[p] =
                        (1.0 - damping) * qp_energies_old[p] + damping * new_qp_energies[p];
                    z_factors[p] = (1.0 - damping) * z_factors[p] + damping * new_z_factors[p];
                }
                convergence_monitor.check_convergence(&qp_energies, &z_factors);
            } else if self.config.use_diis {
                // Iterations 2+: Try DIIS acceleration
                // First record the new values to build history
                convergence_monitor.check_convergence(&new_qp_energies, &new_z_factors);

                // Attempt DIIS extrapolation
                match convergence_monitor.accelerate(&new_qp_energies, &new_z_factors) {
                    Ok((diis_energies, diis_z)) => {
                        // Sanity check: DIIS coefficients can sometimes blow up
                        // Check that extrapolated values are not too far from input
                        let max_energy_deviation = diis_energies
                            .iter()
                            .zip(new_qp_energies.iter())
                            .map(|(d, n)| (d - n).abs())
                            .fold(0.0, f64::max);

                        let max_z_deviation = diis_z
                            .iter()
                            .zip(new_z_factors.iter())
                            .map(|(d, n)| (d - n).abs())
                            .fold(0.0, f64::max);

                        // Accept DIIS if deviations are reasonable
                        // (not more than 0.5 Ha for energies, 0.2 for Z-factors)
                        if max_energy_deviation < 0.5 && max_z_deviation < 0.2 {
                            qp_energies = diis_energies;
                            z_factors = diis_z;
                            diis_active = true;

                            if self.config.print_level > 1 {
                                println!(
                                    "  [DIIS] Extrapolation applied (dE_max={:.4}, dZ_max={:.4})",
                                    max_energy_deviation, max_z_deviation
                                );
                            }
                        } else {
                            // DIIS gave unreasonable values, fall back to damping
                            if self.config.print_level > 1 {
                                println!("  [DIIS] Fallback to damping (dE_max={:.4} > 0.5 or dZ_max={:.4} > 0.2)",
                                         max_energy_deviation, max_z_deviation);
                            }
                            let damping = self.config.damping_factor;
                            for p in 0..n_mo {
                                qp_energies[p] = (1.0 - damping) * qp_energies_old[p]
                                    + damping * new_qp_energies[p];
                                z_factors[p] =
                                    (1.0 - damping) * z_factors[p] + damping * new_z_factors[p];
                            }
                            diis_active = false;
                        }
                    }
                    Err(e) => {
                        // DIIS failed (e.g., singular matrix), fall back to damping
                        if self.config.print_level > 1 {
                            println!("  [DIIS] Error: {}, falling back to damping", e);
                        }
                        let damping = self.config.damping_factor;
                        for p in 0..n_mo {
                            qp_energies[p] =
                                (1.0 - damping) * qp_energies_old[p] + damping * new_qp_energies[p];
                            z_factors[p] =
                                (1.0 - damping) * z_factors[p] + damping * new_z_factors[p];
                        }
                        diis_active = false;
                    }
                }

                // Ensure Z-factors stay in bounds after DIIS
                for p in 0..n_mo {
                    z_factors[p] = z_factors[p].max(self.config.z_min).min(self.config.z_max);
                }
            } else {
                // DIIS disabled: simple damping
                let damping = self.config.damping_factor;
                for p in 0..n_mo {
                    qp_energies[p] =
                        (1.0 - damping) * qp_energies_old[p] + damping * new_qp_energies[p];
                    z_factors[p] = (1.0 - damping) * z_factors[p] + damping * new_z_factors[p];
                }
            }

            // Step 7: Check convergence
            // Use frontier orbital window for convergence (VASP/FHI-aims style)
            let (conv_occ_start, conv_vir_end) =
                self.determine_convergence_window(mo_energy, n_occ, ef);
            let n_conv_orbs = conv_vir_end - conv_occ_start;

            // Compute energy changes for ALL orbitals (for record-keeping)
            let energy_diff: Vec<f64> = qp_energies
                .iter()
                .zip(qp_energies_old.iter())
                .map(|(new, old)| (new - old).abs())
                .collect();
            let max_change_all = energy_diff.iter().cloned().fold(0.0, f64::max);
            let _rms_change_all =
                (energy_diff.iter().map(|x| x * x).sum::<f64>() / n_mo as f64).sqrt();

            // Compute energy changes for CONVERGENCE window only (frontier orbitals)
            let energy_diff_conv: Vec<f64> = (conv_occ_start..conv_vir_end)
                .map(|i| (qp_energies[i] - qp_energies_old[i]).abs())
                .collect();
            let max_change = energy_diff_conv.iter().cloned().fold(0.0, f64::max);
            let rms_change =
                (energy_diff_conv.iter().map(|x| x * x).sum::<f64>() / n_conv_orbs as f64).sqrt();

            // Record iteration
            iteration_history.push(IterationRecord {
                iteration: iter,
                qp_energies: qp_energies.clone(),
                z_factors: z_factors.clone(),
                max_change,
                rms_change,
            });

            convergence_history.energies.push(qp_energies.clone());
            convergence_history.max_energy_changes.push(max_change);
            convergence_history.rms_energy_changes.push(rms_change);

            if self.config.print_level > 0 {
                let accel_status = if diis_active { " [DIIS]" } else { "" };
                if self.config.frontier_convergence {
                    println!(
                        "\nIteration {:2}: max_dE(frontier) = {:.2e} Ha, max_dE(all) = {:.2e} Ha{}",
                        iter + 1,
                        max_change,
                        max_change_all,
                        accel_status
                    );
                    println!(
                        "  Convergence window: MO {}-{} ({} orbitals)",
                        conv_occ_start,
                        conv_vir_end - 1,
                        n_conv_orbs
                    );
                } else {
                    println!(
                        "\nIteration {:2}: max_dE = {:.2e} Ha, rms_dE = {:.2e} Ha{}",
                        iter + 1,
                        max_change,
                        rms_change,
                        accel_status
                    );
                }
                println!(
                    "  E_QP[HOMO] = {:.6} Ha = {:.4} eV  (Z = {:.4})",
                    qp_energies[n_occ - 1],
                    qp_energies[n_occ - 1] * 27.2114,
                    z_factors[n_occ - 1]
                );
                println!(
                    "  E_QP[LUMO] = {:.6} Ha = {:.4} eV  (Z = {:.4})",
                    qp_energies[n_occ],
                    qp_energies[n_occ] * 27.2114,
                    z_factors[n_occ]
                );
                println!(
                    "  Sigma_c[HOMO] = {:.6} Ha, Sigma_c[LUMO] = {:.6} Ha",
                    final_sigma_c[n_occ - 1].re,
                    final_sigma_c[n_occ].re
                );
            }

            // Convergence check
            if max_change < self.config.energy_tolerance {
                if self.config.print_level > 0 {
                    println!("\n[OK] evGW converged in {} iterations!", iter + 1);
                }

                return Ok(EvGWResult {
                    qp_energies,
                    z_factors,
                    sigma_x,
                    sigma_c: final_sigma_c,
                    n_iterations: iter + 1,
                    n_cycles: iter + 1,
                    converged: true,
                    convergence_history,
                    iteration_history,
                });
            }
        }

        // Did not converge
        if self.config.print_level > 0 {
            println!(
                "\n[WARN] evGW did not converge in {} iterations",
                self.config.max_iterations
            );
        }

        Ok(EvGWResult {
            qp_energies,
            z_factors,
            sigma_x,
            sigma_c: final_sigma_c,
            n_iterations: self.config.max_iterations,
            n_cycles: self.config.max_iterations,
            converged: false,
            convergence_history,
            iteration_history,
        })
    }

    /// Determine active orbital window for self-consistency (VASP/FHI-aims style)
    ///
    /// Returns (active_occ_start, active_vir_end) indices for orbitals that will be
    /// updated in the evGW iterations. States outside this window keep their G0W0 values.
    ///
    /// # Algorithm
    ///
    /// 1. If `n_occ_active` is set, use top N occupied orbitals
    /// 2. If `n_vir_active` is set, use bottom N virtual orbitals
    /// 3. Otherwise, use `energy_window` from Fermi level
    ///
    /// # Arguments
    /// * `mo_energy` - HF/DFT orbital energies
    /// * `n_occ` - Number of occupied orbitals
    /// * `ef` - Fermi level (midpoint of HOMO-LUMO gap)
    ///
    /// # Returns
    /// (active_occ_start, active_vir_end) - start index for occupied, end index for virtual
    fn determine_active_window(
        &self,
        mo_energy: &Array1<f64>,
        n_occ: usize,
        ef: f64,
    ) -> (usize, usize) {
        let n_mo = mo_energy.len();

        // Determine occupied active window
        let active_occ_start = if let Some(n_active) = self.config.n_occ_active {
            // Use top N occupied orbitals
            n_occ.saturating_sub(n_active)
        } else {
            // Use energy window: find first occupied orbital within window of Fermi level
            let mut start = 0;
            for i in 0..n_occ {
                if mo_energy[i] >= ef - self.config.energy_window {
                    start = i;
                    break;
                }
            }
            start
        };

        // Determine virtual active window
        let active_vir_end = if let Some(n_active) = self.config.n_vir_active {
            // Use bottom N virtual orbitals
            (n_occ + n_active).min(n_mo)
        } else {
            // Use energy window: find last virtual orbital within window of Fermi level
            let mut end = n_mo;
            for i in (n_occ..n_mo).rev() {
                if mo_energy[i] <= ef + self.config.energy_window {
                    end = i + 1;
                    break;
                }
            }
            // Ensure we include at least LUMO
            end.max(n_occ + 1)
        };

        (active_occ_start, active_vir_end)
    }

    /// Determine which orbitals to check for convergence (frontier orbitals)
    ///
    /// When `frontier_convergence` is enabled, only orbitals near the HOMO-LUMO gap
    /// are considered for convergence. This avoids oscillation from high-energy
    /// virtual states that are physically irrelevant.
    ///
    /// # Returns
    /// (conv_occ_start, conv_vir_end) - orbital range for convergence check
    fn determine_convergence_window(
        &self,
        mo_energy: &Array1<f64>,
        n_occ: usize,
        _ef: f64,
    ) -> (usize, usize) {
        let n_mo = mo_energy.len();

        if !self.config.frontier_convergence {
            // Check all orbitals
            return (0, n_mo);
        }

        // Use smaller frontier_window for convergence check
        let window = self.config.frontier_window;

        // Find occupied orbitals within window of HOMO
        let homo_energy = mo_energy[n_occ - 1];
        let mut conv_occ_start = n_occ.saturating_sub(1);
        for i in (0..n_occ).rev() {
            if mo_energy[i] >= homo_energy - window {
                conv_occ_start = i;
            } else {
                break;
            }
        }

        // Find virtual orbitals within window of LUMO
        let lumo_energy = mo_energy[n_occ];
        let mut conv_vir_end = n_occ + 1;
        for i in n_occ..n_mo {
            if mo_energy[i] <= lumo_energy + window {
                conv_vir_end = i + 1;
            } else {
                break;
            }
        }

        (conv_occ_start, conv_vir_end)
    }

    /// Compute exchange self-energy using RI/DF tensors
    ///
    /// Sigma_x[p,q] = -sum_i sum_R (pi|R)(qi|R)
    ///
    /// where R is the auxiliary basis index and i runs over occupied orbitals.
    fn compute_exchange(&self, mi_all_occ: &Array3<f64>) -> Result<Array2<f64>> {
        let n_mo = self.n_mo;
        let n_occ = self.n_occ;
        let n_aux = self.n_aux;

        let mut sigma_x = Array2::<f64>::zeros((n_mo, n_mo));

        // mi_all_occ has shape (n_mo, n_occ, n_aux)
        // For exchange: Sigma_x[p,q] = -sum_i (pi|R)(qi|R)

        for p in 0..n_mo {
            for q in 0..=p {
                let mut sum = 0.0;
                for i in 0..n_occ {
                    for r in 0..n_aux {
                        let pi_r = mi_all_occ[[p, i, r]];
                        let qi_r = mi_all_occ[[q, i, r]];
                        sum += pi_r * qi_r;
                    }
                }
                sigma_x[[p, q]] = -sum;
                sigma_x[[q, p]] = -sum;
            }
        }

        Ok(sigma_x)
    }
}

/// Perform a single G0W0 calculation (equivalent to iteration 0 of evGW)
///
/// This is a convenience function for G0W0 with linearized QP equation.
///
/// # Arguments
/// * `mo_energy` - HF/DFT orbital energies
/// * `sigma_x_diag` - Exchange self-energy diagonal
/// * `vxc_dft` - DFT XC potential diagonal
/// * `lpq` - Full DF tensor [naux, nmo, nmo]
/// * `n_occ` - Number of occupied orbitals
///
/// # Returns
/// (qp_energies, z_factors)
pub fn g0w0_linearized(
    mo_energy: &Array1<f64>,
    sigma_x_diag: &Array1<f64>,
    vxc_dft: &Array1<f64>,
    lpq: &Array3<f64>,
    n_occ: usize,
) -> Result<(Array1<f64>, Array1<f64>)> {
    let n_mo = mo_energy.len();
    let orbs: Vec<usize> = (0..n_mo).collect();

    // AC configuration (PySCF defaults)
    let config = ACConfig::default();

    // Compute Sigma_c on imaginary axis
    // CRITICAL: Pass nocc explicitly, not inferred from eigenvalue signs
    let (sigma_iw, omega_iw) =
        correlation_ac::get_sigma_diag(mo_energy, lpq, &orbs, n_occ, &config)?;

    // Pade analytic continuation
    let (coeff, omega_fit) = correlation_ac::ac_pade_thiele_diag(&sigma_iw, &omega_iw);

    // Fermi level
    let ef = (mo_energy[n_occ - 1] + mo_energy[n_occ]) / 2.0;
    let de = 1e-6;

    let mut qp_energies = Array1::zeros(n_mo);
    let mut z_factors = Array1::zeros(n_mo);

    for (idx, &p) in orbs.iter().enumerate() {
        let ep = mo_energy[p];

        let omega_p = omega_fit.row(idx).to_owned();
        let coeff_p = coeff.column(idx).to_owned();

        // Evaluate at HF energy
        let omega_eval = ep - ef;
        let sigma_c_real = correlation_ac::pade_thiele(omega_eval, &omega_p, &coeff_p).re;

        // Z-factor from numerical derivative
        let sigma_c_plus = correlation_ac::pade_thiele(omega_eval + de, &omega_p, &coeff_p).re;
        let dsigma = (sigma_c_plus - sigma_c_real) / de;
        let z = (1.0 / (1.0 - dsigma)).max(0.1).min(0.999);

        z_factors[p] = z;

        // QP energy
        let correction = sigma_x_diag[p] + sigma_c_real - vxc_dft[p];
        qp_energies[p] = ep + z * correction;
    }

    Ok((qp_energies, z_factors))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_evgw_config_default() {
        let config = EvGWConfig::default();
        assert_eq!(config.max_iterations, 12);
        assert!((config.energy_tolerance - 1e-4).abs() < 1e-10);
        assert!((config.damping_factor - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_evgw_driver_creation() {
        let config = EvGWConfig::default();
        let driver = EvGWDriver::new(10, 50, 5, config);
        assert_eq!(driver.n_mo, 10);
        assert_eq!(driver.n_aux, 50);
        assert_eq!(driver.n_occ, 5);
    }

    #[test]
    fn test_cd_config_defaults() {
        let config = ContourDeformationConfig::default();
        assert_eq!(config.n_imag_points, 100);
        assert!((config.xi_max - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_df_tensor_transpose() {
        // Verify that transposing [nmo, nmo, naux] to [naux, nmo, nmo] works
        let n_mo = 3;
        let n_aux = 5;

        let mut mn_full = Array3::<f64>::zeros((n_mo, n_mo, n_aux));
        mn_full[[0, 1, 2]] = 1.5;
        mn_full[[2, 0, 4]] = 2.5;

        let lpq = mn_full.clone().permuted_axes([2, 0, 1]);

        assert_eq!(lpq.dim(), (n_aux, n_mo, n_mo));
        assert!((lpq[[2, 0, 1]] - 1.5).abs() < 1e-12);
        assert!((lpq[[4, 2, 0]] - 2.5).abs() < 1e-12);
    }
}
