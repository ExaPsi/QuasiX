//! Newton-Raphson QP Equation Solver
//!
//! This module implements the solved (non-linearized) QP equation using
//! Newton-Raphson iteration with robust convergence strategies.
//!
//! # Physical Background
//!
//! The quasiparticle (QP) equation is:
//!
//! ```text
//! E_QP = ε_HF + Σx + Re[Σc(E_QP)] - Vxc
//! ```
//!
//! Unlike the linearized approximation which evaluates Σc at ε_HF and uses
//! a Z-factor correction, the solved QP equation self-consistently finds
//! the energy where the equation is satisfied.
//!
//! # Algorithm
//!
//! Newton-Raphson iteration:
//!
//! ```text
//! f(E) = E - ε_HF - Σx - Re[Σc(E)] + Vxc = 0
//! f'(E) = 1 - dΣc/dE
//! E_{n+1} = E_n - f(E_n) / f'(E_n)
//! ```
//!
//! With damping to ensure convergence:
//!
//! ```text
//! E_{n+1} = E_n - α × f(E_n) / f'(E_n)
//! ```
//!
//! where α is adaptively adjusted based on convergence behavior.
//!
//! # Advantages over Linearized QP
//!
//! 1. **More accurate for strongly correlated systems**: When Z-factor is small
//!    (strong correlation), linearization error can be ~50-150 meV
//! 2. **Required for converged IPs**: For comparison with high-accuracy codes
//!    like TURBOMOLE, solved QP is essential
//! 3. **Systematic improvability**: Can add higher-order corrections
//!
//! # References
//!
//! - F. Bruneval, GW approximation review, arXiv:1203.4684 (2012)
//! - X. Ren et al., Phys. Rev. B 88, 035120 (2013)

use crate::common::{QuasixError, Result};
use crate::selfenergy::correlation_ac::pade_thiele;
use ndarray::Array1;
use num_complex::Complex64;

/// Configuration for Newton-Raphson QP solver
#[derive(Debug, Clone)]
pub struct NewtonSolverConfig {
    /// Maximum Newton iterations per orbital
    pub max_iterations: usize,
    /// Energy convergence tolerance (Ha)
    pub energy_tolerance: f64,
    /// Initial damping factor (0 < α ≤ 1)
    pub initial_damping: f64,
    /// Minimum damping factor
    pub min_damping: f64,
    /// Maximum damping factor
    pub max_damping: f64,
    /// Maximum step size (Ha) to prevent divergence
    pub max_step: f64,
    /// Finite difference step for numerical derivative (Ha)
    pub derivative_delta: f64,
    /// Use Richardson extrapolation for derivative
    pub use_richardson: bool,
    /// Minimum Z-factor (physical bound)
    pub z_min: f64,
    /// Maximum Z-factor (physical bound)
    pub z_max: f64,
    /// Use bisection fallback if Newton fails
    pub use_bisection_fallback: bool,
    /// Maximum bisection iterations (fallback)
    pub max_bisection_iterations: usize,
    /// Verbosity level
    pub verbose: usize,
}

impl Default for NewtonSolverConfig {
    fn default() -> Self {
        Self {
            max_iterations: 50,
            energy_tolerance: 1e-8, // 0.27 µeV precision
            initial_damping: 1.0,   // Full Newton step initially
            min_damping: 0.1,
            max_damping: 1.0,
            max_step: 1.0,          // Max 1 Ha step
            derivative_delta: 1e-6, // 0.027 meV for numerical derivative
            use_richardson: true,   // Better accuracy for derivative
            z_min: 0.1,
            z_max: 0.999,
            use_bisection_fallback: true,
            max_bisection_iterations: 30,
            verbose: 0,
        }
    }
}

/// Result from Newton-Raphson QP solver
#[derive(Debug, Clone)]
pub struct NewtonSolverResult {
    /// Converged QP energy
    pub qp_energy: f64,
    /// Final residual |f(E)|
    pub residual: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Converged flag
    pub converged: bool,
    /// Final Σc(E_QP)
    pub sigma_c: Complex64,
    /// Effective Z-factor at solution (informational)
    pub z_factor: f64,
    /// Method used ("newton" or "bisection")
    pub method: String,
}

/// Function signature for evaluating Σc(ω) using Padé approximant
pub type SigmaCEvaluator = dyn Fn(f64) -> Complex64;

/// Newton-Raphson solver for single QP equation
pub struct NewtonQPSolver {
    /// Configuration
    config: NewtonSolverConfig,
    /// HF/DFT orbital energy
    epsilon_hf: f64,
    /// Exchange self-energy
    sigma_x: f64,
    /// Exchange-correlation potential
    vxc: f64,
    /// Fermi level
    fermi_level: f64,
}

impl NewtonQPSolver {
    /// Create new Newton-Raphson solver for a single orbital
    ///
    /// # Arguments
    /// * `epsilon_hf` - HF/DFT orbital energy
    /// * `sigma_x` - Exchange self-energy Σx
    /// * `vxc` - DFT XC potential Vxc
    /// * `fermi_level` - Fermi level (midpoint of HOMO-LUMO gap)
    /// * `config` - Solver configuration
    pub fn new(
        epsilon_hf: f64,
        sigma_x: f64,
        vxc: f64,
        fermi_level: f64,
        config: NewtonSolverConfig,
    ) -> Self {
        Self {
            config,
            epsilon_hf,
            sigma_x,
            vxc,
            fermi_level,
        }
    }

    /// Solve the QP equation using Newton-Raphson iteration
    ///
    /// # Arguments
    /// * `sigma_c_eval` - Function to evaluate Σc(ω) at any frequency ω
    /// * `initial_guess` - Initial guess for E_QP (typically linearized result)
    ///
    /// # Returns
    /// `NewtonSolverResult` with converged QP energy
    pub fn solve<F>(&self, sigma_c_eval: &F, initial_guess: f64) -> Result<NewtonSolverResult>
    where
        F: Fn(f64) -> Complex64 + ?Sized,
    {
        let mut energy = initial_guess;
        let mut damping = self.config.initial_damping;
        let mut prev_residual = f64::MAX;

        if self.config.verbose > 0 {
            println!("Newton-Raphson QP solver:");
            println!(
                "  ε_HF = {:.6} Ha, Σx = {:.6} Ha, Vxc = {:.6} Ha",
                self.epsilon_hf, self.sigma_x, self.vxc
            );
            println!("  Initial guess: E = {:.6} Ha", initial_guess);
        }

        for iter in 0..self.config.max_iterations {
            // Evaluate Σc at current energy (relative to Fermi level)
            let omega = energy - self.fermi_level;
            let sigma_c = sigma_c_eval(omega);
            let sigma_c_real = sigma_c.re;

            // Compute residual: f(E) = E - ε_HF - Σx - Re[Σc(E)] + Vxc
            let residual = energy - self.epsilon_hf - self.sigma_x - sigma_c_real + self.vxc;
            let abs_residual = residual.abs();

            if self.config.verbose > 1 {
                println!(
                    "  Iter {:2}: E = {:.8} Ha, Σc = {:.6} Ha, |f| = {:.2e}",
                    iter, energy, sigma_c_real, abs_residual
                );
            }

            // Check convergence
            if abs_residual < self.config.energy_tolerance {
                // Compute Z-factor at solution for informational purposes
                let z = self.compute_z_factor(sigma_c_eval, energy);

                if self.config.verbose > 0 {
                    println!(
                        "  Converged in {} iterations: E_QP = {:.6} Ha (Z = {:.4})",
                        iter + 1,
                        energy,
                        z
                    );
                }

                return Ok(NewtonSolverResult {
                    qp_energy: energy,
                    residual: abs_residual,
                    iterations: iter + 1,
                    converged: true,
                    sigma_c,
                    z_factor: z,
                    method: "newton".to_string(),
                });
            }

            // Compute derivative: f'(E) = 1 - dΣc/dE
            let dsigma_dw = self.compute_derivative(sigma_c_eval, omega);
            let fprime = 1.0 - dsigma_dw;

            // Newton step: ΔE = -f(E) / f'(E)
            let step = if fprime.abs() < 1e-10 {
                // Near zero derivative - use gradient descent
                -residual.signum() * self.config.max_step.min(abs_residual)
            } else {
                -residual / fprime
            };

            // Apply maximum step constraint
            let clipped_step = step.max(-self.config.max_step).min(self.config.max_step);

            // Adaptive damping
            if abs_residual > prev_residual {
                // Diverging - reduce damping
                damping = (damping * 0.5).max(self.config.min_damping);
            } else if abs_residual < 0.5 * prev_residual && damping < self.config.max_damping {
                // Good progress - increase damping
                damping = (damping * 1.2).min(self.config.max_damping);
            }

            // Update energy
            energy += damping * clipped_step;
            prev_residual = abs_residual;
        }

        // Newton failed to converge - try bisection fallback
        if self.config.use_bisection_fallback {
            if self.config.verbose > 0 {
                println!("  Newton did not converge, trying bisection...");
            }
            return self.bisection_fallback(sigma_c_eval, initial_guess);
        }

        // Return best result even if not converged
        let omega = energy - self.fermi_level;
        let _sigma_c = sigma_c_eval(omega);
        let _z = self.compute_z_factor(sigma_c_eval, energy);

        Err(QuasixError::ConvergenceError(format!(
            "Newton-Raphson failed to converge after {} iterations (residual = {:.2e})",
            self.config.max_iterations, prev_residual
        )))
    }

    /// Compute numerical derivative dΣc/dω using central difference
    fn compute_derivative<F>(&self, sigma_c_eval: &F, omega: f64) -> f64
    where
        F: Fn(f64) -> Complex64 + ?Sized,
    {
        let de = self.config.derivative_delta;

        if self.config.use_richardson {
            // Richardson extrapolation for 4th-order accuracy
            // f'(x) ≈ (8f(x+h) - 8f(x-h) - f(x+2h) + f(x-2h)) / (12h)
            let sigma_p1 = sigma_c_eval(omega + de).re;
            let sigma_m1 = sigma_c_eval(omega - de).re;
            let sigma_p2 = sigma_c_eval(omega + 2.0 * de).re;
            let sigma_m2 = sigma_c_eval(omega - 2.0 * de).re;

            (8.0 * (sigma_p1 - sigma_m1) - (sigma_p2 - sigma_m2)) / (12.0 * de)
        } else {
            // Central difference (2nd-order)
            let sigma_plus = sigma_c_eval(omega + de).re;
            let sigma_minus = sigma_c_eval(omega - de).re;
            (sigma_plus - sigma_minus) / (2.0 * de)
        }
    }

    /// Compute Z-factor at given energy
    fn compute_z_factor<F>(&self, sigma_c_eval: &F, energy: f64) -> f64
    where
        F: Fn(f64) -> Complex64 + ?Sized,
    {
        let omega = energy - self.fermi_level;
        let dsigma_dw = self.compute_derivative(sigma_c_eval, omega);
        let z_raw = 1.0 / (1.0 - dsigma_dw);
        z_raw.max(self.config.z_min).min(self.config.z_max)
    }

    /// Bisection fallback for robust convergence
    fn bisection_fallback<F>(
        &self,
        sigma_c_eval: &F,
        initial_guess: f64,
    ) -> Result<NewtonSolverResult>
    where
        F: Fn(f64) -> Complex64 + ?Sized,
    {
        // Define search bracket around initial guess
        let bracket_width = 2.0; // Ha
        let mut e_low = initial_guess - bracket_width;
        let mut e_high = initial_guess + bracket_width;

        // Evaluate residual at bracket bounds
        let f_low = self.residual(sigma_c_eval, e_low);
        let f_high = self.residual(sigma_c_eval, e_high);

        // Check if we have a valid bracket
        if f_low * f_high > 0.0 {
            // No sign change - expand bracket
            for expansion in 1..5 {
                let width = bracket_width * (2.0_f64).powi(expansion);
                e_low = initial_guess - width;
                e_high = initial_guess + width;
                let f_l = self.residual(sigma_c_eval, e_low);
                let f_h = self.residual(sigma_c_eval, e_high);
                if f_l * f_h <= 0.0 {
                    break;
                }
            }
        }

        // Bisection iteration
        for iter in 0..self.config.max_bisection_iterations {
            let e_mid = (e_low + e_high) / 2.0;
            let f_mid = self.residual(sigma_c_eval, e_mid);

            if f_mid.abs() < self.config.energy_tolerance
                || (e_high - e_low) / 2.0 < self.config.energy_tolerance
            {
                let omega = e_mid - self.fermi_level;
                let sigma_c = sigma_c_eval(omega);
                let z = self.compute_z_factor(sigma_c_eval, e_mid);

                if self.config.verbose > 0 {
                    println!(
                        "  Bisection converged in {} iterations: E_QP = {:.6} Ha",
                        iter + 1,
                        e_mid
                    );
                }

                return Ok(NewtonSolverResult {
                    qp_energy: e_mid,
                    residual: f_mid.abs(),
                    iterations: iter + 1,
                    converged: true,
                    sigma_c,
                    z_factor: z,
                    method: "bisection".to_string(),
                });
            }

            let f_low = self.residual(sigma_c_eval, e_low);
            if f_mid * f_low < 0.0 {
                e_high = e_mid;
            } else {
                e_low = e_mid;
            }
        }

        Err(QuasixError::ConvergenceError(
            "Bisection fallback failed to converge".to_string(),
        ))
    }

    /// Compute residual f(E) = E - ε_HF - Σx - Re[Σc(E)] + Vxc
    fn residual<F>(&self, sigma_c_eval: &F, energy: f64) -> f64
    where
        F: Fn(f64) -> Complex64 + ?Sized,
    {
        let omega = energy - self.fermi_level;
        let sigma_c_real = sigma_c_eval(omega).re;
        energy - self.epsilon_hf - self.sigma_x - sigma_c_real + self.vxc
    }
}

/// Solve QP equations for multiple orbitals
///
/// # Arguments
/// * `mo_energy` - HF/DFT orbital energies
/// * `sigma_x_diag` - Exchange self-energy diagonal
/// * `vxc_diag` - DFT XC potential diagonal
/// * `orbs` - Orbital indices to solve
/// * `nocc` - Number of occupied orbitals
/// * `sigma_c_evaluators` - Closure factory returning Σc evaluator for each orbital
/// * `initial_guesses` - Initial QP energy guesses (e.g., from linearized solver)
/// * `config` - Solver configuration
///
/// # Returns
/// Vector of `NewtonSolverResult` for each orbital
pub fn solve_qp_equations<F>(
    mo_energy: &Array1<f64>,
    sigma_x_diag: &Array1<f64>,
    vxc_diag: &Array1<f64>,
    orbs: &[usize],
    nocc: usize,
    sigma_c_evaluators: F,
    initial_guesses: &Array1<f64>,
    config: NewtonSolverConfig,
) -> Vec<Result<NewtonSolverResult>>
where
    F: Fn(usize) -> Box<dyn Fn(f64) -> Complex64>,
{
    let ef = (mo_energy[nocc - 1] + mo_energy[nocc]) / 2.0;

    orbs.iter()
        .enumerate()
        .map(|(idx, &p)| {
            let solver = NewtonQPSolver::new(
                mo_energy[p],
                sigma_x_diag[p],
                vxc_diag[p],
                ef,
                config.clone(),
            );

            let evaluator = sigma_c_evaluators(idx);
            solver.solve(&*evaluator, initial_guesses[p])
        })
        .collect()
}

/// Convenience function to solve a single QP equation
///
/// # Arguments
/// * `epsilon_hf` - HF orbital energy
/// * `sigma_x` - Exchange self-energy
/// * `vxc` - XC potential
/// * `fermi_level` - Fermi level
/// * `sigma_c_eval` - Function to evaluate Σc(ω)
/// * `initial_guess` - Initial guess for E_QP
/// * `config` - Optional solver configuration
///
/// # Returns
/// Converged QP energy or error
pub fn solve_single_qp<F>(
    epsilon_hf: f64,
    sigma_x: f64,
    vxc: f64,
    fermi_level: f64,
    sigma_c_eval: F,
    initial_guess: f64,
    config: Option<NewtonSolverConfig>,
) -> Result<NewtonSolverResult>
where
    F: Fn(f64) -> Complex64,
{
    let cfg = config.unwrap_or_default();
    let solver = NewtonQPSolver::new(epsilon_hf, sigma_x, vxc, fermi_level, cfg);
    solver.solve(&sigma_c_eval, initial_guess)
}

/// Create a Σc evaluator from Padé coefficients
///
/// This function creates a closure that can evaluate Σc(ω) at any real frequency
/// using the Padé approximant coefficients from analytic continuation.
///
/// # Arguments
/// * `pade_coeff` - Thiele continued fraction coefficients
/// * `omega_fit` - Complex frequency points used in fitting
///
/// # Returns
/// A boxed closure `Fn(f64) -> Complex64` that evaluates Σc(ω)
///
/// # Example
/// ```ignore
/// let evaluator = create_sigma_c_evaluator(&coeff_p, &omega_p);
/// let sigma_c = evaluator(0.5); // Evaluate at ω = 0.5 Ha
/// ```
pub fn create_sigma_c_evaluator(
    pade_coeff: Array1<Complex64>,
    omega_fit: Array1<Complex64>,
) -> Box<dyn Fn(f64) -> Complex64 + Send + Sync> {
    Box::new(move |omega: f64| pade_thiele(omega, &omega_fit, &pade_coeff))
}

/// Solve QP equations for all orbitals using Padé coefficients
///
/// This is a high-level function that combines Padé evaluation with Newton-Raphson
/// solving for convenient integration with evGW.
///
/// # Arguments
/// * `mo_energy` - HF/DFT orbital energies
/// * `sigma_x_diag` - Exchange self-energy diagonal
/// * `vxc_diag` - DFT XC potential diagonal
/// * `pade_coeffs` - Padé coefficients [nfit, norbs]
/// * `omega_fits` - Frequency points [norbs, nfit]
/// * `orbs` - Orbital indices
/// * `nocc` - Number of occupied orbitals
/// * `initial_guesses` - Initial QP energy guesses
/// * `config` - Solver configuration
///
/// # Returns
/// Tuple of (qp_energies, z_factors, converged_flags)
pub fn solve_qp_with_pade(
    mo_energy: &Array1<f64>,
    sigma_x_diag: &Array1<f64>,
    vxc_diag: &Array1<f64>,
    pade_coeffs: &ndarray::Array2<Complex64>,
    omega_fits: &ndarray::Array2<Complex64>,
    orbs: &[usize],
    nocc: usize,
    initial_guesses: &Array1<f64>,
    config: NewtonSolverConfig,
) -> (Array1<f64>, Array1<f64>, Vec<bool>) {
    let norbs = orbs.len();
    let ef = (mo_energy[nocc - 1] + mo_energy[nocc]) / 2.0;

    let mut qp_energies = Array1::zeros(norbs);
    let mut z_factors = Array1::zeros(norbs);
    let mut converged = Vec::with_capacity(norbs);

    for (idx, &p) in orbs.iter().enumerate() {
        // Extract Padé coefficients for this orbital
        let coeff_p = pade_coeffs.column(idx).to_owned();
        let omega_p = omega_fits.row(idx).to_owned();

        // Create evaluator
        let evaluator = |omega: f64| pade_thiele(omega, &omega_p, &coeff_p);

        // Create solver
        let solver = NewtonQPSolver::new(
            mo_energy[p],
            sigma_x_diag[p],
            vxc_diag[p],
            ef,
            config.clone(),
        );

        // Solve
        match solver.solve(&evaluator, initial_guesses[p]) {
            Ok(result) => {
                qp_energies[idx] = result.qp_energy;
                z_factors[idx] = result.z_factor;
                converged.push(result.converged);
            }
            Err(_) => {
                // Fallback to linearized result
                qp_energies[idx] = initial_guesses[p];
                z_factors[idx] = 0.8; // Typical Z-factor
                converged.push(false);
            }
        }
    }

    (qp_energies, z_factors, converged)
}

/// Compare linearized vs solved QP energies
///
/// Utility function for benchmarking the accuracy improvement.
///
/// # Returns
/// Statistics about the difference: (mean_diff, max_diff, rms_diff)
pub fn compare_linearized_vs_solved(
    linearized: &Array1<f64>,
    solved: &Array1<f64>,
) -> (f64, f64, f64) {
    let n = linearized.len();
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }

    let diffs: Vec<f64> = linearized
        .iter()
        .zip(solved.iter())
        .map(|(l, s)| (l - s).abs())
        .collect();

    let mean = diffs.iter().sum::<f64>() / n as f64;
    let max = diffs.iter().cloned().fold(0.0, f64::max);
    let rms = (diffs.iter().map(|d| d * d).sum::<f64>() / n as f64).sqrt();

    (mean, max, rms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_newton_solver_config_default() {
        let config = NewtonSolverConfig::default();
        assert_eq!(config.max_iterations, 50);
        assert!((config.energy_tolerance - 1e-8).abs() < 1e-15);
        assert!(config.use_richardson);
    }

    #[test]
    fn test_newton_solver_simple() {
        // Simple case: Σc = constant = -0.1 Ha
        // QP equation: E = ε_HF + Σx + Σc - Vxc = -0.5 + (-0.3) + (-0.1) - (-0.4) = -0.5
        let epsilon_hf = -0.5;
        let sigma_x = -0.3;
        let vxc = -0.4;
        let fermi_level = 0.0;

        let sigma_c_eval = |_omega: f64| Complex64::new(-0.1, 0.0);

        let config = NewtonSolverConfig::default();
        let solver = NewtonQPSolver::new(epsilon_hf, sigma_x, vxc, fermi_level, config);

        let result = solver.solve(&sigma_c_eval, epsilon_hf).unwrap();

        // Expected: E = -0.5 + (-0.3) + (-0.1) - (-0.4) = -0.5
        assert!(result.converged);
        assert_abs_diff_eq!(result.qp_energy, -0.5, epsilon = 1e-7);
    }

    #[test]
    fn test_newton_solver_linear_sigma() {
        // Linear Σc: Σc(ω) = -0.1 - 0.2*ω
        // This tests the derivative computation
        let epsilon_hf = -0.5;
        let sigma_x = -0.3;
        let vxc = -0.4;
        let fermi_level = 0.0;

        let sigma_c_eval = |omega: f64| Complex64::new(-0.1 - 0.2 * omega, 0.0);

        let config = NewtonSolverConfig::default();
        let solver = NewtonQPSolver::new(epsilon_hf, sigma_x, vxc, fermi_level, config);

        let result = solver.solve(&sigma_c_eval, epsilon_hf).unwrap();

        // Verify convergence
        assert!(result.converged);

        // Verify residual is small
        let final_sigma_c = sigma_c_eval(result.qp_energy - fermi_level).re;
        let residual = result.qp_energy - epsilon_hf - sigma_x - final_sigma_c + vxc;
        assert!(residual.abs() < 1e-7);
    }

    #[test]
    fn test_z_factor_computation() {
        // Linear Σc with known derivative: Σc(ω) = -0.3*ω
        // Z = 1 / (1 - dΣc/dω) = 1 / (1 - (-0.3)) = 1/1.3 ≈ 0.769
        let epsilon_hf = -0.5;
        let sigma_x = 0.0;
        let vxc = 0.0;
        let fermi_level = 0.0;

        let sigma_c_eval = |omega: f64| Complex64::new(-0.3 * omega, 0.0);

        let config = NewtonSolverConfig::default();
        let solver = NewtonQPSolver::new(epsilon_hf, sigma_x, vxc, fermi_level, config);

        let z = solver.compute_z_factor(&sigma_c_eval, epsilon_hf);

        // Expected Z = 1 / (1 - (-0.3)) = 1/1.3 ≈ 0.7692
        assert_abs_diff_eq!(z, 1.0 / 1.3, epsilon = 1e-4);
    }

    #[test]
    fn test_richardson_derivative() {
        // Test Richardson extrapolation accuracy
        // Use Σc(ω) = ω³ to test 4th-order convergence
        let epsilon_hf = -0.5;
        let sigma_x = 0.0;
        let vxc = 0.0;
        let fermi_level = 0.0;

        let sigma_c_eval = |omega: f64| Complex64::new(omega.powi(3), 0.0);

        // Analytical derivative: dΣc/dω = 3ω²
        let omega_test: f64 = 0.5;
        let expected_deriv = 3.0 * omega_test.powi(2); // = 0.75

        let mut config = NewtonSolverConfig::default();
        config.use_richardson = true;
        let solver = NewtonQPSolver::new(epsilon_hf, sigma_x, vxc, fermi_level, config);

        let computed_deriv = solver.compute_derivative(&sigma_c_eval, omega_test);

        // Richardson should give ~6th order accuracy for cubic function
        assert_abs_diff_eq!(computed_deriv, expected_deriv, epsilon = 1e-10);
    }

    #[test]
    fn test_bisection_fallback() {
        // Create a case where Newton might struggle but bisection works
        // Σc with a sharp feature
        let epsilon_hf = -0.5;
        let sigma_x = -0.3;
        let vxc = -0.4;
        let fermi_level = 0.0;

        // Σc(ω) = -0.1 / (1 + ω²) - smooth but with nonlinear behavior
        let sigma_c_eval = |omega: f64| Complex64::new(-0.1 / (1.0 + omega * omega), 0.0);

        let mut config = NewtonSolverConfig::default();
        config.use_bisection_fallback = true;
        config.max_iterations = 3; // Force Newton to fail quickly

        let solver = NewtonQPSolver::new(epsilon_hf, sigma_x, vxc, fermi_level, config);

        let result = solver.solve(&sigma_c_eval, epsilon_hf);

        // Should converge (either Newton or bisection)
        assert!(result.is_ok());
        let r = result.unwrap();
        assert!(r.converged);
    }
}
