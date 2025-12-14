//! Quasiparticle calculations (G₀W₀ and evGW)
//!
//! **Sprint 5 Status**: evGW implementation for S5-1
//!
//! This module provides quasiparticle energy solvers for GW calculations.
//! - G₀W₀: One-shot GW with linearized QP equation
//! - evGW: Eigenvalue self-consistent GW (S5-1)
//! - Newton-Raphson: Solved (non-linearized) QP equation

#![allow(clippy::many_single_char_names)] // Mathematical notation

pub mod convergence;
pub mod diagnostics;
pub mod evgw;
pub mod newton_solver;
pub mod numa_alloc;
pub mod solver;

// Re-export main solver function (when implemented)
// pub use solver::solve_quasiparticle_linearized;

use crate::common::Result;
use ndarray::{Array1, Array2};
use num_complex::Complex64;

// Re-export key types
pub use convergence::{
    AccelerationMethod, ConvergenceCriteria, ConvergenceMonitor, ConvergenceStatistics,
    ParallelConvergenceChecker,
};
pub use diagnostics::{ConvergenceHistory, OrbitalAnalysis, QPAnalyzer, QPStatistics};

// Re-export evGW types (S5-1)
pub use evgw::{
    ContourDeformationConfig, EvGWConfig, EvGWConvergenceHistory, EvGWDriver, EvGWResult,
    IterationRecord, QPSolverConfig,
};

// Re-export Newton-Raphson solver types
pub use newton_solver::{
    compare_linearized_vs_solved, create_sigma_c_evaluator, solve_qp_equations, solve_qp_with_pade,
    solve_single_qp, NewtonQPSolver, NewtonSolverConfig, NewtonSolverResult,
};

// ConvergenceCriteria is now imported from the convergence module

/// Quasiparticle equation solver
pub struct QPSolver {
    /// Number of states
    pub nstates: usize,
    /// Convergence criteria
    pub criteria: ConvergenceCriteria,
}

impl QPSolver {
    /// Create a new QP solver
    #[must_use]
    pub fn new(nstates: usize) -> Self {
        Self {
            nstates,
            criteria: ConvergenceCriteria::default(),
        }
    }

    /// Solve linearized QP equation: ε_QP = ε_HF + Z·Re[Σ(ε_HF) - V_xc]
    pub fn solve_linearized(
        &self,
        _e_hf: f64,
        _sigma: Complex64,
        _vxc: f64,
        _z_factor: f64,
    ) -> Result<f64> {
        todo!("Implement linearized QP solver")
    }

    /// Solve full QP equation by Newton-Raphson
    pub fn solve_newton(&self, _e_hf: f64, _sigma_fn: impl Fn(f64) -> Complex64) -> Result<f64> {
        todo!("Implement Newton-Raphson QP solver")
    }

    /// Check convergence
    pub fn check_convergence(&self, _e_old: &Array1<f64>, _e_new: &Array1<f64>) -> Result<bool> {
        todo!("Implement convergence check")
    }
}

/// scGW driver for fully self-consistent GW
pub struct ScGWDriver {
    /// QP solver
    pub solver: QPSolver,
    /// Density mixing parameter
    pub mixing: f64,
}

impl ScGWDriver {
    /// Create a new scGW driver
    #[must_use]
    pub fn new(nstates: usize) -> Self {
        Self {
            solver: QPSolver::new(nstates),
            mixing: 0.3,
        }
    }

    /// Run scGW calculation
    pub fn run(&mut self, _initial_density: &Array2<f64>) -> Result<(Array1<f64>, Array2<f64>)> {
        todo!("Implement scGW iteration loop")
    }

    /// Update density matrix
    pub fn update_density(
        &self,
        _old_density: &Array2<f64>,
        _new_density: &Array2<f64>,
    ) -> Array2<f64> {
        todo!("Implement density mixing")
    }
}

/// Data structure for tracking iterations
#[derive(Debug, Clone)]
pub struct IterationData {
    /// Iteration number
    pub iteration: usize,
    /// QP energies
    pub energies: Array1<f64>,
    /// Z-factors
    pub z_factors: Array1<f64>,
    /// Energy change
    pub energy_change: f64,
    /// Maximum Z-factor
    pub z_max: f64,
    /// Minimum Z-factor
    pub z_min: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convergence_criteria_default() {
        let criteria = ConvergenceCriteria::default();
        assert_eq!(criteria.max_iterations, 30);
        assert_eq!(criteria.energy_tolerance, 1e-6);
    }

    #[test]
    fn test_qp_solver_creation() {
        let solver = QPSolver::new(10);
        assert_eq!(solver.nstates, 10);
    }

    // evGW tests removed - Sprint 4 scope!
    // #[test]
    // fn test_evgw_driver_creation() { ... }

    #[test]
    fn test_scgw_driver_creation() {
        let driver = ScGWDriver::new(10);
        assert_eq!(driver.mixing, 0.3);
    }
}
