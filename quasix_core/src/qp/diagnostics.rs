//! Diagnostics and analysis tools for QP solver
//!
//! This module provides detailed diagnostics for quasiparticle calculations,
//! including convergence analysis, Z-factor validation, and stability metrics.

use crate::common::Result;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Convergence history for QP iterations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceHistory {
    /// Iteration number
    pub iteration: usize,
    /// Energy at this iteration
    pub energy: f64,
    /// Energy change from previous iteration
    pub energy_change: f64,
    /// Function value |F(E)|
    pub residual: f64,
    /// Derivative value F'(E)
    pub derivative: f64,
    /// Step size taken
    pub step_size: f64,
    /// Damping factor used
    pub damping: f64,
}

/// Detailed analysis for a single orbital
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitalAnalysis {
    /// Orbital index
    pub orbital_idx: usize,
    /// Orbital type (occupied/virtual)
    pub orbital_type: String,
    /// Reference energy (Ha)
    pub reference_energy: f64,
    /// QP energy (Ha)
    pub qp_energy: f64,
    /// Energy shift (Ha)
    pub energy_shift: f64,
    /// Z-factor
    pub z_factor: f64,
    /// Physical validity
    pub is_physical: bool,
    /// Convergence history
    pub convergence_history: Vec<ConvergenceHistory>,
    /// Solver method used
    pub solver_method: String,
    /// Total iterations
    pub iterations: usize,
    /// Convergence achieved
    pub converged: bool,
}

/// Statistical analysis of QP results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QPStatistics {
    /// Total number of orbitals
    pub n_orbitals: usize,
    /// Number of occupied orbitals
    pub n_occupied: usize,
    /// Number of virtual orbitals  
    pub n_virtual: usize,
    /// Convergence statistics
    pub convergence_stats: ConvergenceStats,
    /// Z-factor statistics
    pub z_factor_stats: ZFactorStats,
    /// Energy shift statistics
    pub energy_shift_stats: EnergyShiftStats,
    /// Solver performance
    pub performance_stats: PerformanceStats,
}

/// Convergence statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceStats {
    /// Number of converged orbitals
    pub n_converged: usize,
    /// Convergence rate (%)
    pub convergence_rate: f64,
    /// Average iterations to converge
    pub avg_iterations: f64,
    /// Maximum iterations needed
    pub max_iterations: usize,
    /// Minimum iterations needed
    pub min_iterations: usize,
    /// Number using Newton-Raphson
    pub newton_count: usize,
    /// Number using bisection fallback
    pub bisection_count: usize,
}

/// Z-factor statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZFactorStats {
    /// Number with physical Z ∈ (0,1)
    pub n_physical: usize,
    /// Physical rate (%)
    pub physical_rate: f64,
    /// Mean Z-factor
    pub mean: f64,
    /// Standard deviation
    pub std_dev: f64,
    /// Minimum Z-factor
    pub min: f64,
    /// Maximum Z-factor
    pub max: f64,
    /// Z-factor for HOMO
    pub homo_z: Option<f64>,
    /// Z-factor for LUMO
    pub lumo_z: Option<f64>,
    /// Distribution histogram
    pub distribution: HashMap<String, usize>,
}

/// Energy shift statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnergyShiftStats {
    /// Mean energy shift (Ha)
    pub mean: f64,
    /// Standard deviation (Ha)
    pub std_dev: f64,
    /// Maximum shift (Ha)
    pub max_shift: f64,
    /// Minimum shift (Ha)
    pub min_shift: f64,
    /// HOMO energy shift (Ha)
    pub homo_shift: Option<f64>,
    /// LUMO energy shift (Ha)
    pub lumo_shift: Option<f64>,
    /// IP correction (eV)
    pub ip_correction: Option<f64>,
    /// EA correction (eV)
    pub ea_correction: Option<f64>,
}

/// Performance statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceStats {
    /// Total solver time (seconds)
    pub total_time: f64,
    /// Average time per orbital (seconds)
    pub avg_time_per_orbital: f64,
    /// Number of sigma evaluations
    pub n_sigma_evaluations: usize,
    /// Average sigma evaluations per orbital
    pub avg_sigma_per_orbital: f64,
}

/// QP solver analyzer
pub struct QPAnalyzer;

impl QPAnalyzer {
    /// Analyze QP results and generate statistics
    pub fn analyze(
        orbital_solutions: &[super::solver::QPStateSolution],
        mo_energies: &Array1<f64>,
        mo_occ: &Array1<f64>,
        timing_info: Option<f64>,
    ) -> Result<QPStatistics> {
        let n_orbitals = orbital_solutions.len();
        let n_occupied = mo_occ.iter().filter(|&&occ| occ > 0.5).count();
        let n_virtual = n_orbitals - n_occupied;

        // Compute convergence statistics
        let convergence_stats = Self::compute_convergence_stats(orbital_solutions);

        // Compute Z-factor statistics
        let z_factor_stats = Self::compute_z_factor_stats(orbital_solutions, n_occupied, n_virtual);

        // Compute energy shift statistics
        let energy_shift_stats =
            Self::compute_energy_shift_stats(orbital_solutions, mo_energies, n_occupied, n_virtual);

        // Compute performance statistics
        let performance_stats = Self::compute_performance_stats(orbital_solutions, timing_info);

        Ok(QPStatistics {
            n_orbitals,
            n_occupied,
            n_virtual,
            convergence_stats,
            z_factor_stats,
            energy_shift_stats,
            performance_stats,
        })
    }

    /// Compute convergence statistics
    fn compute_convergence_stats(solutions: &[super::solver::QPStateSolution]) -> ConvergenceStats {
        let n_converged = solutions.iter().filter(|s| s.converged).count();
        let convergence_rate = 100.0 * n_converged as f64 / solutions.len() as f64;

        let iterations: Vec<usize> = solutions.iter().map(|s| s.iterations).collect();
        let avg_iterations = iterations.iter().sum::<usize>() as f64 / iterations.len() as f64;
        let max_iterations = *iterations.iter().max().unwrap_or(&0);
        let min_iterations = *iterations.iter().min().unwrap_or(&0);

        let newton_count = solutions.iter().filter(|s| !s.used_bisection).count();
        let bisection_count = solutions.iter().filter(|s| s.used_bisection).count();

        ConvergenceStats {
            n_converged,
            convergence_rate,
            avg_iterations,
            max_iterations,
            min_iterations,
            newton_count,
            bisection_count,
        }
    }

    /// Compute Z-factor statistics
    fn compute_z_factor_stats(
        solutions: &[super::solver::QPStateSolution],
        n_occupied: usize,
        n_virtual: usize,
    ) -> ZFactorStats {
        let z_factors: Vec<f64> = solutions.iter().map(|s| s.z_factor).collect();

        let n_physical = z_factors.iter().filter(|&&z| z > 0.0 && z < 1.0).count();
        let physical_rate = 100.0 * n_physical as f64 / z_factors.len() as f64;

        let mean = z_factors.iter().sum::<f64>() / z_factors.len() as f64;
        let variance =
            z_factors.iter().map(|&z| (z - mean).powi(2)).sum::<f64>() / z_factors.len() as f64;
        let std_dev = variance.sqrt();

        let min = z_factors.iter().copied().fold(f64::INFINITY, f64::min);
        let max = z_factors.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        // HOMO and LUMO Z-factors
        let homo_z = if n_occupied > 0 {
            Some(solutions[n_occupied - 1].z_factor)
        } else {
            None
        };

        let lumo_z = if n_virtual > 0 && n_occupied < solutions.len() {
            Some(solutions[n_occupied].z_factor)
        } else {
            None
        };

        // Distribution histogram
        let mut distribution = HashMap::new();
        for &z in &z_factors {
            let bin = if z < 0.0 {
                "< 0.0"
            } else if z < 0.2 {
                "0.0-0.2"
            } else if z < 0.4 {
                "0.2-0.4"
            } else if z < 0.6 {
                "0.4-0.6"
            } else if z < 0.8 {
                "0.6-0.8"
            } else if z < 1.0 {
                "0.8-1.0"
            } else {
                ">= 1.0"
            };
            *distribution.entry(bin.to_string()).or_insert(0) += 1;
        }

        ZFactorStats {
            n_physical,
            physical_rate,
            mean,
            std_dev,
            min,
            max,
            homo_z,
            lumo_z,
            distribution,
        }
    }

    /// Compute energy shift statistics
    fn compute_energy_shift_stats(
        solutions: &[super::solver::QPStateSolution],
        _mo_energies: &Array1<f64>,
        n_occupied: usize,
        n_virtual: usize,
    ) -> EnergyShiftStats {
        let shifts: Vec<f64> = solutions.iter().map(|s| s.energy_shift).collect();

        let mean = shifts.iter().sum::<f64>() / shifts.len() as f64;
        let variance =
            shifts.iter().map(|&s| (s - mean).powi(2)).sum::<f64>() / shifts.len() as f64;
        let std_dev = variance.sqrt();

        let max_shift = shifts.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min_shift = shifts.iter().copied().fold(f64::INFINITY, f64::min);

        // HOMO and LUMO shifts
        let homo_shift = if n_occupied > 0 {
            Some(solutions[n_occupied - 1].energy_shift)
        } else {
            None
        };

        let lumo_shift = if n_virtual > 0 && n_occupied < solutions.len() {
            Some(solutions[n_occupied].energy_shift)
        } else {
            None
        };

        // IP and EA corrections (in eV)
        const HA_TO_EV: f64 = 27.211_386_245_988;
        let ip_correction = homo_shift.map(|s| -s * HA_TO_EV);
        let ea_correction = lumo_shift.map(|s| -s * HA_TO_EV);

        EnergyShiftStats {
            mean,
            std_dev,
            max_shift,
            min_shift,
            homo_shift,
            lumo_shift,
            ip_correction,
            ea_correction,
        }
    }

    /// Compute performance statistics
    fn compute_performance_stats(
        solutions: &[super::solver::QPStateSolution],
        timing_info: Option<f64>,
    ) -> PerformanceStats {
        let total_time = timing_info.unwrap_or(0.0);
        let avg_time_per_orbital = total_time / solutions.len() as f64;

        // Estimate sigma evaluations (2 per derivative + 1 per iteration)
        let n_sigma_evaluations: usize = solutions.iter()
            .map(|s| 3 * s.iterations) // Approximate
            .sum();
        let avg_sigma_per_orbital = n_sigma_evaluations as f64 / solutions.len() as f64;

        PerformanceStats {
            total_time,
            avg_time_per_orbital,
            n_sigma_evaluations,
            avg_sigma_per_orbital,
        }
    }

    /// Validate Z-factors are physical
    pub fn validate_z_factors(z_factors: &Array1<f64>) -> Result<()> {
        for (i, &z) in z_factors.iter().enumerate() {
            if z <= 0.0 || z >= 1.0 {
                log::error!("Non-physical Z-factor {} for orbital {}", z, i);
                return Err(crate::common::QuasixError::PhysicsError(format!(
                    "Z-factor {} outside (0,1) for orbital {}",
                    z, i
                )));
            }
        }
        Ok(())
    }

    /// Check for convergence issues
    pub fn check_convergence_quality(
        solutions: &[super::solver::QPStateSolution],
        _threshold: f64,
    ) -> Vec<usize> {
        solutions
            .iter()
            .enumerate()
            .filter_map(|(i, s)| {
                if !s.converged || s.iterations > 8 {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Generate summary report
    pub fn generate_report(stats: &QPStatistics) -> String {
        format!(
            "QP Solver Analysis Report\n\
             ========================\n\n\
             System Overview:\n\
             - Total orbitals: {}\n\
             - Occupied: {}, Virtual: {}\n\n\
             Convergence:\n\
             - Converged: {}/{} ({:.1}%)\n\
             - Average iterations: {:.2}\n\
             - Newton-Raphson: {}, Bisection: {}\n\n\
             Z-Factors:\n\
             - Physical (0,1): {}/{} ({:.1}%)\n\
             - Mean: {:.4} ± {:.4}\n\
             - Range: [{:.4}, {:.4}]\n\
             - HOMO Z: {:.4}, LUMO Z: {:.4}\n\n\
             Energy Shifts:\n\
             - Mean: {:.4} ± {:.4} Ha\n\
             - Range: [{:.4}, {:.4}] Ha\n\
             - IP correction: {:.3} eV\n\
             - EA correction: {:.3} eV\n\n\
             Performance:\n\
             - Total time: {:.2} s\n\
             - Time per orbital: {:.3} s\n\
             - Sigma evaluations: {}",
            stats.n_orbitals,
            stats.n_occupied,
            stats.n_virtual,
            stats.convergence_stats.n_converged,
            stats.n_orbitals,
            stats.convergence_stats.convergence_rate,
            stats.convergence_stats.avg_iterations,
            stats.convergence_stats.newton_count,
            stats.convergence_stats.bisection_count,
            stats.z_factor_stats.n_physical,
            stats.n_orbitals,
            stats.z_factor_stats.physical_rate,
            stats.z_factor_stats.mean,
            stats.z_factor_stats.std_dev,
            stats.z_factor_stats.min,
            stats.z_factor_stats.max,
            stats.z_factor_stats.homo_z.unwrap_or(0.0),
            stats.z_factor_stats.lumo_z.unwrap_or(0.0),
            stats.energy_shift_stats.mean,
            stats.energy_shift_stats.std_dev,
            stats.energy_shift_stats.min_shift,
            stats.energy_shift_stats.max_shift,
            stats.energy_shift_stats.ip_correction.unwrap_or(0.0),
            stats.energy_shift_stats.ea_correction.unwrap_or(0.0),
            stats.performance_stats.total_time,
            stats.performance_stats.avg_time_per_orbital,
            stats.performance_stats.n_sigma_evaluations,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_z_factor_validation() {
        let mut z_factors = Array1::from(vec![0.5, 0.7, 0.9]);
        assert!(QPAnalyzer::validate_z_factors(&z_factors).is_ok());

        z_factors[1] = 1.5; // Invalid
        assert!(QPAnalyzer::validate_z_factors(&z_factors).is_err());

        z_factors[1] = -0.1; // Invalid
        assert!(QPAnalyzer::validate_z_factors(&z_factors).is_err());
    }
}
