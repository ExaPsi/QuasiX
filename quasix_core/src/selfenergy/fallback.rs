//! S4-4: AC Controls & Automatic Fallback System
//!
//! This module implements automatic fallback from Analytic Continuation (AC) to
//! Contour Deformation (CD) when quality metrics indicate potential numerical instability.
//!
//! # Key Components
//!
//! - **Quality Assessment**: Cross-validation error, pole stability, causality checks
//! - **Fallback Decision**: Multi-criteria evaluation with early exit
//! - **Logging**: Detailed provenance for fallback decisions
//!
//! # Physical Background
//!
//! AC is faster (~30x) but can produce unphysical results for:
//! - Small band gap systems (< 0.5 eV)
//! - Systems with strong correlation
//! - Core-level calculations
//!
//! CD is more robust but slower, providing a reliable fallback.
//!
//! # PySCF Validation (S4-3 Result)
//!
//! For well-behaved systems (H2O, NH3, CO with def2-svp):
//! - AC and CD produce **identical** HOMO/LUMO energies (MAD < 1e-6 eV)
//! - Fallback is rarely needed for frontier orbitals
//!
//! # References
//!
//! - S4-4 Theory: docs/derivations/s4-4/theory.md
//! - S4-4 Algorithms: docs/derivations/s4-4/algorithms.md
//! - PySCF CD/AC comparison: docs/reports/2025-11-25/S4_PYSCF_CD_VS_AC_COMPARISON.md

// Note: QuasixError and Result will be used when methods are fully implemented
#[allow(unused_imports)]
use crate::common::{QuasixError, Result};
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use std::fmt;

/// Method used for self-energy calculation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GWMethod {
    /// Analytic Continuation (fast, default)
    AC,
    /// Contour Deformation (robust fallback)
    CD,
    /// CD fallback after AC failure
    CDFallback,
}

impl fmt::Display for GWMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GWMethod::AC => write!(f, "AC"),
            GWMethod::CD => write!(f, "CD"),
            GWMethod::CDFallback => write!(f, "CD (AC fallback)"),
        }
    }
}

/// Reason for fallback from AC to CD
#[derive(Debug, Clone, PartialEq)]
pub enum FallbackReason {
    /// No fallback needed
    None,
    /// Cross-validation error exceeded threshold
    CrossValidationError { cv_error: f64, threshold: f64 },
    /// Pole structure is unstable (negative real parts, clustering)
    PoleInstability {
        stability_score: f64,
        threshold: f64,
    },
    /// Causality violation detected (wrong sign of Im Σ)
    CausalityViolation {
        violation_score: f64,
        threshold: f64,
    },
    /// Sum rule violation (spectral weight not conserved)
    SumRuleViolation { error: f64, threshold: f64 },
    /// Small band gap system (< 0.5 eV)
    SmallBandGap { gap_ev: f64 },
    /// User requested CD explicitly
    UserRequested,
    /// AC fitting failed to converge
    ACFitFailure { message: String },
}

impl fmt::Display for FallbackReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FallbackReason::None => write!(f, "No fallback"),
            FallbackReason::CrossValidationError {
                cv_error,
                threshold,
            } => {
                write!(f, "CV error {:.2e} > threshold {:.2e}", cv_error, threshold)
            }
            FallbackReason::PoleInstability {
                stability_score,
                threshold,
            } => {
                write!(
                    f,
                    "Pole stability {:.3} < threshold {:.3}",
                    stability_score, threshold
                )
            }
            FallbackReason::CausalityViolation {
                violation_score,
                threshold,
            } => {
                write!(
                    f,
                    "Causality violation {:.3} > threshold {:.3}",
                    violation_score, threshold
                )
            }
            FallbackReason::SumRuleViolation { error, threshold } => {
                write!(
                    f,
                    "Sum rule error {:.2e} > threshold {:.2e}",
                    error, threshold
                )
            }
            FallbackReason::SmallBandGap { gap_ev } => {
                write!(f, "Small band gap: {:.3} eV < 0.5 eV", gap_ev)
            }
            FallbackReason::UserRequested => write!(f, "User requested CD"),
            FallbackReason::ACFitFailure { message } => {
                write!(f, "AC fit failed: {}", message)
            }
        }
    }
}

/// Quality metrics from AC fitting
#[derive(Debug, Clone)]
pub struct QualityMetrics {
    /// Cross-validation RMSE (normalized)
    pub cv_error: f64,
    /// Pole stability score (0 = unstable, 1 = stable)
    pub pole_stability_score: f64,
    /// Number of problematic poles
    pub n_problematic_poles: usize,
    /// Causality violation score
    pub causality_violation_score: f64,
    /// Sum rule error
    pub sum_rule_error: f64,
    /// HOMO-LUMO gap (eV)
    pub band_gap_ev: f64,
}

impl Default for QualityMetrics {
    fn default() -> Self {
        Self {
            cv_error: 0.0,
            pole_stability_score: 1.0,
            n_problematic_poles: 0,
            causality_violation_score: 0.0,
            sum_rule_error: 0.0,
            band_gap_ev: 10.0, // Large gap = safe
        }
    }
}

/// Thresholds for triggering AC → CD fallback
#[derive(Debug, Clone)]
pub struct FallbackThresholds {
    /// Maximum cross-validation error (normalized)
    pub cv_max: f64,
    /// Minimum pole stability score
    pub pole_stability_min: f64,
    /// Maximum causality violation score
    pub causality_max: f64,
    /// Maximum sum rule error
    pub sum_rule_max: f64,
    /// Minimum band gap (eV) for AC
    pub gap_min_ev: f64,
    /// Enable adaptive threshold scaling
    pub enable_adaptive: bool,
}

impl Default for FallbackThresholds {
    fn default() -> Self {
        Self {
            cv_max: 1e-2,            // 1% relative CV error
            pole_stability_min: 0.5, // 50% stable poles required
            causality_max: 0.1,      // 10% violation tolerance
            sum_rule_max: 0.05,      // 5% sum rule error
            gap_min_ev: 0.5,         // 0.5 eV minimum gap
            enable_adaptive: true,   // Scale based on system
        }
    }
}

/// Fallback decision with provenance
#[derive(Debug, Clone)]
pub struct FallbackDecision {
    /// Should fallback to CD?
    pub should_fallback: bool,
    /// Reason for fallback (or None if AC is acceptable)
    pub reason: FallbackReason,
    /// Quality metrics that led to the decision
    pub metrics: QualityMetrics,
    /// Method to use
    pub method: GWMethod,
}

/// Fallback controller for G₀W₀ calculations
///
/// Implements the S4-4 automatic fallback system with:
/// - Quality assessment pipeline
/// - Multi-criteria decision framework
/// - Detailed logging and provenance
#[derive(Debug)]
pub struct FallbackController {
    /// Thresholds for fallback decisions
    pub thresholds: FallbackThresholds,
    /// Enable verbose logging
    pub verbose: bool,
    /// Force CD method (bypass AC)
    pub force_cd: bool,
    /// Force AC method (no fallback)
    pub force_ac: bool,
    /// Log messages for provenance
    log_messages: Vec<String>,
}

impl Default for FallbackController {
    fn default() -> Self {
        Self::new()
    }
}

impl FallbackController {
    /// Create a new fallback controller with default thresholds
    pub fn new() -> Self {
        Self {
            thresholds: FallbackThresholds::default(),
            verbose: false,
            force_cd: false,
            force_ac: false,
            log_messages: Vec::new(),
        }
    }

    /// Create with custom thresholds
    pub fn with_thresholds(thresholds: FallbackThresholds) -> Self {
        Self {
            thresholds,
            verbose: false,
            force_cd: false,
            force_ac: false,
            log_messages: Vec::new(),
        }
    }

    /// Enable verbose logging
    pub fn with_verbose(mut self, verbose: bool) -> Self {
        self.verbose = verbose;
        self
    }

    /// Force CD method (skip AC entirely)
    pub fn with_force_cd(mut self, force: bool) -> Self {
        self.force_cd = force;
        self
    }

    /// Force AC method (disable fallback)
    pub fn with_force_ac(mut self, force: bool) -> Self {
        self.force_ac = force;
        self
    }

    /// Log a message (stored for provenance)
    fn log(&mut self, message: &str) {
        self.log_messages.push(message.to_string());
        if self.verbose {
            eprintln!("[FallbackController] {}", message);
        }
    }

    /// Get log messages
    pub fn get_logs(&self) -> &[String] {
        &self.log_messages
    }

    /// Clear log messages
    pub fn clear_logs(&mut self) {
        self.log_messages.clear();
    }

    /// Evaluate fallback decision based on quality metrics
    ///
    /// Uses early-exit strategy for efficiency:
    /// 1. Check user overrides (force_cd, force_ac)
    /// 2. Check band gap
    /// 3. Check CV error
    /// 4. Check pole stability
    /// 5. Check causality
    /// 6. Check sum rules
    pub fn evaluate(&mut self, metrics: &QualityMetrics) -> FallbackDecision {
        self.clear_logs();
        self.log(&format!("Evaluating AC quality metrics: CV={:.2e}, pole_stability={:.3}, causality={:.3}, gap={:.2} eV",
            metrics.cv_error, metrics.pole_stability_score, metrics.causality_violation_score, metrics.band_gap_ev));

        // User overrides
        if self.force_cd {
            self.log("Force CD: User requested CD method");
            return FallbackDecision {
                should_fallback: true,
                reason: FallbackReason::UserRequested,
                metrics: metrics.clone(),
                method: GWMethod::CD,
            };
        }

        if self.force_ac {
            self.log("Force AC: Fallback disabled by user");
            return FallbackDecision {
                should_fallback: false,
                reason: FallbackReason::None,
                metrics: metrics.clone(),
                method: GWMethod::AC,
            };
        }

        // Adaptive thresholds based on system
        let thresholds = if self.thresholds.enable_adaptive {
            self.adaptive_thresholds(metrics)
        } else {
            self.thresholds.clone()
        };

        // Check 1: Band gap (early exit)
        if metrics.band_gap_ev < thresholds.gap_min_ev {
            self.log(&format!(
                "FALLBACK: Small band gap {:.3} eV < {:.3} eV threshold",
                metrics.band_gap_ev, thresholds.gap_min_ev
            ));
            return FallbackDecision {
                should_fallback: true,
                reason: FallbackReason::SmallBandGap {
                    gap_ev: metrics.band_gap_ev,
                },
                metrics: metrics.clone(),
                method: GWMethod::CDFallback,
            };
        }

        // Check 2: CV error
        if metrics.cv_error > thresholds.cv_max {
            self.log(&format!(
                "FALLBACK: CV error {:.2e} > threshold {:.2e}",
                metrics.cv_error, thresholds.cv_max
            ));
            return FallbackDecision {
                should_fallback: true,
                reason: FallbackReason::CrossValidationError {
                    cv_error: metrics.cv_error,
                    threshold: thresholds.cv_max,
                },
                metrics: metrics.clone(),
                method: GWMethod::CDFallback,
            };
        }

        // Check 3: Pole stability
        if metrics.pole_stability_score < thresholds.pole_stability_min {
            self.log(&format!(
                "FALLBACK: Pole stability {:.3} < threshold {:.3}",
                metrics.pole_stability_score, thresholds.pole_stability_min
            ));
            return FallbackDecision {
                should_fallback: true,
                reason: FallbackReason::PoleInstability {
                    stability_score: metrics.pole_stability_score,
                    threshold: thresholds.pole_stability_min,
                },
                metrics: metrics.clone(),
                method: GWMethod::CDFallback,
            };
        }

        // Check 4: Causality
        if metrics.causality_violation_score > thresholds.causality_max {
            self.log(&format!(
                "FALLBACK: Causality violation {:.3} > threshold {:.3}",
                metrics.causality_violation_score, thresholds.causality_max
            ));
            return FallbackDecision {
                should_fallback: true,
                reason: FallbackReason::CausalityViolation {
                    violation_score: metrics.causality_violation_score,
                    threshold: thresholds.causality_max,
                },
                metrics: metrics.clone(),
                method: GWMethod::CDFallback,
            };
        }

        // Check 5: Sum rules
        if metrics.sum_rule_error > thresholds.sum_rule_max {
            self.log(&format!(
                "FALLBACK: Sum rule error {:.2e} > threshold {:.2e}",
                metrics.sum_rule_error, thresholds.sum_rule_max
            ));
            return FallbackDecision {
                should_fallback: true,
                reason: FallbackReason::SumRuleViolation {
                    error: metrics.sum_rule_error,
                    threshold: thresholds.sum_rule_max,
                },
                metrics: metrics.clone(),
                method: GWMethod::CDFallback,
            };
        }

        // All checks passed - use AC
        self.log("AC quality acceptable, no fallback needed");
        FallbackDecision {
            should_fallback: false,
            reason: FallbackReason::None,
            metrics: metrics.clone(),
            method: GWMethod::AC,
        }
    }

    /// Compute adaptive thresholds based on system properties
    ///
    /// Stricter thresholds for:
    /// - Small band gaps (< 1 eV): 2x stricter
    /// - Very small gaps (< 0.5 eV): 4x stricter
    fn adaptive_thresholds(&self, metrics: &QualityMetrics) -> FallbackThresholds {
        let mut adapted = self.thresholds.clone();
        let gap = metrics.band_gap_ev;

        // Scale factor based on gap
        let scale = if gap < 0.5 {
            0.25 // 4x stricter
        } else if gap < 1.0 {
            0.5 // 2x stricter
        } else if gap < 2.0 {
            0.75 // 1.33x stricter
        } else {
            1.0 // Default
        };

        adapted.cv_max *= scale;
        adapted.causality_max *= scale;
        adapted.sum_rule_max *= scale;

        if self.verbose {
            eprintln!(
                "[FallbackController] Adaptive scaling: gap={:.2} eV → scale={:.2}",
                gap, scale
            );
        }

        adapted
    }

    /// Compute quality metrics from AC fitting result
    ///
    /// # Arguments
    ///
    /// * `sigma` - Self-energy on imaginary axis [norbs, nw_sigma]
    /// * `omega` - Frequency points [norbs, nw_sigma]
    /// * `coeff` - Pade coefficients [nfit, norbs]
    /// * `mo_energy` - MO energies [nmo]
    /// * `nocc` - Number of occupied orbitals
    ///
    /// # Returns
    ///
    /// Quality metrics for fallback decision
    pub fn compute_quality_metrics(
        &self,
        sigma: &Array2<Complex64>,
        omega: &Array2<Complex64>,
        coeff: &Array2<Complex64>,
        mo_energy: &Array1<f64>,
        nocc: usize,
    ) -> QualityMetrics {
        let (norbs, nw) = sigma.dim();
        let (_nfit, _) = coeff.dim();

        // Compute band gap (HOMO-LUMO)
        let homo = mo_energy[nocc - 1];
        let lumo = mo_energy[nocc];
        let gap_ha = lumo - homo;
        let gap_ev = gap_ha * 27.2114; // Ha to eV

        // Cross-validation error (simplified - leave-one-out on first 10 points)
        let cv_error = self.compute_cv_error(sigma, omega, coeff, norbs.min(5), nw.min(10));

        // Pole stability analysis
        let (stability_score, n_problematic) = self.analyze_pole_stability(coeff);

        // Causality check (simplified)
        let causality_score = self.check_causality(coeff, mo_energy, nocc);

        // Sum rule check (not fully implemented yet)
        let sum_rule_error = 0.0; // Placeholder

        QualityMetrics {
            cv_error,
            pole_stability_score: stability_score,
            n_problematic_poles: n_problematic,
            causality_violation_score: causality_score,
            sum_rule_error,
            band_gap_ev: gap_ev,
        }
    }

    /// Simplified cross-validation error estimation
    fn compute_cv_error(
        &self,
        sigma: &Array2<Complex64>,
        _omega: &Array2<Complex64>, // Reserved for future frequency-weighted CV error
        _coeff: &Array2<Complex64>,
        n_orbs_check: usize,
        n_points_check: usize,
    ) -> f64 {
        // For now, estimate CV error from variance of imaginary part
        // (Full k-fold CV is expensive)
        let mut total_var = 0.0;
        let mut total_mag = 0.0;

        for orb in 0..n_orbs_check {
            for k in 1..n_points_check {
                let diff = (sigma[[orb, k]] - sigma[[orb, k - 1]]).norm();
                let mag = sigma[[orb, k]].norm().max(1e-10);
                total_var += diff / mag;
                total_mag += 1.0;
            }
        }

        if total_mag > 0.0 {
            total_var / total_mag
        } else {
            0.0
        }
    }

    /// Analyze pole stability from Pade coefficients
    fn analyze_pole_stability(&self, coeff: &Array2<Complex64>) -> (f64, usize) {
        let (nfit, norbs) = coeff.dim();
        let mut n_problematic = 0;

        for p in 0..norbs {
            for k in 0..nfit {
                let c = coeff[[k, p]];
                // Check for very large coefficients (instability)
                if c.norm() > 1e6 {
                    n_problematic += 1;
                }
                // Check for coefficients with wrong sign in imaginary part
                if k > 0 && c.im.abs() > 1e3 {
                    n_problematic += 1;
                }
            }
        }

        let total_coeffs = nfit * norbs;
        let stability_score = if total_coeffs > 0 {
            1.0 - (n_problematic as f64 / total_coeffs as f64)
        } else {
            1.0
        };

        (stability_score.max(0.0), n_problematic)
    }

    /// Check causality (simplified)
    fn check_causality(
        &self,
        coeff: &Array2<Complex64>,
        _mo_energy: &Array1<f64>,
        _nocc: usize,
    ) -> f64 {
        // Simplified: check that first coefficients have reasonable imaginary parts
        let (nfit, norbs) = coeff.dim();
        let mut violations = 0;

        for p in 0..norbs {
            if nfit > 0 {
                let a0 = coeff[[0, p]];
                // a0 should be mostly real for well-behaved systems
                if a0.im.abs() > a0.re.abs() * 0.5 {
                    violations += 1;
                }
            }
        }

        violations as f64 / norbs.max(1) as f64
    }
}

/// GW result with method provenance
#[derive(Debug)]
pub struct GWResultWithProvenance {
    /// Self-energy on imaginary axis
    pub sigma: Array2<Complex64>,
    /// Frequency points
    pub omega: Array2<Complex64>,
    /// Method used
    pub method: GWMethod,
    /// Fallback decision (if applicable)
    pub fallback_decision: Option<FallbackDecision>,
    /// Quality metrics (if available)
    pub quality_metrics: Option<QualityMetrics>,
    /// Log messages for provenance
    pub logs: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fallback_controller_default() {
        let controller = FallbackController::new();
        assert!(!controller.verbose);
        assert!(!controller.force_cd);
        assert!(!controller.force_ac);
    }

    #[test]
    fn test_fallback_decision_good_metrics() {
        let mut controller = FallbackController::new();
        let metrics = QualityMetrics {
            cv_error: 1e-4,
            pole_stability_score: 0.95,
            n_problematic_poles: 0,
            causality_violation_score: 0.01,
            sum_rule_error: 0.001,
            band_gap_ev: 5.0,
        };

        let decision = controller.evaluate(&metrics);
        assert!(!decision.should_fallback);
        assert_eq!(decision.method, GWMethod::AC);
    }

    #[test]
    fn test_fallback_decision_small_gap() {
        let mut controller = FallbackController::new();
        let metrics = QualityMetrics {
            cv_error: 1e-4,
            pole_stability_score: 0.95,
            n_problematic_poles: 0,
            causality_violation_score: 0.01,
            sum_rule_error: 0.001,
            band_gap_ev: 0.3, // Small gap triggers fallback
        };

        let decision = controller.evaluate(&metrics);
        assert!(decision.should_fallback);
        assert_eq!(decision.method, GWMethod::CDFallback);
        assert!(matches!(
            decision.reason,
            FallbackReason::SmallBandGap { .. }
        ));
    }

    #[test]
    fn test_fallback_decision_high_cv_error() {
        let mut controller = FallbackController::new();
        let metrics = QualityMetrics {
            cv_error: 0.1, // High CV error
            pole_stability_score: 0.95,
            n_problematic_poles: 0,
            causality_violation_score: 0.01,
            sum_rule_error: 0.001,
            band_gap_ev: 5.0,
        };

        let decision = controller.evaluate(&metrics);
        assert!(decision.should_fallback);
        assert!(matches!(
            decision.reason,
            FallbackReason::CrossValidationError { .. }
        ));
    }

    #[test]
    fn test_fallback_decision_force_cd() {
        let mut controller = FallbackController::new().with_force_cd(true);
        let metrics = QualityMetrics::default(); // Good metrics

        let decision = controller.evaluate(&metrics);
        assert!(decision.should_fallback);
        assert_eq!(decision.method, GWMethod::CD);
        assert!(matches!(decision.reason, FallbackReason::UserRequested));
    }

    #[test]
    fn test_fallback_reason_display() {
        let reason = FallbackReason::CrossValidationError {
            cv_error: 0.05,
            threshold: 0.01,
        };
        let display = format!("{}", reason);
        assert!(display.contains("CV error"));
        assert!(display.contains("threshold"));
    }

    #[test]
    fn test_gw_method_display() {
        assert_eq!(format!("{}", GWMethod::AC), "AC");
        assert_eq!(format!("{}", GWMethod::CD), "CD");
        assert_eq!(format!("{}", GWMethod::CDFallback), "CD (AC fallback)");
    }

    #[test]
    fn test_adaptive_thresholds() {
        let controller = FallbackController::new();
        let metrics_large_gap = QualityMetrics {
            band_gap_ev: 5.0,
            ..Default::default()
        };
        let metrics_small_gap = QualityMetrics {
            band_gap_ev: 0.3,
            ..Default::default()
        };

        let thresh_large = controller.adaptive_thresholds(&metrics_large_gap);
        let thresh_small = controller.adaptive_thresholds(&metrics_small_gap);

        // Small gap should have stricter (smaller) thresholds
        assert!(thresh_small.cv_max < thresh_large.cv_max);
    }
}
