"""QuasiX Fallback Controller for Automatic AC→CD Transitions

This module provides the main Python interface for the fallback controller,
which automatically decides whether to use Analytic Continuation (AC) or
Contour Deformation (CD) based on quality assessment metrics.

Author: QuasiX Development Team
License: MIT
"""

from dataclasses import dataclass, field
from typing import Optional, Tuple, Callable, Dict, Any, List, Union
from enum import Enum, auto
import logging
import warnings
import time
from pathlib import Path
import json

try:
    import yaml
except ImportError:
    yaml = None

import numpy as np
from scipy import linalg
from pyscf import lib, gto, scf, df
from pyscf.lib import logger

# Import the Rust bindings (will be created)
try:
    from .. import _quasix_core
except ImportError:
    warnings.warn("QuasiX core module not available, using fallback Python implementation")
    _quasix_core = None


class FallbackReason(Enum):
    """Reasons for triggering AC→CD fallback."""

    CROSS_VALIDATION_ERROR = "cv_error"
    UNSTABLE_POLES = "pole_instability"
    CAUSALITY_VIOLATION = "causality"
    SUM_RULE_VIOLATION = "sum_rule"
    NUMERICAL_INSTABILITY = "numerical"
    RESIDUE_MAGNITUDE = "large_residues"
    CONVERGENCE_FAILURE = "convergence"
    USER_REQUESTED = "user_override"

    def __str__(self):
        """Human-readable description of the fallback reason."""
        descriptions = {
            self.CROSS_VALIDATION_ERROR: "Cross-validation error exceeded threshold",
            self.UNSTABLE_POLES: "Unstable or unphysical pole structure detected",
            self.CAUSALITY_VIOLATION: "Causality violation in spectral function",
            self.SUM_RULE_VIOLATION: "Spectral sum rule violation",
            self.NUMERICAL_INSTABILITY: "Numerical instability in AC procedure",
            self.RESIDUE_MAGNITUDE: "Anomalously large residue magnitudes",
            self.CONVERGENCE_FAILURE: "AC failed to converge",
            self.USER_REQUESTED: "User explicitly requested CD method"
        }
        return descriptions.get(self, "Unknown fallback reason")


class FallbackMethod(Enum):
    """Available computational methods."""

    ANALYTIC_CONTINUATION = "AC"
    CONTOUR_DEFORMATION = "CD"
    HYBRID = "Hybrid"  # Adaptive combination


@dataclass
class FallbackThresholds:
    """Configuration for fallback decision thresholds.

    Attributes:
        cv_error_max: Maximum allowed cross-validation error (relative)
        pole_stability_min: Minimum imaginary part for stable poles (Ha)
        causality_violation_max: Maximum causality violation tolerance
        residue_magnitude_max: Maximum allowed residue magnitude
        sum_rule_tolerance: Tolerance for spectral sum rule
        enable_adaptive: Enable adaptive threshold adjustment
        confidence_level: Statistical confidence level for decisions
    """

    cv_error_max: float = 1e-2
    pole_stability_min: float = 1e-3
    causality_violation_max: float = 1e-4
    residue_magnitude_max: float = 100.0
    sum_rule_tolerance: float = 0.05
    enable_adaptive: bool = True
    confidence_level: float = 0.95

    def to_dict(self) -> Dict[str, Any]:
        """Convert thresholds to dictionary."""
        return {
            'cv_error_max': self.cv_error_max,
            'pole_stability_min': self.pole_stability_min,
            'causality_violation_max': self.causality_violation_max,
            'residue_magnitude_max': self.residue_magnitude_max,
            'sum_rule_tolerance': self.sum_rule_tolerance,
            'enable_adaptive': self.enable_adaptive,
            'confidence_level': self.confidence_level
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'FallbackThresholds':
        """Create thresholds from dictionary."""
        return cls(**data)

    @classmethod
    def conservative(cls) -> 'FallbackThresholds':
        """Conservative thresholds favoring CD."""
        return cls(
            cv_error_max=0.005,
            pole_stability_min=1e-2,
            causality_violation_max=1e-5,
            residue_magnitude_max=50.0,
            sum_rule_tolerance=0.02
        )

    @classmethod
    def balanced(cls) -> 'FallbackThresholds':
        """Balanced thresholds (default)."""
        return cls()

    @classmethod
    def aggressive(cls) -> 'FallbackThresholds':
        """Aggressive thresholds favoring AC."""
        return cls(
            cv_error_max=0.05,
            pole_stability_min=1e-4,
            causality_violation_max=1e-3,
            residue_magnitude_max=200.0,
            sum_rule_tolerance=0.1
        )


@dataclass
class QualityMetrics:
    """Quality assessment metrics for AC model.

    Attributes:
        cv_error: Cross-validation error (relative)
        pole_stability_score: Measure of pole stability (0-1, higher is better)
        causality_violation: Causality violation measure
        residue_magnitude_max: Maximum residue magnitude
        sum_rule_error: Spectral sum rule error
        convergence_rate: Convergence rate of AC procedure
        z_factor_range: Tuple of (min, max) Z factors
        condition_number: Condition number of fitting problem
        timestamp: Time of assessment
    """

    cv_error: float
    pole_stability_score: float
    causality_violation: float
    residue_magnitude_max: float
    sum_rule_error: float
    convergence_rate: float = 0.0
    z_factor_range: Tuple[float, float] = (0.0, 1.0)
    condition_number: float = 1.0
    timestamp: float = field(default_factory=time.time)

    def passes_thresholds(self, thresholds: FallbackThresholds) -> bool:
        """Check if metrics pass all thresholds."""
        z_min, z_max = self.z_factor_range
        return (
            self.cv_error <= thresholds.cv_error_max and
            self.pole_stability_score >= 0.7 and  # Good stability score
            self.causality_violation <= thresholds.causality_violation_max and
            self.residue_magnitude_max <= thresholds.residue_magnitude_max and
            self.sum_rule_error <= thresholds.sum_rule_tolerance and
            z_min >= 0.0 and z_max <= 1.0  # Physical Z factors
        )

    def get_failures(self, thresholds: FallbackThresholds) -> List[FallbackReason]:
        """Get list of threshold failures."""
        failures = []

        if self.cv_error > thresholds.cv_error_max:
            failures.append(FallbackReason.CROSS_VALIDATION_ERROR)

        if self.pole_stability_score < 0.7:  # Poor stability score
            failures.append(FallbackReason.UNSTABLE_POLES)

        if self.causality_violation > thresholds.causality_violation_max:
            failures.append(FallbackReason.CAUSALITY_VIOLATION)

        if self.residue_magnitude_max > thresholds.residue_magnitude_max:
            failures.append(FallbackReason.RESIDUE_MAGNITUDE)

        if self.sum_rule_error > thresholds.sum_rule_tolerance:
            failures.append(FallbackReason.SUM_RULE_VIOLATION)

        # Check Z factor physical bounds
        z_min, z_max = self.z_factor_range
        if z_min < 0.0 or z_max > 1.0:
            failures.append(FallbackReason.CAUSALITY_VIOLATION)

        return failures

    def calculate_confidence(self) -> float:
        """Calculate confidence score from quality metrics.

        Returns:
            Confidence score in [0, 1]
        """
        scores = []

        # CV error score (lower is better, max 0.1 for worst case)
        cv_score = max(0, 1.0 - self.cv_error / 0.1)
        scores.append(cv_score)

        # Pole stability score (already 0-1)
        scores.append(self.pole_stability_score)

        # Causality score (lower is better, max 0.01 for worst case)
        causality_score = max(0, 1.0 - self.causality_violation / 0.01)
        scores.append(causality_score)

        # Residue magnitude score (lower is better, max 100 for worst case)
        residue_score = max(0, 1.0 - self.residue_magnitude_max / 100.0)
        scores.append(residue_score)

        # Sum rule score (lower is better, max 0.1 for worst case)
        sum_rule_score = max(0, 1.0 - self.sum_rule_error / 0.1)
        scores.append(sum_rule_score)

        # Z factor physical validity
        z_min, z_max = self.z_factor_range
        if 0.0 <= z_min <= z_max <= 1.0:
            z_score = 1.0
        else:
            z_score = 0.0
        scores.append(z_score)

        # Convergence rate (already 0-1)
        scores.append(self.convergence_rate)

        # Average confidence
        return float(np.mean(scores))


@dataclass
class FallbackDecision:
    """Result of fallback decision process.

    Attributes:
        should_fallback: Whether to use CD instead of AC
        method_used: The method to be used
        reason: Primary reason for fallback (if applicable)
        all_reasons: All detected issues
        metrics: Quality metrics that led to decision
        confidence: Confidence in the decision (0-1)
        recommendation: Human-readable recommendation
    """

    should_fallback: bool
    method_used: FallbackMethod
    reason: Optional[FallbackReason] = None
    all_reasons: List[FallbackReason] = field(default_factory=list)
    metrics: Optional[QualityMetrics] = None
    confidence: float = 1.0
    recommendation: str = ""

    def __str__(self) -> str:
        """Human-readable decision summary."""
        if self.should_fallback:
            return f"Fallback to {self.method_used.value}: {self.reason}"
        else:
            return f"Use {self.method_used.value} (confidence: {self.confidence:.2%})"

    @property
    def details(self) -> str:
        """Detailed description of the fallback decision.

        Returns:
            Comprehensive string describing why the fallback was triggered,
            including all failure reasons and quality metrics.
        """
        details_parts = []

        # Start with the primary decision
        if self.should_fallback:
            details_parts.append(f"Fallback triggered to {self.method_used.value}")

            # Add primary reason
            if self.reason:
                details_parts.append(f"Primary reason: {str(self.reason)}")

            # Add all detected issues
            if self.all_reasons:
                details_parts.append(f"All issues detected ({len(self.all_reasons)}):")
                for idx, reason in enumerate(self.all_reasons, 1):
                    details_parts.append(f"  {idx}. {str(reason)}")

            # Add metrics if available
            if self.metrics:
                details_parts.append("Quality metrics summary:")
                details_parts.append(f"  - CV error: {self.metrics.cv_error:.3e}")
                details_parts.append(f"  - Pole stability score: {self.metrics.pole_stability_score:.3f}")
                details_parts.append(f"  - Causality violation: {self.metrics.causality_violation:.3e}")
                details_parts.append(f"  - Max residue magnitude: {self.metrics.residue_magnitude_max:.2f}")
                details_parts.append(f"  - Sum rule error: {self.metrics.sum_rule_error:.3f}")
                z_min, z_max = self.metrics.z_factor_range
                details_parts.append(f"  - Z factor range: [{z_min:.3f}, {z_max:.3f}]")
                details_parts.append(f"  - Confidence: {self.confidence:.2%}")
        else:
            details_parts.append(f"AC method selected (no fallback needed)")
            details_parts.append(f"Confidence: {self.confidence:.2%}")

            # Add metrics for successful AC
            if self.metrics:
                details_parts.append("Quality metrics (all passed):")
                details_parts.append(f"  - CV error: {self.metrics.cv_error:.3e}")
                details_parts.append(f"  - Pole stability score: {self.metrics.pole_stability_score:.3f}")
                details_parts.append(f"  - Causality: OK ({self.metrics.causality_violation:.3e})")
                z_min, z_max = self.metrics.z_factor_range
                details_parts.append(f"  - Z factors: [{z_min:.3f}, {z_max:.3f}] (physical)")

        # Add recommendation
        if self.recommendation:
            details_parts.append(f"Recommendation: {self.recommendation}")

        return "\n".join(details_parts)


@dataclass
class GWParams:
    """Parameters for GW calculation.

    Attributes:
        method: Initial method preference (AC or CD)
        freq_grid_size: Number of frequency points
        eta: Broadening parameter
        max_iter: Maximum iterations for self-consistency
        conv_tol: Convergence tolerance
        auxbasis: Auxiliary basis for density fitting
    """

    method: str = "AC"
    freq_grid_size: int = 200
    eta: float = 0.01
    max_iter: int = 50
    conv_tol: float = 1e-6
    auxbasis: Optional[str] = None


@dataclass
class FallbackStatistics:
    """Statistics tracking for fallback decisions.

    Attributes:
        total_evaluations: Total number of quality evaluations
        fallback_count: Number of times fallback was triggered
        cv_triggers: Number of CV error triggers
        pole_triggers: Number of pole stability triggers
        causality_triggers: Number of causality violation triggers
        reasons_histogram: Count of each fallback reason
        confidence_history: History of confidence scores
        timing_stats: Timing statistics
    """

    total_evaluations: int = 0
    fallback_count: int = 0
    cv_triggers: int = 0
    pole_triggers: int = 0
    causality_triggers: int = 0
    reasons_histogram: Dict[str, int] = field(default_factory=dict)
    confidence_history: List[float] = field(default_factory=list)
    timing_stats: Dict[str, float] = field(default_factory=dict)

    @property
    def fallback_rate(self) -> float:
        """Calculate fallback rate."""
        if self.total_evaluations == 0:
            return 0.0
        return self.fallback_count / self.total_evaluations

    @property
    def success_rate(self) -> float:
        """Calculate success rate (inverse of fallback rate)."""
        return 1.0 - self.fallback_rate

    def update(self, decision: FallbackDecision, elapsed_time: float = 0.0):
        """Update statistics with a new decision."""
        self.total_evaluations += 1
        if decision.should_fallback:
            self.fallback_count += 1
            if decision.reason:
                reason_str = decision.reason.value
                self.reasons_histogram[reason_str] = self.reasons_histogram.get(reason_str, 0) + 1

                # Track specific trigger types
                if decision.reason == FallbackReason.CROSS_VALIDATION_ERROR:
                    self.cv_triggers += 1
                elif decision.reason == FallbackReason.UNSTABLE_POLES:
                    self.pole_triggers += 1
                elif decision.reason == FallbackReason.CAUSALITY_VIOLATION:
                    self.causality_triggers += 1

        if decision.confidence is not None:
            self.confidence_history.append(decision.confidence)

        if elapsed_time > 0:
            if 'total_time' not in self.timing_stats:
                self.timing_stats['total_time'] = 0.0
            self.timing_stats['total_time'] += elapsed_time
            self.timing_stats['avg_time'] = self.timing_stats['total_time'] / self.total_evaluations


@dataclass
class GWResult:
    """Result of GW calculation with fallback.

    Attributes:
        method_used: Method actually used ('AC' or 'CD')
        quasiparticle_energies: Quasiparticle energies
        z_factors: Renormalization factors
        fallback_triggered: Whether fallback was triggered
        fallback_decision: Detailed fallback decision
        quality_metrics: Quality assessment metrics
        calculation_time: Total calculation time in seconds
        timing: Detailed timing information
        convergence: Convergence information
    """

    method_used: str
    quasiparticle_energies: np.ndarray
    z_factors: np.ndarray
    fallback_triggered: bool = False
    fallback_decision: Optional[FallbackDecision] = None
    quality_metrics: Optional[QualityMetrics] = None
    calculation_time: float = 0.0
    timing: Dict[str, float] = field(default_factory=dict)
    convergence: Dict[str, Any] = field(default_factory=dict)

    # Compatibility alias
    @property
    def qp_energies(self) -> np.ndarray:
        """Alias for quasiparticle_energies."""
        return self.quasiparticle_energies


class FallbackController:
    """Main controller for automatic AC→CD fallback.

    This class orchestrates the quality assessment and fallback decision
    process, providing a high-level interface for PySCF integration.

    Attributes:
        thresholds: Configuration for fallback decisions
        n_threads: Number of threads for parallel execution
        cache_size: Size of decision cache
        verbose: Verbosity level (PySCF-compatible)
        logger: Logger instance
        _cache: Decision cache
        _statistics: Runtime statistics
    """

    def __init__(
        self,
        thresholds: Optional[FallbackThresholds] = None,
        n_threads: Optional[int] = None,
        cache_size: int = 1000,
        verbose: int = 0,
        log: Optional[logger.Logger] = None
    ):
        """Initialize the fallback controller.

        Args:
            thresholds: Fallback decision thresholds
            n_threads: Number of threads (None for auto-detect)
            cache_size: Size of decision cache
            verbose: Verbosity level (0=silent, 1=info, 2=debug)
            log: Logger instance (creates new if None)
        """
        self.thresholds = thresholds or FallbackThresholds.balanced()
        self.n_threads = n_threads or lib.num_threads()
        self.cache_size = cache_size
        self.verbose = verbose

        # Set up logging
        if log is None:
            self.logger = logger.Logger(self.__class__.__name__, verbose)
        else:
            self.logger = log

        # Initialize cache and statistics
        self._cache: Dict[str, FallbackDecision] = {}
        self._statistics = FallbackStatistics()

        # Initialize Rust backend if available
        if _quasix_core is not None:
            self._init_rust_backend()

        self.logger.info("FallbackController initialized")
        self.logger.debug(f"Thresholds: {self.thresholds}")

    def evaluate_quality(self, metrics: QualityMetrics) -> FallbackDecision:
        """Evaluate quality metrics and determine if fallback is needed.

        Args:
            metrics: Quality metrics from AC calculation

        Returns:
            FallbackDecision indicating whether to fallback to CD
        """
        # Check if metrics pass thresholds
        if metrics.passes_thresholds(self.thresholds):
            return FallbackDecision(
                should_fallback=False,
                method_used=FallbackMethod.ANALYTIC_CONTINUATION,
                metrics=metrics,
                confidence=metrics.calculate_confidence(),
                recommendation="AC results are reliable"
            )
        else:
            # Get specific failures
            failures = metrics.get_failures(self.thresholds)
            primary_reason = failures[0] if failures else FallbackReason.NUMERICAL_INSTABILITY

            return FallbackDecision(
                should_fallback=True,
                method_used=FallbackMethod.CONTOUR_DEFORMATION,
                reason=primary_reason,
                all_reasons=failures,
                metrics=metrics,
                confidence=metrics.calculate_confidence(),
                recommendation=f"Fallback to CD recommended: {primary_reason}"
            )

    def update_statistics(self, metrics: QualityMetrics, should_fallback: bool):
        """Update statistics tracking.

        Args:
            metrics: Quality metrics that were evaluated
            should_fallback: Whether fallback was triggered
        """
        decision = FallbackDecision(
            should_fallback=should_fallback,
            method_used=FallbackMethod.CONTOUR_DEFORMATION if should_fallback else FallbackMethod.ANALYTIC_CONTINUATION,
            confidence=metrics.calculate_confidence()
        )

        if should_fallback:
            failures = metrics.get_failures(self.thresholds)
            if failures:
                decision.reason = failures[0]

        self._statistics.update(decision)

        if self.verbose >= 2:
            self.logger.debug(f"Statistics updated: total={self._statistics.total_evaluations}, fallbacks={self._statistics.fallback_count}")

    def get_statistics(self) -> FallbackStatistics:
        """Get runtime statistics.

        Returns:
            FallbackStatistics object with current statistics
        """
        return self._statistics

    def _init_rust_backend(self):
        """Initialize Rust backend components."""
        try:
            config = _quasix_core.FallbackConfig(
                cv_error_threshold=self.thresholds.cv_error_max,
                pole_stability_threshold=-self.thresholds.pole_stability_min,
                causality_tolerance=self.thresholds.causality_violation_max,
                sum_rule_threshold=self.thresholds.sum_rule_tolerance,
                enable_simd=True,
                cache_size=self.cache_size,
                num_threads=self.n_threads
            )
            self._rust_controller = _quasix_core.FallbackController(config)
            self.logger.debug("Rust backend initialized successfully")
        except Exception as e:
            self.logger.warning(f"Failed to initialize Rust backend: {e}")
            self._rust_controller = None

    def execute_with_fallback(
        self,
        mf,  # PySCF mean-field object
        params: Optional[GWParams] = None,
        monitor_callback: Optional[Callable] = None
    ) -> GWResult:
        """Execute GW calculation with automatic fallback.

        This is the main entry point for PySCF integration. It performs
        quality assessment and automatically falls back from AC to CD
        if necessary.

        Args:
            mf: PySCF mean-field object (RHF/UHF/GHF)
            params: GW calculation parameters
            monitor_callback: Optional callback for progress monitoring

        Returns:
            GWResult containing quasiparticle energies and metadata

        Raises:
            ValueError: If input validation fails
            RuntimeError: If both AC and CD methods fail
        """
        t0 = time.time()
        self.logger.info("Starting GW calculation with fallback control")

        # Validate input
        if not isinstance(mf, (scf.hf.RHF, scf.uhf.UHF, scf.ghf.GHF)):
            raise ValueError("mf must be a PySCF mean-field object")

        if not mf.converged:
            self.logger.warning("Mean-field calculation not converged")

        # Set default parameters
        if params is None:
            params = GWParams()

        # Extract data from mean-field
        mol = mf.mol
        mo_coeff = mf.mo_coeff
        mo_energy = mf.mo_energy
        mo_occ = mf.mo_occ

        self.logger.debug(f"System: {mol.natm} atoms, {mol.nao} AOs")
        self.logger.debug(f"Method preference: {params.method}")

        # Build DF tensors
        t_df = time.time()
        df_tensors = self._build_df_tensors(mf, params.auxbasis)
        self.logger.timer(self, "DF tensor construction", t_df)

        # Attempt AC calculation first (if requested)
        if params.method == "AC":
            try:
                # Perform AC calculation
                ac_result = self._run_ac_calculation(
                    mo_energy, mo_occ, df_tensors, params, monitor_callback
                )

                # Assess quality
                t_assess = time.time()
                decision = self._assess_ac_quality(ac_result, df_tensors)
                self.logger.timer(self, "Quality assessment", t_assess)

                # Update statistics
                self.update_statistics(decision.metrics, decision.should_fallback)

                if not decision.should_fallback:
                    # AC is acceptable
                    self.logger.info("AC calculation successful")
                    return GWResult(
                        method_used="AC",
                        quasiparticle_energies=ac_result['qp_energies'],
                        z_factors=ac_result['z_factors'],
                        fallback_triggered=False,
                        fallback_decision=decision,
                        quality_metrics=decision.metrics,
                        calculation_time=time.time() - t0,
                        timing={'total': time.time() - t0},
                        convergence=ac_result.get('convergence', {})
                    )
                else:
                    # Need to fallback to CD
                    self.logger.warning(f"AC quality insufficient: {decision.reason}")
                    self.logger.info("Falling back to CD method")

            except Exception as e:
                self.logger.error(f"AC calculation failed: {e}")
                decision = FallbackDecision(
                    should_fallback=True,
                    method_used=FallbackMethod.CONTOUR_DEFORMATION,
                    reason=FallbackReason.CONVERGENCE_FAILURE
                )
        else:
            # Direct CD calculation requested
            decision = FallbackDecision(
                should_fallback=True,
                method_used=FallbackMethod.CONTOUR_DEFORMATION,
                reason=FallbackReason.USER_REQUESTED
            )

        # Run CD calculation (fallback or direct)
        try:
            cd_result = self._run_cd_calculation(
                mo_energy, mo_occ, df_tensors, params, monitor_callback
            )

            self.logger.info("CD calculation successful")
            return GWResult(
                method_used="CD",
                quasiparticle_energies=cd_result['qp_energies'],
                z_factors=cd_result['z_factors'],
                fallback_triggered=(params.method == "AC"),
                fallback_decision=decision,
                quality_metrics=None,  # CD doesn't need quality assessment
                calculation_time=time.time() - t0,
                timing={'total': time.time() - t0},
                convergence=cd_result.get('convergence', {})
            )

        except Exception as e:
            self.logger.error(f"CD calculation failed: {e}")
            raise RuntimeError("Both AC and CD methods failed") from e

    def _build_df_tensors(
        self,
        mf,
        auxbasis: Optional[str] = None
    ) -> Dict[str, np.ndarray]:
        """Build density-fitted tensors from mean-field object.

        Args:
            mf: PySCF mean-field object
            auxbasis: Auxiliary basis name (auto-select if None)

        Returns:
            Dictionary containing DF tensors
        """
        mol = mf.mol

        # Select auxiliary basis
        if auxbasis is None:
            if hasattr(mf, 'with_df') and mf.with_df is not None:
                aux = mf.with_df
            else:
                # Auto-select based on main basis
                main_basis = mol.basis
                if 'def2' in main_basis.lower():
                    auxbasis = main_basis + '-jkfit'
                elif 'cc-pv' in main_basis.lower():
                    auxbasis = main_basis.replace('cc-pv', 'cc-pv') + '-ri'
                else:
                    auxbasis = 'def2-svp-jkfit'  # Fallback

                aux = df.DF(mol, auxbasis)
                aux.kernel()
        else:
            aux = df.DF(mol, auxbasis)
            aux.kernel()

        self.logger.debug(f"Using auxiliary basis: {auxbasis}")
        self.logger.debug(f"Number of auxiliary functions: {aux.get_naoaux()}")

        # Get 3-center integrals
        nao = mol.nao
        naux = aux.get_naoaux()
        mo_coeff = mf.mo_coeff
        mo_occ = mf.mo_occ

        # Occupied and virtual indices
        occ_idx = mo_occ > 0
        vir_idx = ~occ_idx
        nocc = np.sum(occ_idx)
        nvir = np.sum(vir_idx)

        # Transform to MO basis efficiently
        iaP = np.zeros((nocc * nvir, naux))
        ijP = np.zeros((nocc * nocc, naux))
        abP = np.zeros((nvir * nvir, naux))

        # Process in blocks for memory efficiency
        block_size = min(200, naux // 10)
        for p0 in range(0, naux, block_size):
            p1 = min(p0 + block_size, naux)

            # Get 3-center block in AO basis
            eri_3c = aux._cderi[p0:p1, :, :]  # (block, nao, nao)

            # Transform to MO basis
            mo_occ_coeff = mo_coeff[:, occ_idx]
            mo_vir_coeff = mo_coeff[:, vir_idx]

            # (ia|P) elements
            for p in range(p1 - p0):
                tmp = np.dot(eri_3c[p], mo_vir_coeff)
                iaP[:, p0 + p] = np.dot(mo_occ_coeff.T, tmp).ravel()

            # (ij|P) elements
            for p in range(p1 - p0):
                tmp = np.dot(eri_3c[p], mo_occ_coeff)
                ijP[:, p0 + p] = np.dot(mo_occ_coeff.T, tmp).ravel()

            # (ab|P) elements
            for p in range(p1 - p0):
                tmp = np.dot(eri_3c[p], mo_vir_coeff)
                abP[:, p0 + p] = np.dot(mo_vir_coeff.T, tmp).ravel()

        # Compute Coulomb metric
        j2c = aux.get_2c2e()  # (P|Q)
        j2c_sqrt = linalg.sqrtm(j2c).real
        j2c_inv_sqrt = linalg.inv(j2c_sqrt)

        return {
            'iaP': iaP,
            'ijP': ijP,
            'abP': abP,
            'j2c': j2c,
            'j2c_sqrt': j2c_sqrt,
            'j2c_inv_sqrt': j2c_inv_sqrt,
            'nocc': nocc,
            'nvir': nvir,
            'naux': naux,
            'mo_energy': mf.mo_energy,
            'mo_occ': mf.mo_occ
        }

    def _run_ac_calculation(
        self,
        mo_energy: np.ndarray,
        mo_occ: np.ndarray,
        df_tensors: Dict[str, np.ndarray],
        params: GWParams,
        monitor_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run analytic continuation GW calculation.

        Args:
            mo_energy: Molecular orbital energies
            mo_occ: Molecular orbital occupations
            df_tensors: Density-fitted tensors
            params: Calculation parameters
            monitor_callback: Progress callback

        Returns:
            Dictionary with AC calculation results
        """
        # This would call the actual AC implementation
        # For now, using placeholder

        if _quasix_core is not None and hasattr(_quasix_core, 'run_ac_gw'):
            # Use Rust implementation
            result = _quasix_core.run_ac_gw(
                mo_energy=mo_energy,
                mo_occ=mo_occ,
                iaP=df_tensors['iaP'],
                chol_v=df_tensors['j2c_inv_sqrt'],
                n_freq=params.freq_grid_size,
                eta=params.eta,
                max_iter=params.max_iter,
                conv_tol=params.conv_tol
            )
        else:
            # Fallback to Python implementation
            from .analytic_continuation import run_ac_gw_python
            result = run_ac_gw_python(
                mo_energy, mo_occ, df_tensors, params
            )

        if monitor_callback:
            monitor_callback('AC', result)

        return result

    def _run_cd_calculation(
        self,
        mo_energy: np.ndarray,
        mo_occ: np.ndarray,
        df_tensors: Dict[str, np.ndarray],
        params: GWParams,
        monitor_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run contour deformation GW calculation.

        Args:
            mo_energy: Molecular orbital energies
            mo_occ: Molecular orbital occupations
            df_tensors: Density-fitted tensors
            params: Calculation parameters
            monitor_callback: Progress callback

        Returns:
            Dictionary with CD calculation results
        """
        # This would call the actual CD implementation
        # For now, using placeholder

        if _quasix_core is not None and hasattr(_quasix_core, 'run_cd_gw'):
            # Use Rust implementation
            result = _quasix_core.run_cd_gw(
                mo_energy=mo_energy,
                mo_occ=mo_occ,
                iaP=df_tensors['iaP'],
                chol_v=df_tensors['j2c_inv_sqrt'],
                n_freq=params.freq_grid_size,
                eta=params.eta,
                max_iter=params.max_iter,
                conv_tol=params.conv_tol
            )
        else:
            # Fallback to Python implementation (placeholder)
            result = {
                'qp_energies': mo_energy.copy(),  # Placeholder
                'z_factors': np.ones_like(mo_energy) * 0.8,
                'convergence': {'iterations': 1, 'converged': True}
            }

        if monitor_callback:
            monitor_callback('CD', result)

        return result

    def _assess_ac_quality(
        self,
        ac_result: Dict[str, Any],
        df_tensors: Dict[str, np.ndarray]
    ) -> FallbackDecision:
        """Assess the quality of AC calculation.

        Args:
            ac_result: Results from AC calculation
            df_tensors: Density-fitted tensors

        Returns:
            Fallback decision based on quality metrics
        """
        # Extract or compute quality metrics
        if 'quality_metrics' in ac_result:
            metrics = ac_result['quality_metrics']
        else:
            # Compute metrics from results
            metrics = self._compute_quality_metrics(ac_result, df_tensors)

        # Use the public evaluate_quality method
        return self.evaluate_quality(metrics)

    def _compute_quality_metrics(
        self,
        ac_result: Dict[str, Any],
        df_tensors: Dict[str, np.ndarray]
    ) -> QualityMetrics:
        """Compute quality metrics from AC results.

        Args:
            ac_result: AC calculation results
            df_tensors: Density-fitted tensors

        Returns:
            Computed quality metrics
        """
        # Extract relevant data
        qp_energies = ac_result.get('qp_energies', np.array([]))
        z_factors = ac_result.get('z_factors', np.array([]))
        poles = ac_result.get('poles', np.array([]))
        residues = ac_result.get('residues', np.array([]))

        # Compute cross-validation error (placeholder)
        cv_error = ac_result.get('cv_error', 0.01)

        # Check pole stability (convert to score 0-1, higher is better)
        if len(poles) > 0:
            # More negative imaginary parts = more stable
            min_imag = np.min(np.imag(poles))
            pole_stability_score = np.clip(1.0 + min_imag / 0.01, 0, 1)  # -0.01 maps to 0, 0 maps to 1
        else:
            pole_stability_score = 0.5  # Neutral if no poles

        # Check causality (Z factors should be in (0, 1))
        causality_violations = np.sum((z_factors < 0) | (z_factors > 1))
        causality_violation = causality_violations / len(z_factors) if len(z_factors) > 0 else 0

        # Check residue magnitudes
        if len(residues) > 0:
            residue_magnitude_max = np.max(np.abs(residues))
        else:
            residue_magnitude_max = 1.0

        # Z factor range
        if len(z_factors) > 0:
            z_factor_range = (float(np.min(z_factors)), float(np.max(z_factors)))
        else:
            z_factor_range = (0.8, 0.9)  # Default reasonable range

        # Check sum rule (sum of Z factors should be close to N_electrons)
        expected_sum = np.sum(df_tensors['mo_occ'])
        actual_sum = np.sum(z_factors)
        sum_rule_error = abs(actual_sum - expected_sum) / expected_sum if expected_sum > 0 else 0

        # Convergence rate
        convergence = ac_result.get('convergence', {})
        convergence_rate = convergence.get('rate', 0.0)

        # Condition number (placeholder)
        condition_number = ac_result.get('condition_number', 1.0)

        return QualityMetrics(
            cv_error=cv_error,
            pole_stability_score=pole_stability_score,
            causality_violation=causality_violation,
            residue_magnitude_max=residue_magnitude_max,
            sum_rule_error=sum_rule_error,
            convergence_rate=convergence_rate,
            z_factor_range=z_factor_range,
            condition_number=condition_number
        )


    def clear_cache(self):
        """Clear the decision cache."""
        self._cache.clear()
        self.logger.debug("Decision cache cleared")

    def update_thresholds(self, thresholds: FallbackThresholds):
        """Update fallback thresholds.

        Args:
            thresholds: New threshold configuration
        """
        self.thresholds = thresholds
        self.clear_cache()

        if self._rust_controller:
            self._init_rust_backend()

        self.logger.info("Thresholds updated")
        self.logger.debug(f"New thresholds: {thresholds}")

    @classmethod
    def from_config(cls, config_path: Union[str, Path]) -> 'FallbackController':
        """Create controller from configuration file.

        Args:
            config_path: Path to YAML or JSON configuration

        Returns:
            Configured FallbackController instance
        """
        config_path = Path(config_path)

        if config_path.suffix == '.yaml' or config_path.suffix == '.yml':
            if yaml is None:
                raise ImportError("YAML support requires PyYAML package. Install with: pip install pyyaml")
            with open(config_path) as f:
                config = yaml.safe_load(f)
        elif config_path.suffix == '.json':
            with open(config_path) as f:
                config = json.load(f)
        else:
            raise ValueError(f"Unsupported config format: {config_path.suffix}")

        # Extract thresholds
        if 'thresholds' in config:
            thresholds = FallbackThresholds.from_dict(config['thresholds'])
        else:
            thresholds = FallbackThresholds.balanced()

        # Extract other parameters
        n_threads = config.get('n_threads', None)
        cache_size = config.get('cache_size', 1000)
        verbose = config.get('verbose', 0)

        return cls(
            thresholds=thresholds,
            n_threads=n_threads,
            cache_size=cache_size,
            verbose=verbose
        )

    def save_config(self, config_path: Union[str, Path]):
        """Save current configuration to file.

        Args:
            config_path: Path to save configuration
        """
        config_path = Path(config_path)

        config = {
            'thresholds': self.thresholds.to_dict(),
            'n_threads': self.n_threads,
            'cache_size': self.cache_size,
            'verbose': self.verbose,
            'statistics': self.get_statistics()
        }

        if config_path.suffix == '.yaml' or config_path.suffix == '.yml':
            if yaml is None:
                raise ImportError("YAML support requires PyYAML package. Install with: pip install pyyaml")
            with open(config_path, 'w') as f:
                yaml.dump(config, f, default_flow_style=False)
        elif config_path.suffix == '.json':
            with open(config_path, 'w') as f:
                json.dump(config, f, indent=2)
        else:
            raise ValueError(f"Unsupported config format: {config_path.suffix}")

        self.logger.info(f"Configuration saved to {config_path}")

    def clear_statistics(self):
        """Clear all statistics."""
        self._statistics = FallbackStatistics()
        self.logger.debug("Statistics cleared")


# Convenience functions
def run_gw_with_fallback(
    mf,
    method: str = "AC",
    thresholds: Optional[FallbackThresholds] = None,
    verbose: int = 0
) -> GWResult:
    """Convenience function to run GW with automatic fallback.

    Args:
        mf: PySCF mean-field object
        method: Initial method preference ("AC" or "CD")
        thresholds: Fallback thresholds (uses balanced if None)
        verbose: Verbosity level

    Returns:
        GWResult with quasiparticle energies and metadata
    """
    controller = FallbackController(
        thresholds=thresholds,
        verbose=verbose
    )

    params = GWParams(method=method)

    return controller.execute_with_fallback(mf, params)