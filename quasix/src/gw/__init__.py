"""QuasiX GW module - Python interface to high-performance GW/BSE calculations.

This module provides PySCF-compatible interfaces to QuasiX's Rust-based
GW and BSE implementations. Following the QuasiX principle of "tight Python
integration, coarse crossings", all heavy computation is performed in Rust
with Python serving as the interface layer.

Modules:
    analytic_continuation: Analytic continuation from imaginary to real axis
    fallback_controller: Automatic AC→CD fallback for robust GW calculations

Classes:
    AnalyticContinuation: Main interface for analytic continuation
    FallbackController: Controller for automatic AC→CD fallback
    FallbackThresholds: Configuration for fallback thresholds
    QualityMetrics: Quality assessment metrics for AC models
    FallbackDecision: Decision result from fallback evaluation
    FallbackStatistics: Statistics tracking for fallback decisions
    GWResult: Result container for GW calculations with fallback

Functions:
    pade_continuation: Convenience function for Padé continuation
    multipole_continuation: Convenience function for multipole continuation
    auto_continuation: Automatic model selection via cross-validation
    run_gw_with_fallback: Convenience function for GW with automatic fallback

Example:
    >>> from quasix.gw import FallbackController, FallbackThresholds
    >>> controller = FallbackController(thresholds=FallbackThresholds.balanced())
    >>> result = controller.execute_with_fallback(mf)
    >>> print(f"QP energies: {result.qp_energies}")
    >>> print(f"Method used: {result.method_used}")
"""

from .analytic_continuation import (
    AnalyticContinuation,
    perform_analytic_continuation,
)

from .fallback_controller import (
    FallbackController,
    FallbackThresholds,
    FallbackReason,
    FallbackMethod,
    QualityMetrics,
    FallbackDecision,
    FallbackStatistics,
    GWParams,
    GWResult,
    run_gw_with_fallback,
)

__all__ = [
    # Analytic continuation
    'AnalyticContinuation',
    'perform_analytic_continuation',
    # Fallback controller
    'FallbackController',
    'FallbackThresholds',
    'FallbackReason',
    'FallbackMethod',
    'QualityMetrics',
    'FallbackDecision',
    'FallbackStatistics',
    'GWParams',
    'GWResult',
    'run_gw_with_fallback',
]

# Version information
__version__ = '0.1.0'