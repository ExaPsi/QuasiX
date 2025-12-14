"""
Dynamic polarizability denominator updates for evGW iterations.

This module provides high-level Python interfaces for updating P0 denominators
during evGW self-consistency cycles using current quasiparticle energies.
"""

import numpy as np
from typing import Optional, Dict, Any, Tuple
from dataclasses import dataclass
import logging

# Import Rust bindings
try:
    from .quasix import (
        update_polarizability_denominators as _update_denominators_rust,
        PyPolarizabilityBuilder as _PolarizabilityBuilder,
        PyGapStatistics as _GapStatistics,
    )
    RUST_AVAILABLE = True
    # Expose PolarizabilityBuilder publicly
    PolarizabilityBuilder = _PolarizabilityBuilder
except ImportError:
    RUST_AVAILABLE = False
    _update_denominators_rust = None
    _PolarizabilityBuilder = None
    _GapStatistics = None
    PolarizabilityBuilder = None

__all__ = [
    'GapStatistics',
    'PolarizabilityUpdater',
    'PolarizabilityBuilder',
    'update_polarizability_denominators',
    'analyze_gap_evolution',
]

# Configure logging
log = logging.getLogger(__name__)


@dataclass
class GapStatistics:
    """Statistics for gap evolution monitoring during evGW."""
    min_gap: float
    max_gap: float
    mean_gap: float
    n_thresholded: int
    n_negative: int
    n_total: int

    def __str__(self) -> str:
        return (f"Gap Statistics:\n"
                f"  Min: {self.min_gap:.6f} Ha ({self.min_gap*27.2114:.3f} eV)\n"
                f"  Max: {self.max_gap:.6f} Ha ({self.max_gap*27.2114:.3f} eV)\n"
                f"  Mean: {self.mean_gap:.6f} Ha ({self.mean_gap*27.2114:.3f} eV)\n"
                f"  Thresholded: {self.n_thresholded}/{self.n_total}\n"
                f"  Negative: {self.n_negative}")

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'GapStatistics':
        """Create from dictionary returned by Rust."""
        return cls(
            min_gap=data['min_gap'],
            max_gap=data['max_gap'],
            mean_gap=data['mean_gap'],
            n_thresholded=data['n_thresholded'],
            n_negative=data['n_negative'],
            n_total=data['n_total']
        )

    def has_issues(self) -> bool:
        """Check if there are any problematic gaps."""
        return self.n_negative > 0 or self.n_thresholded > 0

    def log_summary(self, level: int = logging.INFO) -> None:
        """Log gap statistics summary."""
        log.log(level, f"Gap range: [{self.min_gap:.6f}, {self.max_gap:.6f}] Ha, "
                       f"mean: {self.mean_gap:.6f} Ha")
        if self.n_thresholded > 0:
            log.warning(f"{self.n_thresholded}/{self.n_total} gaps required thresholding")
        if self.n_negative > 0:
            log.error(f"{self.n_negative} negative gaps detected - check orbital ordering!")


class PolarizabilityUpdater:
    """
    Manager for dynamic P0 denominator updates during evGW.

    This class wraps the Rust PolarizabilityBuilder and provides
    Python-friendly interfaces for updating energy denominators
    during evGW iterations.

    Parameters
    ----------
    nocc : int
        Number of occupied orbitals
    nvirt : int
        Number of virtual orbitals
    naux : int
        Number of auxiliary basis functions
    initial_mo_energies : np.ndarray
        Initial molecular orbital energies (nbasis,)
    df_ia : np.ndarray
        Density-fitted 3-center integrals (ia|P) shape (nocc*nvirt, naux)
    gap_threshold : float, optional
        Minimum gap threshold to avoid singularities (default: 1e-6 Ha)
    eta : float, optional
        Broadening parameter for denominators (default: 1e-4 Ha)
    """

    def __init__(
        self,
        nocc: int,
        nvirt: int,
        naux: int,
        initial_mo_energies: np.ndarray,
        df_ia: np.ndarray,
        gap_threshold: float = 1e-6,
        eta: float = 1e-4
    ):
        self.nocc = nocc
        self.nvirt = nvirt
        self.naux = naux
        self.nbasis = nocc + nvirt

        # Validate inputs
        if initial_mo_energies.shape[0] != self.nbasis:
            raise ValueError(f"MO energies shape {initial_mo_energies.shape} doesn't match "
                           f"nocc+nvirt={self.nbasis}")

        expected_shape = (nocc * nvirt, naux)
        if df_ia.shape != expected_shape:
            raise ValueError(f"DF tensor shape {df_ia.shape} doesn't match "
                           f"expected {expected_shape}")

        # Store initial state
        self.initial_mo_energies = initial_mo_energies.copy()
        self.current_mo_energies = initial_mo_energies.copy()
        self.df_ia = df_ia
        self.gap_threshold = gap_threshold
        self.eta = eta

        # Initialize Rust builder if available
        if RUST_AVAILABLE and _PolarizabilityBuilder is not None:
            self._builder = _PolarizabilityBuilder(
                nocc, nvirt, naux,
                initial_mo_energies,
                df_ia,
                gap_threshold,
                eta
            )
        else:
            self._builder = None
            log.warning("Rust PolarizabilityBuilder not available, using Python fallback")

        # Track gap evolution
        self.gap_history = []

    def update_energies(
        self,
        qp_energies: np.ndarray,
        track_history: bool = True
    ) -> GapStatistics:
        """
        Update energy denominators with current quasiparticle energies.

        Parameters
        ----------
        qp_energies : np.ndarray
            Current quasiparticle energies (nbasis,)
        track_history : bool
            Whether to track gap evolution history

        Returns
        -------
        GapStatistics
            Statistics about the updated gaps
        """
        if qp_energies.shape[0] != self.nbasis:
            raise ValueError(f"QP energies shape {qp_energies.shape} doesn't match "
                           f"nbasis={self.nbasis}")

        # Update stored energies
        self.current_mo_energies = qp_energies.copy()

        # Use Rust implementation if available
        if self._builder is not None:
            stats_dict = self._builder.update_energies(qp_energies)
            stats = GapStatistics.from_dict(stats_dict)
        else:
            # Python fallback
            stats = self._compute_gap_statistics_python(qp_energies)

        # Track history if requested
        if track_history:
            self.gap_history.append(stats)

        # Log statistics
        stats.log_summary()

        return stats

    def _compute_gap_statistics_python(self, qp_energies: np.ndarray) -> GapStatistics:
        """Python fallback for gap statistics computation."""
        occ_energies = qp_energies[:self.nocc]
        virt_energies = qp_energies[self.nocc:]

        # Compute all gaps
        gaps = []
        for i in range(self.nocc):
            for a in range(self.nvirt):
                gap = virt_energies[a] - occ_energies[i]
                gaps.append(gap)

        gaps = np.array(gaps)

        # Compute statistics
        n_negative = np.sum(gaps < 0)
        n_thresholded = np.sum(gaps < self.gap_threshold)

        # Apply thresholding
        gaps_safe = np.maximum(gaps, self.gap_threshold)

        return GapStatistics(
            min_gap=float(np.min(gaps)),
            max_gap=float(np.max(gaps)),
            mean_gap=float(np.mean(gaps)),
            n_thresholded=int(n_thresholded),
            n_negative=int(n_negative),
            n_total=len(gaps)
        )

    def get_gap_evolution(self) -> Dict[str, np.ndarray]:
        """
        Get gap evolution history.

        Returns
        -------
        dict
            Dictionary with arrays of min_gap, max_gap, mean_gap,
            n_thresholded, n_negative over iterations
        """
        if not self.gap_history:
            return {
                'min_gap': np.array([]),
                'max_gap': np.array([]),
                'mean_gap': np.array([]),
                'n_thresholded': np.array([]),
                'n_negative': np.array([])
            }

        return {
            'min_gap': np.array([s.min_gap for s in self.gap_history]),
            'max_gap': np.array([s.max_gap for s in self.gap_history]),
            'mean_gap': np.array([s.mean_gap for s in self.gap_history]),
            'n_thresholded': np.array([s.n_thresholded for s in self.gap_history]),
            'n_negative': np.array([s.n_negative for s in self.gap_history])
        }

    @property
    def occupied_energies(self) -> np.ndarray:
        """Current occupied orbital energies."""
        return self.current_mo_energies[:self.nocc]

    @property
    def virtual_energies(self) -> np.ndarray:
        """Current virtual orbital energies."""
        return self.current_mo_energies[self.nocc:]


def update_polarizability_denominators(
    p0_builder: Optional[PolarizabilityUpdater],
    current_qp_energies: np.ndarray,
    n_occ: int,
    gap_threshold: float = 1e-6
) -> GapStatistics:
    """
    Update P0 energy denominators with current QP energies.

    This is a convenience function for updating denominators without
    creating a full PolarizabilityUpdater object.

    Parameters
    ----------
    p0_builder : PolarizabilityUpdater or None
        Existing builder to update, or None to just compute statistics
    current_qp_energies : np.ndarray
        Current quasiparticle energies
    n_occ : int
        Number of occupied orbitals
    gap_threshold : float
        Minimum gap threshold

    Returns
    -------
    GapStatistics
        Statistics about the gaps
    """
    if p0_builder is not None:
        return p0_builder.update_energies(current_qp_energies)

    # Just compute statistics without builder
    occ_energies = current_qp_energies[:n_occ]
    virt_energies = current_qp_energies[n_occ:]

    gaps = []
    for e_occ in occ_energies:
        for e_virt in virt_energies:
            gaps.append(e_virt - e_occ)

    gaps = np.array(gaps)

    return GapStatistics(
        min_gap=float(np.min(gaps)),
        max_gap=float(np.max(gaps)),
        mean_gap=float(np.mean(gaps)),
        n_thresholded=int(np.sum(gaps < gap_threshold)),
        n_negative=int(np.sum(gaps < 0)),
        n_total=len(gaps)
    )


def analyze_gap_evolution(
    gap_history: list[GapStatistics],
    output_file: Optional[str] = None
) -> Dict[str, Any]:
    """
    Analyze gap evolution over evGW iterations.

    Parameters
    ----------
    gap_history : list[GapStatistics]
        History of gap statistics from evGW iterations
    output_file : str, optional
        If provided, save analysis to file

    Returns
    -------
    dict
        Analysis results including trends, warnings, and recommendations
    """
    if not gap_history:
        return {'status': 'no_data', 'message': 'No gap history available'}

    analysis = {
        'n_iterations': len(gap_history),
        'initial_min_gap': gap_history[0].min_gap,
        'final_min_gap': gap_history[-1].min_gap,
        'initial_mean_gap': gap_history[0].mean_gap,
        'final_mean_gap': gap_history[-1].mean_gap,
        'gap_closure_detected': any(s.n_negative > 0 for s in gap_history),
        'thresholding_required': any(s.n_thresholded > 0 for s in gap_history),
        'warnings': []
    }

    # Check for gap closure
    if analysis['gap_closure_detected']:
        first_negative = next(i for i, s in enumerate(gap_history) if s.n_negative > 0)
        analysis['warnings'].append(
            f"Gap closure detected at iteration {first_negative+1}"
        )

    # Check for excessive thresholding
    max_thresholded = max(s.n_thresholded for s in gap_history)
    if max_thresholded > 0:
        analysis['warnings'].append(
            f"Up to {max_thresholded} gaps required thresholding"
        )

    # Check gap trend
    gap_trend = gap_history[-1].mean_gap - gap_history[0].mean_gap
    if gap_trend < -0.1:  # Gap closing by more than 0.1 Ha
        analysis['warnings'].append(
            f"Significant gap reduction: {gap_trend:.3f} Ha"
        )

    # Recommendations
    analysis['recommendations'] = []
    if analysis['gap_closure_detected']:
        analysis['recommendations'].append(
            "Consider using a larger basis set or checking orbital ordering"
        )
    if max_thresholded > 10:
        analysis['recommendations'].append(
            "Consider increasing gap_threshold or using damping"
        )

    # Save to file if requested
    if output_file:
        import json
        with open(output_file, 'w') as f:
            # Convert GapStatistics to dict for JSON serialization
            json_data = analysis.copy()
            json_data['history'] = [
                {
                    'iteration': i,
                    'min_gap': s.min_gap,
                    'max_gap': s.max_gap,
                    'mean_gap': s.mean_gap,
                    'n_thresholded': s.n_thresholded,
                    'n_negative': s.n_negative
                }
                for i, s in enumerate(gap_history)
            ]
            json.dump(json_data, f, indent=2)
            log.info(f"Gap evolution analysis saved to {output_file}")

    return analysis