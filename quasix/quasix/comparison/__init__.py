"""
CD vs AC comparison harness for QuasiX GW calculations.

This module provides tools for comparing Contour Deformation (CD) and
Analytic Continuation (AC) methods in GW calculations.
"""

from .cd_vs_ac import (
    ComparisonConfig,
    ComparisonResult,
    CDvsACComparator,
    ComparisonReport
)

from .report_generator import HTMLReportGenerator

from .plotting import (
    plot_correlation,
    plot_error_distribution,
    plot_timing_comparison,
    plot_convergence,
    plot_z_factors,
    plot_orbital_comparison,
    create_summary_figure
)

__all__ = [
    # Core classes
    'ComparisonConfig',
    'ComparisonResult',
    'CDvsACComparator',
    'ComparisonReport',

    # Report generation
    'HTMLReportGenerator',

    # Plotting utilities
    'plot_correlation',
    'plot_error_distribution',
    'plot_timing_comparison',
    'plot_convergence',
    'plot_z_factors',
    'plot_orbital_comparison',
    'create_summary_figure'
]