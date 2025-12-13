"""Visualization Tools for Fallback Controller Quality Metrics

This module provides comprehensive visualization utilities for analyzing
fallback decisions, quality metrics, and optimization of thresholds.

Author: QuasiX Development Team
License: MIT
"""

from typing import Optional, List, Dict, Any, Tuple, Union
from pathlib import Path
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

try:
    import seaborn as sns
    # Set style for publication-quality plots
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")
except ImportError:
    # Fallback to basic matplotlib style if seaborn not available
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except:
        plt.style.use('default')
    sns = None

from .fallback_controller import (
    QualityMetrics, FallbackDecision, FallbackThresholds,
    FallbackController, GWResult
)


class FallbackVisualizer:
    """Main visualization class for fallback controller analysis."""

    def __init__(self, figsize: Tuple[int, int] = (12, 8), dpi: int = 150):
        """Initialize visualizer with default settings.

        Args:
            figsize: Default figure size
            dpi: Default DPI for plots
        """
        self.figsize = figsize
        self.dpi = dpi
        self.colors = {
            'AC': '#2E86AB',
            'CD': '#A23B72',
            'threshold': '#F18F01',
            'safe': '#73AB84',
            'unsafe': '#C73E1D'
        }

    def plot_quality_metrics(
        self,
        metrics: QualityMetrics,
        thresholds: Optional[FallbackThresholds] = None,
        save_path: Optional[Path] = None
    ) -> plt.Figure:
        """Create comprehensive quality metrics visualization.

        Args:
            metrics: Quality metrics to visualize
            thresholds: Optional thresholds for comparison
            save_path: Optional path to save figure

        Returns:
            Matplotlib figure
        """
        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.25)

        # Cross-validation error
        ax1 = fig.add_subplot(gs[0, 0])
        self._plot_metric_gauge(
            ax1, metrics.cv_error,
            threshold=thresholds.cv_error_max if thresholds else None,
            title='Cross-Validation Error',
            label='CV Error'
        )

        # Pole stability
        ax2 = fig.add_subplot(gs[0, 1])
        self._plot_pole_stability(
            ax2, metrics.pole_stability_score,
            threshold=0.7 if thresholds else None  # Use fixed threshold for score
        )

        # Causality violation
        ax3 = fig.add_subplot(gs[0, 2])
        self._plot_metric_gauge(
            ax3, metrics.causality_violation,
            threshold=thresholds.causality_violation_max if thresholds else None,
            title='Causality Violation',
            label='Score'
        )

        # Residue magnitude
        ax4 = fig.add_subplot(gs[1, 0])
        self._plot_residue_gauge(
            ax4, metrics.residue_magnitude_max,
            threshold=thresholds.residue_magnitude_max if thresholds else None
        )

        # Sum rule error
        ax5 = fig.add_subplot(gs[1, 1])
        self._plot_metric_gauge(
            ax5, metrics.sum_rule_error,
            threshold=thresholds.sum_rule_tolerance if thresholds else None,
            title='Sum Rule Error',
            label='Error'
        )

        # Overall quality score
        ax6 = fig.add_subplot(gs[1, 2])
        if thresholds:
            self._plot_overall_quality(ax6, metrics, thresholds)
        else:
            ax6.text(0.5, 0.5, 'Thresholds\nNot Provided',
                    ha='center', va='center', fontsize=14)
            ax6.set_xlim(0, 1)
            ax6.set_ylim(0, 1)
            ax6.axis('off')

        fig.suptitle('Quality Metrics Assessment', fontsize=16, fontweight='bold')

        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)

        return fig

    def _plot_metric_gauge(
        self,
        ax: plt.Axes,
        value: float,
        threshold: Optional[float],
        title: str,
        label: str
    ):
        """Create gauge plot for single metric."""
        # Create semi-circular gauge
        theta = np.linspace(np.pi, 0, 100)
        r_outer = 1.0
        r_inner = 0.7

        # Background
        for i, (t1, t2) in enumerate(zip(theta[:-1], theta[1:])):
            color = plt.cm.RdYlGn_r(i / 100)
            vertices = [
                (r_inner * np.cos(t1), r_inner * np.sin(t1)),
                (r_outer * np.cos(t1), r_outer * np.sin(t1)),
                (r_outer * np.cos(t2), r_outer * np.sin(t2)),
                (r_inner * np.cos(t2), r_inner * np.sin(t2))
            ]
            poly = plt.Polygon(vertices, facecolor=color, edgecolor='none')
            ax.add_patch(poly)

        # Value indicator
        if threshold:
            normalized_value = min(value / threshold, 2.0)
        else:
            normalized_value = value

        angle = np.pi * (1 - normalized_value / 2)
        ax.plot([0, 0.85 * np.cos(angle)], [0, 0.85 * np.sin(angle)],
               'k-', linewidth=3)
        ax.plot(0, 0, 'ko', markersize=10)

        # Threshold line
        if threshold:
            threshold_angle = np.pi / 2
            ax.plot([r_inner * np.cos(threshold_angle), r_outer * np.cos(threshold_angle)],
                   [r_inner * np.sin(threshold_angle), r_outer * np.sin(threshold_angle)],
                   'r--', linewidth=2, label='Threshold')

        # Labels
        ax.text(0, -0.3, f'{value:.3f}', ha='center', fontsize=14, fontweight='bold')
        ax.text(0, 1.3, title, ha='center', fontsize=12, fontweight='bold')

        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-0.4, 1.4)
        ax.set_aspect('equal')
        ax.axis('off')

    def _plot_pole_stability(self, ax: plt.Axes, stability_score: float, threshold: Optional[float]):
        """Plot pole stability score indicator."""
        # Create stability bar (score is 0-1, higher is better)
        y_positions = np.array([0])
        bar = ax.barh(y_positions, [stability_score], height=0.3)

        # Color based on threshold
        if threshold and stability_score >= threshold:
            bar[0].set_color(self.colors['safe'])
        else:
            bar[0].set_color(self.colors['unsafe'])

        # Add threshold line
        if threshold:
            ax.axvline(threshold, color='r', linestyle='--', linewidth=2,
                      label=f'Threshold: {threshold:.2f}')

        ax.set_xlim(0, 1.0)
        ax.set_ylim(-0.5, 0.5)
        ax.set_xlabel('Stability Score')
        ax.set_title('Pole Stability', fontweight='bold')
        ax.set_yticks([])
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)

    def _plot_residue_gauge(
        self,
        ax: plt.Axes,
        residue_max: float,
        threshold: Optional[float]
    ):
        """Plot residue magnitude gauge."""
        self._plot_metric_gauge(
            ax, residue_max,
            threshold=threshold,
            title='Max Residue Magnitude',
            label='Magnitude'
        )

    def _plot_overall_quality(
        self,
        ax: plt.Axes,
        metrics: QualityMetrics,
        thresholds: FallbackThresholds
    ):
        """Plot overall quality assessment."""
        # Compute individual pass/fail
        z_min, z_max = metrics.z_factor_range
        checks = {
            'CV Error': metrics.cv_error <= thresholds.cv_error_max,
            'Pole Stability': metrics.pole_stability_score >= 0.7,
            'Causality': metrics.causality_violation <= thresholds.causality_violation_max,
            'Residues': metrics.residue_magnitude_max <= thresholds.residue_magnitude_max,
            'Sum Rule': metrics.sum_rule_error <= thresholds.sum_rule_tolerance,
            'Z Factors': (z_min >= 0.0 and z_max <= 1.0)
        }

        # Create status chart
        y_pos = np.arange(len(checks))
        colors = [self.colors['safe'] if v else self.colors['unsafe']
                 for v in checks.values()]

        bars = ax.barh(y_pos, [1] * len(checks), color=colors)

        # Add check/cross marks
        for i, (name, passed) in enumerate(checks.items()):
            symbol = '✓' if passed else '✗'
            ax.text(0.5, i, symbol, ha='center', va='center',
                   fontsize=20, fontweight='bold', color='white')

        ax.set_yticks(y_pos)
        ax.set_yticklabels(checks.keys())
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_title('Quality Checks', fontweight='bold')

        # Overall verdict
        overall_pass = all(checks.values())
        verdict_color = self.colors['safe'] if overall_pass else self.colors['unsafe']
        verdict_text = 'PASS' if overall_pass else 'FAIL'
        ax.text(1.1, len(checks) / 2 - 0.5, verdict_text,
               fontsize=16, fontweight='bold', color=verdict_color,
               rotation=90, va='center')

    def plot_fallback_decision(
        self,
        decision: FallbackDecision,
        save_path: Optional[Path] = None
    ) -> plt.Figure:
        """Visualize fallback decision details.

        Args:
            decision: Fallback decision to visualize
            save_path: Optional path to save figure

        Returns:
            Matplotlib figure
        """
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=self.dpi)

        # Decision summary
        ax1 = axes[0]
        self._plot_decision_summary(ax1, decision)

        # Confidence meter
        ax2 = axes[1]
        self._plot_confidence_meter(ax2, decision.confidence)

        # Reasons breakdown
        ax3 = axes[2]
        if decision.all_reasons:
            self._plot_reasons_breakdown(ax3, decision.all_reasons)
        else:
            ax3.text(0.5, 0.5, 'No Issues\nDetected', ha='center', va='center',
                    fontsize=14, color=self.colors['safe'])
            ax3.set_xlim(0, 1)
            ax3.set_ylim(0, 1)
            ax3.axis('off')

        fig.suptitle(f'Fallback Decision: {decision.method_used.value}',
                    fontsize=16, fontweight='bold')

        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)

        return fig

    def _plot_decision_summary(self, ax: plt.Axes, decision: FallbackDecision):
        """Plot decision summary panel."""
        ax.axis('off')

        # Title
        method_color = self.colors['CD'] if decision.should_fallback else self.colors['AC']
        ax.text(0.5, 0.9, f"Method: {decision.method_used.value}",
               ha='center', fontsize=16, fontweight='bold', color=method_color)

        # Fallback status
        status_text = "Fallback Triggered" if decision.should_fallback else "AC Accepted"
        status_color = self.colors['unsafe'] if decision.should_fallback else self.colors['safe']
        ax.text(0.5, 0.7, status_text, ha='center', fontsize=14, color=status_color)

        # Primary reason
        if decision.reason:
            ax.text(0.5, 0.5, f"Reason: {decision.reason}", ha='center',
                   fontsize=12, style='italic')

        # Recommendation
        ax.text(0.5, 0.3, decision.recommendation, ha='center',
               fontsize=11, wrap=True)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    def _plot_confidence_meter(self, ax: plt.Axes, confidence: float):
        """Plot confidence meter."""
        # Create circular meter
        theta = np.linspace(0, 2 * np.pi, 100)
        r = 1.0

        # Background circle
        circle = plt.Circle((0, 0), r, fill=False, edgecolor='black', linewidth=2)
        ax.add_patch(circle)

        # Confidence arc
        confidence_theta = theta[:int(confidence * 100)]
        if len(confidence_theta) > 0:
            x = r * np.cos(confidence_theta)
            y = r * np.sin(confidence_theta)
            ax.fill_between(x, 0, y, alpha=0.3, color=self.colors['AC'])
            ax.plot(x, y, linewidth=3, color=self.colors['AC'])

        # Center text
        ax.text(0, 0, f'{confidence:.0%}', ha='center', va='center',
               fontsize=20, fontweight='bold')
        ax.text(0, -0.3, 'Confidence', ha='center', fontsize=12)

        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect('equal')
        ax.axis('off')

    def _plot_reasons_breakdown(self, ax: plt.Axes, reasons: List):
        """Plot breakdown of failure reasons."""
        reason_counts = {}
        for reason in reasons:
            reason_str = str(reason)
            reason_counts[reason_str] = reason_counts.get(reason_str, 0) + 1

        # Pie chart
        labels = list(reason_counts.keys())
        sizes = list(reason_counts.values())
        colors_list = plt.cm.Set3(np.linspace(0, 1, len(labels)))

        wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors_list,
                                          autopct='%1.0f%%', startangle=90)

        # Improve text visibility
        for text in texts:
            text.set_fontsize(10)
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')

        ax.set_title('Issue Breakdown', fontweight='bold')

    def plot_roc_curve(
        self,
        true_failures: np.ndarray,
        predicted_scores: np.ndarray,
        current_threshold: Optional[float] = None,
        save_path: Optional[Path] = None
    ) -> plt.Figure:
        """Plot ROC curve for threshold optimization.

        Args:
            true_failures: True binary labels (1 for should fallback)
            predicted_scores: Predicted scores for fallback
            current_threshold: Current threshold to mark
            save_path: Optional path to save figure

        Returns:
            Matplotlib figure
        """
        from sklearn.metrics import roc_curve, auc

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=self.dpi)

        # ROC Curve
        fpr, tpr, thresholds = roc_curve(true_failures, predicted_scores)
        roc_auc = auc(fpr, tpr)

        ax1.plot(fpr, tpr, linewidth=2, label=f'ROC (AUC = {roc_auc:.3f})')
        ax1.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random')

        # Mark current threshold
        if current_threshold:
            idx = np.argmin(np.abs(thresholds - current_threshold))
            ax1.plot(fpr[idx], tpr[idx], 'ro', markersize=10,
                    label=f'Current ({current_threshold:.3f})')

        ax1.set_xlabel('False Positive Rate')
        ax1.set_ylabel('True Positive Rate')
        ax1.set_title('ROC Curve', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Threshold vs Performance
        ax2.plot(thresholds, tpr, label='Sensitivity (TPR)', linewidth=2)
        ax2.plot(thresholds, 1 - fpr, label='Specificity (1-FPR)', linewidth=2)

        # F1 score
        precision = tpr / (tpr + fpr + 1e-10)
        f1_scores = 2 * (precision * tpr) / (precision + tpr + 1e-10)
        ax2.plot(thresholds, f1_scores, label='F1 Score', linewidth=2)

        if current_threshold:
            ax2.axvline(current_threshold, color='r', linestyle='--',
                       linewidth=2, label='Current Threshold')

        ax2.set_xlabel('Threshold')
        ax2.set_ylabel('Performance Metric')
        ax2.set_title('Threshold Analysis', fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        fig.suptitle('Threshold Optimization Analysis', fontsize=16, fontweight='bold')

        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)

        return fig

    def plot_convergence_monitoring(
        self,
        iterations: List[int],
        metrics_history: List[QualityMetrics],
        save_path: Optional[Path] = None
    ) -> plt.Figure:
        """Plot convergence monitoring during evGW.

        Args:
            iterations: Iteration numbers
            metrics_history: Quality metrics at each iteration
            save_path: Optional path to save figure

        Returns:
            Matplotlib figure
        """
        fig, axes = plt.subplots(2, 2, figsize=self.figsize, dpi=self.dpi)

        # CV Error evolution
        ax1 = axes[0, 0]
        cv_errors = [m.cv_error for m in metrics_history]
        ax1.plot(iterations, cv_errors, 'o-', linewidth=2, markersize=6)
        ax1.set_xlabel('Iteration')
        ax1.set_ylabel('CV Error')
        ax1.set_title('Cross-Validation Error', fontweight='bold')
        ax1.grid(True, alpha=0.3)

        # Pole stability evolution
        ax2 = axes[0, 1]
        stabilities = [m.pole_stability_score for m in metrics_history]
        ax2.plot(iterations, stabilities, 'o-', linewidth=2, markersize=6)
        ax2.set_xlabel('Iteration')
        ax2.set_ylabel('Stability Score')
        ax2.set_title('Pole Stability', fontweight='bold')
        ax2.set_ylim(0, 1)
        ax2.grid(True, alpha=0.3)

        # Causality violation evolution
        ax3 = axes[1, 0]
        causality_violations = [m.causality_violation for m in metrics_history]
        ax3.plot(iterations, causality_violations, 'o-', linewidth=2, markersize=6)
        ax3.set_xlabel('Iteration')
        ax3.set_ylabel('Causality Violation')
        ax3.set_title('Causality Violations', fontweight='bold')
        ax3.grid(True, alpha=0.3)

        # Overall quality indicator
        ax4 = axes[1, 1]
        # Simple quality score (higher is better)
        quality_scores = [
            1.0 - (m.cv_error + m.causality_violation + m.sum_rule_error) / 3
            for m in metrics_history
        ]
        ax4.plot(iterations, quality_scores, 'o-', linewidth=2, markersize=6)
        ax4.set_xlabel('Iteration')
        ax4.set_ylabel('Quality Score')
        ax4.set_title('Overall Quality', fontweight='bold')
        ax4.set_ylim(0, 1)
        ax4.grid(True, alpha=0.3)

        fig.suptitle('Convergence Monitoring', fontsize=16, fontweight='bold')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)

        return fig

    def plot_statistics_summary(
        self,
        controller: FallbackController,
        save_path: Optional[Path] = None
    ) -> plt.Figure:
        """Plot summary statistics from controller.

        Args:
            controller: Fallback controller with statistics
            save_path: Optional path to save figure

        Returns:
            Matplotlib figure
        """
        stats = controller.get_statistics()

        fig, axes = plt.subplots(2, 2, figsize=self.figsize, dpi=self.dpi)

        # Fallback rate pie chart
        ax1 = axes[0, 0]
        if stats['total_assessments'] > 0:
            fallback_rate = stats['fallback_rate']
            sizes = [fallback_rate, 1 - fallback_rate]
            labels = ['Fallback (CD)', 'Direct (AC)']
            colors = [self.colors['CD'], self.colors['AC']]
            ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                   startangle=90)
            ax1.set_title('Method Distribution', fontweight='bold')
        else:
            ax1.text(0.5, 0.5, 'No Data', ha='center', va='center', fontsize=14)
            ax1.set_xlim(0, 1)
            ax1.set_ylim(0, 1)
            ax1.axis('off')

        # Reasons breakdown
        ax2 = axes[0, 1]
        if stats['reasons_count']:
            reasons = list(stats['reasons_count'].keys())
            counts = list(stats['reasons_count'].values())
            y_pos = np.arange(len(reasons))

            ax2.barh(y_pos, counts, color=self.colors['unsafe'])
            ax2.set_yticks(y_pos)
            ax2.set_yticklabels(reasons)
            ax2.set_xlabel('Count')
            ax2.set_title('Fallback Reasons', fontweight='bold')
            ax2.grid(True, alpha=0.3, axis='x')
        else:
            ax2.text(0.5, 0.5, 'No Fallbacks', ha='center', va='center', fontsize=14)
            ax2.set_xlim(0, 1)
            ax2.set_ylim(0, 1)
            ax2.axis('off')

        # Assessment count
        ax3 = axes[1, 0]
        ax3.axis('off')
        metrics_text = f"""
        Total Assessments: {stats.total_evaluations}
        Fallbacks Triggered: {stats.fallback_count}
        Cache Hits: {stats.get('cache_hits', 0)}
        Fallback Rate: {stats.fallback_rate:.1%}
        """
        ax3.text(0.5, 0.5, metrics_text, ha='center', va='center',
                fontsize=12, family='monospace')
        ax3.set_title('Summary Statistics', fontweight='bold')

        # Timing distribution (if available)
        ax4 = axes[1, 1]
        if 'assessment_times' in stats and stats['assessment_times']:
            times = stats['assessment_times']
            ax4.hist(times, bins=20, edgecolor='black', alpha=0.7)
            ax4.set_xlabel('Time (ms)')
            ax4.set_ylabel('Frequency')
            ax4.set_title('Assessment Timing', fontweight='bold')
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(0.5, 0.5, 'No Timing Data', ha='center', va='center', fontsize=14)
            ax4.set_xlim(0, 1)
            ax4.set_ylim(0, 1)
            ax4.axis('off')

        fig.suptitle('Controller Statistics Summary', fontsize=16, fontweight='bold')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)

        return fig


def create_comparison_plot(
    ac_result: GWResult,
    cd_result: GWResult,
    save_path: Optional[Path] = None
) -> plt.Figure:
    """Create comparison plot between AC and CD results.

    Args:
        ac_result: Results from AC calculation
        cd_result: Results from CD calculation
        save_path: Optional path to save figure

    Returns:
        Matplotlib figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # QP energies comparison
    ax1 = axes[0, 0]
    ax1.scatter(ac_result.qp_energies, cd_result.qp_energies, alpha=0.6)

    # Perfect agreement line
    energy_range = [min(ac_result.qp_energies.min(), cd_result.qp_energies.min()),
                   max(ac_result.qp_energies.max(), cd_result.qp_energies.max())]
    ax1.plot(energy_range, energy_range, 'k--', alpha=0.5, label='Perfect Agreement')

    ax1.set_xlabel('AC QP Energies (eV)')
    ax1.set_ylabel('CD QP Energies (eV)')
    ax1.set_title('Quasiparticle Energies', fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Z factors comparison
    ax2 = axes[0, 1]
    ax2.scatter(ac_result.z_factors, cd_result.z_factors, alpha=0.6)
    ax2.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Perfect Agreement')
    ax2.set_xlabel('AC Z Factors')
    ax2.set_ylabel('CD Z Factors')
    ax2.set_title('Renormalization Factors', fontweight='bold')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Energy differences
    ax3 = axes[1, 0]
    energy_diff = cd_result.qp_energies - ac_result.qp_energies
    ax3.hist(energy_diff, bins=30, edgecolor='black', alpha=0.7)
    ax3.axvline(0, color='r', linestyle='--', linewidth=2)
    ax3.set_xlabel('Energy Difference (CD - AC) [eV]')
    ax3.set_ylabel('Frequency')
    ax3.set_title('QP Energy Differences', fontweight='bold')

    # Add statistics
    mean_diff = np.mean(energy_diff)
    std_diff = np.std(energy_diff)
    ax3.text(0.02, 0.98, f'Mean: {mean_diff:.3f} eV\nStd: {std_diff:.3f} eV',
            transform=ax3.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax3.grid(True, alpha=0.3)

    # Z factor differences
    ax4 = axes[1, 1]
    z_diff = cd_result.z_factors - ac_result.z_factors
    ax4.hist(z_diff, bins=30, edgecolor='black', alpha=0.7)
    ax4.axvline(0, color='r', linestyle='--', linewidth=2)
    ax4.set_xlabel('Z Factor Difference (CD - AC)')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Z Factor Differences', fontweight='bold')

    # Add statistics
    mean_z_diff = np.mean(z_diff)
    std_z_diff = np.std(z_diff)
    ax4.text(0.02, 0.98, f'Mean: {mean_z_diff:.3f}\nStd: {std_z_diff:.3f}',
            transform=ax4.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax4.grid(True, alpha=0.3)

    fig.suptitle('AC vs CD Comparison', fontsize=16, fontweight='bold')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)

    return fig


# Convenience functions for quick plotting
def quick_plot_metrics(metrics: QualityMetrics, thresholds: Optional[FallbackThresholds] = None):
    """Quick plot of quality metrics."""
    viz = FallbackVisualizer()
    return viz.plot_quality_metrics(metrics, thresholds)


def quick_plot_decision(decision: FallbackDecision):
    """Quick plot of fallback decision."""
    viz = FallbackVisualizer()
    return viz.plot_fallback_decision(decision)


def quick_plot_comparison(ac_result: GWResult, cd_result: GWResult):
    """Quick comparison plot between AC and CD."""
    return create_comparison_plot(ac_result, cd_result)