"""
Plotting utilities for CD vs AC comparison visualization.

This module provides publication-quality plotting functions for comparing
Contour Deformation and Analytic Continuation GW methods.
"""

from typing import List, Dict, Optional, Tuple
import numpy as np

# Try to import matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: Matplotlib not found. Plotting functionality disabled.")

# Import scipy for statistical functions
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

from .cd_vs_ac import ComparisonResult


def set_publication_style():
    """Set matplotlib style for publication-quality plots"""
    if not HAS_MATPLOTLIB:
        return

    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 11,
        'figure.titlesize': 18,
        'font.family': 'sans-serif',
        'axes.grid': True,
        'grid.alpha': 0.3,
        'axes.spines.top': False,
        'axes.spines.right': False
    })


def plot_correlation(results: List[ComparisonResult],
                    figsize: Tuple[float, float] = (8, 8)) -> Optional['plt.Figure']:
    """
    Create correlation plot between CD and AC quasiparticle energies.

    Parameters
    ----------
    results : List[ComparisonResult]
        List of comparison results
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    matplotlib.figure.Figure or None
        The correlation plot figure
    """
    if not HAS_MATPLOTLIB:
        return None

    set_publication_style()

    fig, ax = plt.subplots(figsize=figsize)

    # Color palette for different molecules
    colors = plt.cm.Set3(np.linspace(0, 1, len(results)))

    # Plot data for each molecule
    all_cd = []
    all_ac = []

    for i, result in enumerate(results):
        # Convert to eV
        cd_energies = result.cd_qp_energies * 27.211
        ac_energies = result.ac_qp_energies * 27.211

        ax.scatter(cd_energies, ac_energies,
                  alpha=0.6, s=50, color=colors[i],
                  label=f"{result.molecule_name}/{result.basis_set}",
                  edgecolors='black', linewidth=0.5)

        all_cd.extend(cd_energies)
        all_ac.extend(ac_energies)

        # Mark outliers with red circles
        if result.outlier_indices:
            outlier_cd = cd_energies[result.outlier_indices]
            outlier_ac = ac_energies[result.outlier_indices]
            ax.scatter(outlier_cd, outlier_ac,
                      s=100, facecolors='none', edgecolors='red',
                      linewidth=2, marker='o')

    # Identity line (perfect agreement)
    all_cd = np.array(all_cd)
    all_ac = np.array(all_ac)
    min_e = min(all_cd.min(), all_ac.min())
    max_e = max(all_cd.max(), all_ac.max())
    ax.plot([min_e, max_e], [min_e, max_e], 'k--',
           alpha=0.5, linewidth=2, label='Perfect agreement')

    # Linear regression
    if HAS_SCIPY and len(all_cd) > 2:
        slope, intercept, r_value, p_value, std_err = stats.linregress(all_cd, all_ac)
        x_fit = np.linspace(min_e, max_e, 100)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, 'r-', alpha=0.5, linewidth=1.5,
               label=f'Linear fit (R² = {r_value**2:.4f})')

    ax.set_xlabel('CD Quasiparticle Energy (eV)')
    ax.set_ylabel('AC Quasiparticle Energy (eV)')
    ax.set_title('CD vs AC Method Comparison')

    # Add legend
    if len(results) <= 5:
        ax.legend(loc='best', framealpha=0.9)

    # Add grid
    ax.grid(True, alpha=0.3)

    # Equal aspect ratio
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    return fig


def plot_error_distribution(results: List[ComparisonResult],
                           figsize: Tuple[float, float] = (10, 6)) -> Optional['plt.Figure']:
    """
    Plot the distribution of energy differences between CD and AC.

    Parameters
    ----------
    results : List[ComparisonResult]
        List of comparison results
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    matplotlib.figure.Figure or None
        The error distribution plot
    """
    if not HAS_MATPLOTLIB:
        return None

    set_publication_style()

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # Collect all errors
    all_errors = []
    molecule_errors = {}

    for result in results:
        errors = (result.cd_qp_energies - result.ac_qp_energies) * 27.211  # Convert to eV
        all_errors.extend(errors)
        molecule_errors[result.molecule_name] = errors

    all_errors = np.array(all_errors)

    # Main histogram
    ax1 = fig.add_subplot(gs[:, 0])
    n, bins, patches = ax1.hist(all_errors, bins=30, alpha=0.7,
                                color='steelblue', edgecolor='black', linewidth=1)

    # Add vertical lines for mean and std
    mean_error = np.mean(all_errors)
    std_error = np.std(all_errors)

    ax1.axvline(x=0, color='green', linestyle='--', linewidth=2,
               label='Zero error', alpha=0.7)
    ax1.axvline(x=mean_error, color='red', linestyle='--', linewidth=2,
               label=f'Mean = {mean_error:.3f} eV', alpha=0.7)
    ax1.axvspan(mean_error - std_error, mean_error + std_error,
               alpha=0.2, color='red', label=f'±1σ = {std_error:.3f} eV')

    ax1.set_xlabel('Energy Difference (CD - AC) [eV]')
    ax1.set_ylabel('Count')
    ax1.set_title('Error Distribution')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)

    # Q-Q plot for normality test
    ax2 = fig.add_subplot(gs[0, 1])
    if HAS_SCIPY:
        stats.probplot(all_errors, dist="norm", plot=ax2)
        ax2.set_title('Q-Q Plot (Normality Test)')
    else:
        # Simple scatter plot as fallback
        sorted_errors = np.sort(all_errors)
        theoretical = np.linspace(all_errors.min(), all_errors.max(), len(all_errors))
        ax2.scatter(theoretical, sorted_errors, alpha=0.5)
        ax2.plot([all_errors.min(), all_errors.max()],
                [all_errors.min(), all_errors.max()], 'r--')
        ax2.set_xlabel('Theoretical Quantiles')
        ax2.set_ylabel('Sample Quantiles')
        ax2.set_title('Q-Q Plot')

    # Box plot by molecule
    ax3 = fig.add_subplot(gs[1, 1])
    if molecule_errors:
        box_data = list(molecule_errors.values())
        box_labels = list(molecule_errors.keys())

        bp = ax3.boxplot(box_data, labels=box_labels, patch_artist=True)

        # Color the box plots
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
            patch.set_alpha(0.7)

        ax3.axhline(y=0, color='green', linestyle='--', alpha=0.5)
        ax3.set_ylabel('Error (eV)')
        ax3.set_title('Error Distribution by Molecule')
        ax3.grid(True, alpha=0.3, axis='y')

        # Rotate x labels if many molecules
        if len(box_labels) > 3:
            ax3.tick_params(axis='x', rotation=45)

    plt.suptitle('CD vs AC Error Analysis', fontsize=16)
    plt.tight_layout()
    return fig


def plot_timing_comparison(results: List[ComparisonResult],
                          figsize: Tuple[float, float] = (10, 6)) -> Optional['plt.Figure']:
    """
    Create timing comparison bar chart for CD vs AC methods.

    Parameters
    ----------
    results : List[ComparisonResult]
        List of comparison results
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    matplotlib.figure.Figure or None
        The timing comparison plot
    """
    if not HAS_MATPLOTLIB:
        return None

    set_publication_style()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    molecules = [r.molecule_name for r in results]
    cd_times = [r.cd_timing for r in results]
    ac_times = [r.ac_timing for r in results]
    speedups = [ac_t / cd_t if cd_t > 0 else 0 for cd_t, ac_t in zip(cd_times, ac_times)]

    x = np.arange(len(molecules))
    width = 0.35

    # Timing comparison
    bars1 = ax1.bar(x - width/2, cd_times, width, label='CD', color='#3498db', alpha=0.8)
    bars2 = ax1.bar(x + width/2, ac_times, width, label='AC', color='#e74c3c', alpha=0.8)

    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}', ha='center', va='bottom', fontsize=9)

    for bar in bars2:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}', ha='center', va='bottom', fontsize=9)

    ax1.set_xlabel('Molecule')
    ax1.set_ylabel('Time (seconds)')
    ax1.set_title('Computation Time Comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels(molecules, rotation=45 if len(molecules) > 4 else 0)
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    # Speedup chart
    colors = ['green' if s > 1 else 'red' for s in speedups]
    bars3 = ax2.bar(x, speedups, color=colors, alpha=0.7, edgecolor='black')

    # Add value labels
    for bar, speedup in zip(bars3, speedups):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{speedup:.2f}×', ha='center', va='bottom', fontsize=9)

    ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='Equal time')
    ax2.set_xlabel('Molecule')
    ax2.set_ylabel('Speedup (AC time / CD time)')
    ax2.set_title('Relative Performance')
    ax2.set_xticks(x)
    ax2.set_xticklabels(molecules, rotation=45 if len(molecules) > 4 else 0)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    plt.suptitle('Performance Comparison: CD vs AC', fontsize=16)
    plt.tight_layout()
    return fig


def plot_convergence(results: List[ComparisonResult],
                    figsize: Tuple[float, float] = (12, 6)) -> Optional['plt.Figure']:
    """
    Plot convergence statistics for CD and AC methods.

    Parameters
    ----------
    results : List[ComparisonResult]
        List of comparison results
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    matplotlib.figure.Figure or None
        The convergence plot
    """
    if not HAS_MATPLOTLIB:
        return None

    set_publication_style()

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    molecules = [r.molecule_name for r in results]
    cd_iters = [r.cd_iterations for r in results]
    ac_iters = [r.ac_iterations for r in results]
    cd_converged = [r.cd_converged for r in results]
    ac_converged = [r.ac_converged for r in results]

    x = np.arange(len(molecules))
    width = 0.35

    # Iteration count comparison
    ax = axes[0]
    bars1 = ax.bar(x - width/2, cd_iters, width, label='CD',
                  color=['#2ecc71' if c else '#e74c3c' for c in cd_converged], alpha=0.8)
    bars2 = ax.bar(x + width/2, ac_iters, width, label='AC',
                  color=['#2ecc71' if c else '#e74c3c' for c in ac_converged], alpha=0.8)

    ax.set_xlabel('Molecule')
    ax.set_ylabel('Iterations')
    ax.set_title('Iteration Count')
    ax.set_xticks(x)
    ax.set_xticklabels(molecules, rotation=45 if len(molecules) > 4 else 0)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Memory usage comparison
    ax = axes[1]
    cd_memory = [r.cd_memory_peak for r in results]
    ac_memory = [r.ac_memory_peak for r in results]

    bars1 = ax.bar(x - width/2, cd_memory, width, label='CD', color='#3498db', alpha=0.8)
    bars2 = ax.bar(x + width/2, ac_memory, width, label='AC', color='#e74c3c', alpha=0.8)

    ax.set_xlabel('Molecule')
    ax.set_ylabel('Peak Memory (MB)')
    ax.set_title('Memory Usage')
    ax.set_xticks(x)
    ax.set_xticklabels(molecules, rotation=45 if len(molecules) > 4 else 0)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Convergence rate
    ax = axes[2]
    cd_rate = sum(cd_converged) / len(cd_converged) * 100
    ac_rate = sum(ac_converged) / len(ac_converged) * 100

    bars = ax.bar(['CD', 'AC'], [cd_rate, ac_rate],
                 color=['#2ecc71' if cd_rate > 90 else '#e74c3c',
                       '#2ecc71' if ac_rate > 90 else '#e74c3c'],
                 alpha=0.8)

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{height:.1f}%', ha='center', va='bottom')

    ax.set_ylabel('Convergence Rate (%)')
    ax.set_title('Overall Convergence')
    ax.set_ylim([0, 105])
    ax.grid(True, alpha=0.3, axis='y')

    plt.suptitle('Convergence Analysis', fontsize=16)
    plt.tight_layout()
    return fig


def plot_z_factors(results: List[ComparisonResult],
                  figsize: Tuple[float, float] = (12, 8)) -> Optional['plt.Figure']:
    """
    Plot Z-factor comparison between CD and AC methods.

    Parameters
    ----------
    results : List[ComparisonResult]
        List of comparison results
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    matplotlib.figure.Figure or None
        The Z-factor comparison plot
    """
    if not HAS_MATPLOTLIB:
        return None

    set_publication_style()

    n_results = len(results)
    n_cols = min(3, n_results)
    n_rows = (n_results + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)

    for idx, result in enumerate(results):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        # Plot Z factors
        n_orbitals = len(result.cd_z_factors)
        orbital_indices = np.arange(n_orbitals)

        ax.scatter(orbital_indices, result.cd_z_factors,
                  alpha=0.6, s=30, color='blue', label='CD', marker='o')
        ax.scatter(orbital_indices, result.ac_z_factors,
                  alpha=0.6, s=30, color='red', label='AC', marker='^')

        # Add physical bounds
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        ax.axhline(y=1, color='black', linestyle='-', alpha=0.3)
        ax.fill_between(orbital_indices, 0, 1, alpha=0.1, color='green',
                       label='Physical range')

        # Mark unphysical values
        unphysical_cd = (result.cd_z_factors < 0) | (result.cd_z_factors > 1)
        unphysical_ac = (result.ac_z_factors < 0) | (result.ac_z_factors > 1)

        if np.any(unphysical_cd):
            ax.scatter(orbital_indices[unphysical_cd],
                      result.cd_z_factors[unphysical_cd],
                      s=100, facecolors='none', edgecolors='red',
                      linewidth=2, marker='o')

        if np.any(unphysical_ac):
            ax.scatter(orbital_indices[unphysical_ac],
                      result.ac_z_factors[unphysical_ac],
                      s=100, facecolors='none', edgecolors='red',
                      linewidth=2, marker='^')

        ax.set_xlabel('Orbital Index')
        ax.set_ylabel('Z Factor')
        ax.set_title(f'{result.molecule_name}/{result.basis_set}')
        ax.set_ylim([-0.1, 1.1])
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)

    # Remove empty subplots
    for idx in range(n_results, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        fig.delaxes(axes[row, col])

    plt.suptitle('Z-Factor Comparison', fontsize=16)
    plt.tight_layout()
    return fig


def plot_orbital_comparison(result: ComparisonResult,
                           orbital_range: Optional[Tuple[int, int]] = None,
                           figsize: Tuple[float, float] = (12, 6)) -> Optional['plt.Figure']:
    """
    Create detailed orbital-by-orbital comparison for a single molecule.

    Parameters
    ----------
    result : ComparisonResult
        Single comparison result
    orbital_range : tuple, optional
        Range of orbitals to plot (start, end)
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    matplotlib.figure.Figure or None
        The orbital comparison plot
    """
    if not HAS_MATPLOTLIB:
        return None

    set_publication_style()

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Determine orbital range
    n_orbitals = len(result.cd_qp_energies)
    if orbital_range:
        start, end = orbital_range
        start = max(0, start)
        end = min(n_orbitals, end)
    else:
        start, end = 0, n_orbitals

    orbital_indices = np.arange(start, end)
    cd_qp = result.cd_qp_energies[start:end] * 27.211  # Convert to eV
    ac_qp = result.ac_qp_energies[start:end] * 27.211
    hf_e = result.hf_energies[start:end] * 27.211

    # QP energies comparison
    ax = axes[0, 0]
    ax.plot(orbital_indices, hf_e, 'k-', label='HF', alpha=0.5, linewidth=1)
    ax.plot(orbital_indices, cd_qp, 'b-', label='CD', linewidth=2, marker='o', markersize=4)
    ax.plot(orbital_indices, ac_qp, 'r--', label='AC', linewidth=2, marker='^', markersize=4)

    ax.set_xlabel('Orbital Index')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('Quasiparticle Energies')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Energy differences
    ax = axes[0, 1]
    cd_diff = cd_qp - hf_e
    ac_diff = ac_qp - hf_e

    ax.bar(orbital_indices - 0.2, cd_diff, 0.4, label='CD - HF', color='blue', alpha=0.7)
    ax.bar(orbital_indices + 0.2, ac_diff, 0.4, label='AC - HF', color='red', alpha=0.7)

    ax.set_xlabel('Orbital Index')
    ax.set_ylabel('QP Correction (eV)')
    ax.set_title('Self-Energy Corrections')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Z factors
    ax = axes[1, 0]
    cd_z = result.cd_z_factors[start:end]
    ac_z = result.ac_z_factors[start:end]

    ax.plot(orbital_indices, cd_z, 'b-', label='CD', linewidth=2, marker='o', markersize=4)
    ax.plot(orbital_indices, ac_z, 'r--', label='AC', linewidth=2, marker='^', markersize=4)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    ax.axhline(y=1, color='black', linestyle='-', alpha=0.3)
    ax.fill_between(orbital_indices, 0, 1, alpha=0.1, color='green')

    ax.set_xlabel('Orbital Index')
    ax.set_ylabel('Z Factor')
    ax.set_title('Quasiparticle Weights')
    ax.set_ylim([-0.1, 1.1])
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Absolute differences
    ax = axes[1, 1]
    energy_diff = np.abs(cd_qp - ac_qp)
    z_diff = np.abs(cd_z - ac_z)

    ax2 = ax.twinx()
    line1 = ax.plot(orbital_indices, energy_diff, 'g-', linewidth=2,
                   marker='s', markersize=4, label='|ΔE| (eV)')
    line2 = ax2.plot(orbital_indices, z_diff, 'm-', linewidth=2,
                    marker='d', markersize=4, label='|ΔZ|')

    ax.set_xlabel('Orbital Index')
    ax.set_ylabel('|ΔE| (eV)', color='g')
    ax2.set_ylabel('|ΔZ|', color='m')
    ax.set_title('Absolute Differences')
    ax.tick_params(axis='y', labelcolor='g')
    ax2.tick_params(axis='y', labelcolor='m')

    # Combine legends
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc='best')
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'Orbital Analysis: {result.molecule_name}/{result.basis_set}',
                fontsize=16)
    plt.tight_layout()
    return fig


def create_summary_figure(results: List[ComparisonResult],
                         figsize: Tuple[float, float] = (16, 10)) -> Optional['plt.Figure']:
    """
    Create a comprehensive summary figure with multiple panels.

    Parameters
    ----------
    results : List[ComparisonResult]
        List of comparison results
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    matplotlib.figure.Figure or None
        The comprehensive summary figure
    """
    if not HAS_MATPLOTLIB:
        return None

    set_publication_style()

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

    # Panel 1: Correlation plot
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    all_cd = []
    all_ac = []
    colors = plt.cm.Set3(np.linspace(0, 1, len(results)))

    for i, result in enumerate(results):
        cd_e = result.cd_qp_energies * 27.211
        ac_e = result.ac_qp_energies * 27.211
        ax1.scatter(cd_e, ac_e, alpha=0.5, s=30, color=colors[i],
                   label=result.molecule_name)
        all_cd.extend(cd_e)
        all_ac.extend(ac_e)

    all_cd = np.array(all_cd)
    all_ac = np.array(all_ac)
    min_e = min(all_cd.min(), all_ac.min())
    max_e = max(all_cd.max(), all_ac.max())
    ax1.plot([min_e, max_e], [min_e, max_e], 'k--', alpha=0.5, linewidth=2)

    ax1.set_xlabel('CD Energy (eV)')
    ax1.set_ylabel('AC Energy (eV)')
    ax1.set_title('Energy Correlation')
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Panel 2: MAD comparison
    ax2 = fig.add_subplot(gs[0, 2])
    mads = [r.mad for r in results]
    molecules = [r.molecule_name for r in results]
    colors_bar = ['green' if m < 0.05 else 'orange' if m < 0.1 else 'red' for m in mads]

    bars = ax2.barh(molecules, mads, color=colors_bar, alpha=0.7)
    ax2.axvline(x=0.05, color='green', linestyle='--', alpha=0.5)
    ax2.set_xlabel('MAD (eV)')
    ax2.set_title('Mean Absolute Deviation')
    ax2.grid(True, alpha=0.3, axis='x')

    # Panel 3: Timing comparison
    ax3 = fig.add_subplot(gs[1, 2])
    speedups = [r.ac_timing / r.cd_timing if r.cd_timing > 0 else 0 for r in results]
    colors_speed = ['green' if s > 1 else 'red' for s in speedups]

    bars = ax3.barh(molecules, speedups, color=colors_speed, alpha=0.7)
    ax3.axvline(x=1, color='black', linestyle='--', alpha=0.5)
    ax3.set_xlabel('Speedup (AC/CD)')
    ax3.set_title('Performance')
    ax3.grid(True, alpha=0.3, axis='x')

    # Panel 4: Error histogram
    ax4 = fig.add_subplot(gs[2, :])
    all_errors = []
    for result in results:
        errors = (result.cd_qp_energies - result.ac_qp_energies) * 27.211
        all_errors.extend(errors)

    ax4.hist(all_errors, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    ax4.axvline(x=0, color='green', linestyle='--', alpha=0.7, linewidth=2)
    ax4.axvline(x=np.mean(all_errors), color='red', linestyle='--', alpha=0.7, linewidth=2)
    ax4.set_xlabel('Energy Difference (CD - AC) [eV]')
    ax4.set_ylabel('Count')
    ax4.set_title('Error Distribution')
    ax4.grid(True, alpha=0.3)

    # Add overall statistics text
    overall_mad = np.mean(mads)
    overall_rmsd = np.sqrt(np.mean([r.rmsd**2 for r in results]))
    avg_speedup = np.mean(speedups)

    stats_text = (f"Overall Statistics:\n"
                 f"MAD: {overall_mad:.3f} eV\n"
                 f"RMSD: {overall_rmsd:.3f} eV\n"
                 f"Avg Speedup: {avg_speedup:.2f}×")

    fig.text(0.98, 0.02, stats_text, transform=fig.transFigure,
            fontsize=10, verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('CD vs AC Comparison Summary', fontsize=18, y=0.98)
    plt.tight_layout()
    return fig