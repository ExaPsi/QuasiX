"""
Visualization utilities for QuasiX GW/BSE calculations.

Provides publication-quality plots for convergence monitoring, spectral functions,
and electronic structure analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Optional, Dict, Any, List, Tuple, Union
import logging

logger = logging.getLogger(__name__)

# Set matplotlib style for publication quality
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.figsize': (8, 6),
    'figure.dpi': 150,
    'lines.linewidth': 2,
    'lines.markersize': 8,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


def plot_evgw_convergence(result: Any,
                          figsize: Tuple[float, float] = (12, 8),
                          save_path: Optional[str] = None) -> Figure:
    """
    Plot evGW convergence history with multiple panels.
    
    Args:
        result: EvGWResult object with convergence_history
        figsize: Figure size
        save_path: Optional path to save figure
        
    Returns:
        matplotlib Figure object
    """
    if result.convergence_history is None:
        raise ValueError("No convergence history available in result")
    
    history = result.convergence_history
    cycles = [h['cycle'] for h in history]
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle('evGW Convergence Analysis', fontsize=16, fontweight='bold')
    
    # Panel 1: Energy convergence
    ax = axes[0, 0]
    energy_diffs = [h['energy_diff'] for h in history]
    ax.semilogy(cycles, energy_diffs, 'b-o', label='Max ΔE')
    ax.axhline(y=1e-4, color='r', linestyle='--', alpha=0.5, label='Conv. threshold')
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Max Energy Difference (Ha)')
    ax.set_title('Energy Convergence')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Z-factor evolution
    ax = axes[0, 1]
    z_mins = [h['z_min'] for h in history]
    z_maxs = [h['z_max'] for h in history]
    z_means = [h['z_mean'] for h in history]
    
    ax.plot(cycles, z_means, 'g-s', label='Mean Z')
    ax.fill_between(cycles, z_mins, z_maxs, alpha=0.3, color='green', label='Min-Max range')
    ax.axhline(y=0.0, color='r', linestyle='--', alpha=0.5)
    ax.axhline(y=1.0, color='r', linestyle='--', alpha=0.5)
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Z-factor')
    ax.set_title('Quasiparticle Weight Evolution')
    ax.set_ylim(-0.1, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 3: Converged orbitals
    ax = axes[1, 0]
    converged_counts = [h['converged_orbitals'] for h in history]
    n_mo = len(result.qp_energies)
    ax.bar(cycles, converged_counts, color='blue', alpha=0.7)
    ax.axhline(y=n_mo, color='g', linestyle='--', alpha=0.5, label=f'Total ({n_mo})')
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Number of Orbitals')
    ax.set_title('Converged Orbitals per Cycle')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Panel 4: QP corrections distribution (final)
    ax = axes[1, 1]
    qp_corrections = (result.qp_energies - result.sigma_x.diagonal() + result.sigma_c_diag) * 27.2114  # to eV
    ax.hist(qp_corrections, bins=20, color='purple', alpha=0.7, edgecolor='black')
    ax.set_xlabel('QP Correction (eV)')
    ax.set_ylabel('Count')
    ax.set_title('Final QP Corrections Distribution')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path)
        logger.info(f"Convergence plot saved to {save_path}")
    
    return fig


def plot_spectral_function(qp_energies: np.ndarray,
                          z_factors: np.ndarray,
                          energy_range: Optional[Tuple[float, float]] = None,
                          broadening: float = 0.1,
                          n_points: int = 1000,
                          figsize: Tuple[float, float] = (10, 6),
                          save_path: Optional[str] = None) -> Figure:
    """
    Plot spectral function A(ω) with QP peaks and satellites.
    
    Args:
        qp_energies: Quasiparticle energies in Ha
        z_factors: Z-factors (quasiparticle weights)
        energy_range: Energy range in eV (auto if None)
        broadening: Gaussian broadening in eV
        n_points: Number of points for energy grid
        figsize: Figure size
        save_path: Optional path to save figure
        
    Returns:
        matplotlib Figure object
    """
    ha2ev = 27.2114
    qp_ev = qp_energies * ha2ev
    
    # Determine energy range
    if energy_range is None:
        e_min = np.min(qp_ev) - 2.0
        e_max = np.max(qp_ev) + 2.0
    else:
        e_min, e_max = energy_range
    
    # Create energy grid
    omega = np.linspace(e_min, e_max, n_points)
    
    # Build spectral function
    spectral = np.zeros_like(omega)
    
    for i, (e_qp, z) in enumerate(zip(qp_ev, z_factors)):
        # QP peak
        spectral += z * np.exp(-(omega - e_qp)**2 / (2 * broadening**2))
        
        # Satellite (simplified model)
        if z < 0.9:  # Significant satellite
            # Place satellite at higher binding energy for occupied states
            if i < len(qp_energies) // 2:  # Occupied
                e_sat = e_qp - 2.0  # Satellite 2 eV below
            else:  # Virtual
                e_sat = e_qp + 2.0  # Satellite 2 eV above
            spectral += (1 - z) * 0.5 * np.exp(-(omega - e_sat)**2 / (2 * (2*broadening)**2))
    
    # Normalize
    spectral /= np.max(spectral)
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot spectral function
    ax.plot(omega, spectral, 'b-', linewidth=2, label='A(ω)')
    ax.fill_between(omega, 0, spectral, alpha=0.3)
    
    # Mark QP positions
    for i, (e_qp, z) in enumerate(zip(qp_ev, z_factors)):
        if i == len(qp_energies) // 2 - 1:  # HOMO
            ax.axvline(x=e_qp, color='r', linestyle='--', alpha=0.7, label='HOMO')
        elif i == len(qp_energies) // 2:  # LUMO
            ax.axvline(x=e_qp, color='g', linestyle='--', alpha=0.7, label='LUMO')
    
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Spectral Function A(ω) (arb. units)')
    ax.set_title('Quasiparticle Spectral Function')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if save_path:
        fig.savefig(save_path)
        logger.info(f"Spectral function plot saved to {save_path}")
    
    return fig


def plot_qp_corrections(mo_energies: np.ndarray,
                        qp_energies: np.ndarray,
                        mo_occ: np.ndarray,
                        orbital_range: Optional[Tuple[int, int]] = None,
                        figsize: Tuple[float, float] = (10, 6),
                        save_path: Optional[str] = None) -> Figure:
    """
    Plot QP corrections vs orbital energy.
    
    Args:
        mo_energies: DFT/HF orbital energies in Ha
        qp_energies: Quasiparticle energies in Ha
        mo_occ: Orbital occupations
        orbital_range: Range of orbitals to plot (auto if None)
        figsize: Figure size
        save_path: Optional path to save figure
        
    Returns:
        matplotlib Figure object
    """
    ha2ev = 27.2114
    mo_ev = mo_energies * ha2ev
    qp_ev = qp_energies * ha2ev
    corrections = qp_ev - mo_ev
    
    # Determine orbital range
    if orbital_range is None:
        nocc = int(np.sum(mo_occ > 0))
        # Show HOMO-5 to LUMO+5
        idx_start = max(0, nocc - 6)
        idx_end = min(len(mo_energies), nocc + 6)
    else:
        idx_start, idx_end = orbital_range
    
    indices = np.arange(idx_start, idx_end)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Panel 1: Energy levels
    ax1.plot(indices, mo_ev[indices], 'b-o', label='DFT/HF', markersize=8)
    ax1.plot(indices, qp_ev[indices], 'r-s', label='GW', markersize=8)
    
    # Mark HOMO/LUMO
    nocc = int(np.sum(mo_occ > 0))
    if nocc - 1 in indices:
        ax1.axvline(x=nocc - 1, color='gray', linestyle='--', alpha=0.5)
        ax1.text(nocc - 1, ax1.get_ylim()[0], 'HOMO', ha='center', va='bottom')
    if nocc in indices:
        ax1.axvline(x=nocc, color='gray', linestyle='--', alpha=0.5)
        ax1.text(nocc, ax1.get_ylim()[0], 'LUMO', ha='center', va='bottom')
    
    ax1.set_xlabel('Orbital Index')
    ax1.set_ylabel('Energy (eV)')
    ax1.set_title('Orbital Energies')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: QP corrections
    colors = ['blue' if occ > 0 else 'red' for occ in mo_occ[indices]]
    ax2.bar(indices, corrections[indices], color=colors, alpha=0.7)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Add legend for colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', alpha=0.7, label='Occupied'),
        Patch(facecolor='red', alpha=0.7, label='Virtual')
    ]
    ax2.legend(handles=legend_elements)
    
    ax2.set_xlabel('Orbital Index')
    ax2.set_ylabel('QP Correction (eV)')
    ax2.set_title('Quasiparticle Corrections')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path)
        logger.info(f"QP corrections plot saved to {save_path}")
    
    return fig


def plot_z_factors(z_factors: np.ndarray,
                  mo_energies: np.ndarray,
                  mo_occ: np.ndarray,
                  figsize: Tuple[float, float] = (10, 6),
                  save_path: Optional[str] = None) -> Figure:
    """
    Plot Z-factors (quasiparticle weights) analysis.
    
    Args:
        z_factors: Z-factors for each orbital
        mo_energies: Orbital energies in Ha
        mo_occ: Orbital occupations
        figsize: Figure size
        save_path: Optional path to save figure
        
    Returns:
        matplotlib Figure object
    """
    ha2ev = 27.2114
    mo_ev = mo_energies * ha2ev
    nocc = int(np.sum(mo_occ > 0))
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Panel 1: Z-factors vs orbital index
    indices = np.arange(len(z_factors))
    colors = ['blue' if i < nocc else 'red' for i in indices]
    
    ax1.scatter(indices, z_factors, c=colors, s=50, alpha=0.7)
    ax1.axhline(y=0.0, color='gray', linestyle='--', alpha=0.5)
    ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax1.axhline(y=0.8, color='green', linestyle=':', alpha=0.5, label='Z=0.8')
    
    # Mark HOMO/LUMO
    ax1.axvline(x=nocc - 0.5, color='black', linestyle='-', alpha=0.3)
    ax1.text(nocc - 0.5, 1.05, 'HOMO|LUMO', ha='center', va='bottom', fontsize=10)
    
    ax1.set_xlabel('Orbital Index')
    ax1.set_ylabel('Z-factor')
    ax1.set_title('Quasiparticle Weights')
    ax1.set_ylim(-0.05, 1.1)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Z-factors vs orbital energy
    ax2.scatter(mo_ev[:nocc], z_factors[:nocc], c='blue', s=50, alpha=0.7, label='Occupied')
    ax2.scatter(mo_ev[nocc:], z_factors[nocc:], c='red', s=50, alpha=0.7, label='Virtual')
    
    # Fit trend lines
    if nocc > 2:
        z_occ = np.polyfit(mo_ev[:nocc], z_factors[:nocc], 1)
        p_occ = np.poly1d(z_occ)
        x_occ = np.linspace(mo_ev[:nocc].min(), mo_ev[:nocc].max(), 100)
        ax2.plot(x_occ, p_occ(x_occ), 'b--', alpha=0.5)
    
    if len(mo_ev) - nocc > 2:
        z_vir = np.polyfit(mo_ev[nocc:], z_factors[nocc:], 1)
        p_vir = np.poly1d(z_vir)
        x_vir = np.linspace(mo_ev[nocc:].min(), mo_ev[nocc:].max(), 100)
        ax2.plot(x_vir, p_vir(x_vir), 'r--', alpha=0.5)
    
    ax2.set_xlabel('Orbital Energy (eV)')
    ax2.set_ylabel('Z-factor')
    ax2.set_title('Z-factor vs Orbital Energy')
    ax2.set_ylim(-0.05, 1.1)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path)
        logger.info(f"Z-factors plot saved to {save_path}")
    
    return fig


def create_evgw_report(result: Any,
                       mo_energies: np.ndarray,
                       mo_occ: np.ndarray,
                       molecule_name: str = "Molecule",
                       save_dir: Optional[str] = None) -> Dict[str, Figure]:
    """
    Create comprehensive evGW analysis report with multiple figures.
    
    Args:
        result: EvGWResult object
        mo_energies: Original MO energies
        mo_occ: MO occupations
        molecule_name: Name of the molecule
        save_dir: Directory to save figures
        
    Returns:
        Dictionary of figure names to Figure objects
    """
    figures = {}
    
    # Convergence plot
    if result.convergence_history:
        fig = plot_evgw_convergence(
            result,
            save_path=f"{save_dir}/convergence.png" if save_dir else None
        )
        figures['convergence'] = fig
    
    # Spectral function
    fig = plot_spectral_function(
        result.qp_energies,
        result.z_factors,
        save_path=f"{save_dir}/spectral.png" if save_dir else None
    )
    figures['spectral'] = fig
    
    # QP corrections
    fig = plot_qp_corrections(
        mo_energies,
        result.qp_energies,
        mo_occ,
        save_path=f"{save_dir}/qp_corrections.png" if save_dir else None
    )
    figures['qp_corrections'] = fig
    
    # Z-factors
    fig = plot_z_factors(
        result.z_factors,
        mo_energies,
        mo_occ,
        save_path=f"{save_dir}/z_factors.png" if save_dir else None
    )
    figures['z_factors'] = fig
    
    # Create summary figure
    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(f'evGW Analysis Report: {molecule_name}', fontsize=18, fontweight='bold')
    
    # Add text summary
    ax = fig.add_subplot(111)
    ax.axis('off')
    
    summary_text = f"""
    Calculation Summary:
    {'='*40}
    Converged: {'Yes' if result.converged else 'No'}
    Cycles: {result.n_cycles}
    Final ΔE: {result.final_energy_diff:.2e} Ha
    Wall time: {result.wall_time:.1f} s
    
    Quasiparticle Properties:
    {'='*40}
    IP: {result.ip_ea()[0]:.2f} eV
    EA: {result.ip_ea()[1]:.2f} eV
    Gap: {result.gap():.2f} eV
    
    Z-factors:
    {'='*40}
    Min: {np.min(result.z_factors):.3f}
    Max: {np.max(result.z_factors):.3f}
    Mean: {np.mean(result.z_factors):.3f}
    Std: {np.std(result.z_factors):.3f}
    """
    
    ax.text(0.5, 0.5, summary_text, transform=ax.transAxes,
            fontsize=14, verticalalignment='center',
            horizontalalignment='center', family='monospace')
    
    if save_dir:
        fig.savefig(f"{save_dir}/summary.png")
    figures['summary'] = fig
    
    logger.info(f"Created evGW report with {len(figures)} figures")
    
    return figures