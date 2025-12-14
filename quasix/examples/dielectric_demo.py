#!/usr/bin/env python3
"""Demonstration of S3-2 dielectric and polarizability functionality.

This example shows how to use the QuasiX Python interface to compute:
1. Independent-particle polarizability P0(ω)
2. Dielectric function ε(ω)
3. Screened Coulomb interaction W
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add parent directory for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from quasix import dielectric


def create_realistic_data(nocc=5, nvirt=10, naux=30):
    """Create realistic mock data for demonstration.
    
    In real usage, this data would come from PySCF calculations.
    """
    # Realistic orbital energies (in Hartree)
    e_occ = np.array([-20.0, -15.0, -10.0, -1.3, -0.7][:nocc])
    e_virt = np.linspace(0.1, 10.0, nvirt)
    
    # Create DF tensors with appropriate magnitude
    # In reality, these come from (ia|P) = ∫ φ_i(r) φ_a(r) χ_P(r) dr
    n_trans = nocc * nvirt
    np.random.seed(42)
    
    # Scale factors to ensure physical results
    df_ia = np.random.randn(n_trans, naux) * 0.01  # Small values
    
    # Create positive definite metric (P|Q)
    A = np.random.randn(naux, naux) * 0.1
    metric = A @ A.T + np.eye(naux) * 2.0  # Ensure positive definite
    
    return df_ia, e_occ, e_virt, metric


def demo_polarizability():
    """Demonstrate polarizability calculation."""
    print("\n" + "="*70)
    print("DEMO 1: Frequency-Dependent Polarizability P0(ω)")
    print("="*70)
    
    # Setup
    nocc, nvirt, naux = 5, 10, 30
    df_ia, e_occ, e_virt, _ = create_realistic_data(nocc, nvirt, naux)
    
    # Create calculator
    calc = dielectric.PolarizabilityCalculator(nocc, nvirt, naux, eta=0.01)
    
    # Compute static polarizability
    p0_static = calc.compute_p0(df_ia, e_occ, e_virt, complex(0.0, 0.01))
    
    print(f"\nStatic Polarizability P0(ω=0):")
    print(f"  Shape: {p0_static.shape}")
    print(f"  Trace: {np.trace(p0_static):.6f}")
    print(f"  Max diagonal: {np.max(np.diag(p0_static)):.6f}")
    print(f"  Min diagonal: {np.min(np.diag(p0_static)):.6f}")
    
    # Frequency scan
    omega_range = np.linspace(-5, 15, 200)
    p0_freq = calc.compute_p0_batch(df_ia, e_occ, e_virt, omega_range, 0.01)
    
    # Plot frequency dependence
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot trace of P0
    trace_p0 = np.array([np.trace(p0_freq[i]) for i in range(len(omega_range))])
    ax1.plot(omega_range, trace_p0.real, 'b-', label='Re[Tr(P0)]')
    ax1.plot(omega_range, trace_p0.imag, 'r--', label='Im[Tr(P0)]')
    ax1.set_xlabel('Frequency ω (a.u.)')
    ax1.set_ylabel('Tr[P0(ω)]')
    ax1.set_title('Polarizability Trace vs Frequency')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot selected diagonal elements
    for i in [0, 5, 10]:
        ax2.plot(omega_range, p0_freq[:, i, i].real, 
                label=f'P0_{{{i},{i}}}')
    
    ax2.set_xlabel('Frequency ω (a.u.)')
    ax2.set_ylabel('Re[P0_ii(ω)]')
    ax2.set_title('Diagonal Elements of P0')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('polarizability_demo.png', dpi=150)
    print(f"\nPlot saved to: polarizability_demo.png")
    
    return p0_static


def demo_dielectric_screening(p0_static):
    """Demonstrate dielectric function and screening."""
    print("\n" + "="*70)
    print("DEMO 2: Dielectric Function and Screening")
    print("="*70)
    
    # Get metric
    _, _, _, metric = create_realistic_data(naux=p0_static.shape[0])
    
    # Create dielectric calculator
    calc = dielectric.DielectricFunction(metric)
    
    print(f"\nMetric properties:")
    print(f"  Condition number: {calc.v_condition:.2e}")
    print(f"  Max eigenvalue: {np.max(np.linalg.eigvalsh(metric)):.4f}")
    print(f"  Min eigenvalue: {np.min(np.linalg.eigvalsh(metric)):.4f}")
    
    # Scale P0 to ensure physical results
    # M = V^(1/2) P0 V^(1/2) must have eigenvalues < 1
    p0_scaled = p0_static * 0.001  # Scale down for stability
    
    # Compute dielectric function
    result = calc.compute_epsilon(p0_scaled)
    epsilon = result['epsilon']
    
    print(f"\nDielectric Matrix ε = 1 - V^(1/2) P0 V^(1/2):")
    print(f"  Shape: {epsilon.shape}")
    print(f"  Condition number: {result['condition_number']:.2e}")
    print(f"  Trace: {np.trace(epsilon):.6f}")
    print(f"  Max diagonal: {np.max(np.diag(epsilon)):.6f}")
    print(f"  Min diagonal: {np.min(np.diag(epsilon)):.6f}")
    
    # Compute inverse for screening
    epsilon_inv = calc.compute_epsilon_inverse(epsilon)
    
    # Verify inversion quality
    identity_check = epsilon @ epsilon_inv
    max_error = np.max(np.abs(identity_check - np.eye(epsilon.shape[0])))
    print(f"\nInverse Dielectric ε^(-1):")
    print(f"  Inversion error: {max_error:.2e}")
    print(f"  Trace[ε^(-1)]: {np.trace(epsilon_inv):.6f}")
    
    # Compute screened Coulomb W
    w = calc.compute_screened_coulomb(epsilon_inv)
    
    print(f"\nScreened Coulomb W = V^(1/2) ε^(-1) V^(1/2):")
    print(f"  Norm: {np.linalg.norm(w):.6f}")
    print(f"  Trace: {np.trace(w):.6f}")
    
    # Visualization
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot matrices
    matrices = [
        (metric, 'Coulomb Metric V', axes[0,0]),
        (p0_scaled, 'Polarizability P0', axes[0,1]),
        (epsilon, 'Dielectric ε', axes[0,2]),
        (epsilon_inv, 'Inverse ε^(-1)', axes[1,0]),
        (w, 'Screened W', axes[1,1]),
        (w - metric, 'W - V', axes[1,2]),
    ]
    
    for mat, title, ax in matrices:
        im = ax.imshow(mat.real, cmap='RdBu_r', aspect='auto')
        ax.set_title(title)
        ax.set_xlabel('Auxiliary Index')
        ax.set_ylabel('Auxiliary Index')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    plt.savefig('dielectric_screening_demo.png', dpi=150)
    print(f"\nPlot saved to: dielectric_screening_demo.png")
    
    return epsilon, w


def demo_physical_interpretation():
    """Demonstrate physical interpretation of results."""
    print("\n" + "="*70)
    print("DEMO 3: Physical Interpretation")
    print("="*70)
    
    # Create data for a "molecule" with clear HOMO-LUMO gap
    nocc = 5
    nvirt = 10
    naux = 20
    
    # Clear HOMO-LUMO gap
    e_homo = -0.3  # HOMO energy
    e_lumo = 0.1   # LUMO energy
    gap = e_lumo - e_homo
    
    e_occ = np.array([-10.0, -5.0, -2.0, -0.5, e_homo])
    e_virt = np.array([e_lumo, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    
    print(f"\nOrbital Energies:")
    print(f"  HOMO: {e_homo:.3f} a.u.")
    print(f"  LUMO: {e_lumo:.3f} a.u.")
    print(f"  HOMO-LUMO gap: {gap:.3f} a.u. = {gap * 27.2114:.2f} eV")
    
    # Create DF tensors
    n_trans = nocc * nvirt
    np.random.seed(42)
    df_ia = np.random.randn(n_trans, naux) * 0.01
    
    # Frequency grid around the gap
    omega_range = np.linspace(-1, 2, 300)
    
    # Create calculator
    calc = dielectric.PolarizabilityCalculator(nocc, nvirt, naux, eta=0.005)
    
    # Compute P0
    p0_freq = calc.compute_p0_batch(df_ia, e_occ, e_virt, omega_range, 0.005)
    
    # Analyze spectral features
    trace_p0 = np.array([np.trace(p0_freq[i]) for i in range(len(omega_range))])
    
    # Find peaks in imaginary part (optical absorption)
    peaks_idx = np.where(np.diff(np.sign(np.diff(trace_p0.imag))) == -2)[0] + 1
    
    print(f"\nSpectral Analysis:")
    print(f"  Number of peaks in Im[P0]: {len(peaks_idx)}")
    if len(peaks_idx) > 0:
        print(f"  First peak at ω = {omega_range[peaks_idx[0]]:.3f} a.u.")
        print(f"  Expected at HOMO-LUMO gap: {gap:.3f} a.u.")
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Top: Real and imaginary parts
    ax1.plot(omega_range, trace_p0.real, 'b-', label='Re[Tr(P0)]', linewidth=2)
    ax1.plot(omega_range, trace_p0.imag, 'r-', label='Im[Tr(P0)]', linewidth=2)
    ax1.axvline(gap, color='g', linestyle='--', alpha=0.5, label='HOMO-LUMO gap')
    ax1.set_xlabel('Frequency ω (a.u.)')
    ax1.set_ylabel('Tr[P0(ω)]')
    ax1.set_title('Polarizability Response Function')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Bottom: Optical absorption (imaginary part)
    ax2.fill_between(omega_range, 0, trace_p0.imag, 
                     where=(trace_p0.imag > 0), alpha=0.3, color='red')
    ax2.plot(omega_range, trace_p0.imag, 'r-', linewidth=2)
    ax2.axvline(gap, color='g', linestyle='--', alpha=0.5, label='HOMO-LUMO gap')
    ax2.set_xlabel('Frequency ω (a.u.)')
    ax2.set_ylabel('Im[Tr(P0(ω))] (arb. units)')
    ax2.set_title('Optical Absorption Spectrum')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, 1.5])
    ax2.set_ylim([0, None])
    
    plt.tight_layout()
    plt.savefig('physical_interpretation_demo.png', dpi=150)
    print(f"\nPlot saved to: physical_interpretation_demo.png")


def main():
    """Run all demonstrations."""
    print("\n" + "="*70)
    print("QuasiX S3-2: Dielectric and Polarizability Demonstration")
    print("="*70)
    
    print("\nThis demonstration shows the key functionality of the S3-2 module:")
    print("  1. Computing frequency-dependent polarizability P0(ω)")
    print("  2. Computing dielectric function ε(ω)")
    print("  3. Computing screened Coulomb interaction W")
    print("  4. Physical interpretation of results")
    
    # Run demonstrations
    p0_static = demo_polarizability()
    epsilon, w = demo_dielectric_screening(p0_static)
    demo_physical_interpretation()
    
    print("\n" + "="*70)
    print("Demonstration Complete!")
    print("="*70)
    print("\nGenerated plots:")
    print("  - polarizability_demo.png: Frequency-dependent polarizability")
    print("  - dielectric_screening_demo.png: Dielectric matrices")
    print("  - physical_interpretation_demo.png: Physical interpretation")
    
    print("\nKey takeaways:")
    print("  ✓ P0(ω) shows resonances at transition energies")
    print("  ✓ Dielectric function ε screens the bare Coulomb interaction")
    print("  ✓ Screened interaction W = V ε^(-1) is used in GW self-energy")
    print("  ✓ Python interface seamlessly wraps high-performance Rust code")


if __name__ == "__main__":
    main()