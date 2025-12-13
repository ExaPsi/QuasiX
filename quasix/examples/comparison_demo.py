#!/usr/bin/env python3
"""
Demonstration of CD vs AC comparison harness for QuasiX.

This script shows various usage patterns for comparing Contour Deformation
and Analytic Continuation methods in GW calculations.
"""

import sys
import os
import numpy as np
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import QuasiX comparison module
try:
    from quasix.comparison import (
        ComparisonConfig,
        CDvsACComparator,
        plot_correlation,
        plot_error_distribution,
        plot_timing_comparison,
        plot_z_factors,
        plot_orbital_comparison,
        create_summary_figure
    )
    COMPARISON_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import comparison module: {e}")
    COMPARISON_AVAILABLE = False

# Try to import PySCF
try:
    from pyscf import gto, scf
    HAS_PYSCF = True
except ImportError:
    HAS_PYSCF = False
    print("Warning: PySCF not found. Using mock data for demonstration.")


def example_basic_comparison():
    """Example 1: Basic comparison on a single molecule"""
    print("\n" + "="*60)
    print("Example 1: Basic Comparison on H2O")
    print("="*60)

    if not COMPARISON_AVAILABLE:
        print("Comparison module not available. Skipping.")
        return

    # Configure comparison
    config = ComparisonConfig(
        n_grid_points=40,
        xi_max=50.0,
        eta=0.01,
        mad_threshold=0.05,  # 50 meV threshold
        parallel_execution=True,
        verbose=1
    )

    # Initialize comparator
    comparator = CDvsACComparator(config)

    # Run comparison on H2O
    print("\nRunning CD vs AC comparison on H2O/cc-pVDZ...")
    result = comparator.compare_molecule('H2O', basis='cc-pVDZ')

    # Display results
    print("\n--- Results ---")
    print(f"Molecule: {result.molecule_name}/{result.basis_set}")
    print(f"MAD: {result.mad:.3f} eV")
    print(f"RMSD: {result.rmsd:.3f} eV")
    print(f"Max deviation: {result.max_deviation:.3f} eV")
    print(f"Correlation (R²): {result.r_squared:.4f}")
    print(f"95% CI for MAD: [{result.mad_ci_lower:.3f}, {result.mad_ci_upper:.3f}] eV")

    # Check if methods agree
    if result.mad < config.mad_threshold:
        print(f"\n✓ Methods agree within threshold ({config.mad_threshold} eV)")
    else:
        print(f"\n✗ Methods differ beyond threshold ({config.mad_threshold} eV)")

    # Performance comparison
    print(f"\n--- Performance ---")
    print(f"CD timing: {result.cd_timing:.2f} s")
    print(f"AC timing: {result.ac_timing:.2f} s")
    speedup = result.ac_timing / result.cd_timing if result.cd_timing > 0 else 0
    print(f"Speedup: {speedup:.2f}×")

    # Convergence information
    print(f"\n--- Convergence ---")
    print(f"CD: {'Converged' if result.cd_converged else 'Not converged'} "
          f"({result.cd_iterations} iterations)")
    print(f"AC: {'Converged' if result.ac_converged else 'Not converged'} "
          f"({result.ac_iterations} iterations)")

    # Outliers
    if result.outlier_indices:
        print(f"\n--- Outliers ---")
        print(f"Found {len(result.outlier_indices)} outlier orbitals:")
        for idx, orbital in zip(result.outlier_indices, result.outlier_orbitals):
            cd_e = result.cd_qp_energies[idx] * 27.211
            ac_e = result.ac_qp_energies[idx] * 27.211
            diff = cd_e - ac_e
            print(f"  {orbital}: CD={cd_e:.3f} eV, AC={ac_e:.3f} eV, Δ={diff:.3f} eV")

    return result


def example_benchmark_suite():
    """Example 2: Run benchmark suite on multiple molecules"""
    print("\n" + "="*60)
    print("Example 2: Benchmark Suite")
    print("="*60)

    if not COMPARISON_AVAILABLE:
        print("Comparison module not available. Skipping.")
        return

    # Configure for benchmark
    config = ComparisonConfig(
        n_grid_points=40,
        xi_max=50.0,
        mad_threshold=0.05,
        bootstrap_samples=1000,
        parallel_execution=True,
        plot_results=True,
        save_intermediates=True,
        output_dir="comparison_results",
        verbose=1
    )

    comparator = CDvsACComparator(config)

    # Define benchmark molecules
    molecules = ['H2O', 'NH3', 'CO']
    basis_sets = ['cc-pVDZ'] * len(molecules)

    print(f"\nRunning benchmark on {len(molecules)} molecules...")
    print(f"Molecules: {', '.join(molecules)}")
    print(f"Basis: {basis_sets[0]}")

    # Run comparisons
    report = comparator.compare_molecules(molecules, basis_sets)

    # Display summary
    print("\n--- Benchmark Summary ---")
    report.print_summary()

    # Generate reports
    print("\n--- Generating Reports ---")

    # HTML report
    html_path = report.generate_html_report()
    print(f"HTML report saved to: {html_path}")

    # JSON export
    json_path = report.to_json()
    print(f"JSON results saved to: {json_path}")

    return report


def example_custom_molecule():
    """Example 3: Custom molecule with PySCF integration"""
    print("\n" + "="*60)
    print("Example 3: Custom Molecule with PySCF")
    print("="*60)

    if not COMPARISON_AVAILABLE:
        print("Comparison module not available. Skipping.")
        return

    if not HAS_PYSCF:
        print("PySCF not available. Skipping.")
        return

    # Define custom molecule
    mol = gto.M(
        atom="""
        O  0.0000  0.0000  0.1173
        H  0.0000  0.7572 -0.4692
        H  0.0000 -0.7572 -0.4692
        """,
        basis='aug-cc-pVDZ',
        symmetry=True,
        verbose=0
    )

    print(f"Custom molecule: {mol.natm} atoms, {mol.nao} basis functions")

    # Run mean-field calculation
    print("Running RHF calculation...")
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-10
    mf.kernel()
    print(f"HF energy: {mf.e_tot:.6f} Ha")

    # Configure comparison with higher accuracy
    config = ComparisonConfig(
        n_grid_points=60,  # More grid points
        xi_max=60.0,       # Larger range
        convergence_tol=1e-7,
        mad_threshold=0.03,  # Stricter threshold
        bootstrap_samples=2000,
        verbose=1
    )

    # Run comparison from mean-field object
    comparator = CDvsACComparator(config)
    print("\nRunning CD vs AC comparison...")
    result = comparator.compare_from_mean_field(mf)

    # Detailed analysis
    print("\n--- Detailed Analysis ---")
    print(f"Number of orbitals: {len(result.cd_qp_energies)}")
    print(f"Occupied orbitals: {np.sum(mf.mo_occ > 0)}")

    # Statistics by orbital type
    nocc = np.sum(mf.mo_occ > 0)
    occ_mad = np.mean(np.abs(
        result.cd_qp_energies[:nocc] - result.ac_qp_energies[:nocc]
    )) * 27.211
    vir_mad = np.mean(np.abs(
        result.cd_qp_energies[nocc:] - result.ac_qp_energies[nocc:]
    )) * 27.211

    print(f"Occupied orbital MAD: {occ_mad:.3f} eV")
    print(f"Virtual orbital MAD: {vir_mad:.3f} eV")

    # Z-factor analysis
    print(f"\n--- Z-Factor Analysis ---")
    print(f"CD Z-factors: mean={np.mean(result.cd_z_factors):.3f}, "
          f"std={np.std(result.cd_z_factors):.3f}")
    print(f"AC Z-factors: mean={np.mean(result.ac_z_factors):.3f}, "
          f"std={np.std(result.ac_z_factors):.3f}")

    # Check for unphysical Z-factors
    unphysical_cd = np.sum((result.cd_z_factors < 0) | (result.cd_z_factors > 1))
    unphysical_ac = np.sum((result.ac_z_factors < 0) | (result.ac_z_factors > 1))

    if unphysical_cd > 0:
        print(f"Warning: {unphysical_cd} unphysical Z-factors in CD")
    if unphysical_ac > 0:
        print(f"Warning: {unphysical_ac} unphysical Z-factors in AC")

    return result, mf


def example_advanced_configuration():
    """Example 4: Advanced configuration options"""
    print("\n" + "="*60)
    print("Example 4: Advanced Configuration")
    print("="*60)

    if not COMPARISON_AVAILABLE:
        print("Comparison module not available. Skipping.")
        return

    # Custom CD and AC parameters
    config = ComparisonConfig(
        n_grid_points=50,
        xi_max=60.0,
        eta=0.005,  # Smaller imaginary shift

        # Method-specific parameters
        cd_params={
            'damping': 0.3,
            'max_iterations': 30
        },
        ac_params={
            'n_poles': 10,
            'pole_threshold': 1e-4
        },

        # Statistical parameters
        mad_threshold=0.04,
        outlier_sigma=2.5,  # Less strict outlier detection
        bootstrap_samples=5000,  # More bootstrap samples
        confidence_level=0.99,  # 99% confidence interval

        # Performance options
        parallel_execution=True,
        n_threads=4,
        cache_intermediates=True,
        memory_limit_gb=32.0,

        # Output options
        verbose=2,  # Detailed logging
        save_intermediates=True,
        plot_results=True,
        output_dir="advanced_comparison"
    )

    print("Configuration summary:")
    print(f"  Grid points: {config.n_grid_points}")
    print(f"  ξ_max: {config.xi_max}")
    print(f"  η: {config.eta}")
    print(f"  MAD threshold: {config.mad_threshold} eV")
    print(f"  Bootstrap samples: {config.bootstrap_samples}")
    print(f"  Confidence level: {config.confidence_level*100}%")
    print(f"  Memory limit: {config.memory_limit_gb} GB")

    # Initialize comparator with custom config
    comparator = CDvsACComparator(config)

    # Run on a challenging system
    print("\nTesting on CH4...")
    result = comparator.compare_molecule('CH4', basis='def2-TZVP')

    print(f"\nResult: MAD = {result.mad:.4f} eV "
          f"({result.mad_ci_lower:.4f}, {result.mad_ci_upper:.4f})")

    return result


def example_visualization():
    """Example 5: Generate publication-quality plots"""
    print("\n" + "="*60)
    print("Example 5: Visualization")
    print("="*60)

    if not COMPARISON_AVAILABLE:
        print("Comparison module not available. Skipping.")
        return

    # Run comparison on multiple molecules for plotting
    config = ComparisonConfig(
        n_grid_points=40,
        verbose=0  # Quiet mode for plotting
    )

    comparator = CDvsACComparator(config)
    molecules = ['H2O', 'NH3', 'HF']

    print(f"Running comparisons for visualization...")
    results = []
    for mol in molecules:
        print(f"  Processing {mol}...")
        result = comparator.compare_molecule(mol, basis='cc-pVDZ')
        results.append(result)

    # Generate various plots
    print("\nGenerating plots...")

    # 1. Correlation plot
    fig1 = plot_correlation(results)
    if fig1:
        fig1.savefig('correlation_plot.png', dpi=150, bbox_inches='tight')
        print("  Saved: correlation_plot.png")

    # 2. Error distribution
    fig2 = plot_error_distribution(results)
    if fig2:
        fig2.savefig('error_distribution.png', dpi=150, bbox_inches='tight')
        print("  Saved: error_distribution.png")

    # 3. Timing comparison
    fig3 = plot_timing_comparison(results)
    if fig3:
        fig3.savefig('timing_comparison.png', dpi=150, bbox_inches='tight')
        print("  Saved: timing_comparison.png")

    # 4. Z-factor comparison
    fig4 = plot_z_factors(results)
    if fig4:
        fig4.savefig('z_factors.png', dpi=150, bbox_inches='tight')
        print("  Saved: z_factors.png")

    # 5. Orbital analysis for first molecule
    fig5 = plot_orbital_comparison(results[0])
    if fig5:
        fig5.savefig('orbital_analysis.png', dpi=150, bbox_inches='tight')
        print("  Saved: orbital_analysis.png")

    # 6. Comprehensive summary
    fig6 = create_summary_figure(results)
    if fig6:
        fig6.savefig('summary_figure.png', dpi=300, bbox_inches='tight')
        print("  Saved: summary_figure.png (high-res)")

    return results


def main():
    """Run all examples"""
    print("\n" + "#"*60)
    print("# CD vs AC Comparison Harness Demonstration")
    print("#"*60)

    # Example 1: Basic comparison
    result1 = example_basic_comparison()

    # Example 2: Benchmark suite
    report = example_benchmark_suite()

    # Example 3: Custom molecule (if PySCF available)
    if HAS_PYSCF:
        result = example_custom_molecule()
        if result is not None:
            result3, mf = result

    # Example 4: Advanced configuration
    result4 = example_advanced_configuration()

    # Example 5: Visualization
    results = example_visualization()

    print("\n" + "#"*60)
    print("# Demonstration Complete")
    print("#"*60)
    print("\nAll examples have been executed successfully!")
    print("Check the generated files for detailed results:")
    print("  - comparison_results/: Intermediate results and reports")
    print("  - *.png: Visualization plots")
    print("  - *.html: Interactive HTML reports")
    print("  - *.json: Machine-readable results")


if __name__ == "__main__":
    main()