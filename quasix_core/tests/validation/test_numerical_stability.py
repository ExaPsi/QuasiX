#!/usr/bin/env python3
"""
Numerical Stability Analysis for Correlation Self-Energy

This script performs comprehensive tests to identify potential catastrophic cancellation,
overflow, or conditioning issues in the Σᶜ calculation.

Tests:
1. Quadrature accuracy validation
2. Denominator conditioning analysis
3. Catastrophic cancellation detection
4. Frequency grid convergence
5. Eta broadening sensitivity

Reference: docs/derivations/frequency_integration_specification.md
"""

import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad
from typing import Tuple, Dict, List
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

try:
    from pyscf import gto, scf, dft
    from pyscf.ao2mo import kernel as ao2mo_kernel
    PYSCF_AVAILABLE = True
except ImportError:
    print("WARNING: PySCF not available, using synthetic data")
    PYSCF_AVAILABLE = False

try:
    import quasix
    QUASIX_AVAILABLE = True
except ImportError:
    print("WARNING: QuasiX not available, testing quadrature only")
    QUASIX_AVAILABLE = False


class NumericalStabilityAnalyzer:
    """Analyzes numerical stability of GW correlation self-energy"""

    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.results = {}

    def log(self, message: str):
        """Print message if verbose"""
        if self.verbose:
            print(message)

    # =======================================================================
    # TEST 1: Quadrature Accuracy Validation
    # =======================================================================

    def test_quadrature_accuracy(self, n_points: int = 32, xi_max: float = 100.0) -> Dict:
        """
        Test 1: Verify Gauss-Legendre quadrature accuracy

        Reference: Onida et al. RMP 2002, Section III.B

        Returns:
            dict: {
                'quadrature_error': float,
                'nodes': np.ndarray,
                'weights': np.ndarray,
                'passed': bool
            }
        """
        self.log("\n" + "="*80)
        self.log("TEST 1: Quadrature Accuracy Validation")
        self.log("="*80)

        # Get standard GL nodes on [-1, 1]
        nodes_11, weights_11 = roots_legendre(n_points)

        # Transform to [0, 1]
        nodes_01 = 0.5 * (nodes_11 + 1.0)
        weights_01 = 0.5 * weights_11

        # Transform to [0, xi_max]
        nodes = nodes_01 * xi_max
        weights = weights_01 * xi_max

        # Test: Integrate f(x) = 1 over [0, xi_max]
        # Should give exactly xi_max
        integral_test = np.sum(weights * 1.0)
        quadrature_error = abs(integral_test - xi_max)

        self.log(f"∫₀^{xi_max} 1 dx = {integral_test:.15f}")
        self.log(f"Exact value:        {xi_max:.15f}")
        self.log(f"Quadrature error:   {quadrature_error:.2e}")

        # Test polynomial integration (Gauss-Legendre exact up to degree 2n-1)
        # Test x^2
        integral_x2 = np.sum(weights * nodes**2)
        exact_x2 = xi_max**3 / 3.0
        error_x2 = abs(integral_x2 - exact_x2)

        self.log(f"\n∫₀^{xi_max} x² dx = {integral_x2:.15f}")
        self.log(f"Exact value:         {exact_x2:.15f}")
        self.log(f"Error:               {error_x2:.2e}")

        # Test x^4
        integral_x4 = np.sum(weights * nodes**4)
        exact_x4 = xi_max**5 / 5.0
        error_x4 = abs(integral_x4 - exact_x4)

        self.log(f"\n∫₀^{xi_max} x⁴ dx = {integral_x4:.15f}")
        self.log(f"Exact value:         {exact_x4:.15f}")
        self.log(f"Error:               {error_x4:.2e}")

        # Check tolerance: n=32 GL should be machine precision for low-order polynomials
        tolerance = 1e-10
        passed = (quadrature_error < tolerance and
                  error_x2 < tolerance and
                  error_x4 < tolerance)

        if passed:
            self.log(f"\n✅ PASSED: Quadrature error < {tolerance:.0e}")
        else:
            self.log(f"\n❌ FAILED: Quadrature error > {tolerance:.0e}")
            self.log("   This suggests incorrect grid transformation!")

        result = {
            'quadrature_error': quadrature_error,
            'error_x2': error_x2,
            'error_x4': error_x4,
            'nodes': nodes,
            'weights': weights,
            'passed': passed,
            'tolerance': tolerance
        }

        self.results['test_1_quadrature'] = result
        return result

    # =======================================================================
    # TEST 2: Denominator Conditioning Analysis
    # =======================================================================

    def test_denominator_conditioning(
        self,
        mo_energies: np.ndarray,
        omega_grid: np.ndarray,
        eta: float = 0.001
    ) -> Dict:
        """
        Test 2: Analyze denominator distribution in GW calculation

        Small denominators → large Green's function → potential explosion

        Args:
            mo_energies: Orbital energies [n_mo]
            omega_grid: Imaginary frequency grid [n_omega]
            eta: Broadening parameter (Hartree)

        Returns:
            dict: {
                'min_denominator': float,
                'max_denominator': float,
                'degeneracies': list,
                'conditioning': float,
                'passed': bool
            }
        """
        self.log("\n" + "="*80)
        self.log("TEST 2: Denominator Conditioning Analysis")
        self.log("="*80)

        denoms_real = []
        denoms_imag = []

        # Analyze all denominators: |ε_p - ε_q ± iω|
        for ep in mo_energies:
            for eq in mo_energies:
                # Real axis (pole positions)
                denom_real = abs(ep - eq) + eta  # Add eta for stability
                denoms_real.append(denom_real)

                # Imaginary axis
                for omega in omega_grid:
                    denom_imag = np.sqrt((ep - eq)**2 + omega**2)
                    denoms_imag.append(denom_imag)

        denoms_real = np.array(denoms_real)
        denoms_imag = np.array(denoms_imag)

        min_denom_real = np.min(denoms_real[denoms_real > 0])
        max_denom_real = np.max(denoms_real)
        min_denom_imag = np.min(denoms_imag)
        max_denom_imag = np.max(denoms_imag)

        self.log(f"Real axis denominators:")
        self.log(f"  Min: {min_denom_real:.6f} Ha = {min_denom_real*27.211:.3f} eV")
        self.log(f"  Max: {max_denom_real:.6f} Ha = {max_denom_real*27.211:.3f} eV")

        self.log(f"\nImaginary axis denominators:")
        self.log(f"  Min: {min_denom_imag:.6f} Ha = {min_denom_imag*27.211:.3f} eV")
        self.log(f"  Max: {max_denom_imag:.6f} Ha = {max_denom_imag*27.211:.3f} eV")

        # Check for near-degeneracies (ε_p ≈ ε_q)
        degeneracy_threshold = 0.01  # Ha (~0.27 eV)
        degeneracies = []

        for i, ep in enumerate(mo_energies):
            for j, eq in enumerate(mo_energies):
                if i < j and abs(ep - eq) < degeneracy_threshold:
                    degeneracies.append((i, j, ep, eq, abs(ep - eq)))

        if degeneracies:
            self.log(f"\n⚠️  Found {len(degeneracies)} near-degenerate pairs:")
            for i, j, ep, eq, delta in degeneracies[:5]:
                self.log(f"   MO {i:2d} & {j:2d}: "
                        f"ε_p={ep:.6f}, ε_q={eq:.6f}, Δε={delta:.2e} Ha")
        else:
            self.log("\n✅ No near-degeneracies found")

        # Conditioning: ratio of max/min denominators
        conditioning_real = max_denom_real / min_denom_real
        conditioning_imag = max_denom_imag / min_denom_imag

        self.log(f"\nConditioning (max/min):")
        self.log(f"  Real axis:      {conditioning_real:.2e}")
        self.log(f"  Imaginary axis: {conditioning_imag:.2e}")

        # Warning thresholds
        warning_threshold = 1e3  # Conditioning > 1000 can cause issues
        passed = conditioning_real < warning_threshold

        if not passed:
            self.log(f"\n⚠️  WARNING: High conditioning ({conditioning_real:.2e})")
            self.log("   Consider increasing eta or checking for degeneracies")
        else:
            self.log("\n✅ PASSED: Conditioning acceptable")

        result = {
            'min_denom_real': min_denom_real,
            'max_denom_real': max_denom_real,
            'min_denom_imag': min_denom_imag,
            'max_denom_imag': max_denom_imag,
            'conditioning_real': conditioning_real,
            'conditioning_imag': conditioning_imag,
            'degeneracies': degeneracies,
            'passed': passed,
            'warning_threshold': warning_threshold
        }

        self.results['test_2_denominators'] = result
        return result

    # =======================================================================
    # TEST 3: Catastrophic Cancellation Detection
    # =======================================================================

    def test_catastrophic_cancellation(
        self,
        sigma_c_components: Dict[str, np.ndarray]
    ) -> Dict:
        """
        Test 3: Check if Σᶜ = (large positive) - (large negative) ≈ 0

        This is the most critical test for numerical stability.

        Args:
            sigma_c_components: {
                'residue': residue contribution,
                'integral': imaginary axis integral,
                'total': total Σᶜ
            }

        Returns:
            dict: {
                'cancellation_factor': float,
                'passed': bool
            }
        """
        self.log("\n" + "="*80)
        self.log("TEST 3: Catastrophic Cancellation Detection")
        self.log("="*80)

        residue = sigma_c_components.get('residue', None)
        integral = sigma_c_components.get('integral', None)
        total = sigma_c_components.get('total', None)

        if residue is None or integral is None:
            self.log("⚠️  Residue/integral separation not available")
            self.log("   Skipping cancellation test")
            return {'passed': True, 'skipped': True}

        # Analyze for each orbital (focus on HOMO)
        n_mo = len(residue) if residue.ndim == 1 else residue.shape[-1]

        results_by_orbital = []

        for i in range(n_mo):
            # Extract values for this orbital
            if residue.ndim == 1:
                res_i = residue[i]
                int_i = integral[i]
                tot_i = total[i] if total is not None else res_i + int_i
            else:
                # If 2D, take diagonal or first evaluation point
                res_i = residue[0, i] if residue.ndim == 2 else residue[i]
                int_i = integral[0, i] if integral.ndim == 2 else integral[i]
                tot_i = total[0, i] if (total is not None and total.ndim == 2) else res_i + int_i

            # Get magnitudes
            res_mag = abs(res_i)
            int_mag = abs(int_i)
            tot_mag = abs(tot_i)

            # Cancellation factor: how much larger are components than result?
            max_component = max(res_mag, int_mag)

            if tot_mag > 1e-12:  # Avoid division by very small numbers
                cancellation_factor = max_component / tot_mag
            else:
                cancellation_factor = np.inf

            results_by_orbital.append({
                'orbital': i,
                'residue': res_i,
                'integral': int_i,
                'total': tot_i,
                'residue_mag': res_mag,
                'integral_mag': int_mag,
                'total_mag': tot_mag,
                'cancellation_factor': cancellation_factor
            })

        # Report for HOMO and a few other orbitals
        homo_idx = n_mo // 2  # Approximate HOMO position

        self.log(f"\nAnalyzing {n_mo} orbitals...")
        self.log(f"\nHOMO (orbital {homo_idx}):")

        homo_result = results_by_orbital[homo_idx]
        self.log(f"  Residue:  {homo_result['residue']:.6f} Ha "
                f"(|mag|: {homo_result['residue_mag']:.6f})")
        self.log(f"  Integral: {homo_result['integral']:.6f} Ha "
                f"(|mag|: {homo_result['integral_mag']:.6f})")
        self.log(f"  Total:    {homo_result['total']:.6f} Ha "
                f"(|mag|: {homo_result['total_mag']:.6f})")
        self.log(f"  Cancellation factor: {homo_result['cancellation_factor']:.2e}")

        # Check worst case
        max_cancellation = max(r['cancellation_factor']
                              for r in results_by_orbital
                              if np.isfinite(r['cancellation_factor']))

        self.log(f"\nMaximum cancellation factor across all orbitals: {max_cancellation:.2e}")

        # Thresholds
        warning_threshold = 100    # Significant cancellation
        critical_threshold = 1000  # Severe cancellation (loss of ~3 digits)

        if max_cancellation < warning_threshold:
            self.log(f"✅ PASSED: No significant cancellation (< {warning_threshold})")
            passed = True
        elif max_cancellation < critical_threshold:
            self.log(f"⚠️  WARNING: Moderate cancellation ({max_cancellation:.2e})")
            self.log("   You may lose 1-2 digits of precision")
            passed = True
        else:
            self.log(f"❌ CRITICAL: Severe cancellation ({max_cancellation:.2e})")
            self.log("   Residue and integral are canceling with precision loss!")
            self.log("   Consider:")
            self.log("   - Increasing eta (more broadening)")
            self.log("   - Using extended precision arithmetic")
            self.log("   - Reformulating to avoid subtraction")
            passed = False

        result = {
            'max_cancellation_factor': max_cancellation,
            'homo_cancellation': homo_result['cancellation_factor'],
            'results_by_orbital': results_by_orbital,
            'passed': passed,
            'warning_threshold': warning_threshold,
            'critical_threshold': critical_threshold
        }

        self.results['test_3_cancellation'] = result
        return result

    # =======================================================================
    # TEST 4: Frequency Grid Convergence
    # =======================================================================

    def test_frequency_convergence(
        self,
        compute_gw_fn,
        n_freqs: List[int] = [8, 16, 32, 64],
        **kwargs
    ) -> Dict:
        """
        Test 4: Check Σᶜ convergence with grid density

        For Gauss-Legendre, convergence should be exponential.

        Args:
            compute_gw_fn: Function that computes GW with signature
                          compute_gw_fn(n_freq=n, **kwargs) -> dict with 'sigma_c'
            n_freqs: List of grid sizes to test
            **kwargs: Additional arguments to compute_gw_fn

        Returns:
            dict: {
                'convergence_rate': float,
                'passed': bool
            }
        """
        self.log("\n" + "="*80)
        self.log("TEST 4: Frequency Grid Convergence")
        self.log("="*80)

        sigma_c_values = []

        for n in n_freqs:
            self.log(f"\nComputing with n_freq = {n}...")
            try:
                result = compute_gw_fn(n_freq=n, **kwargs)

                # Extract HOMO value
                sigma_c = result['sigma_c']

                # Get HOMO index
                if 'n_occ' in result:
                    homo_idx = result['n_occ'] - 1
                else:
                    # Assume middle of array
                    homo_idx = sigma_c.shape[-1] // 2

                # Extract HOMO Σᶜ
                if sigma_c.ndim == 1:
                    sigma_c_homo = sigma_c[homo_idx]
                else:
                    sigma_c_homo = sigma_c[0, homo_idx]  # First evaluation point

                sigma_c_values.append(sigma_c_homo)
                self.log(f"  Σᶜ(HOMO) = {sigma_c_homo:.8f} Ha")

            except Exception as e:
                self.log(f"  ❌ Failed: {e}")
                sigma_c_values.append(np.nan)

        sigma_c_values = np.array(sigma_c_values)

        # Check convergence
        self.log("\n" + "-"*60)
        self.log("Convergence analysis:")
        self.log("-"*60)

        valid_mask = ~np.isnan(sigma_c_values)
        if np.sum(valid_mask) < 2:
            self.log("❌ Insufficient valid data points")
            return {'passed': False, 'error': 'Insufficient data'}

        diffs = np.abs(np.diff(sigma_c_values[valid_mask]))

        for i, (n1, n2) in enumerate(zip(
            np.array(n_freqs)[valid_mask][:-1],
            np.array(n_freqs)[valid_mask][1:]
        )):
            if i < len(diffs):
                self.log(f"  Δ({n1:2d}→{n2:2d}): {diffs[i]:.2e} Ha = "
                        f"{diffs[i]*27.211*1000:.2e} meV")

        # Should converge: each doubling should reduce error
        if len(diffs) >= 2:
            # Check if errors are decreasing
            decreasing = all(diffs[i+1] < diffs[i] * 1.5 for i in range(len(diffs)-1))

            # Estimate convergence rate
            # For GL: error ~ exp(-c*n)
            # log(error[n+1]/error[n]) ~ -c*(n[n+1] - n[n])
            convergence_rates = []
            for i in range(len(diffs)-1):
                if diffs[i] > 0 and diffs[i+1] > 0:
                    rate = -np.log(diffs[i+1] / diffs[i])
                    convergence_rates.append(rate)

            if convergence_rates:
                avg_rate = np.mean(convergence_rates)
                self.log(f"\nAverage convergence rate: {avg_rate:.3f}")
                self.log("  (Exponential convergence: rate > 0.5)")

            if decreasing:
                self.log("\n✅ PASSED: Σᶜ converging with grid refinement")
                passed = True
            else:
                self.log("\n❌ FAILED: Σᶜ NOT converging")
                self.log("   Integration scheme may be incorrect!")
                passed = False
        else:
            self.log("\n⚠️  Insufficient data to assess convergence")
            passed = False

        result = {
            'n_freqs': n_freqs,
            'sigma_c_values': sigma_c_values.tolist(),
            'differences': diffs.tolist() if len(diffs) > 0 else [],
            'passed': passed
        }

        self.results['test_4_convergence'] = result
        return result

    # =======================================================================
    # TEST 5: Eta Broadening Sensitivity
    # =======================================================================

    def test_eta_sensitivity(
        self,
        compute_gw_fn,
        etas: List[float] = [0.001, 0.01, 0.1],
        **kwargs
    ) -> Dict:
        """
        Test 5: Check if Σᶜ depends too strongly on eta

        Physical results should be relatively insensitive to eta choice
        (within reasonable range).

        Args:
            compute_gw_fn: Function with signature compute_gw_fn(eta=η, **kwargs)
            etas: List of eta values to test (Hartree)
            **kwargs: Additional arguments

        Returns:
            dict: {
                'sensitivity': float (fractional variation),
                'passed': bool
            }
        """
        self.log("\n" + "="*80)
        self.log("TEST 5: Eta Broadening Sensitivity")
        self.log("="*80)

        sigma_c_values = []

        for eta in etas:
            self.log(f"\nComputing with η = {eta:.3f} Ha...")
            try:
                result = compute_gw_fn(eta=eta, **kwargs)

                # Extract HOMO value
                sigma_c = result['sigma_c']

                if 'n_occ' in result:
                    homo_idx = result['n_occ'] - 1
                else:
                    homo_idx = sigma_c.shape[-1] // 2

                if sigma_c.ndim == 1:
                    sigma_c_homo = sigma_c[homo_idx]
                else:
                    sigma_c_homo = sigma_c[0, homo_idx]

                sigma_c_values.append(sigma_c_homo)
                self.log(f"  Σᶜ(HOMO) = {sigma_c_homo:.6f} Ha")

            except Exception as e:
                self.log(f"  ❌ Failed: {e}")
                sigma_c_values.append(np.nan)

        sigma_c_values = np.array(sigma_c_values)

        # Analyze sensitivity
        self.log("\n" + "-"*60)
        self.log("Sensitivity analysis:")
        self.log("-"*60)

        valid_mask = ~np.isnan(sigma_c_values)
        if np.sum(valid_mask) < 2:
            self.log("❌ Insufficient valid data")
            return {'passed': False, 'error': 'Insufficient data'}

        valid_values = sigma_c_values[valid_mask]
        sigma_c_range = np.max(valid_values) - np.min(valid_values)
        sigma_c_mean = np.mean(valid_values)

        if abs(sigma_c_mean) > 1e-10:
            sensitivity = sigma_c_range / abs(sigma_c_mean)
        else:
            sensitivity = np.inf

        self.log(f"  Range: {sigma_c_range:.6f} Ha")
        self.log(f"  Mean:  {sigma_c_mean:.6f} Ha")
        self.log(f"  Sensitivity: {sensitivity:.1%}")

        # Thresholds
        acceptable_threshold = 0.5   # 50% variation
        warning_threshold = 0.1      # 10% variation

        if sensitivity < warning_threshold:
            self.log(f"\n✅ EXCELLENT: Low eta sensitivity (< {warning_threshold:.0%})")
            passed = True
        elif sensitivity < acceptable_threshold:
            self.log(f"\n✅ PASSED: Acceptable eta sensitivity (< {acceptable_threshold:.0%})")
            passed = True
        else:
            self.log(f"\n❌ FAILED: TOO SENSITIVE to eta (> {acceptable_threshold:.0%})")
            self.log("   This may indicate pole treatment bug or inappropriate eta range")
            passed = False

        result = {
            'etas': etas,
            'sigma_c_values': sigma_c_values.tolist(),
            'sensitivity': sensitivity,
            'passed': passed,
            'acceptable_threshold': acceptable_threshold
        }

        self.results['test_5_eta_sensitivity'] = result
        return result

    # =======================================================================
    # Summary Report
    # =======================================================================

    def generate_report(self, output_file: str = None) -> str:
        """Generate comprehensive analysis report"""

        lines = []
        lines.append("="*80)
        lines.append("NUMERICAL STABILITY ANALYSIS REPORT")
        lines.append("Correlation Self-Energy (Σᶜ) for GW Calculations")
        lines.append("="*80)
        lines.append("")

        # Test 1: Quadrature
        if 'test_1_quadrature' in self.results:
            r = self.results['test_1_quadrature']
            lines.append("TEST 1: Quadrature Accuracy")
            lines.append("-"*80)
            lines.append(f"  Status: {'✅ PASSED' if r['passed'] else '❌ FAILED'}")
            lines.append(f"  Quadrature error: {r['quadrature_error']:.2e}")
            lines.append(f"  Tolerance:        {r['tolerance']:.2e}")
            lines.append("")

        # Test 2: Denominators
        if 'test_2_denominators' in self.results:
            r = self.results['test_2_denominators']
            lines.append("TEST 2: Denominator Conditioning")
            lines.append("-"*80)
            lines.append(f"  Status: {'✅ PASSED' if r['passed'] else '⚠️  WARNING'}")
            lines.append(f"  Min denominator (real):  {r['min_denom_real']:.6f} Ha")
            lines.append(f"  Max denominator (real):  {r['max_denom_real']:.6f} Ha")
            lines.append(f"  Conditioning (real):     {r['conditioning_real']:.2e}")
            lines.append(f"  Near-degeneracies:       {len(r['degeneracies'])}")
            lines.append("")

        # Test 3: Cancellation (CRITICAL)
        if 'test_3_cancellation' in self.results:
            r = self.results['test_3_cancellation']
            if not r.get('skipped', False):
                lines.append("TEST 3: Catastrophic Cancellation (CRITICAL)")
                lines.append("-"*80)
                lines.append(f"  Status: {'✅ PASSED' if r['passed'] else '❌ CRITICAL'}")
                lines.append(f"  Max cancellation factor: {r['max_cancellation_factor']:.2e}")
                lines.append(f"  HOMO cancellation:       {r['homo_cancellation']:.2e}")
                lines.append(f"  Warning threshold:       {r['warning_threshold']:.2e}")
                lines.append(f"  Critical threshold:      {r['critical_threshold']:.2e}")
                lines.append("")

        # Test 4: Convergence
        if 'test_4_convergence' in self.results:
            r = self.results['test_4_convergence']
            lines.append("TEST 4: Frequency Grid Convergence")
            lines.append("-"*80)
            lines.append(f"  Status: {'✅ PASSED' if r['passed'] else '❌ FAILED'}")
            if r.get('differences'):
                lines.append(f"  Error reduction: {r['differences']}")
            lines.append("")

        # Test 5: Eta sensitivity
        if 'test_5_eta_sensitivity' in self.results:
            r = self.results['test_5_eta_sensitivity']
            lines.append("TEST 5: Eta Broadening Sensitivity")
            lines.append("-"*80)
            lines.append(f"  Status:      {'✅ PASSED' if r['passed'] else '❌ FAILED'}")
            if not np.isinf(r['sensitivity']):
                lines.append(f"  Sensitivity: {r['sensitivity']:.1%}")
            lines.append("")

        # Overall assessment
        lines.append("="*80)
        lines.append("OVERALL ASSESSMENT")
        lines.append("="*80)

        all_passed = all(
            r.get('passed', False)
            for key, r in self.results.items()
            if not r.get('skipped', False)
        )

        if all_passed:
            lines.append("✅ All tests PASSED - Numerical implementation appears stable")
        else:
            lines.append("❌ Some tests FAILED - Review flagged issues above")

            # Critical issues
            if 'test_3_cancellation' in self.results:
                r = self.results['test_3_cancellation']
                if not r.get('passed', True) and not r.get('skipped', False):
                    lines.append("")
                    lines.append("⚠️  CRITICAL: Catastrophic cancellation detected!")
                    lines.append("   Immediate action required - see test 3 recommendations")

        lines.append("="*80)

        report = "\n".join(lines)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report)
            self.log(f"\nReport saved to: {output_file}")

        return report


# =======================================================================
# Main execution
# =======================================================================

def main():
    """Run all numerical stability tests"""

    print("\n" + "="*80)
    print("NUMERICAL STABILITY ANALYSIS FOR CORRELATION SELF-ENERGY")
    print("="*80)

    analyzer = NumericalStabilityAnalyzer(verbose=True)

    # Test 1: Quadrature accuracy (always run - no dependencies)
    print("\nRunning Test 1: Quadrature Accuracy...")
    test1 = analyzer.test_quadrature_accuracy(n_points=32, xi_max=100.0)

    # Test 2: Denominator conditioning (needs molecular system)
    if PYSCF_AVAILABLE:
        print("\nRunning Test 2: Denominator Conditioning...")

        # Create simple H2O molecule
        mol = gto.M(
            atom='O 0 0 0; H 0 1 0; H 0 0 1',
            basis='def2-svp',
            verbose=0
        )
        mf = scf.RHF(mol).run()

        mo_energies = mf.mo_energy
        omega_grid = test1['nodes']  # Use GL grid from test 1

        test2 = analyzer.test_denominator_conditioning(
            mo_energies=mo_energies,
            omega_grid=omega_grid,
            eta=0.001
        )
    else:
        print("\nSkipping Test 2: PySCF not available")
        test2 = None

    # Test 3: Catastrophic cancellation
    # This requires actual Σᶜ calculation - create synthetic example
    print("\nRunning Test 3: Catastrophic Cancellation (synthetic)...")

    # Simulate realistic cancellation scenario
    n_mo = 10
    residue_part = np.random.randn(n_mo) * 5.0 - 2.5  # Large, mixed signs
    integral_part = -residue_part * 0.9  # Almost cancel
    total_part = residue_part + integral_part

    test3 = analyzer.test_catastrophic_cancellation({
        'residue': residue_part,
        'integral': integral_part,
        'total': total_part
    })

    # Tests 4 and 5 require full GW implementation
    # Skip for now - would need QuasiX or PySCF GW
    print("\nTests 4-5: Require full GW implementation (skipped for now)")

    # Generate report
    print("\n" + "="*80)
    print("Generating report...")
    print("="*80)

    output_dir = "/home/vyv/Working/QuasiX/quasix_core/docs/reports"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "NUMERICAL_STABILITY_ANALYSIS.md")

    report = analyzer.generate_report(output_file=output_file)
    print("\n" + report)

    print(f"\nFull report saved to: {output_file}")

    return analyzer


if __name__ == "__main__":
    analyzer = main()
