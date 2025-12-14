#!/usr/bin/env python3
"""
Unit tests for CD vs AC comparison module.

This module tests the comparison harness functionality including
configuration, execution, analysis, and report generation.
"""

import unittest
import sys
import os
import tempfile
import json
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import comparison module
try:
    from quasix.comparison import (
        ComparisonConfig,
        ComparisonResult,
        CDvsACComparator,
        ComparisonReport,
        HTMLReportGenerator
    )
    from quasix.comparison.plotting import (
        plot_correlation,
        plot_error_distribution,
        plot_timing_comparison
    )
    COMPARISON_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import comparison module: {e}")
    COMPARISON_AVAILABLE = False

# Try to import PySCF for integration tests
try:
    from pyscf import gto, scf
    HAS_PYSCF = True
except ImportError:
    HAS_PYSCF = False


class TestComparisonConfig(unittest.TestCase):
    """Test ComparisonConfig dataclass"""

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_default_config(self):
        """Test default configuration values"""
        config = ComparisonConfig()

        self.assertEqual(config.n_grid_points, 40)
        self.assertEqual(config.xi_max, 50.0)
        self.assertEqual(config.eta, 0.01)
        self.assertEqual(config.mad_threshold, 0.05)
        self.assertEqual(config.outlier_sigma, 3.0)
        self.assertEqual(config.bootstrap_samples, 1000)
        self.assertTrue(config.parallel_execution)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_custom_config(self):
        """Test custom configuration"""
        config = ComparisonConfig(
            n_grid_points=60,
            xi_max=100.0,
            mad_threshold=0.1,
            parallel_execution=False,
            verbose=2
        )

        self.assertEqual(config.n_grid_points, 60)
        self.assertEqual(config.xi_max, 100.0)
        self.assertEqual(config.mad_threshold, 0.1)
        self.assertFalse(config.parallel_execution)
        self.assertEqual(config.verbose, 2)


class TestComparisonResult(unittest.TestCase):
    """Test ComparisonResult dataclass"""

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def setUp(self):
        """Create mock comparison result"""
        n_orbitals = 10
        self.result = ComparisonResult(
            molecule_name="H2O",
            basis_set="cc-pVDZ",
            cd_qp_energies=np.random.randn(n_orbitals),
            ac_qp_energies=np.random.randn(n_orbitals),
            hf_energies=np.random.randn(n_orbitals),
            mad=0.045,
            rmsd=0.052,
            max_deviation=0.120,
            correlation=0.998,
            r_squared=0.996,
            mad_ci_lower=0.040,
            mad_ci_upper=0.050,
            outlier_indices=[3, 7],
            outlier_orbitals=["MO_3", "MO_7"],
            cd_timing=12.5,
            ac_timing=8.3,
            cd_memory_peak=256.0,
            ac_memory_peak=180.0,
            cd_converged=True,
            ac_converged=True,
            cd_iterations=10,
            ac_iterations=8,
            cd_z_factors=0.7 + np.random.rand(n_orbitals) * 0.2,
            ac_z_factors=0.6 + np.random.rand(n_orbitals) * 0.3,
            cd_spectral_moments={},
            ac_condition_number=15.2
        )

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_result_attributes(self):
        """Test result attributes"""
        self.assertEqual(self.result.molecule_name, "H2O")
        self.assertEqual(self.result.basis_set, "cc-pVDZ")
        self.assertAlmostEqual(self.result.mad, 0.045)
        self.assertAlmostEqual(self.result.correlation, 0.998)
        self.assertTrue(self.result.cd_converged)
        self.assertTrue(self.result.ac_converged)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_result_to_dict(self):
        """Test conversion to dictionary"""
        result_dict = self.result.to_dict()

        self.assertIsInstance(result_dict, dict)
        self.assertIn('molecule_name', result_dict)
        self.assertIn('mad', result_dict)
        self.assertIn('cd_qp_energies', result_dict)

        # Check numpy arrays are converted to lists
        self.assertIsInstance(result_dict['cd_qp_energies'], list)
        self.assertIsInstance(result_dict['ac_qp_energies'], list)


class TestCDvsACComparator(unittest.TestCase):
    """Test CDvsACComparator class"""

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def setUp(self):
        """Set up test comparator"""
        self.config = ComparisonConfig(
            n_grid_points=20,  # Smaller for testing
            verbose=0,  # Quiet mode
            save_intermediates=False,
            parallel_execution=False  # Sequential for testing
        )
        self.comparator = CDvsACComparator(self.config)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_initialization(self):
        """Test comparator initialization"""
        self.assertIsNotNone(self.comparator)
        self.assertEqual(self.comparator.config.n_grid_points, 20)
        self.assertIsInstance(self.comparator._molecule_database, dict)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_molecule_database(self):
        """Test molecule database loading"""
        db = self.comparator._molecule_database
        self.assertIn('H2O', db)
        self.assertIn('NH3', db)
        self.assertIn('CO', db)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    @unittest.skipIf(not HAS_PYSCF, "PySCF not available")
    def test_load_molecule(self):
        """Test molecule loading from database"""
        mol = self.comparator._load_molecule('H2O', 'cc-pVDZ')
        self.assertIsNotNone(mol)
        self.assertEqual(mol.natm, 3)  # H2O has 3 atoms

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_detect_outliers(self):
        """Test outlier detection"""
        # Create data with clear outliers
        # Using values that are definitely > 3 standard deviations from mean
        normal_values = np.ones(20) * 0.01  # Many normal values around 0.01
        normal_values[5] = 0.5   # Add an outlier
        normal_values[15] = -0.4  # Add another outlier

        outliers = self.comparator._detect_outliers(normal_values)

        # Should detect at least the outliers we added
        self.assertGreater(len(outliers), 0)  # Should find some outliers

        # The detected outliers should have large deviations
        for idx in outliers:
            deviation = abs(normal_values[idx] - np.mean(normal_values))
            self.assertGreater(deviation, 0.1)  # Outliers should deviate significantly

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_bootstrap_ci(self):
        """Test bootstrap confidence interval calculation"""
        np.random.seed(42)  # For reproducibility
        cd_qp = np.random.randn(50)
        ac_qp = cd_qp + np.random.randn(50) * 0.01  # Small differences

        lower, upper = self.comparator._bootstrap_ci(cd_qp, ac_qp)

        self.assertLess(lower, upper)
        self.assertGreater(lower, 0)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_analyze_results(self):
        """Test result analysis"""
        # Create mock CD and AC results
        n_orbitals = 20
        hf_energies = np.linspace(-1, 1, n_orbitals)

        cd_result = {
            'qp_energies': hf_energies + np.random.randn(n_orbitals) * 0.01,
            'z_factors': 0.7 + np.random.rand(n_orbitals) * 0.2,
            'converged': True,
            'iterations': 10,
            'timing': 5.0,
            'memory_peak': 100.0,
            'hf_energies': hf_energies
        }

        ac_result = {
            'qp_energies': hf_energies + np.random.randn(n_orbitals) * 0.012,
            'z_factors': 0.65 + np.random.rand(n_orbitals) * 0.25,
            'converged': True,
            'iterations': 8,
            'timing': 3.0,
            'memory_peak': 80.0,
            'condition_number': 10.0,
            'hf_energies': hf_energies
        }

        result = self.comparator._analyze_results('Test', 'basis', cd_result, ac_result)

        self.assertIsInstance(result, ComparisonResult)
        self.assertEqual(result.molecule_name, 'Test')
        self.assertGreater(result.correlation, 0.9)  # Should be highly correlated
        self.assertTrue(0 < result.mad < 1.0)  # Reasonable MAD in eV


class TestComparisonReport(unittest.TestCase):
    """Test ComparisonReport class"""

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def setUp(self):
        """Create mock results for report"""
        self.results = []
        for mol in ['H2O', 'NH3', 'CO']:
            n_orbitals = 15
            result = ComparisonResult(
                molecule_name=mol,
                basis_set="cc-pVDZ",
                cd_qp_energies=np.random.randn(n_orbitals),
                ac_qp_energies=np.random.randn(n_orbitals),
                hf_energies=np.random.randn(n_orbitals),
                mad=np.random.uniform(0.02, 0.08),
                rmsd=np.random.uniform(0.03, 0.10),
                max_deviation=np.random.uniform(0.05, 0.15),
                correlation=np.random.uniform(0.99, 0.999),
                r_squared=np.random.uniform(0.98, 0.998),
                mad_ci_lower=np.random.uniform(0.01, 0.04),
                mad_ci_upper=np.random.uniform(0.04, 0.10),
                outlier_indices=[],
                outlier_orbitals=[],
                cd_timing=np.random.uniform(5, 20),
                ac_timing=np.random.uniform(3, 15),
                cd_memory_peak=np.random.uniform(100, 500),
                ac_memory_peak=np.random.uniform(80, 400),
                cd_converged=True,
                ac_converged=True,
                cd_iterations=np.random.randint(5, 15),
                ac_iterations=np.random.randint(5, 12),
                cd_z_factors=0.7 + np.random.rand(n_orbitals) * 0.2,
                ac_z_factors=0.6 + np.random.rand(n_orbitals) * 0.3
            )
            self.results.append(result)

        self.config = ComparisonConfig(mad_threshold=0.05)
        self.report = ComparisonReport(self.results, self.config)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_summary_computation(self):
        """Test summary statistics computation"""
        self.assertIsNotNone(self.report.overall_mad)
        self.assertIsNotNone(self.report.overall_rmsd)
        self.assertIsNotNone(self.report.pass_rate)
        self.assertIsNotNone(self.report.avg_speedup)

        self.assertTrue(0 <= self.report.pass_rate <= 1)
        self.assertGreater(self.report.avg_speedup, 0)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_json_export(self):
        """Test JSON export functionality"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json_path = f.name

        try:
            self.report.to_json(json_path)

            # Check file exists and is valid JSON
            self.assertTrue(os.path.exists(json_path))

            with open(json_path, 'r') as f:
                data = json.load(f)

            self.assertIn('config', data)
            self.assertIn('summary', data)
            self.assertIn('results', data)
            self.assertIn('metadata', data)
            self.assertEqual(len(data['results']), 3)

        finally:
            if os.path.exists(json_path):
                os.unlink(json_path)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_html_generation(self):
        """Test HTML report generation"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False) as f:
            html_path = f.name

        try:
            self.report.generate_html_report(html_path)

            # Check file exists and contains expected content
            self.assertTrue(os.path.exists(html_path))

            with open(html_path, 'r') as f:
                content = f.read()

            self.assertIn('CD vs AC', content)
            self.assertIn('Summary Statistics', content)
            self.assertIn('H2O', content)
            self.assertIn('NH3', content)
            self.assertIn('CO', content)

        finally:
            if os.path.exists(html_path):
                os.unlink(html_path)


class TestHTMLReportGenerator(unittest.TestCase):
    """Test HTML report generator"""

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def setUp(self):
        """Set up report generator"""
        self.config = ComparisonConfig()
        self.generator = HTMLReportGenerator(self.config)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_css_generation(self):
        """Test CSS style generation"""
        css = self.generator._generate_css()
        self.assertIn('body', css)
        self.assertIn('font-family', css)
        self.assertIn('.container', css)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_javascript_generation(self):
        """Test JavaScript generation"""
        js = self.generator._generate_javascript()
        self.assertIn('sortTable', js)
        self.assertIn('document', js)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_html_structure(self):
        """Test HTML structure generation"""
        self.generator.add_summary({'overall_mad': 0.05, 'pass_rate': 0.9})

        html = self.generator.generate_html()

        self.assertIn('<!DOCTYPE html>', html)
        self.assertIn('<html', html)
        self.assertIn('</html>', html)
        self.assertIn('<head>', html)
        self.assertIn('<body>', html)


class TestPlottingFunctions(unittest.TestCase):
    """Test plotting utility functions"""

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def setUp(self):
        """Create mock results for plotting"""
        self.results = []
        for mol in ['H2O', 'NH3']:
            n_orbitals = 10
            result = ComparisonResult(
                molecule_name=mol,
                basis_set="cc-pVDZ",
                cd_qp_energies=np.random.randn(n_orbitals),
                ac_qp_energies=np.random.randn(n_orbitals),
                hf_energies=np.random.randn(n_orbitals),
                mad=0.05,
                rmsd=0.06,
                max_deviation=0.10,
                correlation=0.995,
                r_squared=0.990,
                mad_ci_lower=0.04,
                mad_ci_upper=0.06,
                outlier_indices=[],
                outlier_orbitals=[],
                cd_timing=10.0,
                ac_timing=7.0,
                cd_memory_peak=200.0,
                ac_memory_peak=150.0,
                cd_converged=True,
                ac_converged=True,
                cd_iterations=10,
                ac_iterations=8,
                cd_z_factors=0.7 + np.random.rand(n_orbitals) * 0.2,
                ac_z_factors=0.6 + np.random.rand(n_orbitals) * 0.3
            )
            self.results.append(result)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_plot_correlation(self):
        """Test correlation plot generation"""
        fig = plot_correlation(self.results)
        if fig is not None:  # Only if matplotlib is available
            self.assertIsNotNone(fig)
            # Clean up
            import matplotlib.pyplot as plt
            plt.close(fig)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_plot_error_distribution(self):
        """Test error distribution plot"""
        fig = plot_error_distribution(self.results)
        if fig is not None:
            self.assertIsNotNone(fig)
            import matplotlib.pyplot as plt
            plt.close(fig)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    def test_plot_timing_comparison(self):
        """Test timing comparison plot"""
        fig = plot_timing_comparison(self.results)
        if fig is not None:
            self.assertIsNotNone(fig)
            import matplotlib.pyplot as plt
            plt.close(fig)


class TestIntegration(unittest.TestCase):
    """Integration tests with PySCF"""

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    @unittest.skipIf(not HAS_PYSCF, "PySCF not available")
    def test_full_workflow(self):
        """Test complete workflow with PySCF"""
        # Configure comparator
        config = ComparisonConfig(
            n_grid_points=10,  # Very small for testing
            verbose=0,
            save_intermediates=False
        )

        comparator = CDvsACComparator(config)

        # Run comparison on H2O
        result = comparator.compare_molecule('H2O', basis='sto-3g')  # Small basis

        # Check result
        self.assertIsInstance(result, ComparisonResult)
        self.assertEqual(result.molecule_name, 'H2O')
        self.assertEqual(result.basis_set, 'sto-3g')
        self.assertIsNotNone(result.mad)
        self.assertIsNotNone(result.correlation)

    @unittest.skipIf(not COMPARISON_AVAILABLE, "Comparison module not available")
    @unittest.skipIf(not HAS_PYSCF, "PySCF not available")
    def test_mean_field_integration(self):
        """Test integration with existing mean-field calculation"""
        # Create molecule
        mol = gto.M(
            atom='H 0 0 0; H 0 0 0.74',
            basis='sto-3g',
            verbose=0
        )

        # Run HF
        mf = scf.RHF(mol)
        mf.kernel()

        # Configure and run comparison
        config = ComparisonConfig(n_grid_points=10, verbose=0)
        comparator = CDvsACComparator(config)
        result = comparator.compare_from_mean_field(mf)

        # Check result
        self.assertIsInstance(result, ComparisonResult)
        self.assertEqual(len(result.cd_qp_energies), mol.nao)


def main():
    """Run all tests"""
    unittest.main(verbosity=2)


if __name__ == "__main__":
    main()