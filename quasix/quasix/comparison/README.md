# CD vs AC Comparison Module

This module provides a comprehensive comparison harness for evaluating Contour Deformation (CD) and Analytic Continuation (AC) methods in GW calculations within the QuasiX framework.

## Features

- **Automated Comparison**: Compare CD and AC methods on single molecules or benchmark suites
- **Statistical Analysis**: Calculate MAD, RMSD, correlation, bootstrap confidence intervals
- **Outlier Detection**: Identify problematic orbitals using statistical methods
- **Performance Metrics**: Track timing, memory usage, and convergence statistics
- **Publication-Quality Reports**: Generate HTML reports with embedded plots
- **PySCF Integration**: Seamless integration with PySCF molecular calculations
- **Parallel Execution**: Run CD and AC calculations in parallel for efficiency

## Quick Start

### Basic Comparison

```python
from quasix.comparison import CDvsACComparator, ComparisonConfig

# Configure comparison
config = ComparisonConfig(
    n_grid_points=40,
    xi_max=50.0,
    mad_threshold=0.05  # 50 meV threshold
)

# Run comparison
comparator = CDvsACComparator(config)
result = comparator.compare_molecule('H2O', basis='cc-pVDZ')

# Check results
print(f"MAD: {result.mad:.3f} eV")
print(f"R²: {result.r_squared:.4f}")
```

### Benchmark Suite

```python
# Run on multiple molecules
molecules = ['H2O', 'NH3', 'CO']
report = comparator.compare_molecules(molecules)

# Generate HTML report
report.generate_html_report("benchmark_report.html")

# Export to JSON
report.to_json("results.json")
```

### Custom PySCF Integration

```python
from pyscf import gto, scf

# Create and run mean-field calculation
mol = gto.M(atom='H2O.xyz', basis='aug-cc-pVTZ')
mf = scf.RHF(mol).run()

# Compare from mean-field object
result = comparator.compare_from_mean_field(mf)
```

## Configuration Options

### Key Parameters

- `n_grid_points`: Number of frequency grid points (default: 40)
- `xi_max`: Maximum frequency for grid (default: 50.0)
- `eta`: Imaginary shift for CD method (default: 0.01)
- `mad_threshold`: MAD threshold for pass/fail (default: 0.05 eV)
- `parallel_execution`: Run CD/AC in parallel (default: True)
- `bootstrap_samples`: Number of bootstrap samples for CI (default: 1000)

### Advanced Options

```python
config = ComparisonConfig(
    # Method-specific parameters
    cd_params={'damping': 0.3},
    ac_params={'n_poles': 10},

    # Statistical options
    outlier_sigma=3.0,
    confidence_level=0.95,

    # Performance options
    n_threads=4,
    memory_limit_gb=32.0,

    # Output options
    save_intermediates=True,
    plot_results=True,
    verbose=2
)
```

## Analysis Metrics

### Statistical Measures
- **MAD**: Mean Absolute Deviation between CD and AC energies
- **RMSD**: Root Mean Square Deviation
- **R²**: Coefficient of determination (correlation)
- **95% CI**: Bootstrap confidence intervals for MAD

### Performance Metrics
- **Timing**: Wall-clock time for each method
- **Memory**: Peak memory usage during calculation
- **Speedup**: Relative performance (AC/CD timing ratio)
- **Convergence**: Iteration count and convergence status

### Diagnostic Information
- **Z-factors**: Quasiparticle renormalization factors
- **Outliers**: Orbitals with significant CD/AC differences
- **Condition number**: Numerical stability indicator for AC

## Visualization

The module provides publication-quality plotting functions:

```python
from quasix.comparison import (
    plot_correlation,
    plot_error_distribution,
    plot_timing_comparison,
    plot_z_factors,
    create_summary_figure
)

# Create various plots
fig1 = plot_correlation(results)
fig2 = plot_error_distribution(results)
fig3 = plot_timing_comparison(results)
fig4 = create_summary_figure(results)

# Save for publication
fig4.savefig('summary.pdf', dpi=300, bbox_inches='tight')
```

## Report Generation

### HTML Reports

Interactive HTML reports include:
- Summary statistics and pass/fail status
- Detailed per-molecule results
- Convergence and performance metrics
- Embedded visualizations
- Outlier analysis

### JSON Export

Machine-readable JSON output contains:
- Complete numerical results
- Configuration parameters
- Metadata (timestamp, platform, versions)
- Reproducibility information

## Testing

Run the test suite:

```bash
pytest tests/test_comparison.py -v
```

Run the demonstration:

```bash
python examples/comparison_demo.py
```

## Module Structure

```
quasix/comparison/
├── __init__.py           # Module exports
├── cd_vs_ac.py          # Main comparison classes
├── report_generator.py  # HTML report generation
├── plotting.py          # Visualization utilities
└── README.md           # This file
```

## Integration with QuasiX

The comparison module integrates with QuasiX's GW implementations:

- Uses QuasiX DF tensor builders for efficiency
- Interfaces with Rust backend via PyO3 when available
- Falls back to mock implementations for testing
- Compatible with QuasiX evGW classes

## Requirements

- Python ≥ 3.8
- NumPy
- PySCF (for quantum chemistry calculations)
- Matplotlib (optional, for plotting)
- SciPy (optional, for advanced statistics)

## Citation

If you use this module in your research, please cite:

```bibtex
@software{quasix2024,
  title = {QuasiX: High-Performance GW/BSE Implementation},
  author = {QuasiX Development Team},
  year = {2024},
  url = {https://github.com/quasix/quasix}
}
```