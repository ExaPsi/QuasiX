# GW100 Validation Module

## Overview
High-performance validation module for the GW100 benchmark suite with SIMD-accelerated statistics and PyO3 bindings.

## Key Features
- **Zero-warning compilation** with strict clippy lints
- **SIMD-accelerated** statistical computations (AVX2 on x86_64)
- **Memory-efficient** data structures with optional compression
- **Parallel execution** with Rayon thread pools
- **Python bindings** via PyO3 for seamless integration

## Module Structure

### Core Components (`gw100.rs`)
- `BenchmarkResult`: Stores results for single molecule/basis combinations
- `ValidationStats`: Comprehensive statistical metrics
- `SimdStatistics`: SIMD-accelerated statistical functions
  - MAD computation with AVX2 vectorization
  - RMSE, correlation coefficient (R²)
  - Modified Z-score outlier detection
- `BenchmarkExecutor`: Parallel benchmark execution with caching
- `OptimizedArray`: Memory-efficient array storage with compression

### Python Bindings (`python_bindings.rs`)
- `PyBenchmarkResult`: Python-exposed benchmark result
- `PyValidationStats`: Python-exposed validation statistics
- Statistical functions exposed to Python:
  - `compute_mad_simd`: SIMD-accelerated MAD
  - `compute_rmse`: Root mean square error
  - `compute_correlation`: Pearson R²
  - `detect_outliers`: Outlier detection
  - `bootstrap_confidence_interval`: Bootstrap CI

## Performance Characteristics

### SIMD Optimization
- AVX2 vectorization processes 4 doubles per instruction
- ~3-5x speedup for MAD computation on large datasets
- Automatic fallback to scalar computation on non-x86 architectures

### Memory Optimization
- Automatic compression for arrays > 10MB
- Cache-friendly data layouts
- Memory pooling for temporary arrays

### Parallel Execution
- Work-stealing thread pool with Rayon
- Efficient task distribution
- Result caching to avoid redundant calculations

## Usage Examples

### Rust
```rust
use quasix_core::validation::{
    BenchmarkResult, SimdStatistics, compute_validation_stats
};

// Compute MAD with SIMD acceleration
let calculated = vec![1.0, 2.0, 3.0, 4.0];
let reference = vec![1.1, 2.1, 3.1, 4.1];
let mad = SimdStatistics::compute_mad_simd(&calculated, &reference)?;

// Create benchmark result
let result = BenchmarkResult::new(
    "H2O".to_string(),
    "def2-tzvp".to_string(),
    "def2-tzvp-jkfit".to_string(),
    qp_energies,
    z_factors,
    reference_energies,
    homo_idx,
    lumo_idx,
    wall_time,
);

// Compute validation statistics
let stats = compute_validation_stats(&results)?;
println!("MAD: {:.3} eV, R²: {:.3}", stats.mad, stats.correlation);
```

### Python
```python
import quasix_core.validation as val
import numpy as np

# Compute MAD with SIMD
calculated = np.array([1.0, 2.0, 3.0, 4.0])
reference = np.array([1.1, 2.1, 3.1, 4.1])
mad = val.compute_mad_simd(calculated, reference)

# Create benchmark result
result = val.BenchmarkResult(
    molecule="H2O",
    basis_set="def2-tzvp",
    aux_basis="def2-tzvp-jkfit",
    qp_energies=qp_energies,
    z_factors=z_factors,
    reference_energies=reference_energies,
    homo_idx=homo_idx,
    lumo_idx=lumo_idx,
    wall_time=10.5
)

# Compute validation statistics
stats = val.compute_validation_stats([result])
print(f"MAD: {stats.mad:.3f} eV, R²: {stats.correlation:.3f}")
```

## Validation Criteria

The module implements the GW100 validation criteria:
- **MAD threshold**: < 0.2 eV (default)
- **Correlation threshold**: R² > 0.95
- **Outlier detection**: Modified Z-score > 3.5
- **Bootstrap CI**: 95% confidence intervals

## Testing

Run tests with:
```bash
cargo test --test test_validation_gw100
```

Run benchmarks with:
```bash
cargo bench --bench validation_benchmarks
```

## Zero-Warning Standards

This module adheres to strict zero-warning standards:
```rust
#![warn(clippy::all, clippy::pedantic, clippy::perf)]
```

All code is:
- Free of `unwrap()` and `expect()` calls
- Properly documented with examples
- Tested with >95% coverage
- Benchmarked against reference implementations