# QuasiX Core

[![Rust](https://img.shields.io/badge/Rust-1.75%2B-orange)](https://www.rust-lang.org/)
[![Clippy](https://img.shields.io/badge/Clippy-0%20warnings-brightgreen)](https://github.com/rust-lang/rust-clippy)

High-performance Rust implementation of GW/BSE methods for quantum chemistry calculations.

## Overview

`quasix_core` is the computational kernel of the QuasiX project, providing efficient implementations of many-body perturbation theory methods including:

- **G₀W₀/evGW**: Eigenvalue self-consistent GW for quasiparticle calculations
- **BSE-TDA**: Bethe-Salpeter equation with Tamm-Dancoff approximation for optical excitations
- **Parallel Computing**: Rayon-based parallelization with BLAS thread management
- **SIMD Optimization**: AVX2/AVX-512 vectorization for hot paths

## Module Structure

### Core Modules (Production-Ready)

| Module | Description | Status |
|--------|-------------|--------|
| **`common`** | Shared utilities, error types, physical constants | ✅ Production |
| **`df`** | Density fitting with PySCF-validated RI tensors | ✅ Production |
| **`freq`** | Frequency grids (Gauss-Legendre), contour deformation | ✅ Production |
| **`dielectric`** | Polarizability P⁰(iω), screening W(iω) | ✅ Production |
| **`selfenergy`** | Exchange Σˣ, correlation Σᶜ (CD + AC) | ✅ Production |
| **`qp`** | QP solver with Newton-Raphson, DIIS acceleration | ✅ Production |
| **`gw`** | Complete evGW driver with convergence control | ✅ Production |
| **`bse`** | BSE-TDA kernel, Davidson solver, optical properties | ✅ Production |
| **`analysis`** | NTO decomposition, exciton analysis, spectral functions | ✅ Production |
| **`io`** | HDF5/JSON I/O with schema validation | ✅ Production |
| **`linalg`** | BLAS/LAPACK integration, Cholesky, eigensolvers | ✅ Production |
| **`pbc`** | Periodic boundary conditions | ⏳ Sprint 7 |

### Key Features

#### GW Implementation
- **Contour Deformation**: Residues + imaginary axis integral
- **Analytic Continuation**: Padé approximants with configurable order
- **evGW Self-Consistency**: Eigenvalue updates with P⁰ reconstruction
- **Z-Factor Calculation**: Quasiparticle weights with physical bounds

#### BSE-TDA Implementation
- **Matrix-Free Kernel**: O(N_trans × N_aux) memory complexity
- **Exchange Kernel**: K^x via density fitting, singlet/triplet support
- **Direct Kernel**: K^d with W(ω=0) screening
- **Davidson Solver**: Preconditioned with MGS2 orthogonalization
- **Optical Properties**: Oscillator strengths, transition dipoles, spectra

#### Performance Optimizations
- **SIMD Vectorization**: AVX2/FMA for frequency operations (2.5x speedup)
- **Parallel Execution**: Rayon with BlasThreadGuard (85% efficiency on 16 cores)
- **Memory Efficiency**: Blocked algorithms (60% reduction)
- **BLAS Integration**: Direct CBLAS calls for matrix operations
- **Cholesky Inversion**: SPD matrix inversion (1.5-2x vs LU)

## Building

```bash
# Debug build
cargo build

# Release build with optimizations
cargo build --release

# Run tests
cargo test --release

# Run specific module tests
cargo test --release --lib gw
cargo test --release --lib bse

# Check code quality (zero warnings!)
cargo clippy --all-features

# Generate documentation
cargo doc --no-deps --open
```

## Code Quality

**Zero Warnings Policy**: All code passes `cargo clippy` with pedantic lints enabled.

```rust
// Enabled in lib.rs
#![warn(clippy::all, clippy::pedantic, clippy::perf)]
```

### Lint Configuration

Scientific computing-specific allowances:
- `many_single_char_names`: Mathematical notation (i, j, k, p, q)
- `cast_precision_loss`: Array indexing to f64
- `too_many_arguments`: Complex scientific functions

## Performance Metrics

### GW Calculation (40 basis functions)
| Metric | Value |
|--------|-------|
| Single iteration | ~8.5 ms |
| Parallel efficiency (8 cores) | 6.2x speedup |
| Parallel efficiency (16 cores) | 10.8x speedup |
| Memory reduction | 60% via blocking |
| Numerical accuracy | < 1e-8 vs PySCF |

### BSE-TDA (Benzene/STO-3G)
| Metric | Value |
|--------|-------|
| Kernel build | ~50 ms |
| Davidson solve (5 states) | ~150 ms |
| Eigenvalue accuracy | < 1e-12 Ha |
| Hermiticity error | < 1e-12 |

## Dependencies

- `ndarray`: N-dimensional arrays for scientific computing
- `ndarray-linalg`: LAPACK bindings for eigensolvers
- `num-complex`: Complex number support
- `rayon`: Data parallelism
- `serde`: Serialization framework
- `thiserror`/`anyhow`: Error handling
- `hdf5`: HDF5 file I/O

## API Overview

### GW Calculation

```rust
use quasix_core::gw::EvGWDriver;
use quasix_core::qp::ConvergenceCriteria;

// Configure evGW
let criteria = ConvergenceCriteria {
    max_iter: 30,
    energy_tol: 1e-6,
    use_diis: true,
    ..Default::default()
};

// Run calculation
let driver = EvGWDriver::new(mol_data, df_tensors, freq_grid);
let result = driver.run(&criteria)?;

println!("HOMO QP: {:.4} eV", result.qp_energies[nocc - 1] * HARTREE_TO_EV);
```

### BSE-TDA Calculation

```rust
use quasix_core::bse::{BSETDAKernel, DavidsonConfig, SpinType};

// Build BSE kernel
let kernel = BSETDAKernel::new(
    df_tensor,
    qp_energies_occ,
    qp_energies_virt,
    screening_w0,
    SpinType::Singlet,
);

// Solve for excited states
let config = DavidsonConfig {
    nstates: 5,
    max_iter: 100,
    tol: 1e-8,
    ..Default::default()
};

let (eigenvalues, eigenvectors) = kernel.solve_davidson(&config)?;

// Compute oscillator strengths
let fosc = kernel.oscillator_strengths(&eigenvectors, &transition_dipoles);
```

## Testing

```bash
# Run all tests
cargo test --release

# Run specific module tests
cargo test --release --lib df
cargo test --release --lib dielectric
cargo test --release --lib selfenergy
cargo test --release --lib gw
cargo test --release --lib bse

# Run with verbose output
cargo test --release -- --nocapture

# Run benchmarks
cargo bench --bench polarizability_simd_bench
cargo bench --bench screening_bench
```

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 0.6.0 | 2025-12-03 | BSE-TDA complete, zero clippy warnings |
| 0.5.0 | 2025-11 | evGW with GW100 validation |
| 0.4.0 | 2025-09 | G₀W₀ with CD and AC |
| 0.3.0 | 2025-09 | DF/RI infrastructure |
| 0.2.0 | 2025-08 | PyO3 bindings |
| 0.1.0 | 2025-07 | Initial Rust core |

## License

BSD-3-Clause

## Authors

QuasiX Development Team

## Contributing

This is part of the QuasiX project. See the main repository for contribution guidelines.

---

**Current Version**: 0.6.0 - Production-ready GW/BSE for molecules
