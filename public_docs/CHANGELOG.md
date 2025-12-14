# Changelog

All notable changes to QuasiX will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- BSE benzene benchmark validation
- HPC parallelization improvements
- Advanced theory documentation

## [0.2.0] - 2025-12-03

### Added
- **BSE Implementation** (Sprint 6)
  - BSE kernel infrastructure with exchange (Kx) and direct (Kd) components
  - Davidson iterative eigensolver for BSE Hamiltonian
  - Oscillator strengths and optical absorption spectra
  - Exciton analysis: binding energies, densities, spatial characterization
  - Singlet/triplet spin channel support
  - Tamm-Dancoff approximation (TDA) option
  
- **Python BSE Driver**
  - High-level `BSEDriver` class
  - Spectrum generation with Lorentzian broadening
  - Exciton property analysis

- **Validation**
  - BSE GW100 molecule tests
  - BSE vs TD-DFT comparison
  - Reference data for benzene, water, ammonia

### Changed
- Improved SIMD operations in DF module
- Enhanced polarizability and screening calculations
- Optimized correlation self-energy (AC method)

### Fixed
- Various clippy warnings in quasix_core
- Memory handling for large BSE calculations

## [0.1.0] - 2025-11-28

### Added
- **G₀W₀ Implementation** (Sprint 5)
  - Complete G₀W₀ calculation workflow
  - Analytic Continuation (AC) method
  - Contour Deformation (CD) method
  - Linearized and Newton QP solvers
  - SIMD-optimized tensor operations

- **evGW Implementation**
  - Eigenvalue self-consistent GW
  - Convergence monitoring and damping
  - P₀ update scheme

- **GW100 Validation**
  - Full GW100 benchmark suite
  - Tier 1/2/3 molecule categories
  - Statistical validation framework
  - Comparison with TURBOMOLE reference

- **Frequency Integration**
  - Minimax frequency grids
  - Gauss-Legendre grids
  - Adaptive quadrature

- **Python Interface**
  - `G0W0Driver` class
  - `evGWDriver` class
  - PySCF integration

### Performance
- Parallel frequency integration
- RI/DF approximation for efficiency
- Memory-optimized tensor operations

## [0.0.1] - 2025-10-01

### Added
- Initial project structure
- Rust core library (quasix_core)
- Python bindings via PyO3
- Basic density-fitting infrastructure
- Two-electron integral evaluation
- PySCF adapter for molecular data

---

## Version History Summary

| Version | Date | Highlights |
|---------|------|------------|
| 0.2.0 | 2025-12-03 | BSE implementation, optical spectra |
| 0.1.0 | 2025-11-28 | G₀W₀/evGW, GW100 validation |
| 0.0.1 | 2025-10-01 | Initial release, infrastructure |
