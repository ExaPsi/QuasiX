# Changelog

All notable changes to QuasiX will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2025-12 (M1 Milestone)

This is the first public release, marking the M1 milestone.

### Added

**GW Implementation**
- G0W0 with contour deformation method
- evGW (eigenvalue self-consistent GW)
- Linearized and Newton QP solvers
- Minimax and Gauss-Legendre frequency grids
- Full GW100 benchmark validation

**BSE Implementation**
- BSE-TDA for optical excitations
- Singlet and triplet spin channels
- Oscillator strength calculation
- Absorption spectrum generation
- Exciton binding energy analysis

**Performance**
- Rust core with BLAS-3 optimization
- AVX-512 SIMD support
- Parallel frequency integration via Rayon
- W matrix caching for evGW

**Python Interface**
- G0W0Driver class
- evGWDriver class
- BSEDriver class
- Seamless PySCF integration

**Documentation**
- ReadTheDocs integration
- API documentation
- Theory overview
- Benchmark reports

### Validation

- GW100 benchmark: MAD = 1.55 meV vs PySCF
- evGW@PBE0: MAD = 0.29 eV vs experiment
- Performance: 8-40x speedup over PySCF

## [0.5.0] - 2025-11

### Added
- evGW implementation with P0 update scheme
- Newton QP solver for difficult cases
- Convergence monitoring and damping
- Starting point dependence analysis

### Changed
- Improved self-energy interpolation
- Enhanced frequency grid optimization

## [0.4.0] - 2025-11

### Added
- BSE-TDA implementation
- Davidson iterative eigensolver
- Oscillator strength calculation
- Transition dipole moments

### Fixed
- Memory handling for large BSE matrices
- Numerical stability in Davidson solver

## [0.3.0] - 2025-10

### Added
- G0W0 contour deformation method
- Analytic continuation with Pade approximants
- Screened interaction W calculation
- Dielectric function from RPA

### Changed
- Optimized polarizability calculation
- Improved memory layout for tensors

## [0.2.0] - 2025-10

### Added
- Density fitting (RI) infrastructure
- Three-center integral evaluation
- Self-energy evaluation framework
- Exchange self-energy calculation

### Performance
- SIMD-optimized tensor contractions
- Parallel integral evaluation

## [0.1.0] - 2025-09

### Added
- Initial project structure
- Rust core library (quasix_core)
- Python bindings via PyO3
- PySCF adapter for molecular data
- Basic linear algebra utilities

### Infrastructure
- Maturin build system
- pytest test framework
- CI/CD pipeline setup

---

## Version History Summary

| Version | Date | Highlights |
|---------|------|------------|
| **0.6.0** | 2025-12 | **M1 Release**: Full GW+BSE, GW100 validated, ReadTheDocs |
| 0.5.0 | 2025-11 | evGW implementation, Newton solver |
| 0.4.0 | 2025-11 | BSE-TDA, optical spectra |
| 0.3.0 | 2025-10 | G0W0 contour deformation |
| 0.2.0 | 2025-10 | Density fitting, self-energy |
| 0.1.0 | 2025-09 | Initial release |

## Roadmap

Planned features for future releases:

**v0.7.0 (M2)**
- Periodic boundary conditions (k-points)
- scGW (fully self-consistent GW)
- Full BSE (beyond TDA)

**v0.8.0 (M3)**
- GPU acceleration (CUDA/ROCm)
- Out-of-core calculations for large systems
- Spin-orbit coupling

**v1.0.0**
- Stable API
- Comprehensive documentation
- Publication in JCTC
