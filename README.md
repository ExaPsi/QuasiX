# QuasiX

[![License](https://img.shields.io/badge/License-MIT%2FApache--2.0-blue.svg)](LICENSE-MIT)
[![Rust](https://img.shields.io/badge/Rust-1.75%2B-orange)](https://www.rust-lang.org/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![GW100](https://img.shields.io/badge/GW100-Validated-brightgreen)](docs/validation/)

**QuasiX** is a high-performance Rust implementation of the GW approximation and Bethe-Salpeter Equation (BSE) methods for accurate quasiparticle and optical excitation calculations, with seamless Python integration for the PySCF quantum chemistry ecosystem.

## Key Features

- **Production-Ready GW**: G₀W₀ and evGW with contour deformation (CD) and analytic continuation (AC)
- **BSE-TDA**: Complete Bethe-Salpeter equation implementation with Davidson solver and optical properties
- **High Performance**: 8-40x speedup over PySCF's pure Python implementation
- **Validated Accuracy**: Mean absolute deviation (MAD) of 1.55 meV vs PySCF for GW100 benchmark
- **Newton-Raphson Solver**: 77% accuracy improvement over linearized QP solvers (MAD vs TURBOMOLE: 120→28 meV)

## Validation Summary

| Metric | Value | Notes |
|--------|-------|-------|
| GW100 G₀W₀ MAD vs PySCF | 1.55 meV | def2-TZVP basis, CD method |
| evGW@PBE0 MAD vs Expt | 0.47 eV | 22% improvement over G₀W₀@PBE |
| Newton vs Linearized | 77% improvement | MAD: 120→28 meV vs TURBOMOLE |
| BSE-TDA Benzene | < 1e-8 Ha | vs PySCF reference |
| Speedup vs PySCF | 8-40x | System-dependent |

## Installation

### From Source (Recommended)

```bash
# Clone the repository
git clone https://github.com/ExaPsi/QuasiX.git
cd QuasiX

# Set up Python environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install maturin numpy scipy pyscf pytest h5py

# Build and install QuasiX
cd quasix && maturin develop --release && cd ..
```

### Requirements

- **Rust**: 1.75+ (via [rustup](https://rustup.rs/))
- **Python**: 3.10+
- **PySCF**: 2.4+

## Quick Start

### GW Calculation

```python
import quasix
from pyscf import gto, scf

# Run HF calculation with PySCF
mol = gto.M(atom='H 0 0 0; F 0 0 1.1', basis='def2-svp')
mf = scf.RHF(mol).run()

# Perform evGW calculation
result = quasix.gw.run_evgw(mol, mf, auxbasis='def2-svp-ri')

print(f"HOMO QP energy: {result['homo_eV']:.3f} eV")
print(f"LUMO QP energy: {result['lumo_eV']:.3f} eV")
print(f"QP gap: {result['gap_eV']:.3f} eV")
```

### BSE-TDA Calculation

```python
from quasix import run_bse_tda, BSETDAConfig

# Configure BSE-TDA
config = BSETDAConfig(nstates=5, tol=1e-8)

# Run BSE-TDA for optical excitations
result = run_bse_tda(mol, mf, config)

# Print results
print("\nExcited States:")
print("State   Energy (eV)   f_osc")
for i, (e, f) in enumerate(zip(result.eigenvalues_eV, result.oscillator_strengths)):
    print(f"  S{i+1}     {e:8.3f}     {f:.4f}")
```

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     Python Layer (quasix)                       │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐          │
│  │  PySCF Input │  │  QuasiX API  │  │ Visualization │          │
│  │   HF/DFT     │  │  GW & BSE    │  │  & Analysis   │          │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘          │
└─────────┼─────────────────┼─────────────────┼───────────────────┘
          │                 │                 │
          ▼                 ▼                 ▼
┌─────────────────────────────────────────────────────────────────┐
│                   PyO3 Binding Layer                            │
│         Zero-copy arrays │ GIL release │ Error handling         │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Rust Core (quasix_core)                        │
│  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────┐            │
│  │   df    │  │   gw    │  │   bse   │  │  linalg │            │
│  │ DF/RI   │  │ evGW    │  │ BSE-TDA │  │  BLAS   │            │
│  └────┬────┘  └────┬────┘  └────┬────┘  └────┬────┘            │
│       │            │            │            │                  │
│  ┌────┴────┐  ┌────┴────┐  ┌────┴────┐  ┌────┴────┐            │
│  │  freq   │  │selfenergy│  │   qp    │  │   io   │            │
│  │ Grids   │  │ Σˣ, Σᶜ  │  │ Solver  │  │ HDF5   │            │
│  └─────────┘  └─────────┘  └─────────┘  └─────────┘            │
└─────────────────────────────────────────────────────────────────┘
```

## Implementation Status

### Methods

| Method | Molecules | Periodic | Notes |
|--------|-----------|----------|-------|
| G₀W₀ | ✅ | ⏳ Planned | Contour deformation, analytic continuation |
| evGW | ✅ | ⏳ Planned | Eigenvalue self-consistency |
| scGW | ⏳ Planned | ⏳ Planned | Full self-consistency |
| BSE-TDA | ✅ | ⏳ Planned | Davidson solver, NTO analysis |
| BSE-Full | ⏳ Planned | ⏳ Planned | Beyond TDA |

### Frequency Integration

| Method | Description | Default |
|--------|-------------|---------|
| Contour Deformation (CD) | Residues + imaginary axis integral | ✅ Recommended |
| Analytic Continuation (AC) | [15/15] Padé approximants | Available |

## Benchmarks

QuasiX has been validated against the GW100 benchmark set:

- **G₀W₀@HF**: MAD = 1.55 meV vs PySCF (def2-TZVP)
- **evGW@PBE0**: MAD = 0.47 eV vs experiment (vs 0.60 eV for G₀W₀@PBE)
- **Newton-Raphson**: 77% accuracy improvement over linearized solver

### Timing Comparison (vs PySCF)

| Molecule | AO | QuasiX (s) | PySCF (s) | Speedup |
|----------|-----|------------|-----------|---------|
| H₂O | 24 | 0.12 | 1.0 | 8x |
| NH₃ | 30 | 0.14 | 2.4 | 17x |
| CH₄ | 34 | 0.13 | 3.0 | 23x |
| N₂ | 60 | 0.25 | 10.0 | 40x |

## Testing

```bash
# Run Rust tests
cargo test --release

# Run Python tests
pytest tests/ -v

# Run GW100 validation
python tests/benchmarks/gw100/run_tier_benchmark.py --tier 1
```

## Documentation

- [Manuscript Methods Section](docs/manuscripts/drafts/M1/latex/) - JCTC manuscript
- [Theory Derivations](docs/derivations/) - Mathematical foundations
- [Development Guidelines](docs/guidelines/) - Contributing guide
- [Validation Reports](docs/reports/) - Benchmark results

## Theoretical Background

QuasiX implements many-body perturbation theory based on Hedin's equations:

**GW Self-Energy**:
```
Σ(r,r';ω) = i∫ G(r,r';ω+ω')W(r,r';ω') dω'/(2π)
```

**Quasiparticle Equation**:
```
εᵢᵠᴾ = εᵢᴷˢ + ⟨ψᵢ|Σ(εᵢᵠᴾ) - Vxc|ψᵢ⟩
```

**BSE (Tamm-Dancoff)**:
```
A·X = Ω·X
A_ia,jb = (εₐ - εᵢ)δᵢⱼδₐᵦ + Kˣ_ia,jb + Kᵈ_ia,jb
```

## License

QuasiX is dual-licensed under the MIT License and Apache License 2.0. See [LICENSE-MIT](LICENSE-MIT) and [LICENSE-APACHE](LICENSE-APACHE) for details.

## Citation

If you use QuasiX in your research, please cite:

```bibtex
@article{quasix2025,
  title={QuasiX: A High-Performance Rust Implementation of GW and Bethe-Salpeter
         Equation Methods with Seamless Python Integration},
  author={QuasiX Development Team},
  journal={J. Chem. Theory Comput.},
  year={2025},
  note={Submitted}
}
```

## References

1. Hedin, L. (1965). "New method for calculating the one-particle Green's function." *Phys. Rev.* 139, A796.
2. van Setten, M. J., et al. (2015). "GW100: Benchmarking G₀W₀ for molecular systems." *J. Chem. Theory Comput.* 11, 5665.
3. Rohlfing, M. & Louie, S. G. (2000). "Electron-hole excitations and optical spectra from first principles." *Phys. Rev. B* 62, 4927.

## Acknowledgments

- PySCF community for the quantum chemistry framework
- PyO3 developers for Python-Rust bindings
- libcint developers for integral evaluation

---

**QuasiX** - High-performance GW/BSE calculations with Rust efficiency and Python convenience.

https://github.com/ExaPsi/QuasiX
