# QuasiX

[![License](https://img.shields.io/badge/License-MIT%2FApache--2.0-blue.svg)](LICENSE-MIT)
[![Rust](https://img.shields.io/badge/Rust-1.70+-orange.svg)](https://www.rust-lang.org/)
[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)

**QuasiX** is a high-performance implementation of the GW approximation and Bethe-Salpeter Equation (BSE) for calculating quasiparticle energies and optical excitations in molecules and materials.

## Features

- **G₀W₀ and evGW**: One-shot and eigenvalue self-consistent GW calculations
- **BSE Optical Spectra**: Singlet and triplet excitations with oscillator strengths
- **High Performance**: Rust core with SIMD optimizations and parallel execution
- **Python Interface**: Seamless integration with PySCF for molecular calculations
- **GW100 Validated**: Benchmarked against the GW100 test set

## Quick Start

```python
from quasix import G0W0Driver
from pyscf import gto, scf

# Build molecule
mol = gto.M(atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587', basis='def2-svp')
mf = scf.RHF(mol).run()

# Run G0W0
gw = G0W0Driver(mf)
result = gw.kernel()
print(f"HOMO QP energy: {result.homo_qp:.3f} eV")
```

## Installation

See [INSTALL.md](INSTALL.md) for detailed installation instructions.

### Quick Install

```bash
# Clone repository
git clone https://github.com/ExaPsi/QuasiX.git
cd QuasiX

# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install with pip
pip install -e .
```

## Documentation

- [Installation Guide](INSTALL.md)
- [Quick Start Tutorial](QUICKSTART.md)
- [API Reference](API.md)
- [Theory Background](THEORY.md)
- [Benchmarks](BENCHMARKS.md)
- [Contributing](CONTRIBUTING.md)
- [Changelog](CHANGELOG.md)

## Performance

QuasiX achieves significant speedups through:
- Rust implementation with zero-cost abstractions
- SIMD-vectorized tensor operations
- Parallel frequency integration
- Efficient density-fitting (RI) approximation

## Citation

If you use QuasiX in your research, please cite:

```bibtex
@article{quasix2025,
  title={QuasiX: A High-Performance GW/BSE Implementation},
  author={Vchirawongkwin, Viwat and others},
  journal={Journal of Chemical Theory and Computation},
  year={2025},
  note={In preparation}
}
```

## License

QuasiX is dual-licensed under [MIT](LICENSE-MIT) and [Apache 2.0](LICENSE-APACHE) licenses.

## Acknowledgments

QuasiX builds upon:
- [PySCF](https://pyscf.org/) for molecular integrals and SCF
- [ndarray](https://docs.rs/ndarray/) for efficient array operations
- [PyO3](https://pyo3.rs/) for Python-Rust bindings
