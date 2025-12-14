# QuasiX Python Bindings

This directory contains the Python package for QuasiX, providing PyO3 bindings to the high-performance Rust computational kernel.

## Structure

```
quasix/
├── Cargo.toml          # Rust/PyO3 configuration
├── pyproject.toml      # Python package metadata
├── src/
│   ├── lib.rs          # PyO3 module definition
│   └── python.rs       # Python-specific wrappers
├── quasix/
│   ├── __init__.py     # Python module initialization
│   └── core.py         # High-level Python API
└── tests/
    └── test_basic.py   # Python tests
```

## Building

### Development Build

For development, use maturin to build and install the package in development mode:

```bash
# Install maturin if not already installed
pip install maturin

# Build and install in development mode
cd quasix
maturin develop --release

# Or for debug mode (faster compilation, slower execution)
maturin develop
```

### Production Build

To build a wheel for distribution:

```bash
maturin build --release
```

The wheel will be created in `target/wheels/`.

## Testing

Run the Python tests:

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run tests
pytest tests/
```

## Usage

After building, you can use QuasiX in Python:

```python
from pyscf import gto, scf
from quasix import EvGW, BSE_TDA

# Setup molecule
mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='cc-pvdz')
mf = scf.RHF(mol).run()

# Run evGW
gw = EvGW(mf, auxbasis='cc-pvdz-jkfit')
gw.kernel()

print(f"IP: {gw.get_ionization_potential():.3f} eV")
print(f"EA: {gw.get_electron_affinity():.3f} eV")

# Run BSE for optical excitations
bse = BSE_TDA(gw, nroots=10)
bse.kernel()

print(f"First excitation: {bse.omega[0]:.3f} eV")
```

## Dependencies

### Build Dependencies
- Rust (1.70+)
- maturin (1.7+)
- Python (3.10+)

### Runtime Dependencies
- numpy (≥1.24.0)
- scipy (≥1.10.0)
- pyscf (≥2.3.0)
- h5py (≥3.8.0)

### Optional Dependencies
- matplotlib (for visualization)
- jupyter (for notebooks)

## Architecture

The package uses PyO3 to expose Rust functionality to Python:

1. **Rust Core (`src/`)**: PyO3 bindings that wrap the quasix_core crate
2. **Python API (`quasix/`)**: High-level PySCF-compatible classes
3. **Zero-copy Arrays**: NumPy arrays are passed to Rust without copying when possible
4. **GIL Release**: Long computations release the Python GIL for parallelism

## Performance Considerations

- The Python layer handles I/O and PySCF integration
- Heavy computation happens in the Rust kernel
- DF tensor construction is optimized for memory efficiency
- Use HDF5 for out-of-core operations with large basis sets

## Development Notes

When modifying the PyO3 bindings:

1. Update both `lib.rs` and `python.rs` in the `src/` directory
2. Ensure proper error handling with PyResult
3. Use `py.allow_threads()` for GIL release in compute-intensive functions
4. Maintain PySCF API compatibility in the Python layer
5. Add tests for new functionality

## License

MIT OR Apache-2.0