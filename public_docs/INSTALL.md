# Installation Guide

## Requirements

### System Requirements
- **Operating System**: Linux (Ubuntu 20.04+, CentOS 8+), macOS 11+
- **CPU**: x86_64 with AVX2 support (for SIMD optimizations)
- **Memory**: 8 GB RAM minimum, 32 GB+ recommended for large molecules
- **Disk**: 2 GB for installation

### Software Requirements
- **Rust**: 1.70 or later
- **Python**: 3.9 or later
- **C Compiler**: GCC 9+ or Clang 11+ (for libcint)
- **BLAS/LAPACK**: OpenBLAS, MKL, or system BLAS

## Installation Methods

### Method 1: pip install (Recommended)

```bash
# Create virtual environment
python -m venv quasix-env
source quasix-env/bin/activate

# Install from source
git clone https://github.com/ExaPsi/QuasiX.git
cd QuasiX
pip install -e .
```

### Method 2: Development Installation

For development with full build control:

```bash
# Clone repository
git clone https://github.com/ExaPsi/QuasiX.git
cd QuasiX

# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install Python dependencies
pip install numpy scipy pyscf maturin

# Build Rust extension
cd quasix
maturin develop --release
cd ..

# Install Python package
pip install -e .
```

### Method 3: Conda Environment

```bash
# Create conda environment
conda create -n quasix python=3.11 numpy scipy
conda activate quasix

# Install PySCF
pip install pyscf

# Clone and install QuasiX
git clone https://github.com/ExaPsi/QuasiX.git
cd QuasiX
pip install -e .
```

## Verifying Installation

```python
import quasix
print(f"QuasiX version: {quasix.__version__}")

# Quick test
from quasix import G0W0Driver
from pyscf import gto, scf

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
mf = scf.RHF(mol).run()
gw = G0W0Driver(mf)
result = gw.kernel()
print(f"H2 HOMO: {result.homo_qp:.3f} eV")
```

Expected output:
```
QuasiX version: 0.1.0
H2 HOMO: -16.xxx eV
```

## Optional Dependencies

### Intel MKL (for improved BLAS performance)

```bash
pip install mkl mkl-include
export MKLROOT=/path/to/mkl
```

### OpenMP (for multi-threading)

```bash
export OMP_NUM_THREADS=4
export RAYON_NUM_THREADS=4
```

## Troubleshooting

### Rust compilation errors

Ensure Rust is up to date:
```bash
rustup update stable
```

### libcint not found

Install libcint system-wide or let PySCF build it:
```bash
pip install pyscf --no-binary pyscf
```

### Memory errors

For large molecules, increase stack size:
```bash
ulimit -s unlimited
```

## Platform-Specific Notes

### Ubuntu/Debian

```bash
sudo apt update
sudo apt install build-essential libopenblas-dev liblapack-dev
```

### macOS

```bash
brew install openblas
export LDFLAGS="-L/usr/local/opt/openblas/lib"
export CPPFLAGS="-I/usr/local/opt/openblas/include"
```

### HPC Clusters

Load appropriate modules:
```bash
module load gcc/11.2.0
module load python/3.10
module load openblas/0.3.21
```
