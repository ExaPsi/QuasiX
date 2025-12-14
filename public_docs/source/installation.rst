============
Installation
============

This guide covers all installation methods for QuasiX.

System Requirements
===================

Operating Systems
-----------------

* **Linux**: Ubuntu 20.04+, CentOS 8+, Debian 11+
* **macOS**: 11 (Big Sur) or later

Hardware Requirements
---------------------

* **CPU**: x86_64 with AVX2 support (for SIMD optimizations)
* **Memory**: 8 GB RAM minimum, 32 GB+ recommended for large molecules
* **Disk**: 2 GB for installation

Software Requirements
---------------------

* **Python**: 3.9 or later (3.11+ recommended)
* **Rust**: 1.70 or later (for building from source)
* **C Compiler**: GCC 9+ or Clang 11+ (for libcint)
* **BLAS/LAPACK**: OpenBLAS, MKL, or system BLAS

Installation Methods
====================

Method 1: pip install (Recommended)
-----------------------------------

The simplest way to install QuasiX:

.. code-block:: bash

   # Create virtual environment (recommended)
   python -m venv quasix-env
   source quasix-env/bin/activate

   # Clone and install
   git clone https://github.com/ExaPsi/QuasiX.git
   cd QuasiX
   pip install -e .

This will automatically:

1. Install Python dependencies (numpy, scipy, pyscf)
2. Build the Rust extension via maturin
3. Install the quasix package

Method 2: Development Installation
----------------------------------

For development with full build control:

.. code-block:: bash

   # Clone repository
   git clone https://github.com/ExaPsi/QuasiX.git
   cd QuasiX

   # Create virtual environment
   python -m venv .venv
   source .venv/bin/activate

   # Install Python dependencies
   pip install numpy scipy pyscf maturin

   # Build Rust extension in release mode
   cd quasix
   maturin develop --release
   cd ..

   # Install Python package in editable mode
   pip install -e .

Method 3: Conda Environment
---------------------------

Using conda for dependency management:

.. code-block:: bash

   # Create conda environment
   conda create -n quasix python=3.11 numpy scipy
   conda activate quasix

   # Install PySCF
   pip install pyscf

   # Clone and install QuasiX
   git clone https://github.com/ExaPsi/QuasiX.git
   cd QuasiX
   pip install -e .

Verifying Installation
======================

Test that QuasiX is installed correctly:

.. code-block:: python

   import quasix
   print(f"QuasiX version: {quasix.__version__}")

   # Quick functional test
   from quasix import G0W0Driver
   from pyscf import gto, scf

   mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
   mf = scf.RHF(mol).run()
   gw = G0W0Driver(mf)
   result = gw.kernel()
   print(f"H2 HOMO: {result.homo_qp:.3f} eV")

Expected output:

.. code-block:: text

   QuasiX version: 0.6.0
   H2 HOMO: -16.xxx eV

Optional Dependencies
=====================

Intel MKL
---------

For improved BLAS performance on Intel processors:

.. code-block:: bash

   pip install mkl mkl-include
   export MKLROOT=/path/to/mkl

OpenMP Threading
----------------

Configure parallel execution:

.. code-block:: bash

   # Set number of threads for linear algebra
   export OMP_NUM_THREADS=4

   # Set number of threads for Rayon (Rust parallelism)
   export RAYON_NUM_THREADS=4

HDF5 Support
------------

For large calculation checkpointing:

.. code-block:: bash

   pip install h5py

Troubleshooting
===============

Rust Compilation Errors
-----------------------

If Rust compilation fails, ensure Rust is up to date:

.. code-block:: bash

   rustup update stable

   # Verify Rust version
   rustc --version  # Should be 1.70+

libcint Not Found
-----------------

PySCF includes libcint, but if issues occur:

.. code-block:: bash

   # Force rebuild of PySCF with libcint
   pip install pyscf --no-binary pyscf --force-reinstall

Memory Errors
-------------

For large molecules, increase stack size:

.. code-block:: bash

   ulimit -s unlimited

BLAS Threading Conflicts
------------------------

If you experience slow performance or hangs:

.. code-block:: bash

   # Set single-threaded BLAS (let Rayon handle parallelism)
   export OMP_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1

Platform-Specific Notes
=======================

Ubuntu/Debian
-------------

Install system dependencies:

.. code-block:: bash

   sudo apt update
   sudo apt install build-essential libopenblas-dev liblapack-dev

macOS
-----

Using Homebrew:

.. code-block:: bash

   brew install openblas
   export LDFLAGS="-L/usr/local/opt/openblas/lib"
   export CPPFLAGS="-I/usr/local/opt/openblas/include"

For Apple Silicon (M1/M2/M3):

.. code-block:: bash

   brew install openblas
   export LDFLAGS="-L/opt/homebrew/opt/openblas/lib"
   export CPPFLAGS="-I/opt/homebrew/opt/openblas/include"

HPC Clusters
------------

Example module loading for typical HPC environments:

.. code-block:: bash

   module load gcc/11.2.0
   module load python/3.10
   module load openblas/0.3.21

   # Or with Intel MKL
   module load intel/2022.1
   module load python/3.10

Building Documentation
======================

To build the documentation locally:

.. code-block:: bash

   cd public_docs
   pip install -r requirements-docs.txt
   sphinx-build -b html source _build/html

   # View documentation
   open _build/html/index.html  # macOS
   xdg-open _build/html/index.html  # Linux

Next Steps
==========

* :doc:`quickstart` - Your first QuasiX calculation
* :doc:`tutorials/index` - In-depth tutorials
* :doc:`api/index` - API reference
