.. QuasiX documentation master file

==========================================
QuasiX: High-Performance GW/BSE for Python
==========================================

.. image:: https://img.shields.io/badge/License-MIT%2FApache--2.0-blue.svg
   :target: https://github.com/ExaPsi/QuasiX/blob/main/LICENSE-MIT
   :alt: License

.. image:: https://img.shields.io/badge/Rust-1.70+-orange.svg
   :target: https://www.rust-lang.org/
   :alt: Rust

.. image:: https://img.shields.io/badge/Python-3.9+-blue.svg
   :target: https://www.python.org/
   :alt: Python

**QuasiX** is a high-performance implementation of the GW approximation and
Bethe-Salpeter Equation (BSE) for calculating quasiparticle energies and optical
excitations in molecules and materials.

.. note::
   QuasiX v0.6.0 is the current release (M1 milestone). The project is under
   active development with a manuscript in preparation for *Journal of Chemical
   Theory and Computation*.

Key Features
============

* **G0W0 and evGW**: One-shot and eigenvalue self-consistent GW calculations
* **BSE Optical Spectra**: Singlet and triplet excitations with oscillator strengths
* **High Performance**: Rust core with SIMD optimizations (8-40x faster than PySCF)
* **Python Interface**: Seamless integration with PySCF for molecular calculations
* **GW100 Validated**: Benchmarked against the GW100 test set (MAD = 1.55 meV vs PySCF)

Quick Example
=============

.. code-block:: python

   from quasix import G0W0Driver
   from pyscf import gto, scf

   # Build molecule
   mol = gto.M(
       atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587',
       basis='def2-svp'
   )
   mf = scf.RHF(mol).run()

   # Run G0W0
   gw = G0W0Driver(mf)
   result = gw.kernel()
   print(f"HOMO QP energy: {result.homo_qp:.3f} eV")

Performance Highlights
======================

QuasiX achieves significant speedups through:

* Rust implementation with zero-cost abstractions
* SIMD-vectorized tensor operations (AVX-512 support)
* Parallel frequency integration
* Efficient density-fitting (RI) approximation

.. list-table:: Performance Comparison (G₀W₀@PBE/def2-TZVP, 64 threads)
   :header-rows: 1
   :widths: 25 20 20 20

   * - Molecule
     - QuasiX
     - PySCF
     - Speedup
   * - H₂O (43 AO)
     - 1.0 s
     - 8.0 s
     - 8.3×
   * - CH₄ (55 AO)
     - 1.3 s
     - 29.8 s
     - 23.3×
   * - N₂ (62 AO)
     - 1.7 s
     - 67.6 s
     - **40.1×**

Validation
==========

QuasiX has been rigorously validated against the GW100 benchmark set:

* **G₀W₀ vs PySCF**: MAD = 1.55 meV (sub-meV agreement, 11 molecules)
* **evGW@PBE0/def2-TZVP**: MAD = 0.29 eV, MSE = +0.14 eV vs NIST experimental IPs (50 molecules)
* **evGW vs TURBOMOLE**: MAD = 27.7 meV (Newton vs graphical QP solver)
* See :doc:`benchmarks` for detailed validation results with actual data.

Documentation Contents
======================

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   tutorials/index
   benchmarks

.. toctree::
   :maxdepth: 2
   :caption: Theory Background

   theory/gw
   theory/bse

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   changelog

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citation
========

If you use QuasiX in your research, please cite:

.. code-block:: bibtex

   @article{quasix2025,
     title={QuasiX: A High-Performance Rust Implementation of GW and
            Bethe-Salpeter Equation Methods with Seamless Python Integration},
     author={Vchirawongkwin, Viwat},
     journal={J. Chem. Theory Comput.},
     year={2025},
     note={In preparation}
   }

License
=======

QuasiX is dual-licensed under `MIT <https://github.com/ExaPsi/QuasiX/blob/main/LICENSE-MIT>`_
and `Apache 2.0 <https://github.com/ExaPsi/QuasiX/blob/main/LICENSE-APACHE>`_ licenses.

Acknowledgments
===============

QuasiX builds upon:

* `PySCF <https://pyscf.org/>`_ for molecular integrals and SCF
* `ndarray <https://docs.rs/ndarray/>`_ for efficient array operations
* `PyO3 <https://pyo3.rs/>`_ for Python-Rust bindings
