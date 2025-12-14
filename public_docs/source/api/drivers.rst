==============
Driver Modules
==============

This page provides an overview of the high-level driver classes that serve as
the main entry points for QuasiX calculations.

Overview
========

QuasiX provides three main driver classes:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Driver
     - Purpose
   * - :class:`G0W0Driver`
     - One-shot G0W0 quasiparticle calculations
   * - :class:`evGWDriver`
     - Eigenvalue self-consistent GW calculations
   * - :class:`BSEDriver`
     - Bethe-Salpeter equation for optical excitations

Driver Hierarchy
================

.. code-block:: text

   PySCF Mean-Field (RHF/RKS)
           |
           v
   +-------+-------+
   |               |
   v               v
   G0W0Driver    evGWDriver
   |               |
   +-------+-------+
           |
           v
       BSEDriver
           |
           v
   Optical Spectrum

Typical Workflow
================

Complete GW+BSE calculation:

.. code-block:: python

   from pyscf import gto, dft
   from quasix import G0W0Driver, BSEDriver

   # Step 1: Mean-field calculation (PySCF)
   mol = gto.M(atom='...', basis='def2-tzvp')
   mf = dft.RKS(mol)
   mf.xc = 'pbe0'
   mf.kernel()

   # Step 2: GW calculation (QuasiX)
   gw = G0W0Driver(mf)
   gw_result = gw.kernel()

   # Step 3: BSE calculation (QuasiX)
   bse = BSEDriver(mf, gw_result=gw_result, n_states=10)
   bse_result = bse.kernel()

   # Step 4: Analysis
   print(f"QP gap: {gw_result.gap_qp:.3f} eV")
   print(f"Optical gap: {bse_result.excitation_energies[0]:.3f} eV")

Common Options
==============

All drivers share common configuration patterns:

Frequency Integration
---------------------

.. code-block:: python

   # Minimax grid (default, most efficient)
   driver = G0W0Driver(mf, freq_method='minimax', n_freq=32)

   # Gauss-Legendre grid (for comparison/validation)
   driver = G0W0Driver(mf, freq_method='gauss-legendre', n_freq=64)

QP Solver Selection
-------------------

.. code-block:: python

   # Linearized solver (fast, usually sufficient)
   driver = G0W0Driver(mf, qp_solver='linearized')

   # Newton solver (more accurate for difficult cases)
   driver = G0W0Driver(mf, qp_solver='newton', newton_tol=1e-6)

Frozen Core Approximation
-------------------------

.. code-block:: python

   # Freeze core orbitals (speeds up calculation for heavy atoms)
   driver = G0W0Driver(mf, frozen_core=2)  # Freeze 2 lowest orbitals

Memory Management
-----------------

.. code-block:: python

   # For very large systems
   driver = G0W0Driver(mf, incore=False, max_memory=8000)  # 8 GB limit

Best Practices
==============

Choosing Starting Point
-----------------------

For G0W0, the starting point affects results:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Starting Point
     - Recommendation
   * - HF
     - Overestimates gaps, good for benchmarking
   * - PBE
     - Underestimates gaps, good for metals
   * - **PBE0**
     - **Best overall accuracy for molecules**
   * - B3LYP
     - Similar to PBE0, slightly larger gaps

Convergence Parameters
----------------------

.. code-block:: python

   # Default (usually sufficient)
   evgw = evGWDriver(mf, max_iter=30, conv_tol=1e-5)

   # Tight convergence (for publication)
   evgw = evGWDriver(mf, max_iter=50, conv_tol=1e-6)

   # Quick estimate
   evgw = evGWDriver(mf, max_iter=10, conv_tol=1e-3)

Basis Set Guidelines
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Basis
     - Accuracy
     - Use Case
   * - def2-SVP
     - Low (~0.3 eV)
     - Quick tests, large molecules
   * - def2-TZVP
     - Good (~0.1 eV)
     - **Production calculations**
   * - def2-QZVP
     - High (~0.02 eV)
     - Benchmark, basis set limit

Error Handling
==============

All drivers raise specific exceptions:

.. code-block:: python

   from quasix.exceptions import ConvergenceError, BasisError

   try:
       result = driver.kernel()
   except ConvergenceError as e:
       print(f"Calculation did not converge: {e}")
       # Options: increase max_iter, try different starting point
   except BasisError as e:
       print(f"Basis set error: {e}")
       # Options: specify auxiliary basis explicitly

Parallelization
===============

QuasiX uses Rayon for parallelism:

.. code-block:: python

   import os
   # Set BEFORE importing quasix
   os.environ['RAYON_NUM_THREADS'] = '8'

   from quasix import G0W0Driver
   # Now uses 8 threads for parallel operations

Performance Tips
================

1. **Use appropriate basis**: def2-TZVP balances accuracy and cost
2. **Freeze core orbitals**: For heavy atoms, saves time with minimal accuracy loss
3. **Use minimax grid**: Default 32 points usually sufficient
4. **Check convergence**: Monitor evGW convergence history
5. **Active space for BSE**: Restrict n_occ/n_vir for large systems

Detailed API
============

For detailed documentation of each driver:

* :doc:`gw` - G0W0Driver and evGWDriver
* :doc:`bse` - BSEDriver
