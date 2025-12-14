=============
API Reference
=============

This section provides detailed documentation for the QuasiX Python API.

.. note::
   QuasiX provides high-level driver classes for common calculations. For most
   users, these drivers provide the simplest interface.

Overview
========

QuasiX is organized into the following modules:

Core Drivers
------------

The main entry points for calculations:

* :class:`G0W0Driver` - One-shot G0W0 calculations
* :class:`evGWDriver` - Eigenvalue self-consistent GW
* :class:`BSEDriver` - Bethe-Salpeter equation for optical excitations

Result Objects
--------------

Calculation results are returned as structured objects:

* :class:`G0W0Result` - G0W0 calculation results
* :class:`evGWResult` - evGW calculation results (extends G0W0Result)
* :class:`BSEResult` - BSE calculation results

Utility Functions
-----------------

Helper functions for common tasks:

* Basis set handling
* Unit conversion
* Analysis tools

Module Documentation
====================

.. toctree::
   :maxdepth: 2

   gw
   bse
   drivers

Quick Reference
===============

G0W0Driver
----------

.. code-block:: python

   from quasix import G0W0Driver
   from pyscf import gto, scf

   mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='def2-svp')
   mf = scf.RHF(mol).run()

   gw = G0W0Driver(
       mf,                          # PySCF mean-field object
       basis_aux=None,              # Auxiliary basis (auto-selected)
       n_freq=32,                   # Number of frequency points
       freq_method='minimax',       # Frequency integration method
       eta=0.001,                   # Broadening parameter (Ha)
       frozen_core=0,               # Frozen core orbitals
       qp_solver='linearized',      # QP solver type
   )
   result = gw.kernel()

   print(f"HOMO: {result.homo_qp:.3f} eV")
   print(f"Gap: {result.gap_qp:.3f} eV")

evGWDriver
----------

.. code-block:: python

   from quasix import evGWDriver

   gw = evGWDriver(
       mf,                          # PySCF mean-field object
       max_iter=30,                 # Maximum iterations
       conv_tol=1e-5,               # Convergence tolerance (eV)
       mixing=0.5,                  # Damping parameter
   )
   result = gw.kernel()

   print(f"Converged: {result.converged}")
   print(f"Iterations: {result.n_iter}")

BSEDriver
---------

.. code-block:: python

   from quasix import BSEDriver

   bse = BSEDriver(
       mf,                          # PySCF mean-field object
       gw_result=None,              # Optional pre-computed GW result
       n_states=10,                 # Number of excited states
       spin='singlet',              # 'singlet' or 'triplet'
       tda=True,                    # Tamm-Dancoff approximation
       n_occ=None,                  # Active occupied orbitals
       n_vir=None,                  # Active virtual orbitals
   )
   result = bse.kernel()

   print(f"First excitation: {result.excitation_energies[0]:.3f} eV")

Result Attributes
=================

G0W0Result / evGWResult
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Attribute
     - Type
     - Description
   * - qp_energies
     - np.ndarray
     - Quasiparticle energies (eV)
   * - dft_energies
     - np.ndarray
     - DFT/HF orbital energies (eV)
   * - sigma_x
     - np.ndarray
     - Exchange self-energy (eV)
   * - sigma_c
     - np.ndarray
     - Correlation self-energy (eV)
   * - z_factor
     - np.ndarray
     - Renormalization factors
   * - homo_idx
     - int
     - HOMO orbital index
   * - lumo_idx
     - int
     - LUMO orbital index
   * - homo_qp
     - float
     - HOMO quasiparticle energy (eV)
   * - lumo_qp
     - float
     - LUMO quasiparticle energy (eV)
   * - gap_qp
     - float
     - QP band gap (eV)
   * - homo_dft
     - float
     - HOMO DFT energy (eV)
   * - converged
     - bool
     - Convergence status
   * - n_iter
     - int
     - Number of iterations (evGW only)

BSEResult
---------

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Attribute
     - Type
     - Description
   * - excitation_energies
     - np.ndarray
     - Excitation energies (eV)
   * - oscillator_strengths
     - np.ndarray
     - Oscillator strengths
   * - eigenvectors
     - np.ndarray
     - BSE eigenvectors
   * - transition_dipoles
     - np.ndarray
     - Transition dipole moments
   * - n_states
     - int
     - Number of computed states

Error Handling
==============

QuasiX provides specific exception classes:

.. code-block:: python

   from quasix.exceptions import ConvergenceError, BasisError

   try:
       result = gw.kernel()
   except ConvergenceError as e:
       print(f"GW did not converge: {e}")
   except BasisError as e:
       print(f"Basis set error: {e}")

See Also
========

* :doc:`../quickstart` - Getting started with QuasiX
* :doc:`../theory/gw` - GW theory background
* :doc:`../theory/bse` - BSE theory background
* :doc:`../benchmarks` - Performance benchmarks
