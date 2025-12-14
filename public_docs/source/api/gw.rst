======
GW API
======

This module provides classes and functions for GW calculations.

G0W0Driver
==========

One-shot G0W0 quasiparticle calculations.

Class Reference
---------------

.. py:class:: G0W0Driver(mf, **kwargs)

   One-shot G0W0 driver for quasiparticle calculations.

   :param mf: PySCF mean-field object (RHF, RKS, UHF, UKS)
   :type mf: pyscf.scf.hf.SCF
   :param basis_aux: Auxiliary basis for density fitting. If None, automatically
       selected based on the orbital basis.
   :type basis_aux: str or None
   :param n_freq: Number of frequency grid points (default: 32)
   :type n_freq: int
   :param freq_method: Frequency integration method, 'minimax' or 'gauss-legendre'
       (default: 'minimax')
   :type freq_method: str
   :param eta: Broadening parameter in Hartree (default: 0.001)
   :type eta: float
   :param frozen_core: Number of frozen core orbitals (default: 0)
   :type frozen_core: int
   :param qp_solver: QP equation solver, 'linearized' or 'newton' (default: 'linearized')
   :type qp_solver: str
   :param newton_tol: Convergence tolerance for Newton solver (default: 1e-6)
   :type newton_tol: float

   .. py:method:: kernel()

      Run the G0W0 calculation.

      :returns: G0W0Result object containing quasiparticle energies and related quantities
      :rtype: G0W0Result

      :raises ConvergenceError: If the QP equation solver fails to converge
      :raises BasisError: If the auxiliary basis cannot be determined

Example Usage
-------------

Basic G0W0 calculation:

.. code-block:: python

   from quasix import G0W0Driver
   from pyscf import gto, scf

   # Set up molecule and run HF
   mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='def2-svp')
   mf = scf.RHF(mol).run()

   # Run G0W0
   gw = G0W0Driver(mf)
   result = gw.kernel()

   # Access results
   print(f"HOMO QP energy: {result.homo_qp:.4f} eV")
   print(f"LUMO QP energy: {result.lumo_qp:.4f} eV")
   print(f"QP gap: {result.gap_qp:.4f} eV")

With custom parameters:

.. code-block:: python

   gw = G0W0Driver(
       mf,
       n_freq=64,               # More frequency points
       freq_method='minimax',   # Minimax grid
       qp_solver='newton',      # Newton solver
       newton_tol=1e-7,         # Tight convergence
   )
   result = gw.kernel()

evGWDriver
==========

Eigenvalue self-consistent GW calculations.

Class Reference
---------------

.. py:class:: evGWDriver(mf, **kwargs)

   Eigenvalue self-consistent GW driver.

   Inherits all parameters from G0W0Driver plus:

   :param max_iter: Maximum number of self-consistency iterations (default: 30)
   :type max_iter: int
   :param conv_tol: Convergence tolerance in eV (default: 1e-5)
   :type conv_tol: float
   :param mixing: Mixing/damping parameter for eigenvalue updates (default: 0.5)
   :type mixing: float

   .. py:method:: kernel()

      Run the evGW calculation.

      :returns: evGWResult object with converged quasiparticle energies
      :rtype: evGWResult

Example Usage
-------------

.. code-block:: python

   from quasix import evGWDriver
   from pyscf import gto, dft

   # Set up molecule and run DFT
   mol = gto.M(atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587', basis='def2-tzvp')
   mf = dft.RKS(mol)
   mf.xc = 'pbe0'
   mf.run()

   # Run evGW
   gw = evGWDriver(mf, max_iter=20, conv_tol=1e-4)
   result = gw.kernel()

   print(f"Converged in {result.n_iter} iterations")
   print(f"evGW HOMO: {result.homo_qp:.4f} eV")

G0W0Result
==========

Result object for G0W0 calculations.

.. py:class:: G0W0Result

   Container for G0W0 calculation results.

   .. py:attribute:: qp_energies
      :type: numpy.ndarray

      Quasiparticle energies for all orbitals (eV)

   .. py:attribute:: dft_energies
      :type: numpy.ndarray

      Mean-field (DFT/HF) orbital energies (eV)

   .. py:attribute:: sigma_x
      :type: numpy.ndarray

      Exchange self-energy diagonal elements (eV)

   .. py:attribute:: sigma_c
      :type: numpy.ndarray

      Correlation self-energy at QP energy (eV)

   .. py:attribute:: z_factor
      :type: numpy.ndarray

      Renormalization factors (dimensionless)

   .. py:attribute:: homo_idx
      :type: int

      Index of the HOMO orbital

   .. py:attribute:: lumo_idx
      :type: int

      Index of the LUMO orbital

   .. py:attribute:: homo_qp
      :type: float

      HOMO quasiparticle energy (eV)

   .. py:attribute:: lumo_qp
      :type: float

      LUMO quasiparticle energy (eV)

   .. py:attribute:: gap_qp
      :type: float

      Quasiparticle HOMO-LUMO gap (eV)

   .. py:attribute:: homo_dft
      :type: float

      HOMO mean-field energy (eV)

   .. py:attribute:: converged
      :type: bool

      Whether the calculation converged

evGWResult
==========

Result object for evGW calculations.

.. py:class:: evGWResult

   Container for evGW calculation results. Extends G0W0Result.

   .. py:attribute:: n_iter
      :type: int

      Number of self-consistency iterations

   .. py:attribute:: convergence_history
      :type: list

      List of energy changes per iteration (eV)

Accessing Detailed Results
==========================

.. code-block:: python

   result = gw.kernel()

   # All QP energies
   for i, (e_dft, e_qp) in enumerate(zip(result.dft_energies, result.qp_energies)):
       correction = e_qp - e_dft
       print(f"Orbital {i:3d}: DFT={e_dft:8.4f} eV, QP={e_qp:8.4f} eV, "
             f"Correction={correction:+6.4f} eV")

   # Self-energy analysis for HOMO
   homo = result.homo_idx
   print(f"\nHOMO Self-Energy Breakdown:")
   print(f"  Sigma_x: {result.sigma_x[homo]:.4f} eV")
   print(f"  Sigma_c: {result.sigma_c[homo]:.4f} eV")
   print(f"  Z factor: {result.z_factor[homo]:.4f}")

   # QP correction
   vxc_homo = result.dft_energies[homo] - result.qp_energies[homo]  # Approximate
   total_sigma = result.sigma_x[homo] + result.sigma_c[homo]
   print(f"  Total self-energy: {total_sigma:.4f} eV")

See Also
========

* :doc:`../quickstart` - Getting started tutorial
* :doc:`../theory/gw` - GW theory background
* :doc:`bse` - BSE API documentation
