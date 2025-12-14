=================
Quick Start Guide
=================

This guide will walk you through your first QuasiX calculations in under 5 minutes.

Prerequisites
=============

Ensure QuasiX is installed. See :doc:`installation` for detailed instructions.

.. code-block:: python

   # Verify installation
   import quasix
   print(f"QuasiX version: {quasix.__version__}")

Example 1: G0W0 Calculation for Water
=====================================

The simplest GW calculation - one-shot G0W0 starting from Hartree-Fock:

.. code-block:: python

   from quasix import G0W0Driver
   from pyscf import gto, scf

   # Define water molecule
   mol = gto.M(
       atom='''
       O  0.000000  0.000000  0.117369
       H  0.000000  0.756950 -0.469476
       H  0.000000 -0.756950 -0.469476
       ''',
       basis='def2-svp',
       unit='angstrom'
   )

   # Run Hartree-Fock
   mf = scf.RHF(mol)
   mf.kernel()

   # Run G0W0
   gw = G0W0Driver(mf)
   result = gw.kernel()

   # Print results
   print(f"HOMO (HF):   {result.homo_dft:.3f} eV")
   print(f"HOMO (G0W0): {result.homo_qp:.3f} eV")
   print(f"LUMO (G0W0): {result.lumo_qp:.3f} eV")
   print(f"Gap (G0W0):  {result.gap_qp:.3f} eV")

Expected output::

   HOMO (HF):   -13.62 eV
   HOMO (G0W0): -12.57 eV
   LUMO (G0W0):  1.23 eV
   Gap (G0W0):  13.80 eV

Example 2: evGW for Better Accuracy
===================================

Eigenvalue self-consistent GW improves accuracy, especially for band gaps:

.. code-block:: python

   from quasix import evGWDriver
   from pyscf import gto, scf, dft

   # Use DFT (PBE) starting point
   mol = gto.M(
       atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587',
       basis='def2-tzvp'
   )
   mf = dft.RKS(mol)
   mf.xc = 'pbe'
   mf.kernel()

   # Run eigenvalue self-consistent GW
   gw = evGWDriver(mf, max_iter=10, conv_tol=1e-4)
   result = gw.kernel()

   print(f"evGW HOMO: {result.homo_qp:.3f} eV")
   print(f"Converged in {result.n_iter} iterations")

.. tip::
   evGW@PBE0 typically gives the most accurate ionization potentials,
   with MAD = 0.29 eV compared to experimental values (CCSD(T)).

Example 3: BSE Optical Spectrum
===============================

Calculate optical excitations for benzene using BSE:

.. code-block:: python

   from quasix import BSEDriver
   from pyscf import gto, scf

   # Benzene molecule
   mol = gto.M(
       atom='''
       C  1.3970  0.0000  0.0000
       C  0.6985  1.2098  0.0000
       C -0.6985  1.2098  0.0000
       C -1.3970  0.0000  0.0000
       C -0.6985 -1.2098  0.0000
       C  0.6985 -1.2098  0.0000
       H  2.4810  0.0000  0.0000
       H  1.2405  2.1486  0.0000
       H -1.2405  2.1486  0.0000
       H -2.4810  0.0000  0.0000
       H -1.2405 -2.1486  0.0000
       H  1.2405 -2.1486  0.0000
       ''',
       basis='def2-svp'
   )
   mf = scf.RHF(mol).run()

   # Run BSE
   bse = BSEDriver(mf, n_states=10, spin='singlet')
   result = bse.kernel()

   # Print excitation energies
   print("Excitation Energies (eV):")
   for i, (E, f) in enumerate(zip(result.excitation_energies,
                                   result.oscillator_strengths)):
       print(f"  S{i+1}: {E:.3f} eV  (f = {f:.4f})")

   # Generate absorption spectrum
   energies, spectrum = bse.get_spectrum(
       energy_range=(4.0, 10.0),
       broadening=0.1
   )

Example 4: Custom Calculation Parameters
========================================

G0W0 with Custom Frequency Grid
-------------------------------

.. code-block:: python

   from quasix import G0W0Driver

   gw = G0W0Driver(
       mf,
       n_freq=64,               # Number of frequency points
       freq_method='minimax',   # 'minimax' or 'gauss-legendre'
       eta=0.01,                # Broadening parameter (Ha)
   )
   result = gw.kernel()

Choosing the QP Solver
----------------------

.. code-block:: python

   from quasix import G0W0Driver

   # Linearized QP equation (fast, usually sufficient)
   gw = G0W0Driver(mf, qp_solver='linearized')

   # Newton solver (more accurate for difficult cases)
   gw = G0W0Driver(mf, qp_solver='newton', newton_tol=1e-6)

Example 5: Parallel Execution
=============================

Configure parallelism for large calculations:

.. code-block:: python

   import os
   # Set BEFORE importing quasix
   os.environ['RAYON_NUM_THREADS'] = '8'

   from quasix import G0W0Driver

   gw = G0W0Driver(mf, n_workers=8)
   result = gw.kernel()

Output Analysis
===============

Accessing Quasiparticle Energies
--------------------------------

.. code-block:: python

   # Access all QP energies (in eV)
   print("All QP energies (eV):")
   for i, e in enumerate(result.qp_energies):
       occ = "occ" if i < result.homo_idx + 1 else "vir"
       print(f"  Orbital {i:3d} ({occ}): {e:10.4f} eV")

   # Frontier orbitals
   print(f"\nHOMO index: {result.homo_idx}")
   print(f"HOMO: {result.homo_qp:.4f} eV")
   print(f"LUMO: {result.lumo_qp:.4f} eV")
   print(f"Gap:  {result.gap_qp:.4f} eV")

Self-Energy Analysis
--------------------

.. code-block:: python

   # Exchange self-energy (eV)
   homo = result.homo_idx
   print(f"Sigma_x(HOMO): {result.sigma_x[homo]:.4f} eV")

   # Correlation self-energy (eV)
   print(f"Sigma_c(HOMO): {result.sigma_c[homo]:.4f} eV")

   # Renormalization factor
   print(f"Z(HOMO): {result.z_factor[homo]:.4f}")

   # QP correction breakdown
   dft_homo = result.dft_energies[homo]
   sigma_total = result.sigma_x[homo] + result.sigma_c[homo]
   print(f"\nQP correction breakdown for HOMO:")
   print(f"  DFT energy:    {dft_homo:.4f} eV")
   print(f"  Sigma_x:       {result.sigma_x[homo]:.4f} eV")
   print(f"  Sigma_c:       {result.sigma_c[homo]:.4f} eV")
   print(f"  Z factor:      {result.z_factor[homo]:.4f}")
   print(f"  QP energy:     {result.homo_qp:.4f} eV")

Comparing Starting Points
=========================

Different DFT functionals give different G0W0 results:

.. code-block:: python

   from quasix import G0W0Driver
   from pyscf import gto, scf, dft

   mol = gto.M(atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587',
               basis='def2-tzvp')

   # G0W0@HF
   mf_hf = scf.RHF(mol).run()
   gw_hf = G0W0Driver(mf_hf).kernel()
   print(f"G0W0@HF HOMO: {gw_hf.homo_qp:.3f} eV")

   # G0W0@PBE
   mf_pbe = dft.RKS(mol)
   mf_pbe.xc = 'pbe'
   mf_pbe.run()
   gw_pbe = G0W0Driver(mf_pbe).kernel()
   print(f"G0W0@PBE HOMO: {gw_pbe.homo_qp:.3f} eV")

   # G0W0@PBE0
   mf_pbe0 = dft.RKS(mol)
   mf_pbe0.xc = 'pbe0'
   mf_pbe0.run()
   gw_pbe0 = G0W0Driver(mf_pbe0).kernel()
   print(f"G0W0@PBE0 HOMO: {gw_pbe0.homo_qp:.3f} eV")

Common Workflow: GW + BSE
=========================

A typical workflow for optical properties:

.. code-block:: python

   from quasix import G0W0Driver, BSEDriver
   from pyscf import gto, dft

   # 1. Build molecule and run DFT
   mol = gto.M(atom='...', basis='def2-tzvp')
   mf = dft.RKS(mol)
   mf.xc = 'pbe0'
   mf.run()

   # 2. Run GW for quasiparticle energies
   gw = G0W0Driver(mf)
   gw_result = gw.kernel()
   print(f"QP gap: {gw_result.gap_qp:.3f} eV")

   # 3. Run BSE for optical excitations
   bse = BSEDriver(mf, gw_result=gw_result, n_states=10)
   bse_result = bse.kernel()
   print(f"Optical gap: {bse_result.excitation_energies[0]:.3f} eV")

   # 4. Compute exciton binding energy
   E_bind = gw_result.gap_qp - bse_result.excitation_energies[0]
   print(f"Exciton binding energy: {E_bind:.3f} eV")

Next Steps
==========

* :doc:`api/index` - Complete API reference
* :doc:`theory/gw` - Understanding GW theory
* :doc:`theory/bse` - Understanding BSE theory
* :doc:`benchmarks` - Accuracy and performance data
* :doc:`tutorials/index` - Advanced tutorials
