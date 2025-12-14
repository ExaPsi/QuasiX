=========
Tutorials
=========

This section contains in-depth tutorials for using QuasiX.

.. note::
   Tutorials are organized by complexity, from basic to advanced.

Getting Started
===============

Start with these tutorials if you are new to QuasiX:

1. **Basic G0W0 Calculation** - Your first quasiparticle calculation
2. **Comparing Starting Points** - Understanding HF, PBE, PBE0 starting points
3. **evGW for Better Accuracy** - When and how to use self-consistent GW

Available Tutorials
===================

Basic Tutorials
---------------

.. toctree::
   :maxdepth: 1

   basic_g0w0
   starting_points
   evgw_basics

Advanced Topics
---------------

.. toctree::
   :maxdepth: 1

   bse_optical
   basis_convergence
   performance_tuning

Benchmark Reproduction
----------------------

.. toctree::
   :maxdepth: 1

   gw100_reproduction

Planned Tutorials
=================

The following tutorials are under development:

* **GW100 Benchmark Reproduction** - Step-by-step guide to reproducing published results
* **Basis Set Convergence Studies** - How to assess basis set completeness
* **Frequency Integration Parameters** - Optimizing frequency grids
* **Performance Optimization Tips** - Getting the best performance from QuasiX
* **Large System Calculations** - Strategies for molecules with 100+ atoms
* **Periodic Systems** - GW/BSE for crystals (coming in v1.0)

Example Gallery
===============

Quick examples demonstrating specific features:

G0W0 Ionization Potential
-------------------------

.. code-block:: python

   from quasix import G0W0Driver
   from pyscf import gto, scf

   mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='def2-tzvp')
   mf = scf.RHF(mol).run()

   gw = G0W0Driver(mf)
   result = gw.kernel()

   # Ionization potential (negative of HOMO)
   IP = -result.homo_qp
   print(f"Ionization Potential: {IP:.3f} eV")

Optical Gap from BSE
--------------------

.. code-block:: python

   from quasix import G0W0Driver, BSEDriver
   from pyscf import gto, dft

   mol = gto.M(atom='...', basis='def2-svp')
   mf = dft.RKS(mol)
   mf.xc = 'pbe0'
   mf.run()

   gw = G0W0Driver(mf).kernel()
   bse = BSEDriver(mf, gw_result=gw, n_states=5)
   result = bse.kernel()

   print(f"Optical gap: {result.excitation_energies[0]:.3f} eV")

Exciton Binding Energy
----------------------

.. code-block:: python

   # After GW and BSE calculations
   qp_gap = gw_result.gap_qp
   optical_gap = bse_result.excitation_energies[0]
   binding = qp_gap - optical_gap

   print(f"QP gap: {qp_gap:.3f} eV")
   print(f"Optical gap: {optical_gap:.3f} eV")
   print(f"Exciton binding: {binding:.3f} eV")

See Also
========

* :doc:`../quickstart` - Quick start guide
* :doc:`../api/index` - API reference
* :doc:`../benchmarks` - Benchmark data
