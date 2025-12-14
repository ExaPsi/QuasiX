=======
BSE API
=======

This module provides classes and functions for Bethe-Salpeter Equation (BSE)
calculations of optical excitations.

BSEDriver
=========

BSE-TDA calculations for optical excitations.

Class Reference
---------------

.. py:class:: BSEDriver(mf, **kwargs)

   Bethe-Salpeter Equation driver for optical excitations.

   :param mf: PySCF mean-field object (RHF, RKS)
   :type mf: pyscf.scf.hf.SCF
   :param gw_result: Pre-computed GW result for QP energies. If None, runs G0W0 first.
   :type gw_result: G0W0Result or None
   :param n_states: Number of excited states to compute (default: 10)
   :type n_states: int
   :param spin: Spin channel, 'singlet' or 'triplet' (default: 'singlet')
   :type spin: str
   :param tda: Use Tamm-Dancoff approximation (default: True)
   :type tda: bool
   :param n_occ: Number of active occupied orbitals. If None, uses all.
   :type n_occ: int or None
   :param n_vir: Number of active virtual orbitals. If None, uses all.
   :type n_vir: int or None

   .. py:method:: kernel()

      Solve the BSE eigenvalue problem.

      :returns: BSEResult object containing excitation energies and eigenvectors
      :rtype: BSEResult

   .. py:method:: get_spectrum(energy_range=(0, 20), n_points=1000, broadening=0.1)

      Generate absorption spectrum from computed excitations.

      :param energy_range: Energy range in eV (min, max)
      :type energy_range: tuple
      :param n_points: Number of points in spectrum
      :type n_points: int
      :param broadening: Lorentzian broadening in eV
      :type broadening: float
      :returns: Tuple of (energies, spectrum) arrays
      :rtype: tuple[numpy.ndarray, numpy.ndarray]

Example Usage
-------------

Basic BSE calculation:

.. code-block:: python

   from quasix import BSEDriver
   from pyscf import gto, scf

   # Set up molecule and run HF
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
   for i, (E, f) in enumerate(zip(result.excitation_energies,
                                   result.oscillator_strengths)):
       print(f"S{i+1}: {E:.3f} eV  (f = {f:.4f})")

With pre-computed GW:

.. code-block:: python

   from quasix import G0W0Driver, BSEDriver

   # Run GW first
   gw = G0W0Driver(mf)
   gw_result = gw.kernel()

   # Run BSE with GW result
   bse = BSEDriver(mf, gw_result=gw_result, n_states=10)
   bse_result = bse.kernel()

BSEResult
=========

Result object for BSE calculations.

.. py:class:: BSEResult

   Container for BSE calculation results.

   .. py:attribute:: excitation_energies
      :type: numpy.ndarray

      Excitation energies in eV, shape (n_states,)

   .. py:attribute:: oscillator_strengths
      :type: numpy.ndarray

      Oscillator strengths (dimensionless), shape (n_states,)

   .. py:attribute:: eigenvectors
      :type: numpy.ndarray

      BSE eigenvectors, shape (n_states, n_occ * n_vir)

   .. py:attribute:: transition_dipoles
      :type: numpy.ndarray

      Transition dipole moments in atomic units, shape (n_states, 3)

   .. py:attribute:: n_states
      :type: int

      Number of computed excited states

Generating Absorption Spectra
=============================

.. code-block:: python

   # Compute excitations
   bse = BSEDriver(mf, n_states=20)
   result = bse.kernel()

   # Generate spectrum
   energies, spectrum = bse.get_spectrum(
       energy_range=(4.0, 10.0),  # eV
       n_points=1000,
       broadening=0.1  # eV Lorentzian broadening
   )

   # Plot (requires matplotlib)
   import matplotlib.pyplot as plt
   plt.figure(figsize=(8, 5))
   plt.plot(energies, spectrum, 'b-', linewidth=1.5)
   plt.xlabel('Energy (eV)')
   plt.ylabel('Absorption (arb. units)')
   plt.title('BSE Absorption Spectrum')
   plt.show()

Singlet vs Triplet Excitations
==============================

.. code-block:: python

   # Singlet excitations (optically allowed)
   bse_singlet = BSEDriver(mf, n_states=10, spin='singlet')
   singlet_result = bse_singlet.kernel()

   # Triplet excitations (optically forbidden)
   bse_triplet = BSEDriver(mf, n_states=10, spin='triplet')
   triplet_result = bse_triplet.kernel()

   # Compare
   print("Singlet excitations:")
   for i, E in enumerate(singlet_result.excitation_energies[:5]):
       print(f"  S{i+1}: {E:.3f} eV")

   print("\nTriplet excitations:")
   for i, E in enumerate(triplet_result.excitation_energies[:5]):
       print(f"  T{i+1}: {E:.3f} eV")

   # Singlet-triplet splitting
   st_split = singlet_result.excitation_energies[0] - triplet_result.excitation_energies[0]
   print(f"\nS1-T1 splitting: {st_split:.3f} eV")

Active Space Truncation
=======================

For large systems, restrict the active space:

.. code-block:: python

   # Only include HOMO-4 to LUMO+4
   bse = BSEDriver(
       mf,
       n_states=10,
       n_occ=5,    # 5 occupied orbitals (HOMO-4 to HOMO)
       n_vir=5,    # 5 virtual orbitals (LUMO to LUMO+4)
   )
   result = bse.kernel()

Exciton Analysis
================

.. code-block:: python

   result = bse.kernel()

   # Exciton binding energy
   gw = G0W0Driver(mf)
   gw_result = gw.kernel()
   qp_gap = gw_result.gap_qp
   optical_gap = result.excitation_energies[0]
   binding_energy = qp_gap - optical_gap

   print(f"QP gap: {qp_gap:.3f} eV")
   print(f"Optical gap: {optical_gap:.3f} eV")
   print(f"Exciton binding energy: {binding_energy:.3f} eV")

   # Oscillator strength analysis
   total_f = sum(result.oscillator_strengths)
   print(f"\nTotal oscillator strength: {total_f:.3f}")
   print("(Should equal number of electrons for complete basis)")

   # Dominant transitions
   bright_states = [i for i, f in enumerate(result.oscillator_strengths) if f > 0.01]
   print(f"\nBright states (f > 0.01): {bright_states}")

See Also
========

* :doc:`../quickstart` - Getting started tutorial
* :doc:`../theory/bse` - BSE theory background
* :doc:`gw` - GW API documentation
