==========
BSE Theory
==========

The Bethe-Salpeter Equation (BSE) describes optical excitations including
electron-hole interactions (excitons). This page provides a theoretical overview
of the BSE methods implemented in QuasiX.

Introduction
============

While GW calculates charged excitations (adding/removing electrons), the BSE
describes neutral optical excitations where an electron-hole pair is created.

Key concepts:

* **Optical excitation**: Promotion of an electron from occupied to unoccupied state
* **Exciton**: Bound electron-hole pair
* **Exciton binding energy**: Difference between QP gap and optical gap

The BSE naturally captures:

* Electron-hole attraction (screened Coulomb)
* Exchange interactions (for singlet/triplet splitting)
* Oscillator strengths (optical selection rules)

Theoretical Framework
=====================

Two-Particle Green's Function
-----------------------------

The BSE is the equation of motion for the two-particle correlation function:

.. math::

   L(12; 1'2') = G(1,2') G(2,1') + \int d3456 \, G(1,3) G(4,1') \Xi(34;56) L(56;2'2)

where :math:`\Xi` is the BSE kernel containing electron-hole interactions.

BSE Hamiltonian
---------------

In matrix form, the BSE becomes a generalized eigenvalue problem:

.. math::

   \begin{pmatrix} A & B \\ -B^* & -A^* \end{pmatrix}
   \begin{pmatrix} X \\ Y \end{pmatrix}
   = \Omega
   \begin{pmatrix} X \\ Y \end{pmatrix}

where:

* :math:`\Omega` are the excitation energies
* :math:`X, Y` are the electron-hole amplitudes
* :math:`A` is the resonant (particle-hole) block
* :math:`B` is the coupling (particle-hole, hole-particle) block

Matrix Elements
---------------

**Resonant block (A)**:

.. math::

   A_{ia,jb} = (\varepsilon_a^{\text{QP}} - \varepsilon_i^{\text{QP}})\delta_{ij}\delta_{ab}
   + \alpha K^x_{ia,jb} - K^d_{ia,jb}

where:

* :math:`\varepsilon_a^{\text{QP}}, \varepsilon_i^{\text{QP}}` are QP energies from GW
* :math:`K^x` is the exchange (electron-hole exchange) kernel
* :math:`K^d` is the direct (screened Coulomb attraction) kernel
* :math:`\alpha = 2` for singlets, :math:`\alpha = 0` for triplets

**Exchange kernel** (bare Coulomb, electron-hole exchange):

.. math::

   K^x_{ia,jb} = 2 \int d\mathbf{r} d\mathbf{r}' \,
   \psi_i^*(\mathbf{r})\psi_a(\mathbf{r}) \frac{1}{|\mathbf{r}-\mathbf{r}'|}
   \psi_j(\mathbf{r}')\psi_b^*(\mathbf{r}')

.. math::

   K^x_{ia,jb} = 2(ia|jb)

**Direct kernel** (screened Coulomb, electron-hole attraction):

.. math::

   K^d_{ia,jb} = \int d\mathbf{r} d\mathbf{r}' \,
   \psi_i^*(\mathbf{r})\psi_j(\mathbf{r}) W(\mathbf{r},\mathbf{r}')
   \psi_a(\mathbf{r}')\psi_b^*(\mathbf{r}')

.. math::

   K^d_{ia,jb} = (ij|W|ab)

**Coupling block (B)**:

.. math::

   B_{ia,jb} = \alpha K^x_{ia,bj} - K^d_{ia,bj}

Spin Channels
=============

Singlet Excitations
-------------------

For singlet states (:math:`\alpha = 2`), both exchange and direct kernels contribute:

.. math::

   A^{\text{singlet}}_{ia,jb} = \Delta_{ia}\delta_{ij}\delta_{ab} + 2(ia|jb) - (ij|W|ab)

Singlet excitations are **optically bright** (have finite oscillator strength).

Triplet Excitations
-------------------

For triplet states (:math:`\alpha = 0`), only the direct kernel contributes:

.. math::

   A^{\text{triplet}}_{ia,jb} = \Delta_{ia}\delta_{ij}\delta_{ab} - (ij|W|ab)

Triplet excitations are **optically dark** (zero oscillator strength due to spin
selection rules).

.. note::
   The singlet-triplet splitting arises from the exchange kernel :math:`K^x`.
   Typically, triplets lie 0.1-1 eV below singlets.

Tamm-Dancoff Approximation (TDA)
================================

The Tamm-Dancoff approximation sets :math:`B = 0`, reducing the BSE to:

.. math::

   A X = \Omega X

This simplifies the problem to a Hermitian eigenvalue problem with several advantages:

* **Computational**: Smaller matrix, standard diagonalization
* **Numerical**: Always real eigenvalues (no instabilities)
* **Accuracy**: Often sufficient for localized excitations

QuasiX implements BSE-TDA as the default.

**When TDA is accurate:**

* Localized excitations (molecules)
* Large HOMO-LUMO gaps
* Single-reference ground states

**When full BSE is needed:**

* Extended systems with small gaps
* Charge-transfer excitations
* Near-degenerate ground states

Optical Properties
==================

Transition Dipole Moments
-------------------------

The transition dipole moment from the ground state :math:`|0\rangle` to excited
state :math:`|I\rangle`:

.. math::

   \boldsymbol{\mu}_{0I} = \langle 0 | \hat{\boldsymbol{\mu}} | I \rangle
   = \sum_{ia} X_{ia}^I \langle i | \hat{\mathbf{r}} | a \rangle

where :math:`X_{ia}^I` are the BSE eigenvector components for state :math:`I`.

Oscillator Strengths
--------------------

The oscillator strength measures the intensity of optical transitions:

.. math::

   f_I = \frac{2}{3} \Omega_I \sum_{\alpha=x,y,z} |\mu_{0I}^\alpha|^2

**Selection rules:**

* :math:`f > 0`: Optically allowed (bright) transition
* :math:`f = 0`: Optically forbidden (dark) transition

Sum rule:

.. math::

   \sum_I f_I = N_e \quad \text{(Thomas-Reiche-Kuhn sum rule)}

Absorption Spectrum
-------------------

The absorption spectrum is constructed from discrete transitions using broadening:

.. math::

   \sigma(\omega) = \sum_I f_I \, L(\omega - \Omega_I; \gamma)

where :math:`L` is a Lorentzian lineshape function:

.. math::

   L(\omega; \gamma) = \frac{1}{\pi} \frac{\gamma}{(\omega)^2 + \gamma^2}

.. code-block:: python

   # Generate absorption spectrum
   energies, spectrum = bse.get_spectrum(
       energy_range=(4.0, 10.0),  # eV
       n_points=1000,
       broadening=0.1  # eV
   )

Exciton Analysis
================

Exciton Binding Energy
----------------------

The exciton binding energy is the difference between the quasiparticle gap
(fundamental gap) and the optical gap:

.. math::

   E_{\text{bind}} = E_{\text{gap}}^{\text{QP}} - \Omega_1

where :math:`\Omega_1` is the lowest singlet excitation energy.

Typical values:

* Small molecules (H2O): 0.3-0.5 eV
* Conjugated molecules (benzene): 0.8-1.2 eV
* Extended systems: 0.01-0.1 eV

Exciton Wavefunction
--------------------

The exciton wavefunction in real space:

.. math::

   \Psi^I(\mathbf{r}_e, \mathbf{r}_h) = \sum_{ia} X_{ia}^I \psi_a(\mathbf{r}_e) \psi_i^*(\mathbf{r}_h)

This describes the probability amplitude for finding the electron at :math:`\mathbf{r}_e`
and the hole at :math:`\mathbf{r}_h`.

Natural Transition Orbitals
---------------------------

The BSE eigenvector can be analyzed using natural transition orbitals (NTOs),
obtained from SVD of the transition density matrix:

.. math::

   T_{ai} = X_{ia}^I \rightarrow T = U \Sigma V^T

The dominant NTO pairs characterize the excitation:

* **Hole NTO**: :math:`\phi_h = \sum_i V_{i1} \psi_i`
* **Particle NTO**: :math:`\phi_p = \sum_a U_{a1} \psi_a`

Implementation Details
======================

QuasiX Implementation
---------------------

1. **Static Screening**: Uses :math:`W(\omega=0)` from the preceding GW calculation

2. **RI Approximation**: All integrals computed using density fitting

   .. math::

      (ij|W|ab) \approx \sum_{PQ} B^P_{ij} W_{PQ} B^Q_{ab}

3. **Davidson Solver**: Iterative diagonalization for lowest eigenstates

4. **Active Space**: Option to restrict to subset of occupied/virtual orbitals

Computational Workflow
----------------------

1. Run GW calculation to obtain:

   * Quasiparticle energies :math:`\varepsilon_n^{\text{QP}}`
   * Screened interaction :math:`W`

2. Construct BSE matrix elements:

   * Exchange kernel from bare Coulomb integrals
   * Direct kernel from screened Coulomb integrals

3. Solve BSE eigenvalue problem (TDA or full BSE)

4. Compute transition dipoles and oscillator strengths

5. Generate absorption spectrum

.. code-block:: python

   from quasix import G0W0Driver, BSEDriver
   from pyscf import gto, scf

   mol = gto.M(atom='...', basis='def2-svp')
   mf = scf.RHF(mol).run()

   # Step 1: GW calculation
   gw = G0W0Driver(mf)
   gw_result = gw.kernel()

   # Step 2: BSE calculation
   bse = BSEDriver(mf, gw_result=gw_result, n_states=10, spin='singlet')
   bse_result = bse.kernel()

   # Step 3: Analyze results
   print(f"Optical gap: {bse_result.excitation_energies[0]:.3f} eV")
   print(f"Exciton binding: {gw_result.gap_qp - bse_result.excitation_energies[0]:.3f} eV")

References
==========

1. Salpeter, E. E. & Bethe, H. A. (1951). *Phys. Rev.* **84**, 1232.
   "A Relativistic Equation for Bound-State Problems"

2. Rohlfing, M. & Louie, S. G. (2000). *Phys. Rev. B* **62**, 4927.
   "Electron-hole excitations and optical spectra from first principles"

3. Onida, G., Reining, L. & Rubio, A. (2002). *Rev. Mod. Phys.* **74**, 601.
   "Electronic excitations: density-functional versus many-body Green's-function approaches"

4. Blase, X., Duchemin, I. & Jacquemin, D. (2018). *Chem. Soc. Rev.* **47**, 1022.
   "The Bethe-Salpeter equation in chemistry: relations with TD-DFT, applications and challenges"

5. Blase, X., Duchemin, I., Jacquemin, D. & Loos, P.-F. (2020). *J. Phys. Chem. Lett.* **11**, 7371.
   "The Bethe-Salpeter Equation Formalism: From Physics to Chemistry"
