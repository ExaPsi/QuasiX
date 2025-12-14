=========
GW Theory
=========

The GW approximation is a many-body perturbation theory method for calculating
quasiparticle energies. This page provides a theoretical overview of the methods
implemented in QuasiX.

Introduction
============

Quasiparticle energies describe the energies required to add or remove an electron
from a many-body system. Unlike Kohn-Sham DFT eigenvalues (which are auxiliary
quantities), quasiparticle energies have direct physical meaning:

* **Ionization Potential (IP)**: Energy to remove an electron (negative of HOMO QP energy)
* **Electron Affinity (EA)**: Energy released when adding an electron (negative of LUMO QP energy)
* **Fundamental Gap**: IP - EA = HOMO - LUMO (QP gap)

The GW approximation provides a systematic approach to computing these quantities
with typical errors of 0.2-0.5 eV compared to experiment.

Many-Body Perturbation Theory
=============================

Starting Point: Green's Function
--------------------------------

The one-particle Green's function :math:`G(\mathbf{r}, \mathbf{r}'; \omega)` contains
information about the single-particle excitation spectrum. Its poles correspond to
quasiparticle energies.

The Green's function satisfies the Dyson equation:

.. math::

   G = G_0 + G_0 \Sigma G

where :math:`G_0` is the non-interacting Green's function and :math:`\Sigma` is the
self-energy.

The GW Approximation
--------------------

In the GW approximation, the self-energy is:

.. math::

   \Sigma(\mathbf{r}, \mathbf{r}'; \omega) = \frac{i}{2\pi} \int d\omega' \,
   G(\mathbf{r}, \mathbf{r}'; \omega + \omega') W(\mathbf{r}, \mathbf{r}'; \omega')

where:

* :math:`G` is the Green's function
* :math:`W = \varepsilon^{-1} v` is the screened Coulomb interaction
* :math:`v` is the bare Coulomb interaction
* :math:`\varepsilon` is the dielectric function

This is the first-order term in the screened interaction expansion, hence "GW".

Quasiparticle Equation
======================

The central equation of GW calculations is the quasiparticle equation:

.. math::

   \varepsilon_n^{\text{QP}} = \varepsilon_n^{\text{MF}} +
   \langle \psi_n | \Sigma(\varepsilon_n^{\text{QP}}) - V_{xc} | \psi_n \rangle

where:

* :math:`\varepsilon_n^{\text{QP}}` is the quasiparticle energy
* :math:`\varepsilon_n^{\text{MF}}` is the mean-field (DFT/HF) energy
* :math:`\Sigma` is the GW self-energy
* :math:`V_{xc}` is the exchange-correlation potential from DFT (or exchange for HF)
* :math:`\psi_n` are the mean-field orbitals

Self-Energy Components
----------------------

The self-energy naturally splits into exchange and correlation:

.. math::

   \Sigma = \Sigma^x + \Sigma^c

**Exchange Self-Energy** (exact exchange):

.. math::

   \Sigma^x_{nm} = -\sum_i^{\text{occ}} (ni|mi)

This is independent of frequency and gives the Hartree-Fock exchange energy
when :math:`n = m`.

**Correlation Self-Energy**:

.. math::

   \Sigma^c_n(\omega) = \frac{i}{2\pi} \int d\omega' \,
   G_n(\omega + \omega') \left[ W(\omega') - v \right]

The correlation self-energy accounts for screening effects beyond mean-field theory.

G0W0: One-Shot GW
=================

In G0W0 (also called "one-shot GW"), the Green's function and screened interaction
are constructed from the starting mean-field orbitals and evaluated once.

Workflow
--------

1. **Mean-field calculation**: Obtain orbitals :math:`\{\psi_n\}` and energies
   :math:`\{\varepsilon_n\}` from DFT or HF

2. **Polarizability**: Construct the non-interacting polarizability

   .. math::

      P_0(\omega) = -i \int dt \, e^{i\omega t} G_0(t) G_0(-t)

3. **Dielectric function**: Compute the RPA dielectric function

   .. math::

      \varepsilon(\omega) = 1 - v P_0(\omega)

4. **Screened interaction**: Invert to get the screened Coulomb

   .. math::

      W(\omega) = \varepsilon^{-1}(\omega) v

5. **Self-energy**: Evaluate exchange and correlation self-energies

6. **QP energies**: Solve the quasiparticle equation

Linearized QP Equation
----------------------

Direct solution of the QP equation requires iteration since :math:`\Sigma` depends on
:math:`\varepsilon_n^{\text{QP}}`. For efficiency, QuasiX uses the linearized form:

.. math::

   \varepsilon_n^{\text{QP}} \approx \varepsilon_n^{\text{MF}} + Z_n
   \text{Re}\langle \psi_n | \Sigma(\varepsilon_n^{\text{MF}}) - V_{xc} | \psi_n \rangle

where the **renormalization factor** is:

.. math::

   Z_n = \left( 1 - \frac{\partial \text{Re}\Sigma^c_n}{\partial \omega}
   \bigg|_{\omega=\varepsilon_n^{\text{MF}}} \right)^{-1}

The Z-factor is typically 0.7-0.9 and accounts for the energy dependence of the
self-energy.

Newton Solver
-------------

For difficult cases (e.g., strong satellites), QuasiX also provides a Newton solver
that iteratively solves the full QP equation:

.. math::

   \varepsilon_n^{(k+1)} = \varepsilon_n^{(k)} -
   \frac{f(\varepsilon_n^{(k)})}{f'(\varepsilon_n^{(k)})}

where :math:`f(\omega) = \omega - \varepsilon_n^{\text{MF}} -
\text{Re}\langle \psi_n | \Sigma(\omega) - V_{xc} | \psi_n \rangle`.

evGW: Eigenvalue Self-Consistent GW
===================================

In G0W0, results depend on the starting point (HF, PBE, PBE0, etc.). Eigenvalue
self-consistent GW (evGW) reduces this dependence by iterating the QP energies.

Algorithm
---------

1. Initialize: :math:`\varepsilon_n^{(0)} = \varepsilon_n^{\text{MF}}`

2. For iteration :math:`k = 1, 2, \ldots`:

   a. Compute :math:`P_0^{(k)}` using energies :math:`\varepsilon_n^{(k-1)}`
   b. Compute :math:`W^{(k)}` from :math:`P_0^{(k)}`
   c. Evaluate :math:`\Sigma^{(k)}`
   d. Solve QP equation for :math:`\varepsilon_n^{(k)}`
   e. Check convergence: :math:`|\varepsilon_n^{(k)} - \varepsilon_n^{(k-1)}| < \text{tol}`

3. Output converged :math:`\varepsilon_n^{\text{QP}}`

.. note::
   In evGW, only the eigenvalues are updated; the orbitals remain frozen. This
   typically converges in 5-15 iterations with tolerance 10^-4 eV.

Starting Point Dependence
-------------------------

Different starting points give different results:

.. list-table:: Starting Point Effects (H2O/def2-TZVP)
   :header-rows: 1
   :widths: 25 25 25 25

   * - Method
     - Starting Point
     - HOMO (eV)
     - Starting Point Variation
   * - G0W0
     - HF
     - -12.72
     - 0.4 eV
   * - G0W0
     - PBE
     - -12.38
     - (spread)
   * - G0W0
     - PBE0
     - -12.57
     -
   * - evGW
     - HF
     - -12.55
     - 0.05 eV
   * - evGW
     - PBE
     - -12.52
     - (reduced)
   * - evGW
     - PBE0
     - -12.54
     -

evGW significantly reduces starting point dependence.

Density Fitting (RI Approximation)
==================================

QuasiX uses the Resolution of Identity (RI) approximation for computational
efficiency. Orbital products are expanded in an auxiliary basis:

.. math::

   \psi_i(\mathbf{r})\psi_a(\mathbf{r}) \approx \sum_P C_{ia}^P \chi_P(\mathbf{r})

where :math:`\{\chi_P\}` is the auxiliary basis. The fitting coefficients are:

.. math::

   C_{ia}^P = \sum_Q (ia|Q) (Q|P)^{-1/2}

This reduces the scaling of four-center integrals :math:`(ia|jb)` to products of
three-center quantities:

.. math::

   (ia|jb) \approx \sum_{PQ} C_{ia}^P S_{PQ}^{-1} C_{jb}^Q = \sum_P C_{ia}^P C_{jb}^P

**Benefit**: O(N^5) -> O(N^4) scaling with negligible loss of accuracy (~1 meV).

Frequency Integration
=====================

Imaginary Frequency Formulation
-------------------------------

The correlation self-energy is evaluated on the imaginary frequency axis where
functions are smooth:

.. math::

   \Sigma^c(\omega) = \frac{1}{\pi} \int_0^\infty d\omega' \,
   \text{Im}[W(i\omega')] F(\omega, i\omega')

where :math:`F` contains Green's function contributions.

Analytic Continuation
---------------------

Results on the imaginary axis are analytically continued to the real axis using:

* **Pade approximants** (default): Fit rational function to imaginary-axis data
* **Two-pole model**: For difficult cases with strong satellite features

Frequency Grids
---------------

QuasiX supports two frequency integration schemes:

1. **Minimax grids** (default): Optimized for GW, requiring minimal points (32 typical)
2. **Gauss-Legendre**: Standard quadrature, requiring more points (64+)

.. code-block:: python

   # Minimax (recommended)
   gw = G0W0Driver(mf, freq_method='minimax', n_freq=32)

   # Gauss-Legendre
   gw = G0W0Driver(mf, freq_method='gauss-legendre', n_freq=64)

Implementation Details
======================

QuasiX implements GW using the contour deformation approach:

1. **Contour Deformation**: Splits the frequency integral into contributions along
   the real and imaginary axes

2. **RI-GW**: All quantities computed using density-fitted integrals

3. **Parallel Frequency Integration**: Frequency points evaluated in parallel
   using Rayon

4. **SIMD Optimization**: Tensor operations use AVX-512 where available

References
==========

1. Hedin, L. (1965). *Phys. Rev.* **139**, A796.
   "New Method for Calculating the One-Particle Green's Function with Application
   to the Electron-Gas Problem"

2. Hybertsen, M. S. & Louie, S. G. (1986). *Phys. Rev. B* **34**, 5390.
   "Electron correlation in semiconductors and insulators: Band gaps and
   quasiparticle energies"

3. Aryasetiawan, F. & Gunnarsson, O. (1998). *Rep. Prog. Phys.* **61**, 237.
   "The GW method"

4. van Setten, M. J. et al. (2015). *J. Chem. Theory Comput.* **11**, 5665.
   "GW100: Benchmarking G0W0 for Molecular Systems"

5. Golze, D., Dvorak, M. & Rinke, P. (2019). *Front. Chem.* **7**, 377.
   "The GW Compendium: A Practical Guide to Theoretical Photoemission Spectroscopy"
