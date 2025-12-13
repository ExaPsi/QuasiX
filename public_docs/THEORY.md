# Theory Background

This document provides a brief overview of the theoretical methods implemented in QuasiX.

## The GW Approximation

### Many-Body Perturbation Theory

The GW approximation is derived from many-body perturbation theory (MBPT) and provides quasiparticle energies that accurately describe charged excitations (ionization potentials and electron affinities).

The quasiparticle equation reads:

$$
\varepsilon_n^{QP} = \varepsilon_n^{DFT} + \langle \psi_n | \Sigma(\varepsilon_n^{QP}) - V_{xc} | \psi_n \rangle
$$

where $\Sigma$ is the self-energy operator and $V_{xc}$ is the exchange-correlation potential from DFT.

### Self-Energy

The GW self-energy is:

$$
\Sigma(r, r'; \omega) = \frac{i}{2\pi} \int d\omega' \, G(r, r'; \omega + \omega') W(r, r'; \omega')
$$

where:
- $G$ is the Green's function
- $W = \varepsilon^{-1} v$ is the screened Coulomb interaction
- $\varepsilon$ is the dielectric function

### Exchange and Correlation

The self-energy splits into exchange and correlation parts:

$$
\Sigma = \Sigma_x + \Sigma_c
$$

**Exchange self-energy** (exact exchange):
$$
\Sigma_x = -\sum_i^{occ} \frac{\phi_i(r) \phi_i^*(r')}{|r - r'|}
$$

**Correlation self-energy** (from screening):
$$
\Sigma_c(\omega) = \frac{i}{2\pi} \int d\omega' \, G(\omega + \omega') [W(\omega') - v]
$$

---

## G₀W₀: One-Shot GW

In G₀W₀, the Green's function and screened interaction are constructed from DFT orbitals and evaluated once (non-self-consistently).

### Workflow

1. **DFT/HF calculation** → orbitals $\{\psi_n\}$, energies $\{\varepsilon_n\}$
2. **Construct polarizability**: $P_0(\omega) = -i G_0 G_0$
3. **Compute dielectric function**: $\varepsilon = 1 - v P_0$
4. **Screen Coulomb interaction**: $W = \varepsilon^{-1} v$
5. **Evaluate self-energy**: $\Sigma = i G_0 W$
6. **Solve QP equation**: $\varepsilon_n^{QP}$

### Linearized QP Equation

For efficiency, the QP equation is often linearized:

$$
\varepsilon_n^{QP} \approx \varepsilon_n + Z_n \langle \psi_n | \Sigma(\varepsilon_n) - V_{xc} | \psi_n \rangle
$$

where the renormalization factor is:

$$
Z_n = \left( 1 - \frac{\partial \Sigma_c}{\partial \omega} \bigg|_{\omega=\varepsilon_n} \right)^{-1}
$$

---

## evGW: Eigenvalue Self-Consistent GW

In evGW, the quasiparticle energies are updated iteratively while keeping the orbitals fixed:

1. Start with DFT energies: $\varepsilon_n^{(0)} = \varepsilon_n^{DFT}$
2. Iterate until convergence:
   - Compute $W^{(k)}$ using $\varepsilon_n^{(k)}$
   - Evaluate $\Sigma^{(k)}$
   - Update energies: $\varepsilon_n^{(k+1)}$

This improves the description of band gaps and orbital ordering.

---

## Bethe-Salpeter Equation (BSE)

The BSE describes neutral optical excitations (excitons) including electron-hole interactions.

### BSE Hamiltonian

The BSE eigenvalue problem:

$$
\begin{pmatrix} A & B \\ -B^* & -A^* \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix}
= \Omega
\begin{pmatrix} X \\ Y \end{pmatrix}
$$

where $\Omega$ are excitation energies.

### Matrix Elements

**Resonant block**:
$$
A_{ia,jb} = (\varepsilon_a - \varepsilon_i)\delta_{ij}\delta_{ab} + \alpha K^x_{ia,jb} - K^d_{ia,jb}
$$

**Exchange kernel** (electron-hole exchange):
$$
K^x_{ia,jb} = 2 \int dr dr' \, \psi_i^*(r)\psi_a(r) \frac{1}{|r-r'|} \psi_j(r')\psi_b^*(r')
$$

**Direct kernel** (screened Coulomb attraction):
$$
K^d_{ia,jb} = \int dr dr' \, \psi_i^*(r)\psi_j(r) W(r,r') \psi_a(r')\psi_b^*(r')
$$

### Spin Channels

- **Singlet** ($\alpha = 2$): Includes exchange, optically bright
- **Triplet** ($\alpha = 0$): No exchange, optically dark

### Tamm-Dancoff Approximation (TDA)

Setting $B = 0$ gives the TDA, which only solves for $A$:

$$
A X = \Omega X
$$

This is computationally cheaper and often accurate for localized excitations.

---

## Density Fitting (RI Approximation)

QuasiX uses the Resolution of Identity (RI) approximation for efficiency:

$$
\psi_i(r)\psi_a(r) \approx \sum_P C_{ia}^P \chi_P(r)
$$

where $\{\chi_P\}$ is an auxiliary basis. This reduces four-center integrals to three-center:

$$
(ia|jb) \approx \sum_{PQ} C_{ia}^P (P|Q)^{-1} C_{jb}^Q
$$

---

## Frequency Integration

### Imaginary Frequency Formulation

The correlation self-energy is evaluated on the imaginary axis:

$$
\Sigma_c(\omega) = \frac{1}{\pi} \int_0^\infty d\omega' \, \text{Im}[W(i\omega')] \times \text{(Green's function terms)}
$$

### Analytic Continuation

Results are analytically continued to the real axis using:
- **Padé approximants** (default)
- **Two-pole model** (for difficult cases)

### Frequency Grids

- **Minimax grids**: Optimized for GW, minimal number of points
- **Gauss-Legendre**: Standard quadrature, more points needed

---

## References

1. Hedin, L. (1965). *Phys. Rev.* 139, A796. (Original GW paper)
2. Hybertsen, M. S. & Louie, S. G. (1986). *Phys. Rev. B* 34, 5390.
3. Onida, G., Reining, L. & Rubio, A. (2002). *Rev. Mod. Phys.* 74, 601.
4. van Setten, M. J. et al. (2015). *J. Chem. Theory Comput.* 11, 5665. (GW100)
5. Rohlfing, M. & Louie, S. G. (2000). *Phys. Rev. B* 62, 4927. (BSE)
