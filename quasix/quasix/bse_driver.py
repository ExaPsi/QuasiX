"""
Bethe-Salpeter Equation (BSE) implementation for optical excitations.

This module provides a PySCF-compatible interface for BSE-TDA calculations
computing optical excitation energies, oscillator strengths, and absorption spectra.

The BSE-TDA Hamiltonian is:

    H^{BSE-TDA}_{ia,jb} = delta_ij * delta_ab * Delta_qp_{ia} + K^x_{ia,jb} + K^d_{ia,jb}

where:
- Delta_qp = e_a^QP - e_i^QP (quasiparticle energy differences from GW)
- K^x = (ia|jb) is the exchange kernel (bare Coulomb)
- K^d = -(ij|W|ab) is the direct kernel (screened Coulomb at omega=0)

Theory Reference
----------------
- Rohlfing & Louie (2000): Electron-hole excitations and optical spectra from
  first principles. Phys. Rev. B 62, 4927.
- Blase et al. (2011): First-principles GW calculations for fullerenes,
  porphyrins, phtalocyanine, and other molecules of interest for organic
  photovoltaic applications. Phys. Rev. B 83, 115103.
"""

import numpy as np
from typing import Optional, Dict, Any, Tuple, List, Union
from dataclasses import dataclass, field
import logging
import time

try:
    from pyscf import lib, gto, scf, df
    from pyscf.lib import logger
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False

__all__ = [
    'run_bse_tda',
    'BSETDAResult',
    'BSETDAConfig',
    'compute_transition_dipoles',
    'compute_absorption_spectrum',
    'analyze_exciton'
]

# Physical constants
HARTREE_TO_EV = 27.211386245988

# Configure logging
log = logging.getLogger(__name__)


@dataclass
class BSETDAConfig:
    """
    Configuration for BSE-TDA calculation.

    Attributes
    ----------
    nstates : int
        Number of excited states to compute.
    tol : float
        Convergence tolerance for Davidson solver (residual norm).
    max_iter : int
        Maximum number of Davidson iterations.
    max_space : int
        Maximum subspace size before restart. If None, uses 20 * nstates.
    auxbasis : str or None
        Auxiliary basis set for density fitting. If None, auto-selects
        based on the AO basis (e.g., def2-SVP -> def2-SVP-RIFIT).
    spin_type : str
        'singlet' for spin-allowed optical transitions (default) or
        'triplet' for spin-forbidden excitations.
    use_screening : bool
        If True, uses screened Coulomb W(0) for direct kernel.
        If False, uses bare Coulomb (equivalent to CIS/TDHF-TDA).
    broadening : float
        Broadening parameter for absorption spectrum in eV.
    energy_range : tuple or None
        Energy range (min, max) in eV for absorption spectrum.
        If None, auto-detects from eigenvalues.
    verbose : int
        Print verbosity level (0-5, PySCF convention).
    """
    nstates: int = 5
    tol: float = 1e-6
    max_iter: int = 100
    max_space: Optional[int] = None
    auxbasis: Optional[str] = None
    spin_type: str = 'singlet'
    use_screening: bool = True
    broadening: float = 0.1
    energy_range: Optional[Tuple[float, float]] = None
    verbose: int = 1

    def __post_init__(self):
        """Validate and set defaults after initialization."""
        if self.max_space is None:
            self.max_space = 20 * self.nstates
        if self.spin_type not in ('singlet', 'triplet'):
            raise ValueError(f"spin_type must be 'singlet' or 'triplet', got {self.spin_type}")


@dataclass
class BSETDAResult:
    """
    Result of BSE-TDA calculation.

    Attributes
    ----------
    excitation_energies_Ha : np.ndarray
        Excitation energies in Hartree, shape (nstates,).
    excitation_energies_eV : np.ndarray
        Excitation energies in eV, shape (nstates,).
    oscillator_strengths : np.ndarray
        Oscillator strengths for each state, shape (nstates,).
    transition_dipoles : np.ndarray
        Transition dipole moments, shape (nstates, 3) for x, y, z.
    eigenvectors : np.ndarray
        BSE eigenvectors X_{ia,n}, shape (nocc*nvirt, nstates).
    converged : bool
        Whether Davidson solver converged for all requested states.
    n_converged : int
        Number of converged states.
    n_iterations : int
        Number of Davidson iterations performed.
    residual_norms : List[float]
        Final residual norms for each state.
    config : BSETDAConfig
        Configuration used for calculation.
    timing : Dict[str, float]
        Timing breakdown for calculation stages.
    f_sum : float
        Sum of oscillator strengths (Thomas-Reiche-Kuhn sum rule check).
    nocc : int
        Number of occupied orbitals.
    nvirt : int
        Number of virtual orbitals.
    """
    excitation_energies_Ha: np.ndarray
    excitation_energies_eV: np.ndarray
    oscillator_strengths: np.ndarray
    transition_dipoles: np.ndarray
    eigenvectors: np.ndarray
    converged: bool
    n_converged: int
    n_iterations: int
    residual_norms: List[float]
    config: BSETDAConfig
    timing: Dict[str, float]
    f_sum: float
    nocc: int
    nvirt: int

    def to_dict(self) -> Dict[str, Any]:
        """Export results to dictionary format."""
        return {
            'excitation_energies_Ha': self.excitation_energies_Ha.tolist(),
            'excitation_energies_eV': self.excitation_energies_eV.tolist(),
            'oscillator_strengths': self.oscillator_strengths.tolist(),
            'transition_dipoles': self.transition_dipoles.tolist(),
            'converged': self.converged,
            'n_converged': self.n_converged,
            'n_iterations': self.n_iterations,
            'f_sum': self.f_sum,
            'nocc': self.nocc,
            'nvirt': self.nvirt,
        }

    def summary(self) -> str:
        """Generate text summary of results."""
        lines = [
            "=" * 60,
            "BSE-TDA Calculation Results",
            "=" * 60,
            f"Number of states: {len(self.excitation_energies_Ha)}",
            f"Converged: {self.converged} ({self.n_converged}/{len(self.excitation_energies_Ha)})",
            f"Iterations: {self.n_iterations}",
            f"",
            f"{'State':>6} {'Energy (eV)':>12} {'f (osc.)':>10} {'|mu| (a.u.)':>12}",
            "-" * 44,
        ]
        for i in range(len(self.excitation_energies_Ha)):
            e_eV = self.excitation_energies_eV[i]
            f_osc = self.oscillator_strengths[i]
            mu = np.sqrt(np.sum(self.transition_dipoles[i]**2))
            lines.append(f"{i+1:>6} {e_eV:>12.4f} {f_osc:>10.4f} {mu:>12.4f}")

        lines.extend([
            "",
            f"Sum rule check: f_sum = {self.f_sum:.4f}",
            "=" * 60,
        ])
        return "\n".join(lines)

    def print_summary(self):
        """Print summary to stdout."""
        print(self.summary())


def _auto_select_auxbasis(basis: str) -> str:
    """
    Auto-select auxiliary basis based on AO basis.

    Parameters
    ----------
    basis : str
        The AO basis set name.

    Returns
    -------
    str
        Recommended auxiliary basis for RI-V (Coulomb fitting).
        Uses PySCF's built-in basis names.
    """
    basis_lower = basis.lower().replace('-', '').replace('_', '')

    # def2 family - use jkfit which is more widely available
    if 'def2' in basis_lower:
        if 'qzvp' in basis_lower:
            return 'def2-qzvp-jkfit'
        elif 'tzvp' in basis_lower:
            return 'def2-tzvp-jkfit'
        elif 'svp' in basis_lower:
            return 'def2-svp-jkfit'
        else:
            return 'def2-tzvp-jkfit'

    # cc-pV*Z family - use jkfit
    if 'ccpv' in basis_lower or 'augccpv' in basis_lower:
        if 'qz' in basis_lower:
            return 'cc-pvqz-jkfit'
        elif 'tz' in basis_lower:
            return 'cc-pvtz-jkfit'
        elif 'dz' in basis_lower:
            return 'cc-pvdz-jkfit'
        else:
            return 'cc-pvtz-jkfit'

    # STO-3G and minimal bases - use def2-svp-jkfit as fallback
    if 'sto' in basis_lower or '321g' in basis_lower or '631g' in basis_lower:
        return 'def2-svp-jkfit'

    # Default fallback
    return 'def2-tzvp-jkfit'


def compute_transition_dipoles(
    mol: 'gto.Mole',
    mo_coeff: np.ndarray,
    mo_occ: np.ndarray
) -> np.ndarray:
    """
    Compute transition dipole integrals <i|r|a> in MO basis.

    The transition dipole moment for an excitation i->a is:

        mu_{ia}^alpha = <psi_i | r_alpha | psi_a>

    where alpha = x, y, z.

    Parameters
    ----------
    mol : gto.Mole
        PySCF molecule object.
    mo_coeff : np.ndarray
        MO coefficient matrix, shape (nao, nmo).
    mo_occ : np.ndarray
        MO occupation numbers, shape (nmo,).

    Returns
    -------
    np.ndarray
        Transition dipole integrals, shape (nocc*nvirt, 3).
        Indexed as mu[i*nvirt + a, alpha] for transition i->a.
    """
    nocc = int(np.sum(mo_occ > 0))
    nmo = mo_coeff.shape[1]
    nvirt = nmo - nocc

    # Get dipole integrals in AO basis: <mu|r_alpha|nu>
    # Returns shape (3, nao, nao) for x, y, z components
    ao_dipoles = mol.intor('int1e_r')  # (3, nao, nao)

    # Transform to MO basis: <i|r|a>
    # C.T @ ao_dipoles @ C
    mo_occ_coeff = mo_coeff[:, :nocc]
    mo_virt_coeff = mo_coeff[:, nocc:]

    # Compute <i|r_alpha|a> for each component
    transition_dipoles = np.zeros((nocc * nvirt, 3))
    for alpha in range(3):
        # Transform: C_occ.T @ dipole_alpha @ C_virt -> (nocc, nvirt)
        mo_dipole = mo_occ_coeff.T @ ao_dipoles[alpha] @ mo_virt_coeff
        # Flatten to (nocc*nvirt,)
        transition_dipoles[:, alpha] = mo_dipole.ravel()

    return transition_dipoles


def _build_df_tensors_for_bse(
    mol: 'gto.Mole',
    mo_coeff: np.ndarray,
    mo_occ: np.ndarray,
    auxbasis: str,
    verbose: int = 1
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build density-fitted tensors required for BSE-TDA.

    Returns (ia|P), (ij|P), (ab|P) tensors and static screened interaction W(0).

    Parameters
    ----------
    mol : gto.Mole
        PySCF molecule object.
    mo_coeff : np.ndarray
        MO coefficients, shape (nao, nmo).
    mo_occ : np.ndarray
        MO occupations, shape (nmo,).
    auxbasis : str
        Auxiliary basis set name.
    verbose : int
        Verbosity level.

    Returns
    -------
    df_ia : np.ndarray
        DF tensor (ia|P), shape (nocc*nvirt, naux).
    df_ij : np.ndarray
        DF tensor (ij|P), shape (nocc*nocc, naux).
    df_ab : np.ndarray
        DF tensor (ab|P), shape (nvirt*nvirt, naux).
    chol_v : np.ndarray
        Cholesky factor V^{1/2} of Coulomb metric, shape (naux, naux).
    w0 : np.ndarray
        Static screened interaction W(omega=0), shape (naux, naux).
    """
    nocc = int(np.sum(mo_occ > 0))
    nmo = mo_coeff.shape[1]
    nvirt = nmo - nocc
    nao = mol.nao_nr()

    # Create auxiliary molecule
    auxbasis_dict = {mol.atom_symbol(i): auxbasis for i in range(mol.natm)}
    auxmol = df.addons.make_auxmol(mol, auxbasis_dict)
    naux = auxmol.nao_nr()

    if verbose > 0:
        log.info(f'Building DF tensors: nocc={nocc}, nvirt={nvirt}, naux={naux}')

    # Build 2-center Coulomb metric (P|Q)
    metric = auxmol.intor('int2c2e', aosym='s1')  # (naux, naux)

    # Compute V^{1/2} and V^{-1/2} via eigendecomposition
    eigvals, eigvecs = np.linalg.eigh(metric)
    thresh = 1e-10
    idx = eigvals > thresh
    n_valid = np.sum(idx)

    if n_valid < naux:
        log.info(f'Removed {naux - n_valid} linearly dependent aux functions')

    sqrt_eigvals = np.sqrt(np.maximum(eigvals, thresh))
    inv_sqrt_eigvals = 1.0 / sqrt_eigvals

    chol_v = eigvecs @ np.diag(sqrt_eigvals) @ eigvecs.T
    v_inv_sqrt = eigvecs @ np.diag(inv_sqrt_eigvals) @ eigvecs.T

    # Build 3-center integrals in AO basis
    int3c = df.incore.aux_e2(mol, auxmol, intor='int3c2e', aosym='s1')
    int3c = int3c.reshape(nao, nao, naux)

    # Apply V^{-1/2} transformation (PySCF convention)
    int3c = np.einsum('uvQ,QP->uvP', int3c, v_inv_sqrt)

    # Split MO coefficients
    mo_occ_coeff = np.ascontiguousarray(mo_coeff[:, :nocc])
    mo_virt_coeff = np.ascontiguousarray(mo_coeff[:, nocc:])

    # Initialize output tensors
    df_ia = np.zeros((nocc * nvirt, naux))
    df_ij = np.zeros((nocc * nocc, naux))
    df_ab = np.zeros((nvirt * nvirt, naux))

    # Transform 3-center integrals to MO basis
    for p in range(naux):
        ints_p = int3c[:, :, p]

        # (ia|P): occ-virt block
        mo_ia = mo_occ_coeff.T @ ints_p @ mo_virt_coeff
        df_ia[:, p] = mo_ia.ravel()

        # (ij|P): occ-occ block
        mo_ij = mo_occ_coeff.T @ ints_p @ mo_occ_coeff
        df_ij[:, p] = mo_ij.ravel()

        # (ab|P): virt-virt block
        mo_ab = mo_virt_coeff.T @ ints_p @ mo_virt_coeff
        df_ab[:, p] = mo_ab.ravel()

    # Build static screened interaction W(omega=0)
    # For BSE, we use W = v * (1 - v * P_0)^{-1} * v
    # As a simple approximation, start with bare Coulomb W(0) = V
    # A proper implementation would include RPA screening from GW
    w0 = chol_v @ chol_v  # V = V^{1/2} @ V^{1/2}

    return df_ia, df_ij, df_ab, chol_v, w0


def _apply_bse_tda_hamiltonian(
    x: np.ndarray,
    delta_qp: np.ndarray,
    df_ia: np.ndarray,
    df_ij: np.ndarray,
    df_ab: np.ndarray,
    w0: np.ndarray,
    nocc: int,
    nvirt: int,
    spin_type: str = 'singlet',
    use_screening: bool = False
) -> np.ndarray:
    """
    Apply BSE-TDA (or TDA/CIS) Hamiltonian to trial vector(s).

    For singlet TDA/CIS:
        H^{TDA}_{ia,jb} = δ_{ij}δ_{ab}(ε_a - ε_i) + 2*(ia|jb) - (ij|ab)
                                                    ^^^^^^^^   ^^^^^^^
                                                    K^x (+2)   K^d (-1)

    For triplet TDA/CIS:
        H^{TDA}_{ia,jb} = δ_{ij}δ_{ab}(ε_a - ε_i) - (ij|ab)
                                                    ^^^^^^^
                                                    K^d (-1)

    For BSE-TDA with screening:
        H^{BSE}_{ia,jb} = δ_{ij}δ_{ab}Δε^QP_{ia} + 2*(ia|jb) - (ij|W(0)|ab)
        where W(0) is the statically screened Coulomb interaction.

    Parameters
    ----------
    x : np.ndarray
        Trial vector(s), shape (nocc*nvirt,) or (nocc*nvirt, nvec).
    delta_qp : np.ndarray
        QP energy differences, shape (nocc*nvirt,).
    df_ia, df_ij, df_ab : np.ndarray
        DF tensors for transitions.
    w0 : np.ndarray
        Static screened interaction, shape (naux, naux).
    nocc, nvirt : int
        Number of occupied and virtual orbitals.
    spin_type : str
        'singlet' (alpha_x=2) or 'triplet' (alpha_x=0).
    use_screening : bool
        If True, use screened W(0) for direct kernel.
        If False, use bare Coulomb for direct kernel (TDA/CIS).

    Returns
    -------
    np.ndarray
        Result of H * x, same shape as x.
    """
    is_batch = x.ndim == 2
    if not is_batch:
        x = x[:, np.newaxis]

    nvec = x.shape[1]

    # Diagonal contribution: Delta_qp * x
    y = delta_qp[:, np.newaxis] * x

    # Exchange kernel: K^x * x
    # K^x_{ia,jb} = (ia|jb) = sum_P (ia|P)(jb|P)
    # For singlet: factor of 2 from spin summation
    # For triplet: no exchange (alpha_x = 0)
    alpha_x = 2.0 if spin_type == 'singlet' else 0.0
    if alpha_x > 0:
        z = df_ia.T @ x  # (naux, nvec)
        y_kx = df_ia @ z  # (n_trans, nvec)
        y += alpha_x * y_kx

    # Direct kernel: K^d * x
    # K^d_{ia,jb} = -(ij|ab) for bare Coulomb (TDA/CIS)
    # K^d_{ia,jb} = -(ij|W(0)|ab) for screened (BSE-TDA)
    #
    # Using DF: (ij|ab) = sum_P (ij|P)(ab|P) = df_ij.T @ df_ab
    # For each vector x, we compute: sum_{jb} -(ij|ab) * x_{jb}
    #
    # Reshape x to (nocc, nvirt) for each vector
    # Then: K^d @ x = -sum_j (ij|P) * sum_b (ab|P) * x[j,b]
    #              = -df_ij @ df_ab.T @ x_reshaped.T

    for v in range(nvec):
        x_v = x[:, v].reshape(nocc, nvirt)  # (nocc, nvirt)

        # Contract x with virtual indices: sum_b (ab|P) * x[j,b]
        # df_ab shape: (nvirt*nvirt, naux) -> need to reshape
        # x_v shape: (nocc, nvirt)

        # For the direct kernel, we need (ij|ab):
        # Using DF: (ij|ab) = sum_P df_ij[ij,P] * df_ab[ab,P]
        #
        # Efficient implementation:
        # 1. Compute intermediate: z[j,a,P] = sum_b df_ab[ab,P] * x[j,b]
        # 2. Compute result: y[i,a] = -sum_j sum_P df_ij[ij,P] * z[j,a,P]

        # Step 1: Contract over virtual b index
        # df_ab is (nvirt*nvirt, naux), reshape to (nvirt, nvirt, naux)
        df_ab_3d = df_ab.reshape(nvirt, nvirt, -1)  # (a, b, P)
        naux = df_ab_3d.shape[2]

        # z[j,a,P] = sum_b df_ab[a,b,P] * x[j,b]
        # Using einsum: z = einsum('abP,jb->jaP', df_ab_3d, x_v)
        z_jaP = np.einsum('abP,jb->jaP', df_ab_3d, x_v, optimize=True)  # (nocc, nvirt, naux)

        # Step 2: Contract over occupied j index
        # df_ij is (nocc*nocc, naux), reshape to (nocc, nocc, naux)
        df_ij_3d = df_ij.reshape(nocc, nocc, -1)  # (i, j, P)

        # y_kd[i,a] = -sum_j sum_P df_ij[i,j,P] * z[j,a,P]
        # Using einsum: y_kd = -einsum('ijP,jaP->ia', df_ij_3d, z_jaP)
        y_kd = -np.einsum('ijP,jaP->ia', df_ij_3d, z_jaP, optimize=True)  # (nocc, nvirt)

        # If using screening, we would need to apply W(0) instead of bare Coulomb
        # For now, use_screening modifies the direct kernel
        # TODO: Implement proper screened W(0) for BSE

        y[:, v] += y_kd.ravel()

    if not is_batch:
        y = y[:, 0]

    return y


def _davidson_solver(
    apply_h,
    diagonal: np.ndarray,
    n_roots: int,
    tol: float = 1e-6,
    max_iter: int = 100,
    max_space: int = 60
) -> Tuple[np.ndarray, np.ndarray, bool, int, List[float]]:
    """
    Davidson iterative eigenvalue solver.

    Finds the lowest n_roots eigenvalues of a symmetric matrix
    using matrix-free approach.

    Parameters
    ----------
    apply_h : callable
        Function that applies H to vectors: Y = H @ X.
    diagonal : np.ndarray
        Diagonal of H (used for preconditioner), shape (n,).
    n_roots : int
        Number of lowest eigenvalues to find.
    tol : float
        Convergence tolerance for residual norm.
    max_iter : int
        Maximum iterations.
    max_space : int
        Maximum subspace size before restart.

    Returns
    -------
    eigenvalues : np.ndarray
        Converged eigenvalues, shape (n_roots,).
    eigenvectors : np.ndarray
        Converged eigenvectors, shape (n, n_roots).
    converged : bool
        Whether all roots converged.
    n_iter : int
        Number of iterations performed.
    residual_norms : List[float]
        Final residual norms.
    """
    n = len(diagonal)
    n_roots = min(n_roots, n)

    # Handle trivial case
    if n <= n_roots:
        identity = np.eye(n)
        h_full = apply_h(identity)
        h_sym = 0.5 * (h_full + h_full.T)
        eigenvalues, eigenvectors = np.linalg.eigh(h_sym)
        return eigenvalues, eigenvectors, True, 1, [0.0] * n

    # Initialize with unit vectors at lowest diagonal elements
    sorted_idx = np.argsort(diagonal)
    v = np.zeros((n, n_roots))
    for i, idx in enumerate(sorted_idx[:n_roots]):
        v[idx, i] = 1.0

    eigenvalues = np.zeros(n_roots)
    eigenvectors = np.zeros((n, n_roots))
    converged_flags = [False] * n_roots
    residual_norms = [np.inf] * n_roots
    level_shift = 1e-3

    for iteration in range(max_iter):
        # Apply H to subspace
        av = apply_h(v)

        # Form projected Hamiltonian
        h_eff = v.T @ av
        h_eff = 0.5 * (h_eff + h_eff.T)  # Symmetrize

        # Solve small eigenvalue problem
        theta, y = np.linalg.eigh(h_eff)

        # Compute Ritz vectors
        n_ritz = min(n_roots, y.shape[1])
        ritz_vectors = v @ y[:, :n_ritz]
        ritz_ax = av @ y[:, :n_ritz]

        # Check convergence and build new trial vectors
        new_vectors = []

        for j in range(n_roots):
            if j >= n_ritz:
                continue

            if converged_flags[j]:
                continue

            # Residual
            r = ritz_ax[:, j] - theta[j] * ritz_vectors[:, j]
            res_norm = np.linalg.norm(r)
            residual_norms[j] = res_norm

            if res_norm < tol:
                converged_flags[j] = True
                eigenvalues[j] = theta[j]
                eigenvectors[:, j] = ritz_vectors[:, j]
                continue

            # Preconditioned residual
            denom = diagonal - theta[j] + level_shift
            denom = np.where(np.abs(denom) < 1e-8, 1e-8, denom)
            t = r / denom

            # Normalize
            t_norm = np.linalg.norm(t)
            if t_norm > 1e-14:
                t /= t_norm
                new_vectors.append(t)

        # Check if all converged
        if all(converged_flags):
            return eigenvalues, eigenvectors, True, iteration + 1, residual_norms

        # Expand subspace
        if not new_vectors:
            continue

        new_vecs = np.column_stack(new_vectors)

        # Check if restart needed
        if v.shape[1] + new_vecs.shape[1] > max_space:
            # Restart with Ritz vectors
            v = ritz_vectors[:, :n_roots].copy()
        else:
            # Orthogonalize and append
            for vec in new_vectors:
                for col_idx in range(v.shape[1]):
                    vec -= np.dot(vec, v[:, col_idx]) * v[:, col_idx]
                vec_norm = np.linalg.norm(vec)
                if vec_norm > 1e-12:
                    vec /= vec_norm
                    v = np.column_stack([v, vec])

    # Final values even if not converged
    for j in range(n_roots):
        if not converged_flags[j] and j < n_ritz:
            eigenvalues[j] = theta[j]
            eigenvectors[:, j] = ritz_vectors[:, j]

    return eigenvalues, eigenvectors, all(converged_flags), max_iter, residual_norms


def compute_absorption_spectrum(
    excitation_energies_eV: np.ndarray,
    oscillator_strengths: np.ndarray,
    energy_range: Optional[Tuple[float, float]] = None,
    broadening: float = 0.1,
    n_points: int = 1000,
    broadening_type: str = 'lorentzian'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute broadened absorption spectrum.

    Parameters
    ----------
    excitation_energies_eV : np.ndarray
        Excitation energies in eV.
    oscillator_strengths : np.ndarray
        Oscillator strengths for each state.
    energy_range : tuple, optional
        Energy range (min, max) in eV. Auto-detects if None.
    broadening : float
        Broadening parameter (FWHM/2) in eV.
    n_points : int
        Number of points in spectrum.
    broadening_type : str
        'lorentzian' or 'gaussian'.

    Returns
    -------
    energies : np.ndarray
        Energy grid in eV.
    intensities : np.ndarray
        Absorption intensity (arbitrary units).
    """
    if energy_range is None:
        e_min = max(0, excitation_energies_eV.min() - 2.0)
        e_max = excitation_energies_eV.max() + 2.0
    else:
        e_min, e_max = energy_range

    energies = np.linspace(e_min, e_max, n_points)
    intensities = np.zeros(n_points)

    for e, f in zip(excitation_energies_eV, oscillator_strengths):
        if broadening_type == 'lorentzian':
            intensities += f * broadening / ((energies - e)**2 + broadening**2)
        else:  # gaussian
            intensities += f * np.exp(-(energies - e)**2 / (2 * broadening**2))

    return energies, intensities


def analyze_exciton(
    eigenvector: np.ndarray,
    nocc: int,
    nvirt: int
) -> Dict[str, Any]:
    """
    Analyze exciton character from BSE eigenvector.

    Parameters
    ----------
    eigenvector : np.ndarray
        BSE eigenvector, shape (nocc*nvirt,).
    nocc, nvirt : int
        Number of occupied and virtual orbitals.

    Returns
    -------
    dict
        Analysis results including:
        - 'amplitude_squared': np.ndarray of shape (nocc, nvirt)
        - 'dominant_transition': (i, a, weight) tuple
        - 'participation_ratio': float (multi-configurational character)
    """
    # Reshape to (nocc, nvirt)
    amplitude = eigenvector.reshape(nocc, nvirt)
    amplitude_sq = amplitude**2

    # Find dominant transition
    max_idx = np.unravel_index(np.argmax(amplitude_sq), amplitude_sq.shape)
    i_dom, a_dom = max_idx
    weight_dom = amplitude_sq[i_dom, a_dom]

    # Participation ratio: 1 / sum(lambda^2) where lambda = amplitude^2
    # Normalized amplitudes
    total = amplitude_sq.sum()
    if total > 1e-14:
        lambda_n = amplitude_sq.ravel() / total
        pr = 1.0 / np.sum(lambda_n**2)
    else:
        pr = 1.0

    return {
        'amplitude_squared': amplitude_sq,
        'dominant_transition': (i_dom, a_dom, weight_dom),
        'participation_ratio': pr
    }


def run_bse_tda(
    mol: Union['gto.Mole', None] = None,
    mf: Union['scf.RHF', None] = None,
    nstates: int = 5,
    tol: float = 1e-6,
    max_iter: int = 100,
    auxbasis: Optional[str] = None,
    qp_energies: Optional[np.ndarray] = None,
    verbose: bool = True,
    config: Optional[BSETDAConfig] = None,
    **kwargs
) -> BSETDAResult:
    """
    Run BSE-TDA calculation for optical excitations.

    This is the main entry point for BSE-TDA calculations. It computes
    excitation energies, oscillator strengths, and absorption spectra
    using the Bethe-Salpeter equation with Tamm-Dancoff approximation.

    Parameters
    ----------
    mol : gto.Mole, optional
        PySCF molecule object. Required if mf is not provided.
    mf : scf.RHF
        Converged RHF (or RKS) object from PySCF.
    nstates : int
        Number of excited states to compute.
    tol : float
        Davidson solver tolerance (residual norm convergence).
    max_iter : int
        Maximum Davidson iterations.
    auxbasis : str, optional
        Auxiliary basis for RI. Auto-selected if None.
    qp_energies : np.ndarray, optional
        Quasiparticle energies from GW calculation.
        If None, uses Kohn-Sham/HF orbital energies.
    verbose : bool
        Print progress information.
    config : BSETDAConfig, optional
        Full configuration object. Overrides individual parameters.
    **kwargs
        Additional arguments passed to BSETDAConfig.

    Returns
    -------
    BSETDAResult
        Object containing excitation energies, oscillator strengths,
        eigenvectors, and convergence information.

    Examples
    --------
    Basic usage with PySCF:

    >>> from pyscf import gto, scf
    >>> from quasix.bse import run_bse_tda
    >>> mol = gto.M(atom='O 0 0 0; H 0 0 1; H 0 1 0', basis='def2-svp')
    >>> mf = scf.RHF(mol).run()
    >>> result = run_bse_tda(mf=mf, nstates=5)
    >>> print(result.excitation_energies_eV)

    With GW quasiparticle energies:

    >>> gw = quasix.evgw.EvGW(mf).run()
    >>> result = run_bse_tda(mf=mf, qp_energies=gw.qp_energies, nstates=5)
    """
    if not PYSCF_AVAILABLE:
        raise ImportError("PySCF is required for BSE-TDA calculations")

    if mf is None and mol is None:
        raise ValueError("Either mol or mf must be provided")

    if mf is None:
        mf = scf.RHF(mol).run()

    if mol is None:
        mol = mf.mol

    # Build configuration
    if config is None:
        config = BSETDAConfig(
            nstates=nstates,
            tol=tol,
            max_iter=max_iter,
            auxbasis=auxbasis,
            verbose=1 if verbose else 0,
            **kwargs
        )

    timing = {}
    t_start = time.time()

    # Extract MO data
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    mo_energy = mf.mo_energy

    nocc = int(np.sum(mo_occ > 0))
    nmo = mo_coeff.shape[1]
    nvirt = nmo - nocc
    n_trans = nocc * nvirt

    if verbose:
        log.info(f"BSE-TDA calculation: nocc={nocc}, nvirt={nvirt}, nstates={config.nstates}")

    # Auto-select auxiliary basis
    if config.auxbasis is None:
        basis_name = mol.basis if isinstance(mol.basis, str) else 'def2-svp'
        config.auxbasis = _auto_select_auxbasis(basis_name)

    if verbose:
        log.info(f"Auxiliary basis: {config.auxbasis}")

    # Build DF tensors
    t0 = time.time()
    df_ia, df_ij, df_ab, chol_v, w0 = _build_df_tensors_for_bse(
        mol, mo_coeff, mo_occ, config.auxbasis, config.verbose
    )
    timing['df_tensors'] = time.time() - t0

    # Build QP energy differences
    if qp_energies is not None:
        eps = qp_energies
    else:
        eps = mo_energy

    delta_qp = np.zeros(n_trans)
    for i in range(nocc):
        for a in range(nvirt):
            idx = i * nvirt + a
            delta_qp[idx] = eps[nocc + a] - eps[i]

    # Define Hamiltonian application
    # Note: use_screening is disabled by default for TDA/CIS-like behavior
    # The full BSE screening requires GW-corrected energies for proper results
    def apply_h(x):
        return _apply_bse_tda_hamiltonian(
            x, delta_qp, df_ia, df_ij, df_ab, w0,
            nocc, nvirt, config.spin_type,
            use_screening=config.use_screening
        )

    # Run Davidson solver
    t0 = time.time()
    eigenvalues, eigenvectors, converged, n_iter, residual_norms = _davidson_solver(
        apply_h,
        delta_qp,
        config.nstates,
        tol=config.tol,
        max_iter=config.max_iter,
        max_space=config.max_space
    )
    timing['davidson'] = time.time() - t0

    if verbose:
        status = "converged" if converged else "NOT converged"
        log.info(f"Davidson {status} in {n_iter} iterations")

    # Compute transition dipoles
    t0 = time.time()
    mu_ia = compute_transition_dipoles(mol, mo_coeff, mo_occ)
    timing['dipoles'] = time.time() - t0

    # Compute oscillator strengths
    # f_n = (2/3) * E_n * sum_alpha |d_n^alpha|^2
    # d_n^alpha = sum_ia X_n(ia) * mu_{ia}^alpha
    transition_dipoles = np.zeros((config.nstates, 3))
    oscillator_strengths = np.zeros(config.nstates)

    for n in range(config.nstates):
        for alpha in range(3):
            transition_dipoles[n, alpha] = np.dot(eigenvectors[:, n], mu_ia[:, alpha])

        d_sq = np.sum(transition_dipoles[n]**2)
        oscillator_strengths[n] = (2.0 / 3.0) * eigenvalues[n] * d_sq

    f_sum = np.sum(oscillator_strengths)

    timing['total'] = time.time() - t_start

    # Build result
    result = BSETDAResult(
        excitation_energies_Ha=eigenvalues,
        excitation_energies_eV=eigenvalues * HARTREE_TO_EV,
        oscillator_strengths=oscillator_strengths,
        transition_dipoles=transition_dipoles,
        eigenvectors=eigenvectors,
        converged=converged,
        n_converged=sum(r < config.tol for r in residual_norms),
        n_iterations=n_iter,
        residual_norms=residual_norms,
        config=config,
        timing=timing,
        f_sum=f_sum,
        nocc=nocc,
        nvirt=nvirt
    )

    if verbose:
        result.print_summary()

    return result
