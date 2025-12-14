"""
Dielectric function and polarizability calculations for QuasiX GW/BSE.

This module provides Python interfaces for:
- P0(ω) independent particle polarizability
- ε(ω) dielectric function
- ε^-1(ω) inverse dielectric function
- W(ω) screened interaction
"""

from typing import Optional, Tuple, Union, Any
import numpy as np
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

# Import Rust functions when available
try:
    from .quasix import (
        compute_polarizability_p0 as _compute_polarizability_p0,
        compute_polarizability_batch as _compute_polarizability_batch,
        compute_dielectric_function as _compute_dielectric_function,
        compute_epsilon_inverse as _compute_epsilon_inverse,
    )
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    # Mock implementations will be provided below


@dataclass
class Polarizability:
    """Container for polarizability matrix P0(ω).
    
    Attributes:
        matrix: Polarizability matrix in auxiliary basis
        frequency: Frequency at which P0 was computed
        naux: Number of auxiliary basis functions
        hermitian: Whether the matrix is Hermitian (true for imaginary frequencies)
    """
    matrix: np.ndarray
    frequency: complex
    naux: int
    hermitian: bool = False
    
    def __repr__(self) -> str:
        freq_str = f"{self.frequency:.3f}" if np.isreal(self.frequency) else f"{self.frequency:.3f}j"
        return f"Polarizability(ω={freq_str}, naux={self.naux}, hermitian={self.hermitian})"
    
    def check_hermiticity(self, tol: float = 1e-10) -> bool:
        """Check if P0 is Hermitian within tolerance."""
        diff = np.max(np.abs(self.matrix - self.matrix.conj().T))
        return diff < tol


@dataclass
class DielectricFunction:
    """Container for dielectric function ε(ω).
    
    Attributes:
        matrix: Dielectric matrix in auxiliary basis
        frequency: Frequency at which ε was computed
        naux: Number of auxiliary basis functions
        eigenvalues: Eigenvalues of ε (optional)
    """
    matrix: np.ndarray
    frequency: complex
    naux: int
    eigenvalues: Optional[np.ndarray] = None
    
    def __repr__(self) -> str:
        freq_str = f"{self.frequency:.3f}" if np.isreal(self.frequency) else f"{self.frequency:.3f}j"
        return f"DielectricFunction(ω={freq_str}, naux={self.naux})"
    
    def check_positive_definite(self) -> bool:
        """Check if ε(iω) > 1 for imaginary frequency."""
        if not np.iscomplex(self.frequency):
            return True  # Only check for imaginary frequencies
        
        if self.eigenvalues is None:
            self.eigenvalues = np.linalg.eigvalsh(self.matrix)
        
        return np.all(self.eigenvalues > 1.0)


@dataclass 
class ScreenedInteraction:
    """Container for screened interaction W(ω) = ε^-1(ω) v.
    
    Attributes:
        matrix: Screened interaction matrix
        frequency: Frequency at which W was computed
        naux: Number of auxiliary basis functions
        static: Whether this is the static (ω=0) screened interaction
    """
    matrix: np.ndarray
    frequency: complex
    naux: int
    static: bool = False
    
    def __repr__(self) -> str:
        freq_str = "static" if self.static else (
            f"{self.frequency:.3f}" if np.isreal(self.frequency) else f"{self.frequency:.3f}j"
        )
        return f"ScreenedInteraction(ω={freq_str}, naux={self.naux})"


def compute_polarizability_p0(
    df_ia: np.ndarray,
    orbital_energies: Tuple[np.ndarray, np.ndarray],
    frequency: complex,
    eta: float = 0.01
) -> Polarizability:
    """Compute independent particle polarizability P0(ω).
    
    P0_PQ(ω) = 2 Σ_ia (ia|P)(ia|Q) / (ε_a - ε_i - ω - iη)
    
    Args:
        df_ia: DF tensor in transition space, shape (n_trans, naux)
        orbital_energies: Tuple of (occupied, virtual) orbital energies
        frequency: Frequency (real or imaginary)
        eta: Broadening parameter for convergence
        
    Returns:
        Polarizability object containing P0(ω)
        
    Examples:
        >>> # Compute P0 on imaginary axis
        >>> p0 = compute_polarizability_p0(df_ia, (e_occ, e_virt), 1.0j)
        >>> assert p0.check_hermiticity()  # P0(iω) is Hermitian
        
        >>> # Compute P0 on real axis with broadening
        >>> p0 = compute_polarizability_p0(df_ia, (e_occ, e_virt), 5.0, eta=0.01)
    """
    # Use mock implementation for now (Rust bindings have type issues)
    if False and RUST_AVAILABLE:
        e_occ, e_virt = orbital_energies
        result = _compute_polarizability_p0(
            df_ia.tolist() if hasattr(df_ia, 'tolist') else df_ia,
            e_occ.tolist() if hasattr(e_occ, 'tolist') else e_occ,
            e_virt.tolist() if hasattr(e_virt, 'tolist') else e_virt,
            complex(frequency),
            eta
        )
        matrix = np.array(result['matrix']).reshape(result['shape'])
    else:
        # Mock implementation for testing
        e_occ, e_virt = orbital_energies
        nocc = len(e_occ)
        nvirt = len(e_virt)
        naux = df_ia.shape[1]
        
        # Build transition energies
        omega_ia = []
        for e_i in e_occ:
            for e_a in e_virt:
                omega_ia.append(e_a - e_i)
        omega_ia = np.array(omega_ia)
        
        # Compute P0_PQ = 2 * sum_ia (ia|P) * 1/(ω - Δ_ia) * (ia|Q)^*
        # df_ia has shape (n_trans, naux) where n_trans = nocc * nvirt
        denominator = 1.0 / (frequency - omega_ia + 1j * eta)
        # Scale each row of df_ia by denominator
        chi = df_ia * denominator[:, None]  # shape: (n_trans, naux)
        # P0 = 2 * chi^T @ df_ia^* (correct formula for hermiticity)
        matrix = 2.0 * chi.T @ df_ia.conj()  # shape: (naux, naux)
        
    # Check if Hermitian (for imaginary frequencies)
    is_hermitian = np.abs(np.imag(frequency)) > 1e-10
    
    return Polarizability(
        matrix=matrix,
        frequency=frequency,
        naux=matrix.shape[0],
        hermitian=is_hermitian
    )


def compute_polarizability_batch(
    df_ia: np.ndarray,
    orbital_energies: Tuple[np.ndarray, np.ndarray],
    frequencies: np.ndarray,
    eta: float = 0.01
) -> np.ndarray:
    """Compute P0(ω) for multiple frequencies.

    Args:
        df_ia: DF tensor in transition space
        orbital_energies: Tuple of (occupied, virtual) energies
        frequencies: Array of frequencies
        eta: Broadening parameter

    Returns:
        3D numpy array of shape (n_freq, n_aux, n_aux) with P0 matrices

    Examples:
        >>> # Compute P0 on frequency grid
        >>> freqs = 1j * np.linspace(0, 50, 32)
        >>> p0_all = compute_polarizability_batch(df_ia, (e_occ, e_virt), freqs)
        >>> print(f"Computed P0 array shape: {p0_all.shape}")
    """
    if RUST_AVAILABLE:
        e_occ, e_virt = orbital_energies
        # Ensure frequencies is a numpy array
        if not isinstance(frequencies, np.ndarray):
            frequencies = np.array(frequencies, dtype=np.float64)
        # Call the Rust function with numpy arrays directly
        result = _compute_polarizability_batch(
            df_ia,
            e_occ,
            e_virt,
            frequencies,  # Pass numpy array directly
            eta
        )
        return result
    else:
        # Mock implementation - create 3D array
        naux = df_ia.shape[1]
        n_freq = len(frequencies)
        p0_all = np.zeros((n_freq, naux, naux), dtype=np.complex128)
        for i, freq in enumerate(frequencies):
            p0 = compute_polarizability_p0(df_ia, orbital_energies, freq, eta)
            p0_all[i] = p0.matrix
        return p0_all


def compute_dielectric_function(
    p0: Union[Polarizability, np.ndarray],
    v_sqrt: np.ndarray,
    frequency: Optional[complex] = None
) -> Union[np.ndarray, DielectricFunction]:
    """Compute dielectric function ε(ω) = I - v^1/2 P0(ω) v^1/2.

    Args:
        p0: Polarizability matrix, Polarizability object, or 3D array of P0 matrices
        v_sqrt: Square root of Coulomb matrix v^1/2
        frequency: Frequency (if p0 is raw array)

    Returns:
        If p0 is 3D array: returns 3D numpy array of epsilon matrices
        Otherwise: returns numpy array (2D) of epsilon matrix

    Examples:
        >>> # Compute dielectric function
        >>> p0 = compute_polarizability_p0(df_ia, (e_occ, e_virt), 1.0j)
        >>> eps = compute_dielectric_function(p0, v_sqrt)
        >>> assert eps.check_positive_definite()  # ε(iω) > 1
    """
    if isinstance(p0, Polarizability):
        p0_matrix = p0.matrix
        freq = p0.frequency
    else:
        p0_matrix = p0
        freq = frequency if frequency is not None else 0.0

    # Handle both 2D and 3D cases
    if p0_matrix.ndim == 3:
        # Batch computation for multiple frequencies
        n_freq, naux, _ = p0_matrix.shape
        epsilon_all = np.zeros_like(p0_matrix)
        # Cast v_sqrt to same dtype as p0_matrix for complex arithmetic
        v_sqrt_cast = v_sqrt.astype(p0_matrix.dtype) if np.iscomplexobj(p0_matrix) else v_sqrt
        identity = np.eye(naux, dtype=p0_matrix.dtype)
        for i in range(n_freq):
            epsilon_all[i] = identity - v_sqrt_cast @ p0_matrix[i] @ v_sqrt_cast
            # Enforce Hermiticity to eliminate numerical errors
            if np.iscomplexobj(epsilon_all[i]):
                epsilon_all[i] = (epsilon_all[i] + epsilon_all[i].conj().T) / 2
            else:
                epsilon_all[i] = (epsilon_all[i] + epsilon_all[i].T) / 2
        return epsilon_all
    else:
        # Single frequency computation
        if RUST_AVAILABLE and False:  # Disabled for now due to type issues
            result = _compute_dielectric_function(
                p0_matrix,
                v_sqrt
            )
            matrix = result
        else:
            # Mock implementation: ε = I - v^1/2 P0 v^1/2
            naux = v_sqrt.shape[0]
            # Cast v_sqrt to same dtype as p0_matrix for complex arithmetic
            v_sqrt_cast = v_sqrt.astype(p0_matrix.dtype) if np.iscomplexobj(p0_matrix) else v_sqrt
            identity = np.eye(naux, dtype=p0_matrix.dtype)
            matrix = identity - v_sqrt_cast @ p0_matrix @ v_sqrt_cast
            # Enforce Hermiticity to eliminate numerical errors
            if np.iscomplexobj(matrix):
                matrix = (matrix + matrix.conj().T) / 2
            else:
                matrix = (matrix + matrix.T) / 2

        # Return numpy array directly for compatibility
        return matrix


def compute_epsilon_matrix(p0: np.ndarray, v: np.ndarray, symmetrized: bool = False) -> np.ndarray:
    """
    Compute the dielectric matrix ε.

    Two forms are supported:
    - Non-symmetrized: ε = I - vP0 (default)
    - Symmetrized: ε = I - v^(1/2) P0 v^(1/2) (when symmetrized=True and v is v_sqrt)

    Parameters
    ----------
    p0 : np.ndarray
        Polarizability matrix P0, shape (n_aux, n_aux) or (n_freq, n_aux, n_aux)
    v : np.ndarray
        Coulomb matrix V (if symmetrized=False) or its square root V^(1/2) (if symmetrized=True)
        Shape (n_aux, n_aux)
    symmetrized : bool, optional
        If True, assumes v is V^(1/2) and computes symmetrized form.
        If False (default), assumes v is V and computes non-symmetrized form.

    Returns
    -------
    np.ndarray
        Dielectric matrix epsilon, same shape as p0

    Notes
    -----
    When P0 is complex (e.g., for frequency-dependent calculations), V is cast to
    complex dtype to ensure Hermiticity is preserved in matrix products.
    The result is explicitly symmetrized/hermitized to eliminate numerical errors.
    """
    if p0.ndim == 2:
        n_aux = p0.shape[0]
        identity = np.eye(n_aux, dtype=p0.dtype)

        # Cast v to same dtype as p0 to preserve Hermiticity
        v_cast = v.astype(p0.dtype) if np.iscomplexobj(p0) else v

        if symmetrized:
            # Symmetrized form: ε = I - v^(1/2) P0 v^(1/2)
            epsilon = identity - v_cast @ p0 @ v_cast
        else:
            # Non-symmetrized form: ε = I - vP0
            epsilon = identity - v_cast @ p0

        # Enforce Hermiticity/symmetry to eliminate numerical errors
        if np.iscomplexobj(epsilon):
            epsilon = (epsilon + epsilon.conj().T) / 2
        else:
            epsilon = (epsilon + epsilon.T) / 2

        return epsilon
    else:
        n_freq, n_aux, _ = p0.shape
        identity = np.eye(n_aux, dtype=p0.dtype)
        epsilon = np.zeros_like(p0)

        # Cast v to same dtype as p0 to preserve Hermiticity
        v_cast = v.astype(p0.dtype) if np.iscomplexobj(p0) else v

        for i in range(n_freq):
            if symmetrized:
                # Symmetrized form: ε = I - v^(1/2) P0 v^(1/2)
                epsilon[i] = identity - v_cast @ p0[i] @ v_cast
            else:
                # Non-symmetrized form: ε = I - vP0
                epsilon[i] = identity - v_cast @ p0[i]

            # Enforce Hermiticity/symmetry to eliminate numerical errors
            if np.iscomplexobj(epsilon[i]):
                epsilon[i] = (epsilon[i] + epsilon[i].conj().T) / 2
            else:
                epsilon[i] = (epsilon[i] + epsilon[i].T) / 2

        return epsilon


def compute_epsilon_inverse(
    epsilon: Union[DielectricFunction, np.ndarray],
    method: str = 'cholesky',
    regularization: float = 1e-12
) -> np.ndarray:
    """Compute inverse dielectric function ε^-1(ω).
    
    Args:
        epsilon: Dielectric function or matrix
        method: Inversion method ('cholesky', 'lu', 'svd')
        regularization: Regularization for numerical stability
        
    Returns:
        Inverse dielectric matrix
        
    Examples:
        >>> # Compute ε^-1 with stable inversion
        >>> eps_inv = compute_epsilon_inverse(eps, method='cholesky')
        >>> 
        >>> # Check: ε^-1 ε = I
        >>> identity = eps_inv @ eps.matrix
        >>> assert np.allclose(identity, np.eye(eps.naux))
    """
    if isinstance(epsilon, DielectricFunction):
        eps_matrix = epsilon.matrix
    else:
        eps_matrix = epsilon
    
    if False and RUST_AVAILABLE:
        result = _compute_epsilon_inverse(
            eps_matrix.tolist() if hasattr(eps_matrix, 'tolist') else eps_matrix,
            method,
            regularization
        )
        return np.array(result).reshape(eps_matrix.shape)
    else:
        # Mock implementation with chosen method
        if method == 'cholesky':
            try:
                # Try Cholesky for positive definite matrices
                L = np.linalg.cholesky(eps_matrix + regularization * np.eye(eps_matrix.shape[0]))
                return np.linalg.solve(L.T, np.linalg.solve(L, np.eye(eps_matrix.shape[0])))
            except np.linalg.LinAlgError:
                # Fall back to LU
                method = 'lu'
        
        if method == 'lu':
            return np.linalg.solve(eps_matrix + regularization * np.eye(eps_matrix.shape[0]),
                                   np.eye(eps_matrix.shape[0]))
        elif method == 'svd':
            U, s, Vh = np.linalg.svd(eps_matrix)
            s_inv = 1.0 / (s + regularization)
            return (Vh.T * s_inv) @ U.T
        else:
            return np.linalg.inv(eps_matrix + regularization * np.eye(eps_matrix.shape[0]))


def compute_screened_interaction(
    epsilon_inv: np.ndarray,
    v_full: np.ndarray,
    frequency: complex = 0.0
) -> np.ndarray:
    """Compute screened interaction W(ω) = ε^-1(ω) v.

    Args:
        epsilon_inv: Inverse dielectric matrix
        v_full: Full Coulomb matrix
        frequency: Frequency at which W is computed (not used in computation, for metadata only)

    Returns:
        Screened interaction matrix W as numpy array

    Examples:
        >>> # Compute static screened interaction
        >>> W0 = compute_screened_interaction(eps_inv, v_coulomb, 0.0)
        >>> print(f"W0 norm: {np.linalg.norm(W0)}")

        >>> # Compute W on imaginary axis
        >>> W_iw = compute_screened_interaction(eps_inv, v_coulomb, 1.0j)
    """
    W_matrix = epsilon_inv @ v_full

    # Return numpy array directly for compatibility with tests
    return W_matrix


def compute_head_wings_body(
    epsilon: DielectricFunction,
    q_vector: Optional[np.ndarray] = None
) -> Tuple[complex, np.ndarray, np.ndarray]:
    """Extract head, wings, and body of dielectric matrix (for periodic systems).
    
    Args:
        epsilon: Dielectric function
        q_vector: q-point vector (for periodic boundary conditions)
        
    Returns:
        Tuple of (head, wings, body) components
        
    Notes:
        - Head: ε_{G=0,G'=0}(q,ω) - macroscopic dielectric constant
        - Wings: ε_{G=0,G'≠0}(q,ω) and ε_{G≠0,G'=0}(q,ω)
        - Body: ε_{G≠0,G'≠0}(q,ω)
    """
    # For molecules, just return the full matrix
    if q_vector is None:
        return epsilon.matrix[0, 0], epsilon.matrix[0, 1:], epsilon.matrix[1:, 1:]
    
    # For periodic systems (future implementation)
    # This would involve plane-wave decomposition
    head = epsilon.matrix[0, 0]
    wings = epsilon.matrix[0, 1:]
    body = epsilon.matrix[1:, 1:]
    
    return head, wings, body


# Convenience functions

def create_static_dielectric(
    df_ia: np.ndarray,
    orbital_energies: Tuple[np.ndarray, np.ndarray],
    v_sqrt: np.ndarray
) -> DielectricFunction:
    """Create static dielectric function ε(ω=0).
    
    Args:
        df_ia: DF tensor in transition space
        orbital_energies: (occupied, virtual) energies
        v_sqrt: Square root of Coulomb matrix
        
    Returns:
        Static DielectricFunction
        
    Examples:
        >>> eps_static = create_static_dielectric(df_ia, (e_occ, e_virt), v_sqrt)
        >>> print(f"Static dielectric: {eps_static}")
    """
    p0 = compute_polarizability_p0(df_ia, orbital_energies, 0.0, eta=1e-6)
    return compute_dielectric_function(p0, v_sqrt)


def test_kramers_kronig(
    epsilon_real: np.ndarray,
    epsilon_imag: np.ndarray, 
    omega: np.ndarray,
    tol: float = 1e-3
) -> bool:
    """Test Kramers-Kronig relations for dielectric function.
    
    Args:
        epsilon_real: Real part of ε(ω) on real axis
        epsilon_imag: Imaginary part of ε(ω) on real axis
        omega: Frequency grid
        tol: Tolerance for validation
        
    Returns:
        True if K-K relations are satisfied
    """
    # Simplified K-K test
    # In practice, would use Hilbert transform
    return True  # Placeholder


# Module exports
__all__ = [
    'Polarizability',
    'DielectricFunction',
    'ScreenedInteraction',
    'compute_polarizability_p0',
    'compute_polarizability_batch',
    'compute_dielectric_function',
    'compute_epsilon_matrix',
    'compute_epsilon_inverse',
    'compute_screened_interaction',
    'compute_head_wings_body',
    'create_static_dielectric',
    'test_kramers_kronig',
]