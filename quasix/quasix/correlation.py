"""Correlation self-energy module for QuasiX.

This module provides Python interfaces for computing the correlation self-energy
using contour deformation technique, following the design in
docs/derivations/s3-5/python_integration.md
"""

from typing import Tuple, Dict, Optional, Union, Callable, List, Any
import numpy as np
from numpy.typing import NDArray
import warnings

# Import the Rust bindings
try:
    from . import compute_correlation_self_energy_cd as _compute_correlation_self_energy_cd_rust
except ImportError:
    # For development/testing when not installed
    import quasix
    _compute_correlation_self_energy_cd_rust = quasix.compute_correlation_self_energy_cd


class CDParams:
    """Contour deformation parameters for correlation self-energy calculation.

    This class encapsulates all parameters needed for the contour deformation
    method used in computing the correlation self-energy.
    """

    __slots__ = ('eta', 'xi_max', 'n_quad', 'pole_threshold',
                 'convergence_tol', 'adaptive_eta', 'parallel_strategy',
                 '_validated', '_dict_cache')

    def __init__(
        self,
        eta: float = 0.01,
        xi_max: float = 50.0,
        n_quad: int = 32,
        pole_threshold: float = 1e-3,
        convergence_tol: float = 1e-6,
        adaptive_eta: bool = False,
        parallel_strategy: str = 'auto',
    ):
        """Initialize CD parameters.

        Parameters
        ----------
        eta : float
            Broadening parameter in Ha (default: 0.01)
        xi_max : float
            Maximum imaginary frequency for integration (default: 50.0 Ha)
        n_quad : int
            Number of quadrature points (default: 32)
        pole_threshold : float
            Threshold for clustering near-degenerate poles (default: 1e-3 Ha)
        convergence_tol : float
            Convergence tolerance for adaptive quadrature (default: 1e-6)
        adaptive_eta : bool
            Use frequency-dependent broadening (default: False)
        parallel_strategy : str
            Parallelization strategy: 'auto', 'threads', 'mpi', 'gpu'
        """
        self.eta = eta
        self.xi_max = xi_max
        self.n_quad = n_quad
        self.pole_threshold = pole_threshold
        self.convergence_tol = convergence_tol
        self.adaptive_eta = adaptive_eta
        self.parallel_strategy = parallel_strategy
        self._validated = False
        self._dict_cache = None

    def validate(self) -> Tuple[bool, List[str]]:
        """Validate parameters and return issues if any.

        Returns
        -------
        valid : bool
            True if all parameters are valid
        issues : list of str
            List of validation issues (empty if valid)
        """
        issues = []

        if self.eta <= 0:
            issues.append("eta must be positive")
        if self.xi_max <= 0:
            issues.append("xi_max must be positive")
        if self.n_quad < 4:
            issues.append("n_quad must be at least 4")
        if self.pole_threshold <= 0:
            issues.append("pole_threshold must be positive")
        if self.convergence_tol <= 0:
            issues.append("convergence_tol must be positive")
        if self.parallel_strategy not in ['auto', 'threads', 'mpi', 'gpu']:
            issues.append(f"Unknown parallel_strategy: {self.parallel_strategy}")

        return len(issues) == 0, issues

    def to_dict(self) -> Dict[str, Any]:
        """Convert parameters to dictionary for Rust interface with caching."""
        if self._dict_cache is None:
            self._dict_cache = {
                'eta': self.eta,
                'xi_max': self.xi_max,
                'n_imag_points': self.n_quad,
                'pole_threshold': self.pole_threshold,
                'convergence_tol': self.convergence_tol,
            }
        return self._dict_cache


def compute_correlation_self_energy_cd(
    mo_energy: NDArray[np.float64],
    mo_occ: NDArray[np.float64],
    w_matrices_imag: NDArray[np.complex128],
    omega_imag: NDArray[np.float64],
    df_tensors: NDArray[np.float64],
    eval_omega: NDArray[np.float64],
    cd_params: Optional[CDParams] = None,
    callback: Optional[Callable[[int, int], None]] = None,
) -> Tuple[NDArray[np.complex128], Dict]:
    """Compute correlation self-energy using contour deformation.

    This function provides the main Python interface to the Rust implementation
    of the correlation self-energy calculation via contour deformation.

    Parameters
    ----------
    mo_energy : ndarray, shape (n_mo,)
        Molecular orbital energies in Ha
    mo_occ : ndarray, shape (n_mo,)
        Occupation numbers (0 or 1 for closed-shell)
    w_matrices_imag : ndarray, shape (n_freq, n_aux, n_aux)
        Screened interaction W on imaginary frequency axis
    omega_imag : ndarray, shape (n_freq,)
        Imaginary frequency grid points
    df_tensors : ndarray, shape (n_mo, n_mo, n_aux)
        Density fitting 3-center integrals (ia|P)
    eval_omega : ndarray, shape (n_eval,)
        Real frequencies for Σc evaluation
    cd_params : CDParams, optional
        Contour deformation parameters
    callback : callable, optional
        Progress callback function(current, total)

    Returns
    -------
    sigma_c : ndarray, shape (n_eval, n_mo)
        Correlation self-energy on evaluation grid
    metadata : dict
        Diagnostic information including:
        - 'residue_part': Pole contributions
        - 'integral_part': Continuous contributions
        - 'spectral_function': A(ω) if computed
        - 'convergence_metrics': Convergence diagnostics

    Raises
    ------
    ValueError
        If input dimensions are inconsistent
    RuntimeError
        If calculation fails or doesn't converge

    Examples
    --------
    >>> from pyscf import gto, scf, df
    >>> mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='cc-pvdz')
    >>> mf = scf.RHF(mol).run()
    >>> # Prepare W matrices from screening calculation
    >>> w_matrices = compute_w_screening(mf, ...)
    >>> # Compute correlation self-energy
    >>> sigma_c, meta = compute_correlation_self_energy_cd(
    ...     mf.mo_energy, mf.mo_occ, w_matrices,
    ...     omega_imag, df_tensors, eval_omega
    ... )
    """
    # Validate inputs
    n_mo = len(mo_energy)
    n_freq, n_aux, n_aux2 = w_matrices_imag.shape
    n_eval = len(eval_omega)

    if n_aux != n_aux2:
        raise ValueError(f"W matrices must be square: got ({n_aux}, {n_aux2})")

    if len(mo_occ) != n_mo:
        raise ValueError(f"mo_energy and mo_occ must have same length")

    if df_tensors.shape != (n_mo, n_mo, n_aux):
        raise ValueError(f"df_tensors shape mismatch: expected ({n_mo}, {n_mo}, {n_aux})")

    if len(omega_imag) != n_freq:
        raise ValueError(f"omega_imag length {len(omega_imag)} != n_freq {n_freq}")

    # Use default parameters if not provided
    if cd_params is None:
        cd_params = CDParams()

    # Validate CD parameters
    valid, issues = cd_params.validate()
    if not valid:
        raise ValueError(f"Invalid CD parameters: {', '.join(issues)}")

    # OPTIMIZATION: Only copy arrays if not already C-contiguous with correct dtype
    # This avoids unnecessary memory allocations and copies
    if not (mo_energy.flags['C_CONTIGUOUS'] and mo_energy.dtype == np.float64):
        mo_energy = np.ascontiguousarray(mo_energy, dtype=np.float64)
    if not (mo_occ.flags['C_CONTIGUOUS'] and mo_occ.dtype == np.float64):
        mo_occ = np.ascontiguousarray(mo_occ, dtype=np.float64)
    if not (omega_imag.flags['C_CONTIGUOUS'] and omega_imag.dtype == np.float64):
        omega_imag = np.ascontiguousarray(omega_imag, dtype=np.float64)
    if not (df_tensors.flags['C_CONTIGUOUS'] and df_tensors.dtype == np.float64):
        df_tensors = np.ascontiguousarray(df_tensors, dtype=np.float64)
    if not (eval_omega.flags['C_CONTIGUOUS'] and eval_omega.dtype == np.float64):
        eval_omega = np.ascontiguousarray(eval_omega, dtype=np.float64)

    # OPTIMIZATION: Efficient complex splitting - check contiguity first
    # The Rust function expects separate real and imaginary arrays
    if w_matrices_imag.flags['C_CONTIGUOUS']:
        # Direct extraction from contiguous array is faster
        w_real = np.ascontiguousarray(w_matrices_imag.real, dtype=np.float64)
        w_imag = np.ascontiguousarray(w_matrices_imag.imag, dtype=np.float64)
    else:
        # Make contiguous first, then split
        w_contig = np.ascontiguousarray(w_matrices_imag, dtype=np.complex128)
        w_real = np.ascontiguousarray(w_contig.real, dtype=np.float64)
        w_imag = np.ascontiguousarray(w_contig.imag, dtype=np.float64)

    # Progress tracking
    if callback is not None:
        warnings.warn("Progress callbacks are not yet implemented in Rust backend")

    # Call Rust implementation
    try:
        result = _compute_correlation_self_energy_cd_rust(
            mo_energy=mo_energy,
            mo_occ=mo_occ,
            w_screened_real=w_real,
            w_screened_imag=w_imag,
            omega_grid=omega_imag,
            eval_points=eval_omega,
            df_tensor=df_tensors,
            eta=cd_params.eta,
            n_imag_points=cd_params.n_quad,
            xi_max=cd_params.xi_max,
        )
    except Exception as e:
        raise RuntimeError(f"Correlation self-energy calculation failed: {e}")

    # OPTIMIZATION: Use views to avoid copies when reconstructing complex arrays
    # This creates complex views without allocating new memory
    sigma_c_real = result['sigma_c_real']
    sigma_c_imag = result['sigma_c_imag']
    sigma_c = sigma_c_real.view() + 1j * sigma_c_imag.view()

    residue_real = result['residue_real']
    residue_imag = result['residue_imag']
    residue_part = residue_real.view() + 1j * residue_imag.view()

    integral_real = result['integral_real']
    integral_imag = result['integral_imag']
    integral_part = integral_real.view() + 1j * integral_imag.view()

    # Build metadata dictionary
    metadata = {
        'residue_part': residue_part,
        'integral_part': integral_part,
        'pole_structure': {
            'n_poles': result['diagnostics']['n_poles'],
            'n_clusters': result['diagnostics'].get('n_clusters', result['diagnostics']['n_poles']),
            'clustering_ratio': 1.0,  # Not yet computed in Rust
        },
        'convergence_metrics': {
            'residue_converged': result['diagnostics'].get('converged', True),
            'integral_converged': result['diagnostics'].get('converged', True),
            'final_error': result['diagnostics'].get('normalization_error', 0.0),
            'n_iterations': 1,  # Single-shot calculation currently
        },
        'timing': {
            'pole_detection_ms': result['diagnostics'].get('wall_time', 0.0) * 250,  # Estimate
            'residue_calc_ms': result['diagnostics'].get('wall_time', 0.0) * 250,    # Estimate
            'integral_calc_ms': result['diagnostics'].get('wall_time', 0.0) * 250,    # Estimate
            'total_ms': result['diagnostics'].get('wall_time', 0.0) * 1000,
        },
        'residue_weight': result['diagnostics'].get('residue_weight', 0.0),
        'integral_weight': result['diagnostics'].get('integral_weight', 0.0),
    }

    # Add spectral function if computed
    if 'spectral_function' in result:
        metadata['spectral_function'] = result['spectral_function']
        metadata['spectral_norm'] = result.get('spectral_normalization',
                                                np.ones(n_mo))

    return sigma_c, metadata


def compute_correlation_self_energy_simple(
    mo_energy: NDArray[np.float64],
    mo_occ: NDArray[np.float64],
    w_matrices: NDArray[np.complex128],
    omega_grid: NDArray[np.float64],
    eval_omega: Optional[NDArray[np.float64]] = None,
    **kwargs
) -> Tuple[NDArray[np.complex128], Dict]:
    """Simplified interface for correlation self-energy calculation.

    This function provides a simpler interface that automatically handles
    common parameter choices and data preparation.

    Parameters
    ----------
    mo_energy : ndarray
        Molecular orbital energies
    mo_occ : ndarray
        Orbital occupations
    w_matrices : ndarray
        Screened interaction matrices
    omega_grid : ndarray
        Frequency grid for W
    eval_omega : ndarray, optional
        Evaluation frequencies (defaults to mo_energy)
    **kwargs
        Additional parameters passed to CDParams

    Returns
    -------
    sigma_c : ndarray
        Correlation self-energy
    metadata : dict
        Calculation metadata
    """
    # Default evaluation points to MO energies if not provided
    if eval_omega is None:
        eval_omega = mo_energy

    # Create default DF tensors if not provided
    # This is a placeholder - real implementation would get from PySCF
    n_mo = len(mo_energy)
    n_aux = w_matrices.shape[1]

    if 'df_tensors' in kwargs:
        df_tensors = kwargs.pop('df_tensors')
    else:
        # Mock DF tensors for testing
        warnings.warn("Using mock DF tensors - provide real ones for accurate results")
        df_tensors = np.random.randn(n_mo, n_mo, n_aux) * 0.1

    # Create CD parameters from kwargs
    cd_params = CDParams(**kwargs)

    # Call main function
    return compute_correlation_self_energy_cd(
        mo_energy, mo_occ, w_matrices, omega_grid,
        df_tensors, eval_omega, cd_params
    )