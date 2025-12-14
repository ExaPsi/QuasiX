"""
Screened interaction W(ω) calculations for QuasiX GW.

This module implements the screened Coulomb interaction:
W(ω) = ε^-1(ω) v

Where:
- ε(ω) is the dielectric function
- v is the bare Coulomb interaction

Enhanced with S3-3 functionality: advanced dielectric solvers with
stability monitoring and multiple backend options.
"""

from typing import Optional, Tuple, Union, Dict, Any, List
import numpy as np
from dataclasses import dataclass
import logging
import warnings

logger = logging.getLogger(__name__)

# Import Rust functions when available
# Note: These are optional screening accelerators. The main evGW functionality
# uses the run_evgw/evgw functions which are always available when Rust is built.
RUST_AVAILABLE = False
_compute_screened_interaction = None
_compute_screened_interaction_batch = None
_compute_screened_interaction_s33 = None
_compute_polarizability_p0 = None

try:
    from .quasix import compute_polarizability_p0 as _compute_polarizability_p0
    RUST_AVAILABLE = True
except ImportError:
    pass

# These optional batch functions may not be exposed in PyO3 bindings
try:
    from .quasix import (
        compute_screened_interaction as _compute_screened_interaction,
        compute_screened_interaction_batch as _compute_screened_interaction_batch,
        compute_screened_interaction_s33 as _compute_screened_interaction_s33,
    )
except ImportError:
    # Optional batch screening functions not available - this is fine
    pass

# Import Python-only functions from dielectric module
from .dielectric import (
    Polarizability,
    DielectricFunction,
    ScreenedInteraction,
    compute_polarizability_p0,
    compute_dielectric_function,
    compute_epsilon_inverse,
)


@dataclass
class SolverConfig:
    """Configuration for the advanced dielectric solver.

    Attributes:
        solver_type: 'direct', 'regularized', or 'adaptive'
        solver_backend: 'lu', 'cholesky', 'svd', or 'auto'
        regularization: Regularization parameter for ill-conditioned matrices
        condition_threshold: Threshold for detecting ill-conditioning
        self_consistency_tol: Tolerance for self-consistency check
        monitor_condition: Whether to monitor condition number
        svd_threshold: Truncation threshold for SVD
        max_iterations: Maximum iterations for iterative solvers
        block_size: Block size for cache-efficient operations
    """
    solver_type: str = 'direct'
    solver_backend: str = 'auto'
    regularization: float = 1e-10
    condition_threshold: float = 1e12
    self_consistency_tol: float = 1e-8
    monitor_condition: bool = True
    svd_threshold: float = 1e-14
    max_iterations: int = 100
    block_size: int = 1024

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for passing to Rust."""
        return {
            'solver_type': self.solver_type,
            'regularization': self.regularization,
            'condition_threshold': self.condition_threshold,
            'self_consistency_tol': self.self_consistency_tol,
            'monitor_condition': self.monitor_condition,
            'svd_threshold': self.svd_threshold,
            'max_iterations': self.max_iterations,
            'block_size': self.block_size,
        }


def compute_screened_interaction(
    p0: Union[np.ndarray, complex],
    v_sqrt: np.ndarray,
    solver: str = 'auto',
    config: Optional[SolverConfig] = None,
    return_diagnostics: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, Dict[str, Any]]]:
    """
    Compute screened interaction W(ω) using advanced dielectric solver.

    This implements the symmetrized formulation for numerical stability:
    W(ω) = v^{1/2} [1 - v^{1/2} P^0(ω) v^{1/2}]^{-1} v^{1/2}

    Args:
        p0: Polarizability matrix P0(ω), shape (n_aux, n_aux)
            Can be real or complex
        v_sqrt: Square root of Coulomb metric V^(1/2), shape (n_aux, n_aux)
        solver: Solver backend - 'lu', 'cholesky', 'svd', or 'auto'
        config: Optional solver configuration
        return_diagnostics: If True, return (W, diagnostics) tuple

    Returns:
        W: Screened interaction matrix, shape (n_aux, n_aux)
        diagnostics: (optional) Dictionary with solver diagnostics

    Example:
        >>> # Simple usage with default settings
        >>> W = compute_screened_interaction(p0, v_sqrt)

        >>> # Advanced usage with monitoring
        >>> config = SolverConfig(
        ...     solver_type='adaptive',
        ...     monitor_condition=True,
        ...     condition_threshold=1e10
        ... )
        >>> W, diag = compute_screened_interaction(
        ...     p0, v_sqrt, solver='auto',
        ...     config=config, return_diagnostics=True
        ... )
        >>> print(f"Condition number: {diag['condition_number']:.2e}")
    """
    if not RUST_AVAILABLE:
        raise RuntimeError("Rust extensions required for advanced screening calculations")

    # Ensure p0 is complex
    if np.isrealobj(p0):
        p0_complex = p0.astype(np.complex128)
    else:
        p0_complex = p0

    # Split into real and imaginary parts for Rust
    p0_real = np.real(p0_complex)
    p0_imag = np.imag(p0_complex)

    # Use default config if not provided
    if config is None:
        config = SolverConfig()

    # Call Rust implementation
    result = _compute_screened_interaction_s33(
        p0_real, p0_imag, v_sqrt,
        solver, config.to_dict()
    )

    # Reconstruct complex W
    w = result['w_real'] + 1j * result['w_imag']

    if return_diagnostics:
        diagnostics = {
            'condition_number': result['condition_number'],
            'solver_used': result['solver_used'],
            'hermiticity_error': result['hermiticity_error'],
            'self_consistency_error': result['self_consistency_error'],
            'diagnostics': result['diagnostics']
        }
        return w, diagnostics
    else:
        # Issue warnings for potential problems
        if result['condition_number'] > 1e10:
            logger.warning(
                f"Matrix may be ill-conditioned: condition number = {result['condition_number']:.2e}"
            )
        if result['hermiticity_error'] > 1e-8:
            logger.warning(
                f"W not perfectly Hermitian: error = {result['hermiticity_error']:.2e}"
            )
        if result['self_consistency_error'] > 1e-6:
            logger.warning(
                f"Self-consistency check failed: error = {result['self_consistency_error']:.2e}"
            )
        return w


@dataclass
class ScreeningCalculator:
    """Calculator for screened interactions W(ω).

    Enhanced with S3-3 functionality for robust computation.

    Attributes:
        df_ia: DF tensor in transition space
        orbital_energies: (occupied, virtual) orbital energies
        v_sqrt: Square root of Coulomb matrix
        v_full: Full Coulomb matrix
        eta: Broadening parameter
        config: Solver configuration
    """
    df_ia: np.ndarray
    orbital_energies: Tuple[np.ndarray, np.ndarray]
    v_sqrt: np.ndarray
    v_full: Optional[np.ndarray] = None
    eta: float = 0.01
    config: Optional[SolverConfig] = None

    def __post_init__(self):
        """Initialize derived attributes."""
        if self.v_full is None:
            self.v_full = self.v_sqrt @ self.v_sqrt
        if self.config is None:
            self.config = SolverConfig()

    def compute_p0(self, omega: complex) -> np.ndarray:
        """Compute polarizability P0(ω)."""
        if RUST_AVAILABLE:
            return _compute_polarizability_p0(
                self.df_ia,
                self.orbital_energies[0],  # occupied
                self.orbital_energies[1],  # virtual
                omega.real,
                omega.imag,
                self.eta
            )
        else:
            return compute_polarizability_p0(
                self.df_ia,
                self.orbital_energies,
                omega,
                self.eta
            )

    def compute_w_static(self) -> ScreenedInteraction:
        """Compute static screened interaction W(ω=0).

        Returns:
            Static ScreenedInteraction with diagnostics
        """
        # Compute P0(0)
        p0 = self.compute_p0(0.0 + 0j)

        # Compute W using advanced solver
        w, diagnostics = compute_screened_interaction(
            p0, self.v_sqrt,
            solver=self.config.solver_backend,
            config=self.config,
            return_diagnostics=True
        )

        result = ScreenedInteraction(
            matrix=w,
            frequency=0.0,
            naux=w.shape[0],
            static=True
        )
        result.diagnostics = diagnostics  # Attach diagnostics
        return result

    def compute_w_dynamic(
        self,
        frequency: complex,
        monitor: bool = False
    ) -> ScreenedInteraction:
        """Compute dynamic screened interaction W(ω).

        Args:
            frequency: Frequency (real or imaginary)
            monitor: Whether to monitor solver stability

        Returns:
            Dynamic ScreenedInteraction with optional diagnostics
        """
        # Compute P0(ω)
        p0 = self.compute_p0(frequency)

        # Optionally enable monitoring
        config = self.config
        if monitor:
            config = SolverConfig(**self.config.to_dict())
            config.monitor_condition = True

        # Compute W using advanced solver
        if monitor:
            w, diagnostics = compute_screened_interaction(
                p0, self.v_sqrt,
                solver=config.solver_backend,
                config=config,
                return_diagnostics=True
            )
        else:
            w = compute_screened_interaction(
                p0, self.v_sqrt,
                solver=config.solver_backend,
                config=config,
                return_diagnostics=False
            )
            diagnostics = None

        result = ScreenedInteraction(
            matrix=w,
            frequency=frequency,
            naux=w.shape[0],
            static=False
        )
        if diagnostics:
            result.diagnostics = diagnostics
        return result

    def compute_w_batch(
        self,
        frequencies: np.ndarray,
        parallel: bool = True,
        monitor_condition: bool = False
    ) -> List[ScreenedInteraction]:
        """Compute W(ω) for multiple frequencies.

        Args:
            frequencies: Array of frequencies
            parallel: Whether to use parallel computation
            monitor_condition: Whether to monitor condition numbers

        Returns:
            List of ScreenedInteraction objects with diagnostics
        """
        W_list = []
        condition_numbers = []

        for i, freq in enumerate(frequencies):
            W = self.compute_w_dynamic(freq, monitor=monitor_condition)
            W_list.append(W)

            if monitor_condition and hasattr(W, 'diagnostics'):
                cond = W.diagnostics.get('condition_number', np.inf)
                condition_numbers.append(cond)

                # Issue warning for ill-conditioned frequencies
                if cond > 1e12:
                    logger.warning(
                        f"Frequency {i} (ω={freq:.4f}): "
                        f"ill-conditioned matrix (κ={cond:.2e})"
                    )

        # Attach batch statistics if monitoring
        if monitor_condition and condition_numbers:
            batch_stats = {
                'min_condition': np.min(condition_numbers),
                'max_condition': np.max(condition_numbers),
                'mean_condition': np.mean(condition_numbers),
                'problematic_frequencies': [
                    i for i, c in enumerate(condition_numbers) if c > 1e10
                ]
            }
            # Store in first element for convenience
            W_list[0].batch_diagnostics = batch_stats

            logger.info(
                f"Batch statistics: κ ∈ [{batch_stats['min_condition']:.2e}, "
                f"{batch_stats['max_condition']:.2e}], "
                f"{len(batch_stats['problematic_frequencies'])} problematic points"
            )

        return W_list

    def compute_w_minus_v(self, frequency: complex = 0.0) -> np.ndarray:
        """Compute W(ω) - v for correlation self-energy.

        Args:
            frequency: Frequency

        Returns:
            W(ω) - v matrix
        """
        W = self.compute_w_dynamic(frequency)
        return W.matrix - self.v_full

    def analyze_solver_stability(
        self,
        test_frequencies: Optional[np.ndarray] = None
    ) -> Dict[str, Any]:
        """Analyze solver stability across frequency range.

        Args:
            test_frequencies: Frequencies to test (default: logarithmic grid)

        Returns:
            Dictionary with stability analysis
        """
        if test_frequencies is None:
            # Default logarithmic grid from 0 to 100 a.u.
            test_frequencies = np.logspace(-3, 2, 20)

        results = {
            'frequencies': test_frequencies,
            'condition_numbers': [],
            'hermiticity_errors': [],
            'self_consistency_errors': [],
            'solver_recommendations': {}
        }

        for omega in test_frequencies:
            p0 = self.compute_p0(omega + 0j)

            # Test each solver
            for solver in ['lu', 'cholesky', 'svd']:
                try:
                    _, diag = compute_screened_interaction(
                        p0, self.v_sqrt,
                        solver=solver,
                        config=self.config,
                        return_diagnostics=True
                    )

                    if solver == 'lu':  # Store results for default solver
                        results['condition_numbers'].append(diag['condition_number'])
                        results['hermiticity_errors'].append(diag['hermiticity_error'])
                        results['self_consistency_errors'].append(diag['self_consistency_error'])

                    # Track which solvers work
                    if solver not in results['solver_recommendations']:
                        results['solver_recommendations'][solver] = []
                    results['solver_recommendations'][solver].append(omega)

                except Exception as e:
                    logger.debug(f"Solver {solver} failed at ω={omega}: {e}")

        # Convert to arrays
        results['condition_numbers'] = np.array(results['condition_numbers'])
        results['hermiticity_errors'] = np.array(results['hermiticity_errors'])
        results['self_consistency_errors'] = np.array(results['self_consistency_errors'])

        # Add summary statistics
        results['summary'] = {
            'max_condition_number': np.max(results['condition_numbers']),
            'frequencies_ill_conditioned': np.sum(results['condition_numbers'] > 1e10),
            'recommended_solver': self._recommend_solver(results),
            'regularization_needed': np.any(results['condition_numbers'] > 1e12)
        }

        return results

    def _recommend_solver(self, stability_results: Dict) -> str:
        """Recommend best solver based on stability analysis."""
        solver_scores = {}

        for solver, freqs in stability_results['solver_recommendations'].items():
            # Score based on how many frequencies work
            solver_scores[solver] = len(freqs)

        if not solver_scores:
            return 'svd'  # Most robust fallback

        # Return solver with highest score
        return max(solver_scores, key=solver_scores.get)


def verify_self_consistency(
    w: np.ndarray,
    p0: np.ndarray,
    v: np.ndarray,
    tol: float = 1e-8
) -> Tuple[bool, float]:
    """
    Verify self-consistency relation: W = v + vP⁰W

    Args:
        w: Screened interaction W
        p0: Polarizability P⁰
        v: Coulomb interaction
        tol: Tolerance for verification

    Returns:
        (is_consistent, error): Whether consistent and the error norm
    """
    # Compute vP⁰W
    vp0w = v @ p0 @ w

    # Compute residual
    residual = w - v - vp0w

    # Relative error
    error = np.linalg.norm(residual, 'fro') / np.linalg.norm(w, 'fro')

    return error < tol, error


def compute_screened_interaction_symmetric(
    df_ia: np.ndarray,
    orbital_energies: Tuple[np.ndarray, np.ndarray],
    v_sqrt: np.ndarray,
    frequency: complex,
    eta: float = 0.01,
    config: Optional[SolverConfig] = None
) -> np.ndarray:
    """Compute symmetrized screened interaction for numerical stability.

    Uses the symmetrized form:
    W_sym = v^1/2 [I - v^1/2 P0 v^1/2]^-1 v^1/2

    Args:
        df_ia: DF tensor in transition space
        orbital_energies: (occupied, virtual) energies
        v_sqrt: Square root of Coulomb matrix
        frequency: Frequency
        eta: Broadening parameter
        config: Optional solver configuration

    Returns:
        Symmetrized W matrix
    """
    # Create calculator
    calc = ScreeningCalculator(
        df_ia=df_ia,
        orbital_energies=orbital_energies,
        v_sqrt=v_sqrt,
        eta=eta,
        config=config or SolverConfig()
    )

    # Compute P0(ω)
    p0 = calc.compute_p0(frequency)

    # Use advanced solver
    w = compute_screened_interaction(
        p0, v_sqrt,
        solver=config.solver_backend if config else 'auto',
        config=config,
        return_diagnostics=False
    )

    return w


def compute_plasmon_pole_model(
    eps_static: DielectricFunction,
    eps_optical: DielectricFunction,
    plasma_freq: Optional[float] = None
) -> Dict[str, np.ndarray]:
    """Compute plasmon pole model parameters.

    Godby-Needs plasmon pole approximation:
    ε^-1(ω) - I ≈ Ω²_p / (ω² - ω²_p + iηω)

    Args:
        eps_static: Static dielectric ε(0)
        eps_optical: Optical dielectric ε(∞)
        plasma_freq: Plasma frequency (optional)

    Returns:
        Dictionary with plasmon pole parameters
    """
    naux = eps_static.naux

    # Compute ε^-1(0) - I
    eps_inv_static = compute_epsilon_inverse(eps_static)
    chi_static = eps_inv_static - np.eye(naux)

    # Estimate plasma frequency if not provided
    if plasma_freq is None:
        # Simple estimate from f-sum rule
        plasma_freq = np.sqrt(np.trace(chi_static) * 4.0)

    # Compute pole strengths
    omega_p = np.ones(naux) * plasma_freq
    Omega_p = chi_static.diagonal() * omega_p

    return {
        'omega_p': omega_p,  # Pole positions
        'Omega_p': Omega_p,  # Pole strengths
        'plasma_freq': plasma_freq
    }


def apply_plasmon_pole(
    ppm_params: Dict[str, np.ndarray],
    frequency: complex
) -> np.ndarray:
    """Evaluate plasmon pole model at given frequency.

    Args:
        ppm_params: Plasmon pole parameters
        frequency: Evaluation frequency

    Returns:
        ε^-1(ω) from plasmon pole model
    """
    omega_p = ppm_params['omega_p']
    Omega_p = ppm_params['Omega_p']

    # PPM formula
    naux = len(omega_p)
    eps_inv = np.eye(naux)

    for i in range(naux):
        denom = frequency**2 - omega_p[i]**2 + 1j * 0.01 * frequency
        eps_inv[i, i] += Omega_p[i] / denom

    return eps_inv


def compute_head_correction(
    W: ScreenedInteraction,
    q_vector: Optional[np.ndarray] = None,
    cell_volume: Optional[float] = None
) -> float:
    """Compute head correction for q→0 singularity (periodic systems).

    Args:
        W: Screened interaction
        q_vector: q-point vector
        cell_volume: Unit cell volume

    Returns:
        Head correction value
    """
    if q_vector is None or cell_volume is None:
        return 0.0  # No correction for molecules

    # For q→0 limit in periodic systems
    # W_head ~ 4π/q² * 1/ε_macro
    q_norm = np.linalg.norm(q_vector)
    if q_norm < 1e-6:
        # Use macroscopic dielectric constant
        eps_macro = W.matrix[0, 0] / (4 * np.pi)
        correction = 4 * np.pi / (q_norm**2 * eps_macro * cell_volume)
        return correction

    return 0.0


def test_spectral_sum_rule(
    W_frequencies: List[ScreenedInteraction],
    v_coulomb: np.ndarray,
    tol: float = 0.01
) -> Tuple[bool, Dict[str, float]]:
    """Test spectral sum rule for W(ω).

    Sum rule: ∫ Im[W(ω) - v] dω/π = -v

    Args:
        W_frequencies: List of W at different frequencies
        v_coulomb: Bare Coulomb matrix
        tol: Relative tolerance

    Returns:
        (is_satisfied, diagnostics): Whether satisfied and diagnostic info
    """
    if len(W_frequencies) < 2:
        return False, {'error': 'Insufficient frequency points'}

    # Extract imaginary parts
    integrand = []
    frequencies = []
    for W in W_frequencies:
        if np.iscomplex(W.frequency):
            continue  # Skip imaginary frequencies
        integrand.append(np.imag(W.matrix - v_coulomb))
        frequencies.append(np.real(W.frequency))

    if len(integrand) < 2:
        return True, {'error': 'Insufficient real frequency points'}

    # Numerical integration (trapezoidal)
    integrand = np.array(integrand)
    frequencies = np.array(frequencies)

    integral = np.trapz(integrand, frequencies, axis=0) / np.pi

    # Check sum rule
    error = np.linalg.norm(integral + v_coulomb, 'fro') / np.linalg.norm(v_coulomb, 'fro')

    diagnostics = {
        'relative_error': error,
        'max_deviation': np.max(np.abs(integral + v_coulomb)),
        'n_frequencies': len(frequencies),
        'frequency_range': (frequencies[0], frequencies[-1])
    }

    return error < tol, diagnostics


# Module exports
__all__ = [
    'SolverConfig',
    'ScreeningCalculator',
    'compute_screened_interaction',
    'compute_screened_interaction_symmetric',
    'verify_self_consistency',
    'compute_plasmon_pole_model',
    'apply_plasmon_pole',
    'compute_head_correction',
    'test_spectral_sum_rule',
]