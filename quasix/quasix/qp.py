"""
Quasiparticle solver module for QuasiX GW calculations.

This module implements linearized and full self-consistent
quasiparticle equation solvers with PySCF-compatible API.

S3-6 Implementation: Newton-Raphson with bisection fallback
"""

from typing import Optional, Tuple, Union, Dict, Any, List, Callable
import numpy as np
from dataclasses import dataclass, field
import logging
import warnings

logger = logging.getLogger(__name__)

# Import Rust functions when available
try:
    from .quasix import (
        solve_quasiparticle_equations as _solve_quasiparticle_equations,
    )
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    _solve_quasiparticle_equations = None


@dataclass
class QuasiparticleState:
    """Container for a quasiparticle state.
    
    Attributes:
        index: Orbital index
        energy_dft: DFT orbital energy
        energy_qp: Quasiparticle energy
        z_factor: Quasiparticle weight (renormalization factor)
        sigma_x: Exchange self-energy <Σ^x>
        sigma_c: Correlation self-energy <Σ^c>
        vxc: Exchange-correlation potential <v_xc>
        converged: Whether QP equation converged
    """
    index: int
    energy_dft: float
    energy_qp: float
    z_factor: float
    sigma_x: float
    sigma_c: complex
    vxc: float
    converged: bool = True
    
    @property
    def correction(self) -> float:
        """QP correction: E_QP - E_DFT"""
        return self.energy_qp - self.energy_dft
    
    @property
    def self_energy(self) -> complex:
        """Total self-energy: Σ = Σ^x + Σ^c"""
        return self.sigma_x + self.sigma_c
    
    def __repr__(self) -> str:
        return (f"QP(i={self.index}, E_DFT={self.energy_dft:.3f}, "
                f"E_QP={self.energy_qp:.3f}, Z={self.z_factor:.3f})")


@dataclass
class QPSolverConfig:
    """Configuration for quasiparticle solver (S3-6).

    Attributes:
        energy_tolerance: Convergence tolerance for energy (Ha)
        residual_tolerance: Convergence tolerance for residual
        max_newton_iterations: Maximum iterations for Newton-Raphson
        max_bisection_iterations: Maximum iterations for bisection
        initial_damping: Initial damping factor for Newton-Raphson
        min_damping: Minimum damping factor
        max_damping: Maximum damping factor
        max_energy_step: Maximum allowed energy step (Ha)
        derivative_delta: Finite difference step for derivatives
        use_line_search: Use line search for better convergence
        use_richardson: Use Richardson extrapolation for derivatives
        use_bisection_fallback: Use bisection fallback if Newton fails
        z_bounds: Z-factor bounds (min, max) - enforced physical range
        n_threads: Number of threads for parallel execution
        energy_model: Energy dependence model ('lorentzian', 'pole', 'linear')
        energy_dependence_alpha: Energy dependence strength parameter
    """
    # Convergence parameters
    energy_tolerance: float = 1e-6
    residual_tolerance: float = 1e-7
    max_newton_iterations: int = 50
    max_bisection_iterations: int = 100

    # Damping parameters
    initial_damping: float = 1.0
    min_damping: float = 0.01
    max_damping: float = 1.0
    damping_factor: float = 0.5  # Legacy compatibility

    # Numerical parameters
    max_energy_step: float = 0.5
    derivative_delta: float = 1e-5

    # Algorithm options
    use_line_search: bool = True
    use_richardson: bool = True
    use_bisection_fallback: bool = True

    # Physical bounds
    z_bounds: Tuple[float, float] = (0.01, 0.99)

    # Performance
    n_threads: Optional[int] = None

    # Energy dependence model
    energy_model: str = 'lorentzian'
    energy_dependence_alpha: float = 0.5

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for Rust interface."""
        return {
            'tolerance': self.energy_tolerance,
            'max_iterations': self.max_newton_iterations,
            'derivative_delta': self.derivative_delta,
            'damping_factor': self.damping_factor,
            'use_bisection_fallback': self.use_bisection_fallback,
            'energy_model': self.energy_model,
            'energy_dependence_alpha': self.energy_dependence_alpha,
        }


class LinearizedQPSolver:
    """Linearized quasiparticle equation solver.
    
    Solves: E_QP = ε_DFT + Z * Re[Σ(E_QP) - v_xc]
    Where: Z = [1 - ∂Σ/∂ω|_{E_QP}]^{-1}
    """
    
    def __init__(self, config: Optional[QPSolverConfig] = None):
        """Initialize solver with configuration."""
        self.config = config or QPSolverConfig()
    
    def solve(
        self,
        energy_dft: float,
        sigma_x: float,
        sigma_c_func: callable,
        vxc: float,
        orbital_index: int = 0
    ) -> QuasiparticleState:
        """Solve linearized QP equation for single orbital.
        
        Args:
            energy_dft: DFT orbital energy
            sigma_x: Exchange self-energy (real, frequency-independent)
            sigma_c_func: Function that returns Σ^c(ω) for given ω
            vxc: Exchange-correlation potential
            orbital_index: Orbital index for tracking
            
        Returns:
            QuasiparticleState with QP energy and Z-factor
        """
        # Initial guess at DFT energy
        omega = energy_dft
        
        for iteration in range(self.config.max_iter):
            # Evaluate correlation self-energy
            sigma_c = sigma_c_func(omega)
            
            # Compute derivative for Z-factor
            delta = 1e-4
            sigma_c_plus = sigma_c_func(omega + delta)
            sigma_c_minus = sigma_c_func(omega - delta)
            dsigma_dw = (sigma_c_plus - sigma_c_minus) / (2 * delta)
            
            # Z-factor
            z_factor = 1.0 / (1.0 - np.real(dsigma_dw))
            
            # Clip Z to physical range
            z_factor = np.clip(z_factor, 0.01, 0.99)
            
            # QP energy update
            sigma_total = sigma_x + np.real(sigma_c)
            omega_new = energy_dft + z_factor * (sigma_total - vxc)
            
            # Check convergence
            if abs(omega_new - omega) < self.config.tol_energy:
                return QuasiparticleState(
                    index=orbital_index,
                    energy_dft=energy_dft,
                    energy_qp=omega_new,
                    z_factor=z_factor,
                    sigma_x=sigma_x,
                    sigma_c=sigma_c,
                    vxc=vxc,
                    converged=True
                )
            
            omega = omega_new
        
        # Not converged
        logger.warning(f"QP solver did not converge for orbital {orbital_index}")
        return QuasiparticleState(
            index=orbital_index,
            energy_dft=energy_dft,
            energy_qp=omega,
            z_factor=z_factor,
            sigma_x=sigma_x,
            sigma_c=sigma_c,
            vxc=vxc,
            converged=False
        )
    
    def solve_batch(
        self,
        energies_dft: np.ndarray,
        sigma_x: np.ndarray,
        sigma_c_func: callable,
        vxc: np.ndarray
    ) -> List[QuasiparticleState]:
        """Solve QP equations for multiple orbitals.
        
        Args:
            energies_dft: Array of DFT energies
            sigma_x: Array of exchange self-energies
            sigma_c_func: Function that returns Σ^c(i, ω) for orbital i
            vxc: Array of xc potentials
            
        Returns:
            List of QuasiparticleState objects
        """
        qp_states = []
        for i, e_dft in enumerate(energies_dft):
            # Create orbital-specific Σ^c function
            def sigma_c_i(omega):
                return sigma_c_func(i, omega)
            
            qp = self.solve(
                e_dft,
                sigma_x[i],
                sigma_c_i,
                vxc[i],
                orbital_index=i
            )
            qp_states.append(qp)
        
        return qp_states


class QPSolver:
    """Main quasiparticle equation solver with PySCF-compatible API.

    This class provides the primary interface for solving quasiparticle
    equations using the high-performance Rust backend (S3-6).

    Attributes:
        config: Solver configuration
        verbose: Verbosity level (0-5)
        converged: Whether the last solve converged
        qp_energies: Quasiparticle energies from last solve
        z_factors: Z-factors from last solve
        diagnostics: Detailed diagnostics from last solve
        orbital_solutions: Per-orbital solution details
    """

    def __init__(self, config: Optional[QPSolverConfig] = None, verbose: int = 1):
        """Initialize QP solver.

        Args:
            config: Solver configuration (uses defaults if None)
            verbose: Verbosity level for logging
        """
        self.config = config or QPSolverConfig()
        self.verbose = verbose

        # Results from last solve
        self.converged = False
        self.qp_energies = None
        self.z_factors = None
        self.diagnostics = None
        self.orbital_solutions = None
        self._last_result = None

        # Progress callback for Jupyter notebooks
        self._progress_callback = None

    def kernel(
        self,
        mo_energies: np.ndarray,
        mo_occ: np.ndarray,
        sigma_x: np.ndarray,
        sigma_c: Union[np.ndarray, Callable],
        vxc_diagonal: np.ndarray,
        initial_guess: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """Run QP solver (PySCF-style kernel method).

        Args:
            mo_energies: MO energies from DFT/HF (n_mo,)
            mo_occ: MO occupations (n_mo,)
            sigma_x: Exchange self-energy matrix (n_mo, n_mo) or diagonal (n_mo,)
            sigma_c: Correlation self-energy - either:
                - Static diagonal array (n_mo,) for simplified calculation
                - Callable(orbital_idx, energy) -> complex for full energy dependence
            vxc_diagonal: DFT exchange-correlation diagonal (n_mo,)
            initial_guess: Initial guess for QP energies (default: mo_energies)

        Returns:
            Quasiparticle energies array

        Raises:
            RuntimeError: If Rust extension is not available
            ValueError: If input dimensions are inconsistent
        """
        if not RUST_AVAILABLE:
            raise RuntimeError("QuasiX Rust extension not available. Please rebuild.")

        # Input validation
        n_mo = len(mo_energies)
        if len(mo_occ) != n_mo:
            raise ValueError(f"mo_occ length {len(mo_occ)} != mo_energies length {n_mo}")
        if len(vxc_diagonal) != n_mo:
            raise ValueError(f"vxc_diagonal length {len(vxc_diagonal)} != n_mo")

        # Process sigma_x
        if sigma_x.ndim == 2:
            if sigma_x.shape != (n_mo, n_mo):
                raise ValueError(f"sigma_x shape {sigma_x.shape} != expected ({n_mo}, {n_mo})")
            sigma_x_matrix = sigma_x
        elif sigma_x.ndim == 1:
            if len(sigma_x) != n_mo:
                raise ValueError(f"sigma_x length {len(sigma_x)} != n_mo {n_mo}")
            # Convert diagonal to matrix
            sigma_x_matrix = np.diag(sigma_x)
        else:
            raise ValueError(f"sigma_x must be 1D or 2D, got shape {sigma_x.shape}")

        # Process sigma_c
        if callable(sigma_c):
            # Dynamic correlation self-energy - use simplified diagonal for Rust interface
            # The Rust code will apply energy dependence internally
            sigma_c_diagonal = np.zeros(n_mo)
            for i in range(n_mo):
                # Evaluate at reference energy
                sigma_c_val = sigma_c(i, mo_energies[i])
                sigma_c_diagonal[i] = np.real(sigma_c_val)
        else:
            # Static correlation self-energy
            if len(sigma_c) != n_mo:
                raise ValueError(f"sigma_c length {len(sigma_c)} != n_mo {n_mo}")
            sigma_c_diagonal = np.asarray(sigma_c)

        # Log solver start
        if self.verbose > 0:
            logger.info("="*60)
            logger.info("QuasiX QP Solver (S3-6)")
            logger.info("="*60)
            logger.info(f"Number of orbitals: {n_mo}")
            logger.info(f"Number of occupied: {np.sum(mo_occ > 0.5)}")
            logger.info(f"Energy tolerance: {self.config.energy_tolerance:.2e} Ha")
            logger.info(f"Max Newton iterations: {self.config.max_newton_iterations}")
            logger.info(f"Energy model: {self.config.energy_model}")

        # Call Rust solver
        result = _solve_quasiparticle_equations(
            mo_energies,
            mo_occ,
            sigma_x_matrix,
            sigma_c_diagonal,
            vxc_diagonal,
            self.config.to_dict()
        )

        # Extract results
        self.qp_energies = result['qp_energies']
        self.z_factors = result['z_factors']
        self.diagnostics = result.get('diagnostics', {})
        self.orbital_solutions = result.get('orbital_solutions', [])
        self._last_result = result

        # Check convergence
        convergence_flags = result.get('convergence_flags', np.ones(n_mo, dtype=bool))
        self.converged = np.all(convergence_flags)

        # Log results
        if self.verbose > 0:
            self._log_results()

        # Call progress callback if set
        if self._progress_callback:
            self._progress_callback(self)

        return self.qp_energies

    def _log_results(self):
        """Log solver results."""
        if not self.diagnostics:
            return

        logger.info("-"*60)
        logger.info("QP Solver Results:")
        logger.info(f"  Converged: {self.converged}")
        logger.info(f"  Converged orbitals: {self.diagnostics.get('converged_count', 0)}/{len(self.qp_energies)}")
        logger.info(f"  Average iterations: {self.diagnostics.get('avg_iterations', 0):.1f}")
        logger.info(f"  Z-factor range: ({self.diagnostics.get('z_factor_min', 0):.3f}, "
                   f"{self.diagnostics.get('z_factor_max', 0):.3f})")
        logger.info(f"  Average energy shift: {self.diagnostics.get('avg_energy_shift', 0):.4f} Ha")

        if self.diagnostics.get('bisection_count', 0) > 0:
            logger.info(f"  Bisection fallbacks: {self.diagnostics['bisection_count']}")

        # Warn about problematic Z-factors
        z_min = self.diagnostics.get('z_factor_min', 1.0)
        z_max = self.diagnostics.get('z_factor_max', 0.0)
        if z_min < 0.1 or z_max > 0.9:
            logger.warning(f"  ⚠ Problematic Z-factors detected (outside 0.1-0.9 range)")
            logger.warning(f"    Consider vertex corrections or higher-level theory")

    def get_qp_corrections(self, mo_energies: Optional[np.ndarray] = None) -> np.ndarray:
        """Get QP corrections ΔE = E_QP - E_KS.

        Args:
            mo_energies: Reference energies (uses last solve if None)

        Returns:
            Array of QP corrections
        """
        if self.qp_energies is None:
            raise RuntimeError("No QP solution available. Run kernel() first.")

        if mo_energies is None:
            # Try to get from last result
            if self._last_result and 'mo_energies' in self._last_result:
                mo_energies = self._last_result['mo_energies']
            else:
                raise ValueError("mo_energies required when not available from last solve")

        return self.qp_energies - mo_energies

    def get_spectral_function(
        self,
        omega_grid: Optional[np.ndarray] = None,
        broadening: float = 0.1,
        units: str = 'eV'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute spectral function A(ω).

        Args:
            omega_grid: Frequency grid (default: auto around QP energies)
            broadening: Lorentzian broadening (eV)
            units: Energy units ('Ha' or 'eV')

        Returns:
            (omega_grid, spectral_function) tuple
        """
        if self.qp_energies is None or self.z_factors is None:
            raise RuntimeError("No QP solution available. Run kernel() first.")

        # Unit conversion
        if units == 'eV':
            factor = 27.211  # Ha to eV
        else:
            factor = 1.0

        qp_energies_scaled = self.qp_energies * factor

        if omega_grid is None:
            # Auto-generate grid
            e_min = np.min(qp_energies_scaled) - 5.0
            e_max = np.max(qp_energies_scaled) + 5.0
            omega_grid = np.linspace(e_min, e_max, 1000)

        spectral = np.zeros_like(omega_grid)

        # Add quasiparticle peaks
        for e_qp, z in zip(qp_energies_scaled, self.z_factors):
            spectral += z * broadening / ((omega_grid - e_qp)**2 + broadening**2)

        # Add satellite contributions (simplified)
        for e_qp, z in zip(qp_energies_scaled, self.z_factors):
            if z < 0.95:  # Only if significant correlation
                satellite_weight = 1 - z
                e_sat = e_qp - 2.0 * factor  # Simplified satellite position
                spectral += satellite_weight * (2*broadening) / ((omega_grid - e_sat)**2 + (2*broadening)**2)

        # Normalize
        spectral *= 1.0 / np.pi

        return omega_grid, spectral

    def set_progress_callback(self, callback: Callable):
        """Set callback for progress monitoring (e.g., in Jupyter).

        Args:
            callback: Function called with solver instance after each solve
        """
        self._progress_callback = callback

    # PySCF-style method chaining
    def set_config(self, **kwargs) -> 'QPSolver':
        """Update configuration parameters.

        Args:
            **kwargs: Configuration parameters to update

        Returns:
            Self for method chaining
        """
        for key, value in kwargs.items():
            if hasattr(self.config, key):
                setattr(self.config, key, value)
            else:
                raise ValueError(f"Unknown configuration parameter: {key}")
        return self

    def set_tolerance(self, tol: float) -> 'QPSolver':
        """Set energy tolerance.

        Args:
            tol: Energy tolerance in Ha

        Returns:
            Self for method chaining
        """
        self.config.energy_tolerance = tol
        self.config.residual_tolerance = tol * 0.1
        return self

    def set_max_iterations(self, max_iter: int) -> 'QPSolver':
        """Set maximum iterations.

        Args:
            max_iter: Maximum Newton-Raphson iterations

        Returns:
            Self for method chaining
        """
        self.config.max_newton_iterations = max_iter
        return self


class GraphicalQPSolver:
    """Graphical solution of QP equation.
    
    Finds intersection: ω - ε_DFT - Re[Σ(ω) - v_xc] = 0
    """
    
    def __init__(self, config: Optional[QPSolverConfig] = None):
        """Initialize solver."""
        self.config = config or QPSolverConfig(method='graphical')
    
    def solve(
        self,
        energy_dft: float,
        sigma_x: float,
        sigma_c_func: callable,
        vxc: float,
        omega_range: Optional[Tuple[float, float]] = None,
        orbital_index: int = 0
    ) -> QuasiparticleState:
        """Find QP energy by graphical solution.
        
        Args:
            energy_dft: DFT orbital energy
            sigma_x: Exchange self-energy
            sigma_c_func: Correlation self-energy function
            vxc: XC potential
            omega_range: Search range for QP energy
            orbital_index: Orbital index
            
        Returns:
            QuasiparticleState
        """
        if omega_range is None:
            # Default search range around DFT energy
            omega_range = (energy_dft - 5.0, energy_dft + 5.0)
        
        # Define function to find root
        def qp_func(omega):
            sigma_c = sigma_c_func(omega)
            return omega - energy_dft - (sigma_x + np.real(sigma_c) - vxc)
        
        # Bisection method for robustness
        omega_left, omega_right = omega_range
        f_left = qp_func(omega_left)
        f_right = qp_func(omega_right)
        
        # Check if root is bracketed
        if f_left * f_right > 0:
            # Try to find by scanning
            omega_scan = np.linspace(omega_left, omega_right, 100)
            f_scan = [qp_func(w) for w in omega_scan]
            
            # Find sign changes
            sign_changes = np.where(np.diff(np.sign(f_scan)))[0]
            if len(sign_changes) > 0:
                # Take first sign change
                idx = sign_changes[0]
                omega_left = omega_scan[idx]
                omega_right = omega_scan[idx + 1]
            else:
                # No solution found, return DFT energy
                logger.warning(f"No QP solution found for orbital {orbital_index}")
                return QuasiparticleState(
                    index=orbital_index,
                    energy_dft=energy_dft,
                    energy_qp=energy_dft,
                    z_factor=1.0,
                    sigma_x=sigma_x,
                    sigma_c=sigma_c_func(energy_dft),
                    vxc=vxc,
                    converged=False
                )
        
        # Bisection iteration
        for _ in range(self.config.max_iter):
            omega_mid = 0.5 * (omega_left + omega_right)
            f_mid = qp_func(omega_mid)
            
            if abs(f_mid) < self.config.tol_energy:
                # Converged, compute Z-factor
                delta = 1e-4
                sigma_c_plus = sigma_c_func(omega_mid + delta)
                sigma_c_minus = sigma_c_func(omega_mid - delta)
                dsigma_dw = (sigma_c_plus - sigma_c_minus) / (2 * delta)
                z_factor = 1.0 / (1.0 - np.real(dsigma_dw))
                z_factor = np.clip(z_factor, 0.01, 0.99)
                
                return QuasiparticleState(
                    index=orbital_index,
                    energy_dft=energy_dft,
                    energy_qp=omega_mid,
                    z_factor=z_factor,
                    sigma_x=sigma_x,
                    sigma_c=sigma_c_func(omega_mid),
                    vxc=vxc,
                    converged=True
                )
            
            # Update bracket
            if f_mid * f_left < 0:
                omega_right = omega_mid
            else:
                omega_left = omega_mid
        
        # Not converged
        omega_final = 0.5 * (omega_left + omega_right)
        return QuasiparticleState(
            index=orbital_index,
            energy_dft=energy_dft,
            energy_qp=omega_final,
            z_factor=1.0,
            sigma_x=sigma_x,
            sigma_c=sigma_c_func(omega_final),
            vxc=vxc,
            converged=False
        )


def solve_quasiparticle_equations(
    energies_dft: np.ndarray,
    sigma_x: np.ndarray,
    sigma_c_func: callable,
    vxc: np.ndarray,
    config: Optional[QPSolverConfig] = None
) -> List[QuasiparticleState]:
    """Main interface for solving quasiparticle equations.
    
    Args:
        energies_dft: DFT orbital energies
        sigma_x: Exchange self-energies
        sigma_c_func: Function returning Σ^c(i, ω)
        vxc: XC potentials
        config: Solver configuration
        
    Returns:
        List of QuasiparticleState objects
    """
    config = config or QPSolverConfig()
    
    if RUST_AVAILABLE and config.method == 'linearized':
        # Use Rust implementation if available
        try:
            results = _solve_quasiparticle_equations(
                energies_dft.tolist(),
                sigma_x.tolist(),
                sigma_c_func,  # Need wrapper for Rust
                vxc.tolist(),
                config.__dict__
            )
            return [QuasiparticleState(**r) for r in results]
        except:
            pass  # Fall back to Python
    
    # Use Python implementation
    if config.method == 'linearized':
        solver = LinearizedQPSolver(config)
    elif config.method == 'graphical':
        solver = GraphicalQPSolver(config)
    else:
        raise ValueError(f"Unknown solver method: {config.method}")
    
    return solver.solve_batch(energies_dft, sigma_x, sigma_c_func, vxc)


def compute_qp_gap(
    qp_states: List[QuasiparticleState],
    homo_index: int
) -> float:
    """Compute quasiparticle gap.
    
    Args:
        qp_states: List of QP states
        homo_index: Index of HOMO orbital
        
    Returns:
        QP gap (LUMO - HOMO)
    """
    if homo_index >= len(qp_states) - 1:
        raise ValueError("Invalid HOMO index")
    
    homo_qp = qp_states[homo_index].energy_qp
    lumo_qp = qp_states[homo_index + 1].energy_qp
    
    return lumo_qp - homo_qp


def compute_ionization_potential(
    qp_states: List[QuasiparticleState],
    homo_index: int
) -> float:
    """Compute ionization potential (negative of HOMO QP energy).
    
    Args:
        qp_states: List of QP states
        homo_index: HOMO index
        
    Returns:
        Ionization potential
    """
    return -qp_states[homo_index].energy_qp


def compute_electron_affinity(
    qp_states: List[QuasiparticleState],
    homo_index: int
) -> float:
    """Compute electron affinity (negative of LUMO QP energy).
    
    Args:
        qp_states: List of QP states
        homo_index: HOMO index
        
    Returns:
        Electron affinity
    """
    return -qp_states[homo_index + 1].energy_qp


def validate_z_factors(
    qp_states: List[QuasiparticleState],
    tol: float = 0.01
) -> bool:
    """Validate that all Z-factors are in physical range.
    
    Args:
        qp_states: List of QP states
        tol: Tolerance for Z ≈ 0 or Z ≈ 1
        
    Returns:
        True if all Z-factors are valid
    """
    for qp in qp_states:
        if qp.z_factor < tol or qp.z_factor > 1.0 - tol:
            logger.warning(f"Unphysical Z-factor {qp.z_factor:.3f} for orbital {qp.index}")
            return False
    return True


# Module exports
__all__ = [
    'QuasiparticleState',
    'QPSolverConfig',
    'LinearizedQPSolver',
    'GraphicalQPSolver',
    'solve_quasiparticle_equations',
    'compute_qp_gap',
    'compute_ionization_potential',
    'compute_electron_affinity',
    'validate_z_factors',
]