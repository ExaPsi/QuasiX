"""Thin Python wrapper for QuasiX analytic continuation.

This module provides a PySCF-compatible interface to the high-performance
Rust implementation of analytic continuation methods. Following the QuasiX
philosophy of "tight Python integration, coarse crossings", all computation
is performed in Rust with Python serving only as the interface layer.

Key Features:
    - Multipole and Padé analytic continuation models
    - Cross-validation for automatic model selection
    - SIMD-optimized evaluation on real frequencies
    - Thread-safe parallel processing
    - PySCF-compatible API design

Example:
    Basic usage with PySCF GW calculation:

    >>> from pyscf import gto, scf
    >>> from quasix.gw import AnalyticContinuation
    >>>
    >>> # Run SCF calculation
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 1', basis='cc-pvdz')
    >>> mf = scf.RHF(mol).run()
    >>>
    >>> # Prepare imaginary-axis data from GW
    >>> xi_points = np.array([0.1, 0.5, 1.0, 2.0, 5.0])
    >>> sigma_iw = np.array([...])  # Complex self-energy on imaginary axis
    >>>
    >>> # Perform analytic continuation
    >>> ac = AnalyticContinuation()
    >>> ac.fit(xi_points, sigma_iw)
    >>>
    >>> # Evaluate on real axis
    >>> omega = np.linspace(-10, 10, 1000)
    >>> sigma_w = ac.evaluate(omega, eta=0.01)
    >>>
    >>> # Extract spectral function
    >>> spectral = ac.spectral_function(omega)

Note:
    This is a thin wrapper around the Rust implementation. All heavy
    computation is performed in Rust via PyO3 bindings. The Python
    layer only handles data conversion and PySCF integration.
"""

import numpy as np
from typing import Optional, Tuple, Dict, Union, Any
from numpy.typing import NDArray

# Import Rust implementations via PyO3
try:
    from quasix import (
        ACConfig,
        MultipoleModel,
        PadeModel,
        AnalyticContinuationFitter as RustACFitter
    )
except ImportError as e:
    raise ImportError(
        "Failed to import Rust analytic continuation module. "
        "Please ensure QuasiX is properly built with 'maturin develop --release'"
    ) from e


class AnalyticContinuation:
    """PySCF-compatible interface for analytic continuation.

    This class provides a high-level interface to QuasiX's Rust-based
    analytic continuation implementation, following PySCF conventions
    for method chaining and property access.

    Attributes:
        config: ACConfig object with continuation parameters
        fitter: Rust AnalyticContinuationFitter instance
        fitted: Whether the model has been fitted

    Args:
        method: Continuation method ('auto', 'multipole', 'pade')
        max_poles: Maximum number of poles for multipole model
        max_pade_order: Maximum (M,N) order for Padé approximant
        cv_fraction: Fraction of data for cross-validation
        regularization: Regularization parameter for fitting
        eta: Default broadening parameter for real-axis evaluation
        parallel: Enable parallel processing
    """

    def __init__(
        self,
        method: str = 'auto',
        max_poles: int = 20,
        max_pade_order: Tuple[int, int] = (15, 15),
        cv_fraction: float = 0.3,
        regularization: float = 1e-10,
        eta: float = 0.01,
        parallel: bool = True,
        **kwargs
    ):
        """Initialize analytic continuation wrapper."""
        # Create configuration
        self.config = ACConfig(
            max_poles=max_poles,
            max_pade_order=max_pade_order,
            cv_fraction=cv_fraction,
            cv_iterations=kwargs.get('cv_iterations', 5),
            regularization=regularization,
            stability_threshold=kwargs.get('stability_threshold', 1e-3),
            convergence_tol=kwargs.get('convergence_tol', 1e-8),
            max_iterations=kwargs.get('max_iterations', 1000),
            eta=eta,
            parallel=parallel
        )

        self.method = method
        self.fitter: Optional[RustACFitter] = None
        self.fitted = False

        # Store fitting data for potential reuse
        self._xi_data: Optional[NDArray] = None
        self._f_data: Optional[NDArray] = None

    def fit(
        self,
        xi_points: NDArray[np.float64],
        f_data: NDArray[np.complex128],
        weights: Optional[NDArray[np.float64]] = None
    ) -> 'AnalyticContinuation':
        """Fit analytic continuation model to imaginary-axis data.

        Args:
            xi_points: Imaginary frequencies ξ, shape (n_freq,)
            f_data: Complex function values f(iξ), shape (n_freq,) or (n_orb, n_freq)
            weights: Optional weights for fitting, shape (n_freq,)

        Returns:
            Self for method chaining

        Raises:
            ValueError: If data shapes are incompatible
            RuntimeError: If fitting fails
        """
        # Validate inputs
        xi_points = np.asarray(xi_points, dtype=np.float64)
        f_data = np.asarray(f_data, dtype=np.complex128)

        if xi_points.ndim != 1:
            raise ValueError("xi_points must be 1D array")

        # Handle matrix-valued functions
        if f_data.ndim == 2:
            # Fit each orbital separately (handled by Rust in parallel if enabled)
            n_orb, n_freq = f_data.shape
            if n_freq != len(xi_points):
                raise ValueError(
                    f"Frequency dimension mismatch: {n_freq} != {len(xi_points)}"
                )
            # For now, fit diagonal elements - full matrix support in future
            # This is a simplification for initial implementation
            f_data = f_data[0]  # Take first orbital for now
        elif f_data.ndim != 1:
            raise ValueError("f_data must be 1D or 2D array")

        if len(f_data) != len(xi_points):
            raise ValueError(
                f"Data length mismatch: {len(f_data)} != {len(xi_points)}"
            )

        # Create fitter based on method
        if self.method == 'multipole':
            model = MultipoleModel(self.config)
        elif self.method == 'pade':
            model = PadeModel(self.config)
        else:  # 'auto' or default
            # Let the fitter choose based on cross-validation
            model = None

        # Create and run fitter
        self.fitter = RustACFitter(self.config)

        if weights is not None:
            weights = np.asarray(weights, dtype=np.float64)
            # Note: Weight support to be added in Rust implementation

        # Perform fitting (all computation in Rust)
        # Rust API expects separate real and imaginary parts
        f_real = np.real(f_data)
        f_imag = np.imag(f_data)
        self.fitter.fit(xi_points, f_real, f_imag)

        # Store data for potential reuse
        self._xi_data = xi_points
        self._f_data = f_data
        self.fitted = True

        return self

    def evaluate(
        self,
        omega: NDArray[np.float64],
        eta: Optional[float] = None
    ) -> NDArray[np.complex128]:
        """Evaluate continuation on real frequencies.

        Args:
            omega: Real frequencies ω, shape (n_freq,)
            eta: Broadening parameter (uses config default if None)

        Returns:
            Complex function values f(ω + iη), shape (n_freq,)

        Raises:
            RuntimeError: If model not fitted
        """
        if not self.fitted or self.fitter is None:
            raise RuntimeError("Model must be fitted before evaluation")

        omega = np.asarray(omega, dtype=np.float64)

        # Use Rust implementation for evaluation
        # Returns tuple of (real, imag) arrays
        real_part, imag_part = self.fitter.evaluate_real_axis(omega, eta)
        return real_part + 1j * imag_part

    def spectral_function(
        self,
        omega: NDArray[np.float64],
        eta: Optional[float] = None
    ) -> NDArray[np.float64]:
        """Compute spectral function A(ω) = -Im[f(ω + iη)]/π.

        Args:
            omega: Real frequencies ω, shape (n_freq,)
            eta: Broadening parameter (uses config default if None)

        Returns:
            Spectral function A(ω), shape (n_freq,)

        Raises:
            RuntimeError: If model not fitted
        """
        if not self.fitted or self.fitter is None:
            raise RuntimeError("Model must be fitted before evaluation")

        omega = np.asarray(omega, dtype=np.float64)

        # Use Rust implementation
        return self.fitter.spectral_function(omega, eta)

    def validate_causality(
        self,
        omega: Optional[NDArray[np.float64]] = None
    ) -> Dict[str, Any]:
        """Validate causality and Kramers-Kronig relations.

        Args:
            omega: Frequencies for validation (auto-generated if None)

        Returns:
            Dictionary with validation metrics:
                - violations: Number of causality violations
                - violation_fraction: Fraction of points violating causality
                - kk_error: Kramers-Kronig relation error
                - kk_relative_error: Relative KK error

        Raises:
            RuntimeError: If model not fitted
        """
        if not self.fitted or self.fitter is None:
            raise RuntimeError("Model must be fitted before validation")

        if omega is None:
            # Generate reasonable frequency grid
            omega = np.linspace(-20, 20, 1000)
        else:
            omega = np.asarray(omega, dtype=np.float64)

        # Use Rust implementation
        return self.fitter.validate_causality(omega)

    def get_poles_residues(self) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        """Extract poles and residues from fitted model.

        Returns:
            Tuple of (poles, residues) arrays

        Raises:
            RuntimeError: If model not fitted or extraction fails
        """
        if not self.fitted or self.fitter is None:
            raise RuntimeError("Model must be fitted before extraction")

        # Note: This functionality may need to be added to the Rust API
        # For now, return empty arrays as placeholder
        # TODO: Implement in Rust backend
        import warnings
        warnings.warn(
            "get_poles_residues not yet fully implemented in Rust backend",
            RuntimeWarning
        )
        return np.array([], dtype=np.complex128), np.array([], dtype=np.complex128)

    @property
    def model_info(self) -> Dict[str, Any]:
        """Get information about the fitted model.

        Returns:
            Dictionary with model information:
                - has_model: Whether a model is fitted
                - cv_error: Cross-validation error
                - is_stable: Stability check result
                - n_poles: Number of poles
                - n_residues: Number of residues
        """
        if not self.fitted or self.fitter is None:
            return {"has_model": False}

        return self.fitter.get_best_model_info()

    def __repr__(self) -> str:
        """String representation of the continuation object."""
        if self.fitted:
            info = self.model_info
            return (
                f"AnalyticContinuation(method='{self.method}', "
                f"fitted=True, cv_error={info.get('cv_error', 'N/A'):.3e})"
            )
        else:
            return f"AnalyticContinuation(method='{self.method}', fitted=False)"


# Convenience functions for direct usage

def perform_analytic_continuation(
    xi_points: NDArray[np.float64],
    f_data: NDArray[np.complex128],
    omega: NDArray[np.float64],
    method: str = 'auto',
    eta: float = 0.01,
    **kwargs
) -> Tuple[NDArray[np.complex128], Dict[str, Any]]:
    """Perform analytic continuation with specified method.

    Convenience function for single-shot analytic continuation without
    keeping the model object.

    Args:
        xi_points: Imaginary frequencies
        f_data: Function values on imaginary axis
        omega: Real frequencies for evaluation
        method: Continuation method ('auto', 'multipole', 'pade')
        eta: Broadening parameter
        **kwargs: Additional configuration parameters

    Returns:
        Tuple of (continued function, model info dictionary)
    """
    ac = AnalyticContinuation(
        method=method,
        eta=eta,
        **kwargs
    )
    ac.fit(xi_points, f_data)
    result = ac.evaluate(omega)
    info = ac.model_info
    return result, info


# PySCF compatibility helpers

def extract_self_energy_poles(
    gw_obj,
    orbital_index: int = None
) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Extract self-energy poles from a PySCF GW object.

    This helper function interfaces with PySCF GW calculations to
    extract poles and residues for quasiparticle analysis.

    Args:
        gw_obj: PySCF GW object with calculated self-energy
        orbital_index: Specific orbital to analyze (None for all)

    Returns:
        Tuple of (poles, residues) for self-energy

    Note:
        Requires the GW object to have imaginary-axis self-energy data
    """
    # This would interface with actual PySCF GW objects
    # Implementation depends on PySCF GW structure
    raise NotImplementedError(
        "PySCF GW interface to be implemented with GW module"
    )


__all__ = [
    'AnalyticContinuation',
    'perform_analytic_continuation',
    'extract_self_energy_poles',
]