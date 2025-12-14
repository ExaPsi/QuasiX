"""
Frequency grid utilities for QuasiX GW/BSE calculations.

This module provides frequency grid generation and transformation utilities
for contour deformation (CD) and analytical continuation (AC) methods.
"""

from typing import Optional, Tuple, Dict, Any, Union
import numpy as np
from dataclasses import dataclass

# Temporarily use mock implementations until Rust bindings are ready
# from .quasix import (
#     create_frequency_grid as _create_frequency_grid,
#     gauss_legendre_quadrature as _gauss_legendre_quadrature,
#     create_contour_deformation as _create_contour_deformation,
#     create_ac_fitter as _create_ac_fitter,
# )

# Mock implementations for testing
def _create_frequency_grid(grid_type, n_points, bounds=None, omega_max=None):
    """Mock frequency grid creation with proper quadrature rules."""
    if grid_type == 'gauss_legendre':
        # Use proper GL quadrature for Gauss-Legendre grids
        from numpy.polynomial import legendre
        
        # Get GL nodes and weights on [0, omega_max]
        nodes_std, weights_std = legendre.leggauss(n_points)
        omega_max = omega_max or 100.0
        
        # Transform from [-1, 1] to [0, omega_max]
        points = 0.5 * omega_max * (nodes_std + 1.0)
        weights = 0.5 * omega_max * weights_std
    else:
        # For other grids, use uniform spacing (temporary)
        if bounds:
            points = np.linspace(bounds[0], bounds[1], n_points)
        else:
            points = np.linspace(-10, 10, n_points)
        weights = np.ones(n_points) * ((bounds[1] - bounds[0]) / n_points if bounds else 20.0 / n_points)
    
    result = {
        'points': points,
        'weights': weights,
        'grid_type': grid_type,
        'n_points': n_points
    }
    if omega_max:
        result['omega_max'] = omega_max
    if bounds:
        result['bounds'] = bounds
    return result

def _gauss_legendre_quadrature(n_points, a=-1.0, b=1.0):
    """Generate proper Gauss-Legendre quadrature nodes and weights."""
    # Use numpy's polynomial module for actual GL quadrature
    from numpy.polynomial import legendre
    
    # Get nodes and weights on [-1, 1]
    nodes_std, weights_std = legendre.leggauss(n_points)
    
    # Transform to interval [a, b]
    nodes = 0.5 * (b - a) * nodes_std + 0.5 * (b + a)
    weights = 0.5 * (b - a) * weights_std
    
    return {'nodes': nodes, 'weights': weights, 'interval': [a, b]}

def _create_contour_deformation(energy_min, energy_max, n_points):
    """Mock contour deformation."""
    return {
        'energy_min': energy_min,
        'energy_max': energy_max,
        'n_points': n_points,
        'method': 'contour_deformation'
    }

def _create_ac_fitter(n_poles, regularization=1e-8):
    """Mock AC fitter."""
    return {
        'n_poles': n_poles,
        'regularization': regularization,
        'method': 'analytical_continuation'
    }


@dataclass
class FrequencyGrid:
    """Frequency grid for GW calculations.
    
    Attributes:
        points: Frequency points (real or imaginary)
        weights: Quadrature weights
        grid_type: Type of grid ('gauss_legendre', 'contour_deformation', 'minimax')
        n_points: Number of frequency points
        omega_max: Maximum frequency for imaginary axis grids
        bounds: Energy bounds for contour deformation
    """
    points: np.ndarray
    weights: np.ndarray
    grid_type: str
    n_points: int
    omega_max: Optional[float] = None
    bounds: Optional[Tuple[float, float]] = None
    
    # Provide frequencies as an alias for points (for compatibility)
    @property
    def frequencies(self) -> np.ndarray:
        """Alias for points for compatibility."""
        return self.points
    
    def __repr__(self) -> str:
        """String representation of frequency grid."""
        info = f"FrequencyGrid(type='{self.grid_type}', n_points={self.n_points}"
        if self.omega_max is not None:
            info += f", omega_max={self.omega_max:.1f}"
        if self.bounds is not None:
            info += f", bounds=({self.bounds[0]:.1f}, {self.bounds[1]:.1f})"
        info += ")"
        return info
    
    def integrate(self, func: np.ndarray) -> float:
        """Integrate function values using quadrature weights.
        
        Args:
            func: Function values at grid points
            
        Returns:
            Integral value
        """
        if len(func) != self.n_points:
            raise ValueError(f"Function array size {len(func)} doesn't match grid size {self.n_points}")
        return np.sum(self.weights * func)
    
    def to_imaginary(self) -> 'FrequencyGrid':
        """Transform grid to imaginary axis.
        
        Returns:
            New FrequencyGrid on imaginary axis
        """
        if self.grid_type != 'gauss_legendre':
            raise ValueError("Only Gauss-Legendre grids can be transformed to imaginary axis")
        
        # Transform points to imaginary axis: ω → iω
        imag_points = 1j * self.points
        
        return FrequencyGrid(
            points=imag_points,
            weights=self.weights.copy(),
            grid_type='imaginary_axis',
            n_points=self.n_points,
            omega_max=self.omega_max,
            bounds=None
        )
    
    @classmethod
    def gauss_legendre(cls, n_points: int, omega_max: float = 100.0) -> 'FrequencyGrid':
        """Create Gauss-Legendre frequency grid.
        
        Args:
            n_points: Number of quadrature points
            omega_max: Maximum frequency for grid
            
        Returns:
            FrequencyGrid with GL quadrature
        """
        return create_frequency_grid('gauss_legendre', n_points, omega_max=omega_max)
    
    @classmethod
    def contour_deformation(cls, n_points: int, bounds: Tuple[float, float]) -> 'FrequencyGrid':
        """Create contour deformation frequency grid.
        
        Args:
            n_points: Number of contour points
            bounds: Energy bounds (min, max)
            
        Returns:
            FrequencyGrid for contour deformation
        """
        return create_frequency_grid('contour_deformation', n_points, bounds=bounds)
    
    @classmethod
    def imaginary_axis(cls, n_points: int, omega_max: float = None, xi_max: float = None) -> 'FrequencyGrid':
        """Create imaginary axis frequency grid.
        
        Args:
            n_points: Number of points
            omega_max: Maximum imaginary frequency (or use xi_max)
            xi_max: Alternative name for omega_max
            
        Returns:
            FrequencyGrid on imaginary axis
        """
        # Support both omega_max and xi_max for compatibility
        if xi_max is not None:
            omega_max = xi_max
        elif omega_max is None:
            omega_max = 100.0
        return create_frequency_grid('imaginary_axis', n_points, omega_max=omega_max)


def create_frequency_grid(
    grid_type: str,
    n_points: int,
    bounds: Optional[Tuple[float, float]] = None,
    omega_max: Optional[float] = None
) -> FrequencyGrid:
    """Create a frequency grid for GW calculations.
    
    Args:
        grid_type: Type of grid to create
            - 'gauss_legendre': Gauss-Legendre quadrature on imaginary axis
            - 'contour_deformation': Contour deformation in complex plane
            - 'minimax': Minimax grid for analytical continuation
            - 'imaginary_axis': Grid on imaginary axis (returns real ξ values)
        n_points: Number of frequency points
        bounds: Optional bounds (min, max) for contour deformation
        omega_max: Optional maximum frequency for imaginary axis grids
        
    Returns:
        FrequencyGrid object containing points and weights
        
    Examples:
        >>> # Create Gauss-Legendre grid for imaginary axis
        >>> grid = create_frequency_grid('gauss_legendre', 32, omega_max=100.0)
        >>> print(f"Grid has {grid.n_points} points from 0 to {grid.omega_max}")
        
        >>> # Create contour deformation grid
        >>> grid = create_frequency_grid('contour_deformation', 64, bounds=(-20.0, 20.0))
        >>> print(f"Contour spans energy range {grid.bounds}")
        
        >>> # Integrate a function
        >>> func_values = np.exp(-grid.points**2)
        >>> integral = grid.integrate(func_values)
    """
    # Validate grid type
    valid_types = ['gauss_legendre', 'contour_deformation', 'minimax', 'imaginary_axis']
    if grid_type not in valid_types:
        raise ValueError(f"Invalid grid_type: {grid_type}. Must be one of {valid_types}")
    
    # Call Rust implementation
    result = _create_frequency_grid(grid_type, n_points, bounds, omega_max)
    
    # Extract results
    points = np.array(result['points'])
    weights = np.array(result['weights'])
    
    # Handle imaginary axis transformation
    if grid_type == 'imaginary_axis':
        # For imaginary axis, return real ξ values (where ω = iξ)
        # The points are already real from the GL grid
        base_grid = _create_frequency_grid('gauss_legendre', n_points, bounds, omega_max)
        points = np.array(base_grid['points'])  # Real ξ values
        weights = np.array(base_grid['weights'])
    
    return FrequencyGrid(
        points=points,
        weights=weights,
        grid_type=grid_type,
        n_points=n_points,
        omega_max=result.get('omega_max'),
        bounds=tuple(result['bounds']) if 'bounds' in result else None
    )


def gauss_legendre_quadrature(
    n_points: int,
    interval: Tuple[float, float] = (-1.0, 1.0)
) -> Dict[str, np.ndarray]:
    """Generate Gauss-Legendre quadrature nodes and weights.
    
    Args:
        n_points: Number of quadrature points
        interval: Integration interval (a, b)
        
    Returns:
        Dictionary with 'nodes' and 'weights' arrays
        
    Examples:
        >>> # Standard GL quadrature on [-1, 1]
        >>> quad = gauss_legendre_quadrature(16)
        >>> nodes, weights = quad['nodes'], quad['weights']
        
        >>> # GL quadrature on custom interval
        >>> quad = gauss_legendre_quadrature(32, interval=(0.0, 100.0))
        >>> integral = np.sum(weights * np.exp(-nodes))
    """
    a, b = interval
    result = _gauss_legendre_quadrature(n_points, a, b)
    
    return {
        'nodes': np.array(result['nodes']),
        'weights': np.array(result['weights']),
        'interval': interval
    }


@dataclass
class ContourDeformation:
    """Contour deformation calculator for correlation self-energy.
    
    Attributes:
        energy_range: Energy range (min, max) for contour
        n_points: Number of contour points
        method: Integration method ('contour_deformation')
    """
    energy_range: Tuple[float, float]
    n_points: int
    method: str = 'contour_deformation'
    
    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> 'ContourDeformation':
        """Create from configuration dictionary."""
        return cls(
            energy_range=(config['energy_min'], config['energy_max']),
            n_points=config['n_points'],
            method=config.get('method', 'contour_deformation')
        )
    
    def get_contour_points(self, eta: float = 0.01) -> np.ndarray:
        """Get complex contour points.
        
        Args:
            eta: Imaginary shift for convergence
            
        Returns:
            Array of complex contour points
        """
        e_min, e_max = self.energy_range
        real_points = np.linspace(e_min, e_max, self.n_points)
        return real_points + 1j * eta


def create_contour_deformation(
    energy_min: float,
    energy_max: float,
    n_points: int
) -> ContourDeformation:
    """Create contour deformation calculator.
    
    Args:
        energy_min: Minimum energy for contour
        energy_max: Maximum energy for contour
        n_points: Number of contour points
        
    Returns:
        ContourDeformation object
        
    Examples:
        >>> # Create contour for valence states
        >>> cd = create_contour_deformation(-20.0, 0.0, 100)
        >>> contour_points = cd.get_contour_points(eta=0.01)
        >>> print(f"Contour has {len(contour_points)} points")
    """
    config = _create_contour_deformation(energy_min, energy_max, n_points)
    return ContourDeformation.from_config(config)


@dataclass
class ACFitter:
    """Analytical continuation fitter using Padé approximants.
    
    Attributes:
        n_poles: Number of poles in rational approximation
        regularization: Regularization parameter for fitting
        method: Fitting method ('analytical_continuation')
    """
    n_poles: int
    regularization: float = 1e-8
    method: str = 'analytical_continuation'
    
    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> 'ACFitter':
        """Create from configuration dictionary."""
        return cls(
            n_poles=config['n_poles'],
            regularization=config['regularization'],
            method=config.get('method', 'analytical_continuation')
        )
    
    def fit(self, iw_data: np.ndarray, iw_points: np.ndarray) -> np.ndarray:
        """Fit data on imaginary axis.
        
        Args:
            iw_data: Function values on imaginary axis
            iw_points: Imaginary frequency points
            
        Returns:
            Fitting parameters
        """
        # Placeholder for actual Padé fitting
        # Real implementation would solve linear system for Padé coefficients
        params = np.zeros((self.n_poles, 2))  # Poles and residues
        return params
    
    def evaluate(self, omega: Union[float, np.ndarray], params: np.ndarray) -> np.ndarray:
        """Evaluate continuation on real axis.
        
        Args:
            omega: Real frequency points
            params: Fitting parameters from fit()
            
        Returns:
            Function values on real axis
        """
        # Placeholder for Padé evaluation
        if isinstance(omega, float):
            omega = np.array([omega])
        return np.zeros_like(omega)


def create_ac_fitter(
    n_poles: int,
    regularization: float = 1e-8
) -> ACFitter:
    """Create analytical continuation fitter.
    
    Args:
        n_poles: Number of poles for rational approximation
        regularization: Regularization parameter (default: 1e-8)
        
    Returns:
        ACFitter object
        
    Examples:
        >>> # Create fitter with 16 poles
        >>> fitter = create_ac_fitter(16, regularization=1e-6)
        >>> 
        >>> # Fit imaginary axis data
        >>> iw_grid = create_frequency_grid('imaginary_axis', 32, omega_max=50.0)
        >>> iw_data = some_function(iw_grid.points)
        >>> params = fitter.fit(iw_data, iw_grid.points)
        >>> 
        >>> # Evaluate on real axis
        >>> w_real = np.linspace(-10, 10, 200)
        >>> real_data = fitter.evaluate(w_real, params)
    """
    config = _create_ac_fitter(n_poles, regularization)
    return ACFitter.from_config(config)


# Convenience functions for common use cases

def create_gw_frequency_grid(
    method: str = 'cd',
    n_freq: int = 32,
    omega_max: float = 100.0,
    energy_range: Optional[Tuple[float, float]] = None
) -> FrequencyGrid:
    """Create frequency grid optimized for GW calculations.
    
    Args:
        method: Frequency integration method
            - 'cd': Contour deformation (default)
            - 'ac': Analytical continuation
        n_freq: Number of frequency points
        omega_max: Maximum frequency for AC method
        energy_range: Energy range for CD method (default: (-50, 50))
        
    Returns:
        FrequencyGrid configured for GW
        
    Examples:
        >>> # Contour deformation grid
        >>> grid = create_gw_frequency_grid('cd', n_freq=64)
        >>> 
        >>> # Analytical continuation grid
        >>> grid = create_gw_frequency_grid('ac', n_freq=32, omega_max=50.0)
    """
    if method == 'cd':
        if energy_range is None:
            energy_range = (-50.0, 50.0)
        return create_frequency_grid(
            'contour_deformation',
            n_freq,
            bounds=energy_range
        )
    elif method == 'ac':
        return create_frequency_grid(
            'gauss_legendre',
            n_freq,
            omega_max=omega_max
        )
    else:
        raise ValueError(f"Unknown method: {method}. Use 'cd' or 'ac'")


def transform_to_imaginary(omega_real: np.ndarray) -> np.ndarray:
    """Transform real frequencies to imaginary axis.
    
    Args:
        omega_real: Real frequency points
        
    Returns:
        Imaginary frequency points (iω)
        
    Examples:
        >>> omega_real = np.array([1.0, 2.0, 3.0])
        >>> omega_imag = transform_to_imaginary(omega_real)
        >>> print(omega_imag)  # [1j, 2j, 3j]
    """
    return 1j * omega_real


def transform_to_matsubara(
    n_points: int,
    temperature: float,
    omega_max: Optional[float] = None
) -> FrequencyGrid:
    """Create Matsubara frequency grid for finite temperature.
    
    Args:
        n_points: Number of Matsubara frequencies
        temperature: Temperature in Hartree
        omega_max: Maximum frequency cutoff
        
    Returns:
        FrequencyGrid on Matsubara frequencies
        
    Examples:
        >>> # Room temperature grid
        >>> T = 300 / 27211.4  # 300K in Hartree
        >>> grid = transform_to_matsubara(64, T)
    """
    # Matsubara frequencies: ω_n = (2n+1)πT for fermions
    n_arr = np.arange(n_points)
    points = (2 * n_arr + 1) * np.pi * temperature
    
    if omega_max is not None:
        mask = points <= omega_max
        points = points[mask]
        n_points = len(points)
    
    # Weights for Matsubara sum (uniform)
    weights = np.ones(n_points) * 2 * temperature
    
    return FrequencyGrid(
        points=1j * points,  # Imaginary frequencies
        weights=weights,
        grid_type='matsubara',
        n_points=n_points,
        omega_max=omega_max,
        bounds=None
    )


# Module exports
__all__ = [
    'FrequencyGrid',
    'ContourDeformation',
    'ACFitter',
    'create_frequency_grid',
    'gauss_legendre_quadrature',
    'create_contour_deformation',
    'create_ac_fitter',
    'create_gw_frequency_grid',
    'transform_to_imaginary',
    'transform_to_matsubara',
]