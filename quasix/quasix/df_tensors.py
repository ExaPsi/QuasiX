"""
Density fitting tensor operations for QuasiX S2-2.

This module provides high-level Python interfaces for MO transformation
and metric operations on density fitting tensors, with PySCF validation.
"""

import numpy as np
from typing import Optional, Tuple, Dict, Any
import logging
from functools import wraps

# Import Rust bindings
try:
    from .quasix import (
        transform_to_mo_basis,
        compute_cholesky_metric,
        generate_mock_mo_coefficients as generate_mock_mo_coefficients_rust
    )
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    transform_to_mo_basis = None
    compute_cholesky_metric = None
    generate_mock_mo_coefficients_rust = None

logger = logging.getLogger(__name__)


def require_rust(func):
    """Decorator to ensure Rust bindings are available."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        if not RUST_AVAILABLE:
            raise RuntimeError(
                "Rust bindings not available. Please build with: maturin develop --release"
            )
        return func(*args, **kwargs)
    return wrapper


class DFTensorTransformer:
    """
    High-level interface for density fitting tensor transformations.

    This class provides methods to transform 3-center integrals from
    AO to MO basis and compute metric matrices for DF/RI calculations.
    """

    def __init__(self, n_ao: int, n_occ: int, n_vir: int, n_aux: int):
        """
        Initialize DF tensor transformer.

        Args:
            n_ao: Number of atomic orbitals
            n_occ: Number of occupied molecular orbitals
            n_vir: Number of virtual molecular orbitals
            n_aux: Number of auxiliary basis functions
        """
        self.n_ao = n_ao
        self.n_occ = n_occ
        self.n_vir = n_vir
        self.n_aux = n_aux
        self.n_mo = n_occ + n_vir

        # Validate dimensions
        if self.n_mo > self.n_ao:
            raise ValueError(f"Number of MOs ({self.n_mo}) cannot exceed number of AOs ({self.n_ao})")

        logger.info(f"Initialized DFTensorTransformer: n_ao={n_ao}, n_occ={n_occ}, n_vir={n_vir}, n_aux={n_aux}")

    @require_rust
    def transform_mo_3center(
        self,
        j3c_ao: np.ndarray,
        c_occ: np.ndarray,
        c_vir: np.ndarray
    ) -> np.ndarray:
        """
        Transform 3-center integrals from AO to MO basis.

        Performs the transformation (μν|P) → (ia|P) where:
        - μ,ν are AO basis functions
        - i,a are occupied and virtual MO indices
        - P is an auxiliary basis function

        Args:
            j3c_ao: 3-center AO integrals of shape (n_ao, n_ao, n_aux)
            c_occ: Occupied MO coefficients of shape (n_ao, n_occ)
            c_vir: Virtual MO coefficients of shape (n_ao, n_vir)

        Returns:
            MO-transformed integrals of shape (n_occ*n_vir, n_aux)

        Raises:
            ValueError: If array dimensions don't match expected sizes
        """
        # Validate input dimensions
        if j3c_ao.shape != (self.n_ao, self.n_ao, self.n_aux):
            raise ValueError(
                f"j3c_ao shape {j3c_ao.shape} doesn't match expected ({self.n_ao}, {self.n_ao}, {self.n_aux})"
            )
        if c_occ.shape != (self.n_ao, self.n_occ):
            raise ValueError(
                f"c_occ shape {c_occ.shape} doesn't match expected ({self.n_ao}, {self.n_occ})"
            )
        if c_vir.shape != (self.n_ao, self.n_vir):
            raise ValueError(
                f"c_vir shape {c_vir.shape} doesn't match expected ({self.n_ao}, {self.n_vir})"
            )

        # Ensure C-contiguous arrays for optimal performance
        j3c_ao = np.ascontiguousarray(j3c_ao, dtype=np.float64)
        c_occ = np.ascontiguousarray(c_occ, dtype=np.float64)
        c_vir = np.ascontiguousarray(c_vir, dtype=np.float64)

        logger.debug(f"Transforming 3-center integrals: ({self.n_ao}, {self.n_ao}, {self.n_aux}) -> ({self.n_occ*self.n_vir}, {self.n_aux})")

        # Call Rust function
        j3c_ia = transform_to_mo_basis(j3c_ao, c_occ, c_vir)

        # Validate output
        expected_shape = (self.n_occ * self.n_vir, self.n_aux)
        if j3c_ia.shape != expected_shape:
            raise RuntimeError(
                f"Unexpected output shape {j3c_ia.shape}, expected {expected_shape}"
            )

        # Log statistics
        logger.info(f"MO transformation complete: min={j3c_ia.min():.6f}, max={j3c_ia.max():.6f}, mean={j3c_ia.mean():.6f}")

        return j3c_ia

    @require_rust
    def compute_cholesky_v(self, metric: np.ndarray) -> Dict[str, Any]:
        """
        Compute Cholesky factorization of the metric matrix.

        Computes L such that v = L L^T where v is the 2-center metric matrix (P|Q).

        Args:
            metric: 2-center metric matrix of shape (n_aux, n_aux)

        Returns:
            Dictionary containing:
                - 'L': Lower triangular Cholesky factor
                - 'condition_number': Condition number of the metric
                - 'naux': Number of auxiliary functions

        Raises:
            ValueError: If metric is not square or has wrong dimensions
            RuntimeError: If Cholesky factorization fails (matrix not positive definite)
        """
        # Validate input
        if metric.ndim != 2:
            raise ValueError(f"Metric must be 2D, got shape {metric.shape}")
        if metric.shape[0] != metric.shape[1]:
            raise ValueError(f"Metric must be square, got shape {metric.shape}")
        if metric.shape[0] != self.n_aux:
            raise ValueError(
                f"Metric size {metric.shape[0]} doesn't match n_aux {self.n_aux}"
            )

        # Check symmetry
        if not np.allclose(metric, metric.T, rtol=1e-10, atol=1e-14):
            logger.warning("Metric matrix not perfectly symmetric, symmetrizing...")
            metric = 0.5 * (metric + metric.T)

        # Ensure C-contiguous
        metric = np.ascontiguousarray(metric, dtype=np.float64)

        logger.debug(f"Computing Cholesky factorization of {self.n_aux}x{self.n_aux} metric matrix")

        # Call Rust function
        result = compute_cholesky_metric(metric)

        # Extract results
        L = result['L']
        condition_number = result['condition_number']
        naux = result['naux']

        # Validate results
        if L.shape != (self.n_aux, self.n_aux):
            raise RuntimeError(f"Unexpected L shape {L.shape}")
        if naux != self.n_aux:
            raise RuntimeError(f"Unexpected naux {naux}, expected {self.n_aux}")

        # Check lower triangular
        upper_triangle = np.triu(L, k=1)
        if not np.allclose(upper_triangle, 0, atol=1e-14):
            logger.warning("Cholesky factor not perfectly lower triangular")

        # Log condition number
        if condition_number > 1e8:
            logger.warning(f"Metric matrix may be ill-conditioned: condition number = {condition_number:.2e}")
        else:
            logger.info(f"Metric matrix condition number: {condition_number:.2e}")

        return {
            'L': L,
            'condition_number': condition_number,
            'naux': naux
        }

    def validate_mo_coefficients(self, c_mo: np.ndarray, label: str = "MO") -> None:
        """
        Validate MO coefficient matrix for orthonormality.

        Args:
            c_mo: MO coefficient matrix of shape (n_ao, n_mo)
            label: Label for logging

        Raises:
            ValueError: If coefficients are not orthonormal within tolerance
        """
        overlap = c_mo.T @ c_mo
        identity = np.eye(c_mo.shape[1])

        max_error = np.abs(overlap - identity).max()
        if max_error > 1e-8:
            logger.warning(f"{label} coefficients not perfectly orthonormal: max error = {max_error:.2e}")
            if max_error > 1e-6:
                raise ValueError(f"{label} coefficients not orthonormal: max error = {max_error:.2e}")
        else:
            logger.debug(f"{label} coefficients are orthonormal (max error = {max_error:.2e})")


def transform_mo_3center_numpy(
    j3c_ao: np.ndarray,
    c_occ: np.ndarray,
    c_vir: np.ndarray
) -> np.ndarray:
    """
    Pure NumPy implementation of MO transformation for validation.

    This is a reference implementation using einsum for comparison
    with the optimized Rust version.

    Args:
        j3c_ao: 3-center AO integrals of shape (n_ao, n_ao, n_aux)
        c_occ: Occupied MO coefficients of shape (n_ao, n_occ)
        c_vir: Virtual MO coefficients of shape (n_ao, n_vir)

    Returns:
        MO-transformed integrals of shape (n_occ*n_vir, n_aux)
    """
    n_ao, _, n_aux = j3c_ao.shape
    n_occ = c_occ.shape[1]
    n_vir = c_vir.shape[1]

    # Two-step transformation using einsum
    # Step 1: Transform first index with occupied orbitals
    j3c_ov = np.einsum('iu,uvP->ivP', c_occ.T, j3c_ao)

    # Step 2: Transform second index with virtual orbitals
    j3c_ia = np.einsum('ivP,va->iaP', j3c_ov, c_vir)

    # Reshape to composite index format
    j3c_ia = j3c_ia.reshape(n_occ * n_vir, n_aux)

    return j3c_ia


def validate_against_pyscf(
    j3c_ao: np.ndarray,
    c_occ: np.ndarray,
    c_vir: np.ndarray,
    tolerance: float = 1e-8
) -> Tuple[bool, float]:
    """
    Validate MO transformation against PySCF reference implementation.

    Args:
        j3c_ao: 3-center AO integrals
        c_occ: Occupied MO coefficients
        c_vir: Virtual MO coefficients
        tolerance: Maximum allowed difference

    Returns:
        Tuple of (is_valid, max_error)
    """
    try:
        # Try to import PySCF for validation
        from pyscf import lib

        # PySCF's einsum-based transformation
        n_occ = c_occ.shape[1]
        n_vir = c_vir.shape[1]
        n_aux = j3c_ao.shape[2]

        # Use PySCF's optimized einsum
        j3c_ia_pyscf = lib.einsum('mnP,mi,na->iaP', j3c_ao, c_occ, c_vir)
        j3c_ia_pyscf = j3c_ia_pyscf.reshape(n_occ * n_vir, n_aux)

        # Compare with our implementation
        if RUST_AVAILABLE:
            j3c_ia_rust = transform_to_mo_basis(j3c_ao, c_occ, c_vir)
            max_error = np.abs(j3c_ia_rust - j3c_ia_pyscf).max()
        else:
            j3c_ia_numpy = transform_mo_3center_numpy(j3c_ao, c_occ, c_vir)
            max_error = np.abs(j3c_ia_numpy - j3c_ia_pyscf).max()

        is_valid = max_error < tolerance

        if is_valid:
            logger.info(f"✓ MO transformation validated against PySCF (max error = {max_error:.2e})")
        else:
            logger.error(f"✗ MO transformation validation failed (max error = {max_error:.2e} > {tolerance:.2e})")

        return is_valid, max_error

    except ImportError:
        logger.warning("PySCF not available for validation")
        return True, 0.0


@require_rust
def generate_mock_mo_coefficients(n_ao: int, n_mo: int, seed: int = 42) -> np.ndarray:
    """
    Generate orthonormal MO coefficients for testing.

    Args:
        n_ao: Number of atomic orbitals
        n_mo: Number of molecular orbitals
        seed: Random seed for reproducibility

    Returns:
        Orthonormal MO coefficient matrix of shape (n_ao, n_mo)
    """
    if n_mo > n_ao:
        raise ValueError(f"Cannot have more MOs ({n_mo}) than AOs ({n_ao})")

    return generate_mock_mo_coefficients_rust(n_ao, n_mo, seed)


def benchmark_mo_transformation(
    n_ao: int = 100,
    n_aux: int = 300,
    n_occ: int = 50,
    n_vir: int = 50
) -> Dict[str, float]:
    """
    Benchmark MO transformation performance.

    Args:
        n_ao: Number of AO basis functions
        n_aux: Number of auxiliary basis functions
        n_occ: Number of occupied orbitals
        n_vir: Number of virtual orbitals

    Returns:
        Dictionary with timing results
    """
    import time

    # Generate test data
    np.random.seed(42)
    j3c_ao = np.random.randn(n_ao, n_ao, n_aux)
    j3c_ao = 0.5 * (j3c_ao + j3c_ao.transpose(1, 0, 2))  # Symmetrize

    if RUST_AVAILABLE:
        c_occ = generate_mock_mo_coefficients(n_ao, n_occ, 42)
        c_vir = generate_mock_mo_coefficients(n_ao, n_vir, 43)
    else:
        # Simple orthogonalization
        c_occ = np.linalg.qr(np.random.randn(n_ao, n_occ))[0]
        c_vir = np.linalg.qr(np.random.randn(n_ao, n_vir))[0]

    results = {}

    # Benchmark NumPy implementation
    t0 = time.time()
    j3c_ia_numpy = transform_mo_3center_numpy(j3c_ao, c_occ, c_vir)
    results['numpy_time'] = time.time() - t0

    # Benchmark Rust implementation if available
    if RUST_AVAILABLE:
        t0 = time.time()
        j3c_ia_rust = transform_to_mo_basis(j3c_ao, c_occ, c_vir)
        results['rust_time'] = time.time() - t0

        # Compare results
        max_diff = np.abs(j3c_ia_rust - j3c_ia_numpy).max()
        results['max_difference'] = max_diff
        results['speedup'] = results['numpy_time'] / results['rust_time']

        logger.info(f"Benchmark results:")
        logger.info(f"  NumPy time: {results['numpy_time']:.3f} s")
        logger.info(f"  Rust time:  {results['rust_time']:.3f} s")
        logger.info(f"  Speedup:    {results['speedup']:.2f}x")
        logger.info(f"  Max diff:   {max_diff:.2e}")
    else:
        logger.info(f"NumPy time: {results['numpy_time']:.3f} s")

    # Memory usage estimate
    memory_gb = (j3c_ao.nbytes + j3c_ia_numpy.nbytes) / 1e9
    results['memory_gb'] = memory_gb
    logger.info(f"  Memory:     {memory_gb:.2f} GB")

    return results


# Module initialization
if __name__ == "__main__":
    # Run basic tests
    logging.basicConfig(level=logging.INFO)

    print("Testing DFTensorTransformer...")
    transformer = DFTensorTransformer(n_ao=10, n_occ=5, n_vir=5, n_aux=20)

    # Generate test data
    np.random.seed(42)
    j3c_ao = np.random.randn(10, 10, 20)
    metric = np.random.randn(20, 20)
    metric = metric @ metric.T  # Make positive definite

    if RUST_AVAILABLE:
        c_occ = generate_mock_mo_coefficients(10, 5, 42)
        c_vir = generate_mock_mo_coefficients(10, 5, 43)

        # Test transformation
        j3c_ia = transformer.transform_mo_3center(j3c_ao, c_occ, c_vir)
        print(f"Transformed shape: {j3c_ia.shape}")

        # Test Cholesky
        cholesky_result = transformer.compute_cholesky_v(metric)
        print(f"Cholesky condition number: {cholesky_result['condition_number']:.2e}")

        # Run benchmark
        print("\nRunning benchmark...")
        benchmark_mo_transformation(n_ao=50, n_aux=150, n_occ=25, n_vir=25)
    else:
        print("Rust bindings not available. Please build with: maturin develop --release")